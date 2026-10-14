#!/usr/bin/env python3
"""Assign whole tomograms to RELION half-sets and write a sibling STAR file.

Usage:
  python3 regroup_halfsets_by_tomogram.py Refine3D/job012/run_data.star
  python3 regroup_halfsets_by_tomogram.py Refine3D/job012/run_data.star --dry-run

Only rlnRandomSubset is changed (or added). Other tokens and comments are
preserved. Existing input/output files are never overwritten. This prepares
particle metadata for a new reconstruction/refinement; it does not change or
validate existing half maps. No third-party Python packages are required.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import random
import sys
import tempfile


class Token:
    __slots__ = ('value', 'start', 'end', 'quoted')

    def __init__(self, value, start, end, quoted=False):
        self.value = value
        self.start = start
        self.end = end
        self.quoted = quoted


class Loop:
    __slots__ = ('block', 'labels', 'values')

    def __init__(self, block, labels, values):
        self.block = block
        self.labels = labels
        self.values = values


def tokenize(text):
    """Read STAR tokens while retaining their exact source spans."""
    result = []
    i, size = 0, len(text)
    while i < size:
        if text[i].isspace() or (i == 0 and text[i] == '\ufeff'):
            i += 1
            continue
        if text[i] == '#':
            end = text.find('\n', i)
            i = size if end < 0 else end + 1
            continue
        start = i
        if text[i] == ';' and (i == 0 or text[i - 1] == '\n'):
            closing = text.find('\n;', i + 1)
            if closing < 0:
                raise ValueError('Unterminated STAR semicolon text field')
            i = closing + 2
            if i < size and not text[i].isspace() and text[i] != '#':
                raise ValueError('Invalid STAR semicolon field terminator')
            result.append(Token(text[start + 1:closing], start, i, True))
            continue
        if text[i] in "\"'":
            quote = text[i]
            i += 1
            content_start = i
            while i < size:
                if text[i] == quote and (i + 1 == size or text[i + 1].isspace()):
                    result.append(Token(text[content_start:i], start, i + 1, True))
                    i += 1
                    break
                i += 1
            else:
                raise ValueError('Unterminated STAR quoted field')
            continue
        while i < size and not text[i].isspace():
            i += 1
        result.append(Token(text[start:i], start, i))
    return result


def is_control(token):
    if token.quoted:
        return False
    value = token.value.lower()
    return (value in ('loop_', 'stop_', 'global_') or value.startswith('data_')
            or value.startswith('save_') or value.startswith('_'))


def read_loops(text):
    tokens = tokenize(text)
    loops = []
    block = ''
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if not token.quoted and token.value.lower().startswith('data_'):
            block = token.value[5:]
        if token.quoted or token.value.lower() != 'loop_':
            i += 1
            continue
        i += 1
        labels = []
        while i < len(tokens) and not tokens[i].quoted and tokens[i].value.startswith('_'):
            labels.append(tokens[i])
            i += 1
        if not labels:
            raise ValueError('STAR loop without column labels')
        values = []
        while i < len(tokens):
            if is_control(tokens[i]):
                break
            values.append(tokens[i])
            i += 1
        if len(values) % len(labels):
            raise ValueError('Incomplete row in data_' + block)
        loops.append(Loop(block, labels, values))
    return loops


def select_loop(text, column=None, block=None):
    candidates = []
    for loop in read_loops(text):
        names = [label.value for label in loop.labels]
        source = column or next((name for name in ('_rlnTomoName', '_rlnMicrographName')
                                 if name in names), None)
        if source in names and (block is None or loop.block == block):
            candidates.append((loop, source))
    if len(candidates) != 1:
        raise ValueError('Expected one particle loop with rlnTomoName or rlnMicrographName; '
                         'found {}. Use --block/--tomo-column if needed.'.format(len(candidates)))
    loop, source = candidates[0]
    names = [label.value for label in loop.labels]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate column labels in the selected particle loop')
    if not loop.values:
        raise ValueError('Particle table is empty')
    return loop, source


def regroup(text, seed, column=None, block=None):
    loop, source = select_loop(text, column, block)
    names = [label.value for label in loop.labels]
    width = len(names)
    source_index = names.index(source)
    half_index = names.index('_rlnRandomSubset') if '_rlnRandomSubset' in names else None
    rows = [loop.values[i:i + width] for i in range(0, len(loop.values), width)]
    groups = defaultdict(list)
    old_counts = Counter()
    old_groups = defaultdict(set)
    for index, row in enumerate(rows):
        key = row[source_index].value
        if key in ('', '.', '?') or '\n' in key or '\r' in key:
            raise ValueError('Missing or invalid tomogram identity at particle row {}'.format(index + 1))
        groups[key].append(index)
        if half_index is not None:
            half = row[half_index].value
            if half not in ('0', '1', '2'):
                raise ValueError('Invalid rlnRandomSubset at particle row {}'.format(index + 1))
            old_counts[half] += 1
            old_groups[key].add(half)
    if len(groups) < 2:
        raise ValueError('At least two tomograms are required for two nonempty half-sets')

    # Largest groups first minimizes imbalance; the seed controls ties.
    rng = random.Random(seed)
    keys = sorted(groups)
    rng.shuffle(keys)
    keys.sort(key=lambda key: len(groups[key]), reverse=True)
    counts = [0, 0]
    mapping = {}
    for key in keys:
        side = rng.randrange(2) if counts[0] == counts[1] else int(counts[0] > counts[1])
        mapping[key] = side + 1
        counts[side] += len(groups[key])

    edits = []
    changed = 0
    newline = '\r\n' if '\r\n' in text else '\n'
    if half_index is None:
        last_label = loop.labels[-1]
        header_end = text.find('\n', last_label.end)
        if header_end < 0 or loop.values[0].start < header_end:
            raise ValueError('Cannot add a column to this nonstandard STAR header layout')
        edits.append((header_end + 1, header_end + 1,
                      '_rlnRandomSubset #{}{}'.format(width + 1, newline)))
    for row in rows:
        half = str(mapping[row[source_index].value])
        if half_index is None:
            edits.append((row[-1].end, row[-1].end, ' ' + half))
        else:
            token = row[half_index]
            changed += token.value != half
            edits.append((token.start, token.end, half))
    pieces = []
    previous = 0
    for start, end, value in sorted(edits, key=lambda edit: edit[0]):
        pieces.extend((text[previous:start], value))
        previous = end
    pieces.append(text[previous:])
    output = ''.join(pieces)

    # Check the assembled file before publishing it; also protect all other fields.
    check_loop, check_source = select_loop(output, source, loop.block)
    check_names = [label.value for label in check_loop.labels]
    check_width = len(check_names)
    check_rows = [check_loop.values[i:i + check_width]
                  for i in range(0, len(check_loop.values), check_width)]
    if len(check_rows) != len(rows):
        raise ValueError('Output particle count changed')
    check_index = {name: i for i, name in enumerate(check_names)}
    check_groups = defaultdict(set)
    check_counts = Counter()
    for before, after in zip(rows, check_rows):
        for i, name in enumerate(names):
            if name != '_rlnRandomSubset' and before[i].value != after[check_index[name]].value:
                raise ValueError('A non-half-set particle field changed')
        key = after[check_index[check_source]].value
        half = after[check_index['_rlnRandomSubset']].value
        check_groups[key].add(half)
        check_counts[half] += 1
    if any(len(halves) != 1 for halves in check_groups.values()):
        raise ValueError('Output contains a tomogram split between half-sets')
    if check_counts != Counter({'1': counts[0], '2': counts[1]}):
        raise ValueError('Output half-set counts do not match the assignment')
    summary = {
        'group_by': source, 'block': loop.block, 'seed': seed,
        'particles': len(rows), 'tomograms': len(groups),
        'old_half_counts': dict(old_counts),
        'old_mixed_half_tomograms': sum(len(halves) > 1 for halves in old_groups.values())
        if half_index is not None else None,
        'new_half_counts': {'1': counts[0], '2': counts[1]},
        'new_tomogram_counts': dict(Counter(str(half) for half in mapping.values())),
        'new_mixed_half_tomograms': 0,
        'changed_half_labels': changed if half_index is not None else len(rows),
        'added_half_column': half_index is None,
        'largest_tomogram_particles': max(map(len, groups.values())),
        'other_particle_fields_preserved': True,
    }
    return output, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', type=Path, help='RELION run_data.star or its refine job directory')
    parser.add_argument('--seed', type=int, default=0, help='Reproducible grouping seed (default: 0)')
    parser.add_argument('--output', type=Path, help='New STAR path; must be beside the input file')
    parser.add_argument('--tomo-column', help='Override source column, e.g. rlnTomoName')
    parser.add_argument('--block', help='Select a data block if more than one particle loop exists')
    parser.add_argument('--dry-run', action='store_true', help='Print the proposed grouping without writing')
    args = parser.parse_args()
    try:
        requested = args.input.expanduser()
        if requested.is_dir():
            requested = requested / 'run_data.star'
        input_path = requested.resolve(strict=True)
        if not input_path.is_file():
            raise ValueError('Input is not a regular file')
        if requested.absolute() != input_path:
            # Resolve ordinary relative paths, but reject a symlink that silently
            # moves the output away from the requested job's directory.
            if requested.is_symlink() or requested.parent.resolve() != input_path.parent:
                raise ValueError('Input symlink changes the requested job location; use its real path explicitly')
        raw = input_path.read_bytes()
        text = raw.decode('utf-8')
        column = args.tomo_column
        if column and not column.startswith('_'):
            column = '_' + column
        block = args.block
        if block and block.startswith('data_'):
            block = block[5:]
        output, summary = regroup(text, args.seed, column, block)
        if args.output is None:
            output_path = input_path.with_name(input_path.stem + '_halfset_tomogram.star')
        else:
            proposed = args.output.expanduser()
            output_path = (Path.cwd() / proposed if not proposed.is_absolute() else proposed)
            output_path = output_path.parent.resolve(strict=True) / output_path.name
        if output_path.parent != input_path.parent:
            raise ValueError('Output must be in the same directory as the input STAR')
        if output_path.exists() or output_path.is_symlink():
            raise ValueError('Output already exists; choose another --output name: ' + str(output_path))
        payload = output.encode('utf-8')
        summary.update(input=str(input_path), output=str(output_path),
                       input_sha256=hashlib.sha256(raw).hexdigest(),
                       output_sha256=hashlib.sha256(payload).hexdigest(),
                       dry_run=args.dry_run)
        if not args.dry_run:
            # Publish by an exclusive hard link: an existing output cannot be
            # overwritten, even if another process creates it during this run.
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(prefix='.' + output_path.name + '.',
                                                 suffix='.tmp', dir=str(input_path.parent),
                                                 delete=False) as handle:
                    temporary = Path(handle.name)
                    handle.write(payload)
                    handle.flush()
                    os.fsync(handle.fileno())
                if input_path.read_bytes() != raw:
                    raise ValueError('Input changed during processing; refusing to publish')
                os.link(str(temporary), str(output_path))
            finally:
                if temporary is not None:
                    temporary.unlink()
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return 0
    except (OSError, UnicodeError, ValueError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
