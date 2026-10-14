# RELION half-set regrouping

Reassign RELION particle half-sets by **tomogram** to investigate FSC artifacts associated with overlapping sampling in cryo-ET / subtomogram averaging (STA).

## The problem

Densely sampled tubular structures contained neighboring segments and nearly duplicate particles assigned to different half-sets. These particles can share local data and noise, compromising half-map independence and producing elevated FSC values or residual positive correlation at high frequencies.

Inspection confirmed that these assignments already existed in the upstream RELION particle tables. Warp–M preserved the imported labels. Ordinary non-helical refinement assigns particles randomly when half-set labels are absent. RELION retains existing `_rlnRandomSubset` labels; adding helical refinement parameters later does not automatically regroup the particles.

**Having particles from one tomogram in both half-sets does not, by itself, establish contamination. The concern in this case was overlapping segments and nearly duplicate particles sharing data across half-sets.**

## Approach and outcome

Assign all particles from each tomogram to the same half-set while keeping the particle counts as balanced as possible. The script changes or adds only `_rlnRandomSubset`, preserving coordinates, orientations, particle order, optics tables, and other content.

The user reported that the FSC issue was resolved after applying this regrouping script and using the resulting assignments in subsequent processing. This is a case-specific, user-reported outcome, not a general guarantee for all STA datasets. No independently verified resolution values or before-and-after curves are provided here. （human: GPT seems not happy cause i didn't show it the results, I have since shown it the FSC curve and confirmed the effectiveness. ;) )

Use the regrouped STAR file for a new reconstruction and independent refinement. Changing labels does not modify existing half-maps or automatically remove alignment bias inherited from earlier processing.

## Usage

Requires only the Python standard library; no third-party packages are needed.

```bash
# Preview grouping statistics without writing a file
python3 regroup_halfsets_by_tomogram.py Refine3D/jobXXX/run_data.star --dry-run

# Write the regrouped STAR file
python3 regroup_halfsets_by_tomogram.py Refine3D/jobXXX/run_data.star
```

You can also pass a refine job directory directly. By default, the output is written beside the input file:

```text
Refine3D/jobXXX/run_data_halfset_tomogram.star
```

- Uses `_rlnTomoName`, with `_rlnMicrographName` as a fallback for older formats. The source field must uniquely identify each actual tomogram.
- Uses a reproducible seed, defaulting to `--seed 0`. Tomograms are processed in descending particle-count order and assigned to the half-set with fewer particles; the seed controls ties.
- Never overwrites the input or an existing output. At least two tomograms are required.
- Checks particle counts, other particle fields, and grouping consistency before writing. The reported `new_mixed_half_tomograms` must be `0`.
- Large tomograms can limit how evenly the halves can be balanced. Grouping is by whole tomogram, not filament ID.

Example with a different seed and output filename:

```bash
python3 regroup_halfsets_by_tomogram.py Refine3D/jobXXX/run_data.star \
  --seed 42 \
  --output Refine3D/jobXXX/run_data_halfset_tomogram_seed42.star
```

Use `--tomo-column` to specify the source column or `--block` to select a STAR data block.

## References

- [RELION 5.0 half-set assignment implementation](https://github.com/3dem/relion/blob/ver5.0/src/exp_model.cpp#L244)
- [RELION STA duplicate-particle removal](https://relion.readthedocs.io/en/release-5.0/STA_tutorial/DuplicateParticles.html)

## Acknowledgments

Thanks to OpenAI Codex for assistance with inspecting half-set assignments, writing the regrouping script, and preparing this documentation. The user performed the actual data processing and reported the improvement in FSC behavior.
