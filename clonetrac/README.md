# Clonetrac analysis

Open `demo.ipynb` to run the seven original NB01–NB07 tree/profile pairs in this
folder. The completed run uses balanced peak ranking with weight `1.0`,
`max_k_limit=0.70`, and the DP's `prefer_fewer=False`. No explicit outlier size
threshold is supplied. Outputs are in `phytclust_outputs_balanced_weight1`.
The executed notebook is also saved as `demo.balanced_weight1.ipynb`; earlier
notebooks and trial outputs are preserved.

To reproduce the current run from the repository folder:

```sh
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/phytclust-mpl \
  .venv/bin/python clonetrac/run_analysis.py \
  --max-k-limit 0.70 --peak-partition-preference balanced \
  --peak-partition-weight 1.0 \
  --output-dir clonetrac/phytclust_outputs_balanced_weight1
```

Use `--sample NB03` to run one sample or `--top-n 1` to save only the first
ranked partition. The maximum k is `ceil(0.70 * tumour_cells)`. Weight `1.0`
ranks detected peaks by partition quality alone, penalizing uneven cluster
sizes and singleton cells. Score values and detected candidates stay unchanged.
NB01 has no detected peaks, so its saved k = 2 partition uses the highest-score
fallback and is marked in `summary.tsv`.

Only top-level `NB*_final_tree.new` files are included. Files in `trees` and
`old_results`, previously pruned copies, and the unpaired `polytomy.nwk` example
are excluded. Original inputs are preserved. `demo.before_batch.ipynb` contains
the previous notebook.

Normal-cell identifiers are read from the supplied TSV and text lists, which
must agree. Names are matched after removing `CTR_` and `.final.bam`. Normal
cells are removed from the tree and profile table. Internal-node profiles stay
in the filtered table but are not plotted. The `diploid` reference is retained
in filtered inputs and excluded from clustering and plotted tumour rows.

For NB01, NB02, and NB05, trees labelled `with_ecdna` are paired with the supplied
profiles labelled `without_ecdna`. The notebook uses the `major` column as total
copy number, matching the previous demo. Segments are ordered by chromosome and
position and drawn in proportion to their measured lengths; unmeasured gaps are
omitted. A missing cell or segment stops that sample instead of silently omitting
rows. Tree, cluster strip, and profile rows use identical leaf coordinates.

Results are in `phytclust_outputs/<sample>/`:

- Filtered Newick tree, filtered profile TSV, and removed normal-cell list.
- `scores.tsv`, `scores.png`, and `scores.pdf`.
- `run_summary.json` and `all_selected_assignments.tsv`.
- A `rank_<rank>_k<count>` folder per selected partition, with PhytClust assignment
  and alpha TSVs, cell assignments, a cluster tree PNG, and aligned tree/CNP PNG
  and PDF plots. `plot_row_order.tsv` records the exact cell order in the figures.

`phytclust_outputs/summary.tsv` summarizes the completed samples. Cluster IDs are
kept consistent between assignment tables and plots; singleton clusters are not
relabelled as outliers during export.

Many input branches have zero length. PhytClust preserves these lengths and does
not split zero-length branches. If no peak is found, the highest finite score is
saved and explicitly marked as a fallback. NB01 used this fallback at k = 2 in
the saved run. Invalid scores are retained in the TSV; the score plotting function
omits values that cannot be displayed on its logarithmic axis.

## Fewer-outliers trial

`demo.fewer_outliers.ipynb` contains a separate completed run with
`OutlierConfig(size_threshold=3, prefer_fewer=True)`. At each k, the algorithm
prioritizes fewer clusters with fewer than three cells, then partition cost.
Peak selection uses the resulting scores. This does not force equal cluster
sizes and does not guarantee fewer tiny clusters when comparing different k.
The original run remains in `demo.original_settings.ipynb` and `phytclust_outputs`.

To repeat this trial:

```sh
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/phytclust-mpl \
  .venv/bin/python clonetrac/run_analysis.py \
  --prefer-fewer-outliers --outlier-size-threshold 3 \
  --output-dir clonetrac/phytclust_outputs_fewer_outliers
```

The trial's output folder contains `comparison.tsv` for the first ranked
partition and `comparison_all_ranks.tsv` for every selected partition.
Cluster-size variation is measured as standard deviation divided by mean
cluster size; lower values indicate more even sizes. Comparisons include the
selected k because partitions with different k have different size constraints.

## Earlier reduced-limit trial

`demo.fewer_outliers_70pct_maxk.ipynb` contains the completed
seven-sample run using `prefer_fewer=True`, `size_threshold=3`, and 70% of each
sample's previous maximum k, rounded down. The previous limits are read from
`phytclust_outputs_fewer_outliers/<sample>/run_summary.json`.

| Sample | Previous maximum k | New maximum k | Selected k, in rank order |
| --- | ---: | ---: | --- |
| NB01 | 164 | 114 | 2 (fallback) |
| NB02 | 164 | 114 | 9, 13, 6 |
| NB03 | 142 | 99 | 3, 43, 16 |
| NB04 | 117 | 81 | 27, 30, 66 |
| NB05 | 131 | 91 | 2, 4 |
| NB06 | 95 | 66 | 12, 28, 9 |
| NB07 | 94 | 65 | 45, 23, 52 |

Results, plots, the summary, and comparison tables are in
`phytclust_outputs_fewer_outliers_70pct_maxk`. Run all cells in either accepted
notebook to reproduce the analysis.

## Earlier singleton-outlier run

`demo.singleton_outliers_maxk70.ipynb` contains the corrected
completed run with `OutlierConfig(prefer_fewer=True)` and `max_k_limit=0.70`.
`size_threshold` remains unset. Singleton clusters count as outliers when
`prefer_fewer=True`; an explicit threshold still overrides that default.

The maximum k is `ceil(0.70 * tumour_cells)`. It is not a percentage of the
previous maximum. To reproduce this run from the repository folder:

```sh
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/phytclust-mpl \
  .venv/bin/python clonetrac/run_analysis.py \
  --prefer-fewer-outliers --max-k-limit 0.70 \
  --output-dir clonetrac/phytclust_outputs_singletons_maxk70
```

Results are in `phytclust_outputs_singletons_maxk70`; the previous trial folders
are preserved. `comparison.tsv` compares first-ranked partitions against the
original run, counting only singleton clusters as small clusters.

## NB06: balanced peak-ranking trial

This separate trial keeps the DP's `prefer_fewer=False`, uses `max_k_limit=0.70`,
and ranks detected peaks with `PeakConfig(partition_preference="balanced",
partition_weight=0.5)`. Its outputs are in `phytclust_outputs_balanced_peaks/NB06`.
Unlike the DP preference, this option leaves score values and candidate peaks
unchanged. `peak_ranking.tsv` records the original peak metrics and partition
quality. `weight_comparison.tsv` compares ranking weights 0.5, 0.75, and 1.0.

```sh
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/phytclust-mpl \
  .venv/bin/python clonetrac/run_analysis.py --sample NB06 \
  --max-k-limit 0.70 --peak-partition-preference balanced \
  --peak-partition-weight 0.5 \
  --output-dir clonetrac/phytclust_outputs_balanced_peaks
```
