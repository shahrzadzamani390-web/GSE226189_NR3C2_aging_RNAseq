# Quality control summary

## PCA

- PCA was performed using variance-stabilized gene counts from all 82 samples.
- PC1 explained 57% of the variance.
- PC2 explained 13% of the variance.
- Samples showed clear structure, but PC1 and PC2 were not strongly correlated with continuous age.
- No single extreme outlier sample was obvious from the PCA.

## Library size

- Total gene counts ranged from about 73 million to 371 million counts per sample.
- Library size varied across samples, but the lowest-depth samples still had substantial counts.
- DESeq2 size-factor normalization was used to account for library size differences.

## Sample distance heatmap

- The sample distance heatmap showed several expression clusters.
- Clustering was not explained completely by age group or sex.
- No sample was removed based on QC.

## Interpretation

The dataset shows strong sample-level structure, likely reflecting biological donor variation and/or hidden technical factors. Because available metadata does not fully explain this structure, downstream results should be interpreted carefully and ideally validated in an independent dataset.
