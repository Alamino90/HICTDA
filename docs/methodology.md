# HiCTDA Methodological Specification

## 1. Input

The primary input is a `.hic` file containing observed intra-chromosomal Hi-C
contact data. The user specifies:

- chromosome;
- genomic start coordinate;
- genomic end coordinate;
- bin resolution;
- normalization method.

If start/end are omitted, the complete selected chromosome is analyzed.

## 2. Normalization

Normalization is mandatory at the high-level API. The current implementation
passes the requested normalization to `hic-straw`. Knight-Ruiz balancing (`KR`)
is the default.

## 3. Contact matrix

For a selected chromosome and genomic interval, the normalized observed
contact matrix is retrieved at the requested base-pair resolution.

NaN values are replaced with zero and negative numerical artifacts are clipped
to zero.

## 4. Contact-to-distance transformation

The default transformation is:

d_ij = 1 / (C_ij + epsilon)

where C_ij is the normalized contact frequency and epsilon prevents division by
zero.

This defines the dissimilarity matrix used by the Vietoris-Rips calculation.

## 5. Persistent homology

The distance matrix is supplied to Ripser as a precomputed distance matrix.
Persistent homology is calculated through H1 by default.

The H1 persistence lifetime of feature i is:

p_i = death_i - birth_i.

Features with infinite death time are excluded from the finite H1 summary
statistics.

## 6. Summary statistics

The package reports:

- number of finite H1 features;
- mean H1 persistence;
- median H1 persistence;
- standard deviation of H1 persistence.

## 7. Visualization

The standard output consists of:

1. normalized/log-transformed Hi-C contact map;
2. persistence diagram;
3. persistence barcode;
4. H1 persistence distribution;
5. H1 persistence box plot;
6. combined summary figure.

## 8. Reproducibility

All numerical parameters and analysis metadata are written to `summary.json`.
Numerical arrays can additionally be saved as `.npy` files.
