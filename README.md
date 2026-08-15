# HICTDA

**HICTDA** is a Python package for persistent-homology analysis of Hi-C
chromatin contact maps stored in `.hic` format and identifying Topologically Associating Domains
(TADs).

## Author

**Alamin Mustafa**, MBBS
Faculty of Medicine, Al-Neelain University 
Khartoum, Sudan. ORCID: https://orcid.org/0000-0003-1129-6284

## Citation

If you use HICTDA in your research, please cite: Alamin Mustafa. (2026). Alamino90/HICTDA: HICTDA v0.1.1 (Version v0.1.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21959829

### BibTeX:
bibtex
@software{mustafa2026hictda,
  author       = {Mustafa, Alamin},
  title        = {Alamino90/HICTDA: HICTDA v0.1.1},
  year         = {2026},
  version      = {v0.1.1},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21959829},
  url          = {https://doi.org/10.5281/zenodo.21959829}
}

### The package provides:

- chromosome and genomic-region selection;
- configurable Hi-C bin resolution;
- mandatory matrix normalization;
- conversion of contact frequencies to dissimilarities;
- Vietoris-Rips persistent homology through dimension H1;
- H1 feature count, mean, median, and standard deviation;
- contact-map visualization;
- persistence diagram;
- persistence barcode;
- H1 persistence distribution;
- H1 persistence box plot;
- machine-readable CSV/JSON output;
- a single high-level analysis function and a command-line interface.

## Installation

Create a virtual environment and install:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .
```

For development:

```bash
pip install -e ".[dev]"
```

## Basic Python usage

```python
from hictda import analyze_hic

results = analyze_hic(
    hic_file="ENCFF718AWL.hic",
    chromosome="3",
    start=140_000_000,
    end=190_000_000,
    resolution=50_000,
    normalization="KR",
    output_dir="results"
)

print(results.statistics)
```

## Entire chromosome

If `start` and `end` are omitted, the package analyzes the complete selected
chromosome, subject to the size and memory limits of the machine.

```python
results = analyze_hic(
    hic_file="sample.hic",
    chromosome="3",
    resolution=50_000,
    normalization="KR"
)
```

## Single-chromosome `.hic` files

A `.hic` file does not need to contain an entire genome. The selected
chromosome is validated against the chromosomes actually present in the file.

## CLI

```bash
hictda analyze sample.hic \
  --chromosome 3 \
  --start 140000000 \
  --end 190000000 \
  --resolution 50000 \
  --normalization KR \
  --output results
```

## Outputs

The output directory contains:

```text
contact_map.png
persistence_diagram.png
persistence_barcode.png
h1_distribution.png
h1_boxplot.png
hictda_summary.png
h1_persistence.csv
summary.json
```

## Scientific note

The default contact-to-dissimilarity transformation is:

\[
d_{ij} = \frac{1}{C_{ij}+\epsilon}.
\]

This maps larger contact frequency to smaller dissimilarity. The
transformation is exposed as a parameter so that alternative transformations
can be added without changing the rest of the pipeline.

The package currently uses the supplied Hi-C normalization method through
`hic-straw`; the default and recommended method for this implementation is
Knight-Ruiz (`KR`).

## Reproducibility

The package records the input file, chromosome, genomic interval, resolution,
normalization, distance method, epsilon, and persistent-homology parameters
in `summary.json`.

### Associated Publication

A manuscript describing HICTDA: A Persistent Homology Framework for Multiscale Identification of Topologically Associating Domains (TADs) in Hi-C Data is currently in preparation and will be published soon.

### License
MIT
