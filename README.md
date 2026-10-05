# Arnav

**Arnav** (Analysis of RNA variants) is a lightweight statistical method for
detecting mutations from RNA sequencing data, implemented in C++.

> Ramu A. *Arnav: Site specific error models to identify variants in RNA.*
> bioRxiv (2018). [doi:10.1101/397539](https://www.biorxiv.org/content/10.1101/397539v1)

## Background

RNA-seq data are an attractive source for finding somatic and post-zygotic
mutations because large public cohorts such as
[GTEx](https://gtexportal.org/) sequence many tissues from many people. Calling
variants from RNA is hard, though. Allele-specific expression means the true
variant allele fraction (VAF) of a mutation is unknown, and RNA-specific
artifacts (RNA editing, mapping errors near splice junctions, library
preparation errors) create recurrent false positives at particular positions.

Arnav handles this by learning a **site-specific error model** from a panel of
samples. For every position, the background rate of non-reference reads is
estimated from all *other* samples in the panel. A sample's alternate-allele
count is then tested against that background with a one-sided binomial test,
and p-values are corrected for multiple testing. Sites that are noisy across the panel are penalized automatically,
without assuming a particular VAF.

Arnav was used to identify post-zygotic mutations from GTEx RNA-seq data in
Rockweiler *et al.*, *The origins and functional effects of postzygotic
mutations throughout the human life span*,
[Science (2023)](https://www.science.org/doi/10.1126/science.abn7113).

## Install

Requires a C++11 compiler and CMake. The [fdr](https://github.com/gatoravi/fdr)
dependency is a git submodule, so clone recursively:

```bash
git clone --recursive https://github.com/gatoravi/arnav
cd arnav
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=release
make
```

Arnav bundles the standalone
[Rmath](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library)
library (for `pbinom`) and [cereal](https://uscilab.github.io/cereal/) (to
serialize the site-specific priors to disk), so there are no other
dependencies to install.

## Workflow

### 1. Generate per-sample read counts

Arnav works on per-position read counts produced by
[mpileup2readcounts](https://github.com/gatoravi/mpileup2readcounts), compressed
with `gzip` or `bgzip`:

```bash
samtools mpileup -f ref.fa -l regions.bed sample1.bam \
    | mpileup2readcounts sample1 | bgzip > sample1_readcounts.gz
```

Then list every sample in a two-column, tab-separated file of sample name and
path to its read-count file:

```
sample1    sample1_readcounts.gz
sample2    sample2_readcounts.gz
```

### 2. Estimate site-specific priors

Build the background error model for a set of sites given in a BED file:

```bash
./arnav prior-dump-fixed samples.tsv merged_priors.dump regions.bed.gz
```

For large cohorts, priors can be computed in batches with `prior-dump` and
combined with `prior-merge`, which takes a two-column file of dump name and dump
path.

### 3. Call variants

```bash
./arnav call-using-merged samples.tsv merged_priors.dump
```

Each candidate needs at least two reads supporting the alternate allele. The
output contains the read counts for the site, the background reference and
alternate counts from the rest of the panel, the binomial p-value and the
sample name. Only sites passing the multiple-testing cutoff are reported (the
current code uses a fixed Bonferroni threshold; a Benjamini-Hochberg FDR option
is in `src/arnav.cc`).

A longer worked example, detecting putative somatic mutations in RNA-seq data,
is in
[this blog post](http://gatoravi.github.io/blog/genetics/2017/10/12/mutations-rnaseq.html).

### All commands

Run `./arnav` with no arguments to list every subcommand:

| Command | Purpose |
| --- | --- |
| `prior-and-call` | Estimate priors and call variants in one step |
| `prior-dump` | Estimate priors at all sites and save them to disk |
| `prior-dump-fixed` | Estimate priors only at sites in a BED file |
| `prior-merge` | Merge several prior dumps into one |
| `prior-print` | Print the contents of a prior dump |
| `call-using-merged` | Call variants using a merged prior dump |

## License

GPL-3.0. See [LICENSE](LICENSE).
