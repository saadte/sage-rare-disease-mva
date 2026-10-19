# mvapipe

mvapipe uses the GRCh38 whole-genome VCF of one patient to rank variants as possible causes
of a rare genetic disease. A candidate can contain one variant, or two heterozygous
variants in the same gene, since a recessive disease may involve different variants on the
two copies of the gene. The VCF is normalised with bcftools and annotated using VEP. From
these annotations, mvapipe calculates a score for each candidate and writes the candidates
in decreasing score order to the predictions file for Track 1 of the MVA Hackathon 2026.
The submitted predictions, methods report and run commands are kept in
[submission/](submission/README.md).

Team nadia: Sadat Bashir (Hugging Face: `saadte`) and Athanasia Karagiannopoulou
(Hugging Face: `nasia2612`).

## Stages

In a run with `--vcf`, bcftools 1.21 selects the patient sample and normalises the variants
against GRCh38. Chromosomes are renamed to 1 through 22, X, Y and MT before normalisation,
where multi-allelic sites are split and indels are left-aligned against the reference.
Sites with an alternate allele in the sample are retained in the normalised VCF. The REF
alleles are also checked against GRCh38, and the run stops at the first mismatch.

Ensembl VEP 116.2 annotates the normalised variants in chromosome groups, with chromosomes
having few variants sharing one VEP run. A consequence is selected for each variant and
gene on the MANE Select transcript, or according to the default transcript order of VEP
where a MANE Select transcript is absent. VEP also adds the gnomAD allele frequencies and
AlphaMissense scores. After VEP annotation, bcftools adds the ClinVar classification by
exact allele match, and `bcftools +split-vep` writes `variants.tsv.gz`, which contains one
row for each variant and gene. The candidates are scored from this table. The predictions
file is then checked against the Track 1 submission rules.

An interrupted run can reuse completed stages because each stage writes a marker in the
output directory's `.stages` folder after producing its output files. Markers are written
for normalisation, each VEP chromosome group, the ClinVar merge and table extraction.
The VCF path, size and modification time are recorded in the marker, together with the
sample name, mvapipe version and contents of `refs/manifest.json`. Before reusing a stage,
mvapipe compares these values and checks that its outputs exist, so a changed input or a
missing output causes the stage to run again. Setup markers record step completion
without checking the resource and environment definitions; a new installation directory
is required after changing `config/resources.env` or `environment.lock.yml`.

## Score

For each candidate, the score is calculated under a recessive, dominant or X-linked
inheritance model as

`score = genotype evidence * gene prior * phenotype factor`

The genotype evidence is calculated from the candidate's variants, while the gene prior
is obtained from the GenCC classification of the gene for that inheritance model. The
phenotype factor depends on the similarity between the patient's HPO terms and the terms
of diseases annotated to the gene. The constants used in these calculations are given in
[mvapipe/scoring_config.py](mvapipe/scoring_config.py) and apply to each gene.

### Genotype evidence

The impact assigned to a variant starts with its VEP consequence and lies between 0 and
1, with the values below applying before the ClinVar adjustment. For a missense variant,
the impact increases with the AlphaMissense score. Stop gained and frameshift variants
in the last exon have their loss-of-function impact multiplied by 0.8, which reduces the
impact from 0.95 to 0.76.

| Consequence | Impact |
|---|---|
| Stop gained, frameshift, splice donor, splice acceptor, transcript ablation | 0.95 |
| Stop gained or frameshift in the last exon | 0.95 * 0.8 = 0.76 |
| Start lost | 0.70 |
| Stop lost | 0.60 |
| In-frame insertion or deletion, protein altering | 0.50 |
| Missense | 0.15 + 0.65 * AlphaMissense score; 0.35 if the score is absent |
| Splice region | 0.25 |
| Synonymous | 0.03 |
| Any other consequence | 0.01 |

For a ClinVar pathogenic or likely pathogenic allele, the impact is raised to at least
1.00 if the review status has two or more stars, 0.95 with one star and 0.85 with no
stars. A benign or likely benign allele with at least one star has its impact capped at
0.02. Conflicting classifications use the counts of pathogenic, benign and uncertain
submissions. If a pathogenic submission is present, the minimum impact is 0.60; among
the remaining conflicting records, benign submissions outnumbering uncertain submissions
cap the impact at 0.05. A conflicting record with some benign submissions and at least
as many uncertain submissions instead has its impact multiplied by 0.5. The consequence
impact is retained for other conflicting records.

The weight of a variant is the product of its impact, rarity weight and call-quality
factor. The allele frequency for the rarity calculation is the largest of the gnomAD
exome frequency, gnomAD genome frequency and maximum population frequency reported by VEP.
Under a recessive model the rarity weight remains 1 up to a frequency of 0.001, then
decreases as a linear function of the logarithm of frequency until it reaches 0 at 0.02;
for a dominant model the corresponding frequencies are 0.00001 and 0.001. Missing
frequencies give a rarity weight of 1. There is also an exception for ClinVar pathogenic
alleles with at least two review stars, which retain a recessive rarity weight of 1 while
the frequency is below 0.05.

The call-quality factor starts from 1 and is multiplied by 0.5 for a failed variant-caller
filter, by 0.5 for genotype quality below 20 and by 0.5 for read depth below 8. These
factors accumulate. For a heterozygous call, a further factor of 0.7 applies if the fraction
of reads carrying the alternate allele is below 0.25 or above 0.85. Read depth can be
obtained by summing the allele depths if the depth field is missing, but a check is
skipped when the values needed for that check are absent.

For compound candidates, mvapipe takes at most the 12 heterozygous variants with the
highest recessive weights in a gene, excluding variants whose weight is below 0.05.
The phase of a pair is read from the VCF. A pair placed on the same chromosome copy is
removed, whereas a pair with unknown phase is retained with its genotype evidence
multiplied by 0.9. The genotype evidence for each kind of candidate is given below.

| Candidate | Model | Genotype evidence |
|---|---|---|
| Homozygous variant | `AR_hom` | Weight of the variant |
| Pair of heterozygous variants of the same gene | `AR_comphet` | Geometric mean of the two weights |
| Heterozygous variant with an impact of at least 0.85, in a gene whose recessive classification in GenCC is supportive or stronger | `AR_single` | 0.25 * weight of the variant |
| Heterozygous variant with an impact of at least 0.15 | `AD` | 0.8 * weight of the variant |

For a gene on chromosome X, the heterozygous variants are scored under `XL` with the
genotype evidence of an `AD` candidate, and `AR_comphet` and `AR_single` candidates are
not formed. Homozygous and hemizygous variants use their full recessive weight; the model
is `XL` if the sample is inferred as male or the call is hemizygous, and `AR_hom`
otherwise. To infer sex, mvapipe uses the retained X-chromosome calls outside the
pseudo-autosomal regions and takes the sample as male if at least 20 calls are present
and fewer than 15% are heterozygous. Below 20 calls the result is unknown.

### Gene prior

GenCC contains the gene-disease classifications submitted by different organisations,
with an inheritance mode for each classification. The strongest classification of a gene
under recessive, dominant or X-linked inheritance gives the prior for candidates of the
corresponding model. Disputed, refuted, no known disease relationship and animal model
only classifications are treated as no evidence. For a gene lacking a classification
under the candidate's inheritance, a lower prior is used, and this prior also depends on
whether GenCC has usable evidence for the gene under another inheritance.

| GenCC classification for the candidate's inheritance | Gene prior |
|---|---|
| Definitive or strong | 1.0 |
| Moderate | 0.9 |
| Supportive | 0.8 |
| Limited | 0.55 |
| No classification for this inheritance; usable evidence for another inheritance | 0.50 recessive, 0.15 dominant, 0.40 X-linked |
| No usable classification for the gene | 0.45 recessive, 0.30 dominant, 0.40 X-linked |

### Phenotype factor

The patient's HPO terms are compared with the terms annotated to each disease of a gene,
giving a similarity between 0 and 1. The phenotype factor is calculated as

`phenotype factor = 0.30 + 0.70 * similarity`

The information content of an HPO term is the negative logarithm of the fraction of
annotated genes having that term or one of its descendants. For each patient term,
mvapipe finds the most informative ancestor shared with a disease term and divides the
information content of this ancestor by that of the patient term. These ratios are
averaged over the patient terms to give the disease similarity, with a term of zero
information content contributing zero. If several diseases are annotated to the same
gene, the highest disease similarity is used, while a gene without an annotated disease
has similarity 0 and a phenotype factor of 0.30. In a blind run the HPO terms are not read,
and the phenotype factor is 1 for each gene.

### Ranking and EPCR

The same set of alleles may produce candidates under different inheritance models, in
which case the candidate with the highest score is retained. Candidates are ordered by
decreasing score, with ties ordered by gene, alleles and model.

Track 1 requires an EPCR value above 0 and up to 1 for each row. The EPCR of the first
candidate is its score capped at 0.99. For the following candidates, EPCR is scaled by
the square of the candidate's score relative to the first score, and a value that would
be within 0.0001 of the preceding EPCR is lowered. This gives strictly decreasing EPCR
over the first 50 candidates. EPCR has not been calibrated as a clinical probability.

## Running mvapipe

The tools and reference data need about 40 GB of free disk space on Linux. They are
installed by `setup/setup_env.sh` in the directory given by `MVA_ROOT`, which must be set
before running the script. Setup took 10 to 14 minutes on Google Colab.

```bash
export MVA_ROOT=/path/to/installation
bash setup/setup_env.sh
```

For a patient run, give the VCF and a file of HPO terms to `mvapipe run` from the
repository root. The HPO file contains one ID per line, for example `HP:0004322`, and
text after the ID or after a `#` is ignored. Each ID is checked against HPO; an obsolete
ID or a term that does not describe a patient feature, such as an inheritance mode, stops
the run. Whole-genome annotation took about 71 minutes on a Colab machine with 8 vCPU.
The ranking took less than a minute.

```bash
"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe run \
    --vcf sample.vcf.gz --hpo terms.txt --refs "$MVA_ROOT/refs" --out runs/sample
```

The same variants can be ranked without HPO terms by replacing `--hpo terms.txt` with
`--blind`. With the same VCF, reference data and output directory, the completed stages
are reused and the ranking is repeated with a phenotype factor of 1.

```bash
"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe run \
    --vcf sample.vcf.gz --blind --refs "$MVA_ROOT/refs" --out runs/sample
```

The predictions file contains at most ten candidates, as required by Track 1. The run
options include:

| Option | Use |
|---|---|
| `--min-relative-score` | Minimum fraction of the best score for additional candidates (default 0; 0.8 for the submitted file) |
| `--calibration-scores` | File with the true candidate's score for each benchmark case |
| `--coverage` | Target conformal coverage for the score threshold (default 0.9); used with `--calibration-scores` |
| `--sample` | Sample name; required for a VCF with more than one sample |
| `--proband-id` | Proband ID in the predictions file (default `PROBAND01`) |
| `--threads` | Number of threads (default: the machine's CPU count) |

### Conformal candidate sets

With `--calibration-scores`, mvapipe calculates a split conformal threshold from benchmark
cases whose causal variants are known, and applies this threshold to the candidate
scores of the patient. `--min-relative-score` cannot be supplied in the same command.
The calibration file is prepared from the benchmark runs as follows:

1. Rank each benchmark case with the mvapipe version and ranking settings used for the
   patient, including the choice between `--hpo` and `--blind`.
2. Copy the true candidate's score from each case's `ranked.<label>.tsv` to a separate line
   of the calibration file.
3. Use 0 if the true candidate is absent from that table, which contains the first 50
   candidates.

The file contains one finite, nonnegative score for each benchmark case and must have
at least one score. Blank lines and text after `#` are ignored. Scores copied from the
ranked table have four decimal places, so the threshold calculated from them is also
rounded. This rounding can affect coverage for candidates close to the threshold.

The target coverage `c` is given by `--coverage` and defaults to 0.9, with allowed values
between 0 and 1 excluding both endpoints. For `n` calibration scores, mvapipe calculates
`k = ceil((n + 1) * c)`. The k-th highest score gives the threshold, except when `k`
exceeds `n`, for which the threshold is 0. With unrounded scores and exchangeable
calibration and patient cases scored by the same procedure, the threshold is at or below
the true candidate's score with probability at least `c`.

The best candidate is written to the CSV even when its score is below the threshold.
Other candidates reaching the threshold follow in decreasing score order, but Track 1
allows at most ten rows, so a true candidate reaching the threshold can still be left
out. A true candidate missing from the generated candidates is also absent from the CSV;
these exclusions prevent the score threshold's coverage guarantee from extending to the
predictions file. The run manifest records the threshold, target coverage and number of
calibration cases in its `conformal_set` field.

### Output files

In the output file names, `<label>` is `phenotype` for a run with `--hpo` and `blind` for
a run with `--blind`; a different label can be given with `--label`.

| File | Content |
|---|---|
| `<label>.csv` | Track 1 predictions |
| `ranked.<label>.tsv` | First 50 candidates, with scores, EPCR, genotype evidence, gene priors, phenotype similarity and variant annotations |
| `manifest.<label>.json` | mvapipe, tool and reference versions; stage counts and timings; scoring constants and candidate selection settings |
| `variants.tsv.gz` | Annotated table with one row per variant and gene |
| `normalised.vcf.gz`, `chunks/`, `annotated.vcf.gz` | VCF after normalisation, VEP annotation of each chromosome group, and the ClinVar merge |
| `logs/` | bcftools and VEP output for each stage |

The ranking can also be repeated from `variants.tsv.gz` by using
`--table variants.tsv.gz` in place of `--vcf sample.vcf.gz`, without reading or annotating
the VCF again.

The predictions file is checked at the end of a run, and `mvapipe validate` applies the
same checks to an existing file. These include the required column names and allele
fields, a maximum of ten rows and the same proband ID in each row. The EPCR values must
be above 0 and up to 1, with a single row having the highest value.

```bash
"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe validate \
    runs/sample/phenotype.csv
```

## Scope

mvapipe ranks single-nucleotide variants and short insertions and deletions from one
sample on chromosomes 1 through 22, X and Y. Variants with a VEP impact of HIGH, MODERATE
or LOW in a gene are read for ranking, including protein-changing, splice-region and
synonymous variants. Variants classified as pathogenic, likely pathogenic or conflicting
in ClinVar are also read irrespective of their VEP consequence. The remaining non-coding
variants, structural and copy-number variants, repeat expansions and mitochondrial
variants are excluded, and mosaicism is not modelled. For compound pairs the phase comes
from the VCF, with unknown phase retained at the penalty described above. A gene without
a disease annotation in HPO receives the lowest phenotype factor.

## Code

| Directory | Content |
|---|---|
| `mvapipe/` | VCF processing, candidate scoring and ranking, CSV output and validation |
| `setup/` | Environment setup and generation of `refs/manifest.json` |
| `config/` | Reference versions and download addresses, with the GenCC snapshot for this release |
| `submission/` | Track 1 submission files |

Functions and variables use `snake_case`; classes and enum members use
`Capitalized_snake_case`, and constants use `UPPER_SNAKE_CASE`. These conventions follow
Jane Street's [notes on naming](https://blog.janestreet.com/notes-on-naming/).

## Reference data

The reference versions and download addresses are fixed in
[config/resources.env](config/resources.env), and the source of the GenCC snapshot is
recorded in [config/gencc_snapshot.json](config/gencc_snapshot.json). The GenCC file
itself is kept in `config/`. At the end of setup, `refs/manifest.json` records the tool
versions, reference versions and sizes of downloaded files, and the contents of this
manifest are copied into `manifest.<label>.json` for each run.

| Resource | Version | Licence |
|---|---|---|
| Ensembl VEP and cache | 116.2; cache 116, GRCh38 | Apache 2.0 |
| GRCh38 primary assembly | Ensembl release 116 | open |
| gnomAD allele frequencies | VEP cache 116 | open |
| ClinVar | VCF file date 2026-10-04 | public domain |
| AlphaMissense | hg38 scores, 2023 | CC BY 4.0 |
| Human Phenotype Ontology | v2026-09-01 | free with attribution |
| GenCC | export of 2026-10-13, three columns | CC0 1.0 |

## Data in this repository

This repository contains the code, GenCC snapshot and submission files. The submission
includes `submission/nadia_mvapipe.csv` with the single submitted candidate and
`submission/hpo.txt` with eight HPO codes, both derived from the patient's data. The
patient VCF, derived variant tables and clinical document are kept outside the repository.
VCFs, variant tables, run outputs, DOCX documents, term files and credentials are excluded
by `.gitignore`, with exceptions for the two submission files above.

## Licence

The code and submission documents are licensed under [CC BY 4.0](LICENSE), with
attribution to nadia. Reference data remain subject to their providers' terms. The
resources are cited in the [methods report](submission/nadia_track1_report.md), which also
includes the acknowledgement required by the organisers.
