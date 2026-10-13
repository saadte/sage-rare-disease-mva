# mvapipe

mvapipe takes the GRCh38 whole-genome VCF of a single patient and ranks the variants as the
possible causes of a rare genetic disease. The variants are first normalised using bcftools
and annotated using VEP, after which a score is calculated for each variant and also for
each pair of heterozygous variants lying in the same gene, because in a recessive disease
the two copies of the gene can carry different variants. The candidates are then written in
decreasing order of this score in the predictions file required for Track 1 of the MVA
Hackathon 2026. The predictions file, the methods report and the commands of the submitted
run are kept in [submission/](submission/README.md).

**Team:** nadia  
**Members:** Sadat Bashir (Hugging Face: `saadte`) and Athanasia Karagiannopoulou (Hugging Face: `nasia2612`)

## Stages

When `python -m mvapipe run` is given a VCF, the sample is first selected and its variants
are normalised using bcftools 1.21, where the chromosomes are renamed to 1..22, X, Y and
MT, the multi-allelic sites are split and the indels are left-aligned against GRCh38. In
this step only the sites having an alternate allele in the sample are kept, and the run
stops at the first record whose REF allele differs from the reference. The normalised
variants are then annotated using Ensembl VEP 116.2, which is run separately for each
chromosome, except that the chromosomes having few variants are annotated together. VEP
reports one consequence for each variant and gene, which is calculated on the MANE Select
transcript, and where a MANE Select transcript does not exist the default transcript order
of VEP is followed. The gnomAD allele frequencies and the AlphaMissense score are also
added by VEP. After annotation, the ClinVar classification is added for each exact allele
using bcftools, and `bcftools +split-vep` is used for writing `variants.tsv.gz`, which is a
table having one row for each variant and gene. The candidates are scored and ordered from
this table, and finally the predictions file is checked against the submission rules of
Track 1.

During a run, a marker is written in the `.stages` directory of the output when
normalisation, the annotation of a chromosome, the ClinVar merge or the table extraction is
completed. When the same command is given again, the stages having a marker are skipped, so
a run which was interrupted continues from the stage at which it stopped. The marker stores
the path, size and modification time of the VCF along with the sample name, the version of
mvapipe and the contents of `refs/manifest.json`. If any of these has changed, or if an
output file is missing, the marker is removed and the stage is run again, therefore a
marker is present only for a stage completed with the current input. The completed steps of
`setup/setup_env.sh` are skipped in the same way. However, a new installation directory is
required when `config/resources.env` or `environment.lock.yml` is changed.

## Score

A candidate is either a single variant or a pair of heterozygous variants of the same gene,
and each candidate is scored under a recessive, dominant or X-linked model of inheritance.
The score of a candidate is calculated as

`score = genotype evidence × gene prior × phenotype factor`

Here the genotype evidence is calculated from the variants of the candidate. The gene prior
is taken from the GenCC classification of the gene for the inheritance model, while the
phenotype factor depends on the similarity of the HPO terms of the patient with the HPO
terms of the diseases of the gene. The values used in this calculation are given in
[mvapipe/scoring_config.py](mvapipe/scoring_config.py), and they are the same for every
gene.

### Genotype evidence

For calculating the genotype evidence, an impact between 0 and 1 is first assigned to each
variant from the consequence reported by VEP, and the values of this impact are given in
the table below. For a missense variant, the impact increases with its AlphaMissense score,
whereas a stop gained or frameshift variant lying in the last exon of the transcript gets a
lower impact than the other loss-of-function variants.

| Consequence | Impact |
|---|---|
| Stop gained, frameshift, splice donor, splice acceptor, transcript ablation | 0.95 |
| Stop gained or frameshift in the last exon | 0.95 × 0.8 = 0.76 |
| Start lost | 0.70 |
| Stop lost | 0.60 |
| In-frame insertion or deletion, protein altering | 0.50 |
| Missense | 0.15 + 0.65 × AlphaMissense score, and 0.35 if the score is absent |
| Splice region | 0.25 |
| Synonymous | 0.03 |
| Any other consequence | 0.01 |

This impact is then changed according to the ClinVar classification of the allele. For an
allele classified as pathogenic or likely pathogenic, the impact is raised to at least 1.00
if the review status has two or more stars, to at least 0.95 if it has one star and to at
least 0.85 if it has no star, whereas for an allele classified as benign or likely benign
with at least one star, the impact is limited to 0.02. In the case of a conflicting
classification, the pathogenic, benign and uncertain submissions are counted.
If a pathogenic submission is present, the impact is raised to at least 0.60. For the
remaining conflicting records, the impact is limited to 0.05 when the benign submissions
outnumber the uncertain submissions, and otherwise it is multiplied by 0.5 if at least one
benign submission is present.

The weight of a variant is the product of its impact, a rarity weight and a call-quality
factor. For the rarity weight, the allele frequency is taken as the highest of the gnomAD
exome frequency, the gnomAD genome frequency and the maximum population frequency reported
by VEP. Under a recessive model, the rarity weight is 1 up to an allele frequency of 0.001,
after which it decreases on a logarithmic scale of the frequency and becomes 0 at 0.02,
while under a dominant model the corresponding frequencies are 0.00001 and 0.001. An allele
which is absent from gnomAD is also given the rarity weight of 1. In addition, an allele
which is pathogenic in ClinVar with two or more stars keeps the rarity weight of 1 under a
recessive model as long as its frequency is below 0.05.

The call-quality factor starts from 1, and it is multiplied by 0.5 for a call which has
failed a filter of the variant caller, by 0.5 for a genotype quality below 20 and by 0.5
for a read depth below 8. For a heterozygous call, the factor is also multiplied by 0.7
when the fraction of the reads carrying the alternate allele is below 0.25 or above 0.85.
If the genotype quality, the depth or the allele depths are missing in the VCF, the
corresponding check is skipped.

The candidates of a gene are then formed from the weights of its variants, and the genotype
evidence of each kind of candidate is given in the table below. For a pair of heterozygous
variants, the phase is read from the VCF. If the VCF places both variants on the same copy
of the chromosome, the pair is removed, whereas a pair whose phase is unknown is kept and
its genotype evidence is multiplied by 0.9. For forming the pairs of a gene, at most 12
heterozygous variants having the highest weights are used, and a variant with a weight
below 0.05 is left out.

| Candidate | Model | Genotype evidence |
|---|---|---|
| Homozygous variant | `AR_hom` | Weight of the variant |
| Pair of heterozygous variants of the same gene | `AR_comphet` | Geometric mean of the two weights |
| Heterozygous variant with an impact of at least 0.85, in a gene whose recessive classification in GenCC is supportive or stronger | `AR_single` | 0.25 × weight of the variant |
| Heterozygous variant with an impact of at least 0.15 | `AD` | 0.8 × weight of the variant |

For a gene on chromosome X, the `AR_comphet` and `AR_single` candidates are not formed. A
heterozygous variant on chromosome X is scored under the X-linked model `XL` with the same
genotype evidence as an `AD` candidate. A homozygous variant on chromosome X gets its full
weight, where the model is `XL` if the sample is male or the call is hemizygous, and
`AR_hom` otherwise. The sex is inferred from the calls on chromosome X lying outside the
pseudo-autosomal regions, where the sample is taken as male if less than 15% of these calls
are heterozygous. With fewer than 20 such calls, the sex remains unknown.

### Gene prior

The gene prior is obtained from GenCC, which contains the gene-disease classifications
submitted by different organisations along with the inheritance of the disease. For each
gene, the strongest classification submitted with a recessive, dominant or X-linked
inheritance is used for the candidates of the corresponding model, and the gene prior for
each classification is given in the table below. The classifications disputed, refuted, no
known disease relationship and animal model only are not used. When the classification for
the inheritance of the model is missing, a lower gene prior is used, and its value depends
on whether GenCC classifies the gene under another inheritance.

| GenCC classification of the gene for the inheritance of the model | Gene prior |
|---|---|
| Definitive or strong | 1.0 |
| Moderate | 0.9 |
| Supportive | 0.8 |
| Limited | 0.55 |
| No classification for this inheritance, while the gene is classified under another inheritance | 0.50 for a recessive model, 0.15 for a dominant model and 0.40 for the X-linked model |
| No classification of the gene in GenCC | 0.45 for a recessive model, 0.30 for a dominant model and 0.40 for the X-linked model |

### Phenotype factor

The phenotype factor is calculated as

`phenotype factor = 0.30 + 0.70 × similarity`

where the similarity lies between 0 and 1 and is obtained by comparing the HPO terms of the
patient with the HPO terms annotated to each disease of the gene. For this comparison, the
information content of an HPO term is calculated as the negative logarithm of the fraction
of the annotated genes which have this term or one of its descendants. For each patient
term, the most informative ancestor shared with a term of the disease is found, and the
information content of this ancestor is divided by the information content of the patient
term. These ratios are averaged over the patient terms, which gives the similarity of the
patient with the disease. Where a gene has several diseases, the highest of these
similarities is used. For a gene without any annotated disease, the similarity is 0 and the
phenotype factor is therefore 0.30. In a run with `--blind`, the HPO terms are not read and
the phenotype factor is 1 for every gene.

### Ranking and EPCR

When the same set of alleles is obtained under more than one model, only the candidate
having the highest score is kept. The candidates are then ordered by decreasing score, and
equal scores are ordered by gene, alleles and model. The predictions file of Track 1 also
requires a value named EPCR, between 0 and 1, in each row. The first candidate gets its
score as EPCR, limited to 0.99, and the EPCR of each following candidate is in proportion
to the square of its score. Where the value obtained for a candidate is within 0.0001 of
the previous value, it is lowered, so that EPCR decreases strictly over the first 50
candidates. EPCR has not been calibrated as a clinical probability.

## Running mvapipe

mvapipe is run on Linux, where about 40 GB of free disk space is required for the tools and
the reference data. These are installed by `setup/setup_env.sh` in the directory given by
`MVA_ROOT`, and the installation took 10 to 14 minutes on Google Colab.

```bash
export MVA_ROOT=/path/to/installation
bash setup/setup_env.sh
```

For ranking the variants of a patient, the VCF and a file of HPO terms are given to
`mvapipe run`, from the root of this repository and inside the installed environment. The
file of HPO terms contains one HPO ID in each line, for example `HP:0004322`, where the
text after the ID and the text after a `#` are ignored. The run stops at the first ID which
is obsolete in HPO or which does not describe a feature of a patient, for example a mode of
inheritance. Annotation of a whole genome took about 71 minutes on a Colab machine with 8
vCPU, while the ranking is completed within a minute.

```bash
"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe run \
    --vcf sample.vcf.gz --hpo terms.txt --refs "$MVA_ROOT/refs" --out runs/sample
```

The same VCF can also be ranked without the HPO terms by giving `--blind` in place of
`--hpo`. With the same output directory, the completed stages are skipped and only the
ranking is repeated.

```bash
"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe run \
    --vcf sample.vcf.gz --blind --refs "$MVA_ROOT/refs" --out runs/sample
```

At most ten candidates are written in the predictions file, which is the limit of Track 1.
The other options of `mvapipe run` are given below.

| Option | Use |
|---|---|
| `--min-relative-score` | A candidate after the first one is written only if its score is at least this fraction of the first score (0.8 for the submitted file) |
| `--sample` | The sample to be analysed, when the VCF contains more than one sample |
| `--proband-id` | The proband ID written in the predictions file (default `PROBAND01`) |
| `--threads` | The number of threads (default: the number of CPUs of the machine) |

`mvapipe run` writes the following files in the output directory. In the file names,
`<label>` is `phenotype` for a run with `--hpo` and `blind` for a run with `--blind`, and
another label can be given with `--label`.

| File | Content |
|---|---|
| `<label>.csv` | The predictions file in the format of Track 1 |
| `ranked.<label>.tsv` | The first 50 candidates with the score, EPCR, genotype evidence, gene prior and phenotype similarity, and the annotations of each variant |
| `manifest.<label>.json` | The version of mvapipe, the versions of the tools and the reference data, the number of records and the time taken in each stage, and the values used for scoring |
| `variants.tsv.gz` | The annotated table having one row for each variant and gene |
| `normalised.vcf.gz`, `chunks/`, `annotated.vcf.gz` | The VCF after normalisation, after the annotation of each chromosome and after the ClinVar merge |
| `logs/` | The output of bcftools and VEP for each stage |

The ranking can be repeated without the VCF by giving `variants.tsv.gz` with `--table` in
place of `--vcf`.

The predictions file is checked at the end of each run, and the same checks can be run
separately using `mvapipe validate`. The checks follow the submission rules of Track 1,
which include the required columns, a maximum of ten rows, the same proband ID in every
row, EPCR above 0 and up to 1, and a single row having the highest EPCR.

```bash
"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe validate \
    runs/sample/phenotype.csv
```

## Scope

mvapipe ranks the single-nucleotide variants and the short insertions and deletions present
in the VCF of one sample, on chromosomes 1 to 22, X and Y. A variant is read for ranking
when VEP gives it an impact of HIGH, MODERATE or LOW in a gene, which covers the variants
changing the protein, the variants of a splice region and the synonymous variants. In
addition, a variant which is pathogenic, likely pathogenic or conflicting in ClinVar is
read irrespective of its consequence. The remaining non-coding variants, the structural and
copy-number variants, the repeat expansions and the variants of the mitochondrial genome
are not ranked, and mosaicism is not modelled. Since only one sample is analysed, the phase
of two heterozygous variants is available only where the variant caller has written it in
the VCF. A gene whose disease is not annotated in HPO gets the lowest phenotype factor.

## Code

| Directory | Content |
|---|---|
| `mvapipe/` | The stages run on the VCF, the ranking, the predictions file and its checks |
| `setup/` | `setup_env.sh`, and the program which writes `refs/manifest.json` |
| `config/` | The versions and download addresses of the reference data, and the GenCC file used for this release |
| `submission/` | The files submitted to Track 1 |

In the code, functions and variables are written in `snake_case`, classes and enum members
in `Capitalized_snake_case` and constants in `UPPER_SNAKE_CASE`, following the
[notes on naming](https://blog.janestreet.com/notes-on-naming/) of Jane Street.

## Reference data

The versions and download addresses of the reference data are fixed in
[config/resources.env](config/resources.env), while the GenCC file is kept in `config/` and
its source is given in [config/gencc_snapshot.json](config/gencc_snapshot.json). At the end
of the installation, `refs/manifest.json` is written with the versions of the tools, the
versions of the reference data and the sizes of the downloaded files. The same contents are
written in `manifest.<label>.json` of each run.

| Resource | Version | Licence |
|---|---|---|
| Ensembl VEP and cache | 116.2; cache 116, GRCh38 | Apache 2.0 |
| GRCh38 primary assembly | Ensembl release 116 | open |
| gnomAD allele frequencies | as given in the VEP cache 116 | open |
| ClinVar | VCF file date 2026-10-04 | public domain |
| AlphaMissense | hg38 scores, 2023 | CC BY 4.0 |
| Human Phenotype Ontology | v2026-09-01 | free with attribution |
| GenCC | export of 2026-10-13, three columns | CC0 1.0 |

## Data in this repository

This repository contains the code, the GenCC file and the submission files. Among the
submission files, `submission/nadia_mvapipe.csv` contains the one submitted candidate and
`submission/hpo.txt` contains the eight HPO codes, and both are derived from the data of
the patient. The VCF of the patient, the variant tables obtained from it and the clinical
document are kept outside the repository. `.gitignore` excludes the VCF files, the variant
tables, the run outputs, the DOCX documents, the term files and the credentials, where the
two submission files given above are the exceptions.

## Licence

The code and the submission documents are licensed under [CC BY 4.0](LICENSE), with
attribution to **nadia**. The reference data remain under the terms of their providers. The
resources are cited in the [methods report](submission/nadia_track1_report.md), which also
contains the acknowledgement required by the organisers.
