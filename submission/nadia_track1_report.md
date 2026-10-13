# mvapipe: variant prioritisation for Track 1

**Team:** nadia  
**Members:** Sadat Bashir (Hugging Face: saadte) and Athanasia Karagiannopoulou (Hugging Face: nasia2612)  
**Model:** mvapipe phenotype, version 0.4.0  
**Prepared:** 13 October 2026  
**Challenge:** Rare Disease, Real Kid: MVA Hackathon 2026, Track A / Track 1

## Abstract

Interpretation of a single-sample VCF requires prioritisation of individual
variants as well as possible compound-heterozygous pairs, depending on the
inheritance model considered for the disease. In this work, we developed mvapipe
for a GRCh38 VCF, where the variants are first normalised using bcftools and
annotated using VEP, after which a score is calculated under recessive, dominant,
and X-linked models. For calculating this score, the sequence consequence,
ClinVar classification, gnomAD frequency, and quality of the variant call are
combined with GenCC gene-disease evidence and the similarity of the supplied HPO
terms with the disease phenotype; in the case of a missense variant, the
AlphaMissense score is also used for calculating its impact. Using the eight HPO
terms supplied for the analysed sample, the *BUB1B* variants chr15:40209701 T>G
and chr15:40220612 T>G are obtained as the first-ranked compound-heterozygous
pair, and the same pair remains at the first rank when the phenotype factor is
removed. The first variant is stop gained and has a reviewed pathogenic/likely
pathogenic assertion in ClinVar, whereas the second variant is missense with an
AlphaMissense score of 0.923. Both variants are heterozygous, and the allele
depths and genotype quality support the two calls. However, their phase is not
known, therefore segregation evidence and further assessment of the missense
variant are required for interpreting this pair under a recessive model.

To evaluate the ranking, 120 synthetic disease genotypes were inserted separately
into each of the two GIAB genomes, and the ClinVar records of the inserted alleles
were masked before running mvapipe. For HG001, the exact disease genotype is
ranked first in 97.5% of the cases using either exact or noisy HPO terms, while it
is ranked first in 89.2% of the cases when the phenotype terms are removed. The
corresponding values for HG002 are 94.2%, 90.0%, and 73.3% with exact, noisy, and
no phenotype terms, respectively. Both GIAB genomes were used during development
of the ranker, and the allele-balance threshold was changed after the first run
on the challenge genome. Therefore, these results are obtained for the selected
synthetic disease genotypes having clean variant calls. Clinical diagnostic
accuracy of mvapipe has not been evaluated.

## Data and annotation

The organisers provided a single-sample GRCh38 VCF and a clinical document in
which the HPO codes were already given. We extracted these codes from the
document and reviewed the extracted codes before using them in mvapipe; the
eight codes used for the submitted run are given in [hpo.txt](hpo.txt). Recurrent
spontaneous abortion was also present among the supplied codes. Since this term
can describe reproductive or family history, its interpretation as a phenotype
of the proband is not direct. However, we kept the term in the submitted run
because it was present in the clinical document. We also ran mvapipe after
removing this term, where the same *BUB1B* pair was obtained at the first rank
and the output again contained one row.

The public resources used for annotation, ranking and evaluation are given below.

| Resource | Version or snapshot | Use |
|---|---|---|
| GRCh38 primary assembly | Ensembl release 116 | Allele normalisation |
| Ensembl VEP and cache | 116.2; cache 116, GRCh38 | Consequences and gnomAD frequencies |
| ClinVar | VCF file date 2026-10-04 | Allele assertions and review status |
| AlphaMissense | hg38 scores, 2023 | Missense effect prediction |
| Human Phenotype Ontology | 2026-09-01 | Term ancestry and disease annotations |
| GenCC | 2026-10-13 export; required columns retained | Gene-disease validity and inheritance |
| GIAB HG001 and HG002 | GRCh38 benchmark calls v4.2.1, autosomes | Synthetic benchmark backgrounds |

For GenCC, the export downloaded on 13 October 2026 was used. mvapipe reads
`gene_symbol`, `classification_title` and `moi_title` from this export, therefore
these three columns were retained for all 30,345 rows and the gene validity
records obtained from the retained file were compared with the records from the
original export. The records were identical for 6,109 genes. The source, retrieval
date and file sizes are given in
[gencc_snapshot.json](../config/gencc_snapshot.json). The VEP plugin is also fixed
to a commit of release 116.

In mvapipe, the sample is first selected and its variants are normalised using
bcftools 1.21, where the contigs are renamed, the multiallelic records are split
and the indels are left-aligned against GRCh38. The normalisation stops at the
first record whose reference allele differs from the assembly. Since splitting
changes the representation of a multiallelic site, the original record and the
index of the alternate allele are retained with each split record for
interpreting its genotype. The normalised variants are then annotated using VEP
116.2 with the offline cache and the AlphaMissense plugin, one chromosome at a
time. For each variant and gene, VEP reports one consequence, which is
calculated on the MANE Select transcript, and where a MANE Select transcript
does not exist the default transcript order of VEP is followed. After
annotation, the ClinVar record of each exact allele is added, and the table used
for ranking is extracted using bcftools `+split-vep`. Each ClinVar record is
used as evidence for the gene which it names.

All the reference data and the two genomes used for evaluation are public,
whereas the VCF of the patient requires authorised access to the challenge data.
The weights used in mvapipe are fixed values selected by us, and AlphaMissense
enters the ranking through its precomputed scores. The ranking is calculated
using the Python standard library, while the extraction of the HPO codes from
the document and the benchmark evaluation are separate programs.

## Ranking model

For each hypothesis, that is, a variant or a pair of variants considered under
an inheritance model, the score is calculated as

`score = genotype evidence × gene prior × phenotype factor`.

For calculating the genotype evidence, an impact is first assigned to each
variant from its sequence consequence. A loss-of-function variant gets an impact
of 0.95, which is multiplied by 0.8 for a truncating variant lying in the last
exon. For a missense variant, the impact is calculated as `0.15 + 0.65 ×
AlphaMissense`, and 0.35 is used when the AlphaMissense score is absent. For a
pathogenic/likely pathogenic assertion in ClinVar, the impact is then raised to
a minimum value depending on the review status, whereas a reviewed benign
assertion limits the impact. For a conflicting record, the impact is raised or
lowered according to its numbers of pathogenic and benign submissions. All these
values, including the impact of the other consequences, are given in
[scoring_config.py](../mvapipe/scoring_config.py).

The impact is multiplied by a rarity weight, which decreases on a logarithmic
scale of the allele frequency. For a recessive hypothesis, the rarity weight is
1 up to a frequency of 0.001 and becomes 0 at 0.02, while for a dominant
hypothesis these two frequencies are 0.00001 and 0.001. A reviewed pathogenic
allele having at least two stars keeps the rarity weight of 1 under a recessive
hypothesis as long as its frequency is below 0.05, and an allele which is absent
from gnomAD also gets the rarity weight of 1.

The quality of the variant call enters as another factor of the variant weight.
This factor is multiplied by 0.5 for a non-PASS filter, by 0.5 for a depth below
8 and by 0.5 for a genotype quality below 20. For a heterozygous call, the
allele balance is also used, where a value below 0.25 or above 0.85 gives a
factor of 0.7. If the genotype quality, the depth or the allele depths are
missing in the VCF, the corresponding factor remains 1.

From the weights of the variants of a gene, homozygous, compound-pair and
single-hit hypotheses are formed under the recessive model, along with dominant
and X-linked hypotheses. For a compound pair, the genotype evidence is the
geometric mean of the two variant weights. This value is multiplied by 0.9 when
the phase of the pair is unknown, while a pair which is phased in cis in the VCF
is discarded. For forming the pairs of a gene, at most 12 heterozygous variants
having the highest weights are used, where a weight of at least 0.05 is
required. For the X-linked hypotheses, the sex of the sample is inferred from
the heterozygosity of at least 20 X calls lying outside the pseudoautosomal
regions.

The gene prior is obtained from GenCC for the inheritance of the hypothesis,
where a strong or definitive classification gives the full weight, and lower
weights are used for a weaker classification or when the gene is classified only
under another inheritance. For the phenotype factor, the HPO terms of the
patient are compared with the terms of each disease associated with the gene.
The information content of a term is calculated as the negative logarithm of the
fraction of annotated genes which are associated with the term or its
descendants. Each patient term is compared with a disease through their most
informative common ancestor, which gives the Resnik similarity. After division
by the information content of the patient term, these similarities are averaged.
Since a gene can have several associated diseases, this average is calculated
for each disease and the highest value is used. The phenotype factor is then
obtained as `0.30 + 0.70 × similarity`, while in the blind ranking this factor
is set to 1.

When the same set of alleles is obtained under more than one hypothesis, only
the highest score is retained, and equal scores are ordered by gene, alleles and
inheritance model. For the challenge output, EPCR is calculated from the scores
such that it decreases strictly over the first 50 candidates, and the later
candidates get zero. EPCR has not been calibrated as a clinical probability.

## Patient result

For the submitted row, the following *BUB1B* pair was obtained, where the
positions are given on GRCh38.

| Position | Alleles | Consequence | Evidence in the patient run |
|---|---|---|---|
| chr15:40209701 | T>G | BUB1B p.Leu737*, NM_001211.6 | ClinVar pathogenic/likely pathogenic, two stars; heterozygous; AD 21,25; GQ 99; PASS |
| chr15:40220612 | T>G | BUB1B p.Asn1002Lys, NM_001211.6 | AlphaMissense 0.923; absent from the run's ClinVar and gnomAD annotations; heterozygous; AD 15,13; GQ 99; PASS |

This pair was obtained at the first rank in the phenotype run as well as in the
blind run. In the phenotype run, its score is 0.565983 with an HPO similarity of
approximately 0.61, while the score of the next candidate is approximately
0.441. For the submission, we used `--min-relative-score 0.8`, where the first
candidate is always written and another candidate is written only if its score
is at least 80% of the first score, up to a maximum of ten rows. Only the
*BUB1B* pair satisfied this condition, therefore the submitted file contains one
row. The alleles of this row and their order were written by the code, and we
then reviewed the annotations, the call quality and the available phase
information of the pair.

The two positions are separated by 10,911 bases. Both calls are heterozygous,
however their phase is not given in the VCF, and the recessive interpretation
requires the two variants to be in trans. Therefore, the unknown-phase factor is
applied to this pair in the ranking. Biallelic disruption of *BUB1B* is
associated with mosaic variegated aneuploidy [5], which supports the
interpretation at the level of the gene. The AlphaMissense score of 0.923 for
p.Asn1002Lys is a computational prediction, which does not establish the
pathogenicity of this missense allele. Therefore, we report the pair as a
candidate recessive explanation, for which independent assessment of the
missense allele and segregation or other phase evidence are required.

The ranking in mvapipe is automatic and is followed by our review of the
candidates. Structural and copy-number variants, repeats, mosaic mechanisms and
most of the non-coding variation were outside the scope of this analysis, and
other possible second alleles have not been excluded. The public submission
contains the selected primary pair, while the broad ranking tables of the
patient and the unrelated incidental alleles remain private.

## Public-genome evaluation

For the evaluation, 120 cases were generated using HG002 and the seed 20261012,
namely 40 compound-pair cases, 40 homozygous cases and 40 dominant cases. A gene
was selected only if it has a strong or definitive classification in GenCC and a
disease with at least five HPO terms. The inserted variants are autosomal
ClinVar alleles which are pathogenic/likely pathogenic with at least two review
stars and are not longer than 50 bases. A depth of 30, a genotype quality of 99
and balanced allele depths were assigned to the inserted calls, and the same
cases were also evaluated in HG001.

Each case was ranked with the ClinVar records of the inserted variants present,
and again after masking these records. The phenotype input was also changed. In
the exact condition, the terms were sampled from the HPO annotation of the
disease, whereas in the noisy condition each term was replaced by a parent term
with a probability of 0.5 and two terms of other diseases were added. In the
blind condition, the phenotype terms were not used. The following results are
for the exact inserted genotype with the ClinVar records masked.

| Background | Phenotype | First / 120 | Top 3 / 120 | Top 5 / 120 |
|---|---|---|---|---|
| HG002 | Exact | 113 (94.2%) | 120 (100%) | 120 (100%) |
| HG002 | Noisy | 108 (90.0%) | 120 (100%) | 120 (100%) |
| HG002 | Blind | 88 (73.3%) | 109 (90.8%) | 119 (99.2%) |
| HG001 | Exact | 117 (97.5%) | 120 (100%) | 120 (100%) |
| HG001 | Noisy | 117 (97.5%) | 120 (100%) | 120 (100%) |
| HG001 | Blind | 107 (89.2%) | 118 (98.3%) | 119 (99.2%) |

With the ClinVar records present, all 120 genotypes were ranked first in both
genomes under each phenotype condition.

The ranker was first developed using HG002 alone, with HG001 as a holdout. After
the first evaluation of HG001, its results were used for changing the gene
assignment of ClinVar evidence. Later, the minimum allele balance of a
heterozygous call was raised from 0.15 to 0.25 after inspecting the first run on
the patient, which moved the submitted pair from the fifth rank to the first
rank in the blind ranking.

The synthetic phenotype terms were sampled from the same disease annotations
which are used in the similarity calculation, therefore the exact condition is
favourable for the ranking. The selection from ClinVar also favours the
well-studied genes and the established alleles. Masking removes the direct
ClinVar information of the inserted variants, however AlphaMissense retains
indirect label information from its calibration. The inserted calls have clean
quality values and all the cases are autosomal. Also, the same 120 disease
genotypes were inserted in both genomes, so the two evaluations differ only in
the background variants. An independent clinical cohort and the ground truth of
the organisers were not included in this evaluation.

## Reproduction and running time

The repository contains the Linux environment lock, the setup script, the
sources of the reference data, the selected CSV and the HPO codes, and the
commands are given in [submission/README.md](README.md). For repeating the full
run, authorised access to the challenge data is required, since the VCF of the
patient and the genotype-scale tables derived from it are kept outside the
repository.

The patient run was measured on a Colab CPU runtime having 8 vCPU and
approximately 50 GB of RAM. The annotation took 4,273 seconds, while
normalisation, merge and extraction took 15.4, 54.0 and 35.3 seconds,
respectively, and another 40.9 seconds were taken for loading and ranking. In
this run, 4,962,060 normalised records were obtained, from which 31,997
variant/gene rows were retained for ranking. The setup took 10 to 14 minutes in
the later runs. The first patient session used approximately 0.48 Colab compute
units over 107.6 minutes including the session overhead, and the cost in
currency was not measured. The genome of the patient was annotated using version
0.2.2 of mvapipe, and the submitted ranking was calculated using version 0.4.0
from the annotated variant table.

## AI use

We used OpenAI ChatGPT Plus with training disabled and no feedback given. The assistant
was used for code assistance, the phenotype extraction and the evidence of the candidates. 
The annotation and the ranking are calculated by the code from the published constants, 
without any LLM inference step.

## References and acknowledgement

1. McLaren et al. The Ensembl Variant Effect Predictor. *Genome Biology* 17, 122
   (2016). [doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4).
2. Cheng et al. Accurate proteome-wide missense variant effect prediction with
   AlphaMissense. *Science* 381, eadg7492 (2023).
   [doi:10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492).
3. Gargano et al. The Human Phenotype Ontology in 2024: phenotypes around the world.
   *Nucleic Acids Research* 52, D1333-D1346 (2024).
   [doi:10.1093/nar/gkad1005](https://doi.org/10.1093/nar/gkad1005).
4. DiStefano et al. The Gene Curation Coalition: A global effort to harmonize
   gene-disease evidence resources. *Genetics in Medicine* 24, 1732-1744 (2022).
   [PubMed:35507016](https://pubmed.ncbi.nlm.nih.gov/35507016/).
5. Hanks et al. Constitutional aneuploidy and cancer predisposition caused by
   biallelic mutations in BUB1B. *Nature Genetics* 36, 1159-1161 (2004).
   [doi:10.1038/ng1449](https://doi.org/10.1038/ng1449).
6. [ClinVar data access and attribution](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).
7. [GIAB benchmark releases](https://github.com/genome-in-a-bottle/giab_latest_release).
8. [Rare Disease, Real Kid: MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026).

This work was made possible through the Hackathon, organized by Sage Bionetworks
in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking,
Evaluation, and Assessment Consortium for Science), with prize sponsorship from
AWS and Anthropic. We are deeply grateful to the child and their family who
generously contributed their data and their story to advance research into this
rare disease. We acknowledge their trust in making this Hackathon possible.

Copyright 2026 nadia. Original submission materials are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Third-party resources
retain their respective terms. This method uses the Human Phenotype Ontology,
release 2026-09-01.
