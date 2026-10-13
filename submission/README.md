# Track 1 submission

**Team:** nadia  
**Members:** Sadat Bashir (Hugging Face: `saadte`) and Athanasia Karagiannopoulou (Hugging Face: `nasia2612`)  
**Model:** mvapipe phenotype, version 0.4.0

| File | Content |
|---|---|
| [nadia_mvapipe.csv](nadia_mvapipe.csv) | The predictions file, having one row |
| [nadia_track1_report.pdf](nadia_track1_report.pdf) | The methods report, which is also given as [Markdown](nadia_track1_report.md) |
| [nadia_track1_report.tex](nadia_track1_report.tex) | The LaTeX source of the PDF |
| [hpo.txt](hpo.txt) | The eight HPO codes of the patient used for the submitted run |

## Reproducing the prediction

The commands are run on Linux from the root of the repository, with `VCF` set to the path
of the challenge VCF.

```bash
export MVA_ROOT=/path/to/installation
export VCF=/path/to/authorised/challenge.vcf.gz
bash setup/setup_env.sh

"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe run \
    --vcf "$VCF" --hpo submission/hpo.txt --refs "$MVA_ROOT/refs" \
    --out runs/submission --proband-id PROBAND01 --label nadia_mvapipe \
    --min-relative-score 0.8 --threads 8

"$MVA_ROOT/bin/micromamba" run -p "$MVA_ROOT/env" python -m mvapipe validate \
    runs/submission/nadia_mvapipe.csv
cmp submission/nadia_mvapipe.csv runs/submission/nadia_mvapipe.csv
```
