# Design Assignment 2 — v3

Pitch detection: **autocorrelation** (time domain) vs **cepstrum** (frequency domain). MATLAB R2025b + Signal Processing Toolbox.

| File | What it is |
|---|---|
| `D2_report_2026EEY7565.pdf` | **the report to submit** |
| `D2_report_2026EEY7565.md` | report source (opens in Obsidian, figures and maths render) |
| `D2_code_2026EEY7565.m` | **the code to submit** — one script, both detectors, all tests, all figures |
| `p3.mat`, `p4.mat` | pitch tracks for s3 and s4 (`p3_acf`, `p3_cep`, `p4_acf`, `p4_cep`, 300 × 1 Hz, 0 = unvoiced) |
| `assets/` | 15 figures, `results_log.txt`, `metrics_s1_s2.csv`, `pitch_tracks.csv` |
| `build_pdf.py` | report `.md` → `.pdf` (needs Edge and internet for MathJax) |

Data is read from `../speech_samples/`.

## 1. Regenerate results and figures

Git Bash, from `ADSP-2/v3/`. Takes about a minute.

```bash
matlab -batch "D2_code_2026EEY7565"
```

Verify: last line printed is `Figures written to assets`, and `assets/results_log.txt` has a fresh timestamp.

```bash
ls -la assets/results_log.txt
```

## 2. Rebuild the PDF

Only needed after editing the `.md`.

```bash
PYTHONUTF8=1 uv run --no-project --with markdown-it-py --with mdit-py-plugins python build_pdf.py
```

Verify: prints `wrote ...D2_report_2026EEY7565.pdf`.

## 3. If the numbers change

Every number in the report comes from `assets/results_log.txt`. After changing a parameter in Section 2 of the script, rerun step 1, compare the log with the tables in the `.md`, edit, then rerun step 2.
