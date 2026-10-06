# Design Assignment 2 — v4 (integrated)

One Live Script with **both** methods, run on the same frames:

| Method | Domain | Code |
|---|---|---|
| Autocorrelation | time | partner's functions, copied line for line: `track_pitch`, `frame_correlation`, `best_period`, `show_method` |
| Cepstrum | frequency | `track_pitch_cepstrum`, `frame_cepstrum`: reuse his frame, quiet-frame rule and peak picker; only the searched curve differs |
| Shared | both | `smooth_track` (5-frame median), `score_track`, noise test, plots |

No toolbox needed (base MATLAB only).

| File | What it is |
|---|---|
| `D2_code_2026EEY7565.mlx` | **the Live Script**, explanations, code and outputs embedded |
| `D2_report_2026EEY7565.pdf` | **the report**, built from `D2_report_2026EEY7565.md` with `build_pdf.py` |
| `STUDY-PATH.md` | how to learn the assignment: which animation, report section and code function go together |
| `results/` | `s1..s4_pitch.csv`, `p3.mat`, `p4.mat`, `report_plots/01..08_*.png` |

Data is read from `../speech_samples/`.

## 1. Run

Open `D2_code_2026EEY7565.mlx` in MATLAB and press **Run**. About one minute.

Or from Git Bash in `ADSP-2/v4/`:

```bash
matlab -batch "D2_code_2026EEY7565"
```

Verify: the pooled rows of the score table read

```
s1+s2  autocorrelation  median     307  15  24  254      95.3%   92.7%      5.00     302/331         7.3%
s1+s2  cepstrum         median     322  10   9  259      97.0%   97.3%      3.07     320/331         3.5%
```

and the partner's own numbers reappear in the `none` rows (s1: TP=147 FP=13 FN=9 TN=131; s2: TP=160 FP=2 FN=15 TN=123).

## 2. Do not keep a `.m` of the same name in this folder

MATLAB then cannot tell the two apart and the Live Script fails with "Line numbers may not exceed the last expression". If a plain-text copy is needed, use **Save As** to another name or folder.
