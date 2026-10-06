# Study path — understand Design Assignment 2 in about 90 minutes

Three things carry the same content at three depths. Use them together, in this order.

| Layer | What | Where |
|---|---|---|
| **Watch** | 16 short animations, picture first, formula last (10.5 min in total) | `C:\Users\NAVY\manimations\DA-2_anima\media\480p15\` |
| **Read** | the report, same story with the numbers | `D2_report_2026EEY7565.pdf` in this folder |
| **Run** | the Live Script, the actual code with its outputs | `D2_code_2026EEY7565.mlx` in this folder |

Rule for every phase: **watch, then read, then find it in the code, then answer the question without looking.** If the question beats you, watch the scene again, not the next one.

Short on time: phases 2, 3 and 7 are the assignment. About 35 minutes.

## Phase 0 — What is being asked (5 min)

| Do | Item |
|---|---|
| Watch | `S01_Assignment.mp4` |
| Read | report Section 1 |

> [!question]- What does a pitch detector output, and how often?
> One number every 10 ms: the pitch in Hz if the frame is voiced, 0 if it is unvoiced or silent. So two decisions per frame: *is there a voice*, and *how fast does it repeat*. 300 frames per recording.

## Phase 1 — Where pitch comes from, and frames (10 min)

| Do | Item |
|---|---|
| Watch | `S02_SpeechModel.mp4`, `S03_Framing.mp4` |
| Read | report Section 2, the table |
| Code | Live Script Section 4, first six lines of `frame_correlation` (cut the frame, remove the mean, window it) |

> [!question]- Why does 8000 / k give the pitch in Hz?
> k is the period counted in samples and there are 8000 samples in a second, so the pattern repeats 8000 / k times per second. Period 60 samples → 133 Hz.

> [!question]- Why a 40 ms frame and not 10 ms?
> The frame must hold a few repeats of the pattern, or there is nothing to compare. The lowest voice searched is 70 Hz, a period of 14 ms; 40 ms holds almost three of those.

## Phase 2 — Autocorrelation (20 min) ★

| Do | Item |
|---|---|
| Watch | `S04_SlideMultiply.mp4` (twice if needed), then `S05_Taper.mp4` |
| Read | report Section 3.1, Figures 1 and 2 |
| Code | Live Script Section 4, the `for j` loop in `frame_correlation` |

Symbol to picture to code:

| In the formula | In the animation | In the code |
|---|---|---|
| $y(n)$ | the white frame | `a = frame(1:end-k)` |
| $y(n+k)$ | the blue copy, slid by k | `b = frame(k+1:end)` |
| $\sum$ of the products | the green and red area, added up | `sum(a.*b)` |
| divide by the energies | S05: a fair score that does not shrink with lag | `/sqrt(sum(a.^2)*sum(b.^2))` |
| $r(k)$ | one dot on the bottom curve | `r(j)` |

One difference to know: S05 corrects the shrinking overlap by multiplying with N/(N−k). The Live Script gets the same effect by dividing by the energies of the two overlapping parts. Same purpose, different arithmetic.

`S06_CentreClip.mp4` is **optional**. Centre clipping is not in our code. Watch it for one thing only: it shows that the plain autocorrelation has extra bumps caused by the mouth. Those bumps are the reason for the doubling error in phase 7.

> [!question]- In Figure 1, bottom panel, there is a second tall peak near lag 74. What is it?
> Twice the period. If the frame repeats every 37 samples it also repeats every 74. That is why the peak picker takes the *earliest* peak that is nearly as tall as the tallest (the 82 % rule in `best_period`).

> [!question]- Why does a quiet frame get 0 Hz even if its correlation curve had a peak?
> `track_pitch` checks the RMS first (`energy < .008`) and skips the frame. The correlation is normalised, so loudness does not show in it; silence has to be caught separately.

## Phase 3 — Cepstrum (25 min) ★

Watch these five in one sitting. Each one is a single step.

| Scene | The one idea |
|---|---|
| `S07_TimeFrequency.mp4` | pulses every P in time = harmonics every 1/P in frequency |
| `S08_CombTimesEnvelope.mp4` | a voiced spectrum is a comb (pitch) × an envelope (mouth) |
| `S09_LogSplits.mp4` | the log turns × into +: slow envelope + fast ripple |
| `S10_CepstrumProbe.mp4` | lay a cosine along the frequency axis, multiply, add: a spike when its crests sit on the harmonics |
| `S11_ReadTheCepstrum.mp4` | low quefrency = mouth, spike = period, no spike = unvoiced |

| Do | Item |
|---|---|
| Read | report Section 4.1 (its paragraphs follow S08, S09, S10 in that order), Figure 3 |
| Code | Live Script Section 6, `frame_cepstrum`: three lines |

Symbol to picture to code:

| In the formula | In the animation | In the code |
|---|---|---|
| $\lvert X(f)\rvert$ | S08: the spectrum of the frame | `abs(fft(frame,1024))` |
| $\log\lvert X(f)\rvert$ | S09 and S10: the wavy white line | `log(...+eps)` |
| $\cos(2\pi f q)$ and the sum | S10: the pink cosine, the green and red area | `ifft(...)` does it for every q at once |
| $c(q)$ | S10: one dot on the bottom curve | `c(lags+1)` |
| $q$ | how tightly the cosine's crests are packed; a time | the same `lags` as the autocorrelation |

S10 is S04 again on a different curve. In S04 you slide a copy along **time**. In S10 you stretch a cosine along **frequency**. Multiply, add, look for the peak: the same three moves.

> [!question]- Quefrency is measured in which unit, and why?
> Seconds (here, samples). The log spectrum is a curve along frequency in Hz. "How often does it ripple per Hz" has the unit 1/Hz, which is seconds. A ripple every 133 Hz is a quefrency of 1/133 s = 7.5 ms, the pitch period.

> [!question]- Why does an unvoiced frame have no cepstral peak?
> No regular pulses, so no harmonics, so no regular ripple on the log spectrum. Nothing for any cosine to line up with.

> [!question]- Why is the search started at 2.5 ms and not at 0?
> Small quefrencies hold the slow shape of the mouth, which is large (S11, the orange band). Skipping them is called liftering. 2.5 ms is also 400 Hz, the highest pitch searched.

## Phase 4 — Why the two are relatives (5 min)

| Do | Item |
|---|---|
| Watch | `S12_OneRecipe.mp4` |
| Read | report Section 4.3 |

> [!question]- Both methods are "an inverse FFT of something". Of what?
> Autocorrelation: of the squared spectrum. Cepstrum: of the log spectrum. Squaring lets the strongest parts dominate (sturdy in noise, but the mouth's resonances leak in). The log levels everything (clean sharp peak, but noisy parts count as much as strong ones). Every result in phase 7 follows from this.

## Phase 5 — From 300 guesses to a track (10 min)

| Do | Item |
|---|---|
| Watch | `S13_DecideAndSmooth.mp4` |
| Read | report Section 2: peak picker, quiet-frame test, smoothing |
| Code | `track_pitch`, `track_pitch_cepstrum`, `best_period`, `smooth_track` |

The animation shows an earlier version of the detector. The three steps are the same in the Live Script: quiet-frame test (RMS below 0.008), peak-height threshold (0.50 autocorrelation, 0.12 cepstrum), 5-frame median.

Put `track_pitch` and `track_pitch_cepstrum` side by side. They differ in **two things only**: the curve handed to `best_period`, and the threshold. That is the whole design of the comparison.

> [!question]- What does the median filter fix, and what can it not fix?
> One wrong frame between right ones. It cannot fix a run of wrong frames: four doubled frames in a row survive a 5-frame median (s2, phase 7).

## Phase 6 — Scoring (5 min)

| Do | Item |
|---|---|
| Watch | `S14_Scoring.mp4` |
| Read | report Section 2, error measures |
| Code | `score_track` |

The animation and the report use different names for the same things:

| Animation | Report and code |
|---|---|
| V→UV, voice missed | FN |
| UV→V, voice invented | FP |
| gross error | a TP frame more than 20 % off |
| FFE | frame error |

## Phase 7 — Results and why (15 min) ★

| Do | Item |
|---|---|
| Watch | `S15_Verdict.mp4` |
| Read | report Sections 5.2, 5.3, 7 and 8; Figures 4, 5 and 8 |
| Run | the Live Script; find the score table and the noise table in its output |

Four numbers to know: clean speech **7.3 %** (autocorrelation) against **3.5 %** (cepstrum); at 10 dB SNR **11.9 %** against **14.2 %**.

Two autocorrelation errors, each explained by a scene you have watched:

| Error | Where | Why | Scene |
|---|---|---|---|
| Misses the low-pitched ending of s2 | Figure 5, after 2.2 s | long lag through a tapered window: even a perfect 85 Hz tone only scores 0.69, close to the 0.50 threshold | S05 |
| Doubles the pitch, 324 Hz instead of 165 Hz | Figure 5, near 0.26 s | a strong second harmonic gives a peak at half the lag; the "earliest strong peak" rule takes it | S06, S12 |

> [!question]- Why does the cepstrum win on clean speech but lose in noise?
> Same cause both times: the log. It levels the spectrum, so the mouth's resonances cannot mislead it (clean: fewer errors). But it also gives the noisy gaps between harmonics full weight (noise: the peak sinks below the threshold).

## Phase 8 — s3 and s4 (5 min)

| Do | Item |
|---|---|
| Watch | `S16_S3S4.mp4` for the idea only |
| Read | report Section 6, Figures 6 and 7 |

The animation draws the tracks of the earlier detector, so its numbers and its circled errors differ from the report. The idea carries over unchanged: with no answer key, two methods that fail in different ways check each other.

> [!question]- s4: the autocorrelation ends near 385 Hz, the cepstrum near 105 Hz. Which is right, and how do you know without a reference?
> The cepstrum. The pitch has been falling smoothly to about 105 Hz just before; a jump to 385 Hz in one frame is not something a voice does. 385 Hz is also the short end of the search range, where a detector lands when it has lost the period.

## Phase 9 — Say it aloud (5 min)

Answer each in two sentences, no notes.

1. What does autocorrelation measure, and where is the pitch in it?
2. What is a cepstrum, in one line of operations?
3. What is the single difference between the two methods in our code?
4. Which is better on clean speech, by how much, and why?
5. Which is better in noise, and why?
6. What are the two weaknesses of our comparison? (thresholds tuned on the scored files; only 600 frames from two speakers)

If all six come out clean, the assignment is understood.
