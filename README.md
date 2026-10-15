# Non-Contact Respiratory-Rate Detection in Dairy Goats Using FMCW Radar

MATLAB implementation of the signal-processing framework described in:

> *Precision Non-Contact Detection of Respiratory Rate in Dairy Goats Using FMCW Radar and Adaptive Signal Processing.*

The framework estimates the respiratory rate (RR) of dairy goats from 77-GHz frequency-modulated continuous-wave (FMCW) radar data, using a single receive channel (RX1) and a three-stage processing pipeline that mitigates motion interference and automates respiratory-component extraction.

## Pipeline overview

1. **Preprocessing** — raw DCA1000 samples are de-interleaved, RX1 is retained, and the slow-time sampling rate is decimated by a factor of two. A Hamming window and a range FFT are applied, static clutter is removed by slow-time mean subtraction, the first five range bins are zeroed, and the target range bin is selected by the maximum-energy criterion. The phase is then extracted and differentiated (centered difference).
2. **Body-motion detection and restoration** — a short-time Fourier transform (STFT) of the phase-difference signal yields seven spectral-morphology features (spectral flatness, spectral entropy, spectral flux, spectral bandwidth, spectral crest, 95% roll-off frequency, and delta log-energy), which are robustly standardized and combined into a composite score. A dual-threshold hysteresis rule produces a binary motion mask. Contaminated intervals are reconstructed by interpolating the instantaneous frequency estimated from adjacent valid segments.
3. **Adaptive narrowband filtering** — a smoothing window is chosen by maximizing the local SNR, the dominant respiratory frequency is estimated within [0.1, 2.0] Hz, and a fourth-order Butterworth band-pass filter with a fixed ±0.25 Hz half-bandwidth is applied (zero-phase).
4. **AIMF-EMD** — empirical mode decomposition (EMD) is followed by an adaptive IMF selection based on three criteria (energy consistency, spectral alignment, and autocorrelation-based periodicity stability) combined multiplicatively. The RR is obtained from the periodogram peak of the selected IMF.

## File structure

| File | Description |
| --- | --- |
| `main.m` | Entry script: runs the full pipeline on one `.bin` file and prints the RR. |
| `readDCA1000.m` | Reads a DCA1000 raw binary, de-interleaves, keeps RX1, decimates slow time. |
| `rangePhaseExtract.m` | Range FFT, clutter removal, maximum-energy range-bin selection, phase extraction. |
| `detectMotion.m` | Seven-feature STFT motion detection with robust standardization and hysteresis. |
| `restoreMotionPhase.m` | Frequency-phase trajectory reconstruction of motion-corrupted segments. |
| `adaptiveNarrowband.m` | Adaptive narrowband filtering (smoothing-window selection, ±0.25 Hz band). |
| `aimfEmd.m` | EMD decomposition and AIMF-based IMF selection, periodogram RR estimation. |

## Parameter settings

The parameter values below correspond to those reported in the manuscript (Section 3.3.4, Table 4).

| Parameter | Value | Description |
| --- | --- | --- |
| `Q_hi` | 3 | high hysteresis coefficient |
| `Q_lo` | 1.5 | low hysteresis coefficient |
| `delta` | 0.25 Hz | adaptive-filter half-bandwidth |
| `eps_E` | 0.01 | energy-consistency epsilon |
| `eps_f` | 0.001 | spectral-alignment epsilon |
| `eps_p` | 1e-3 | periodicity epsilon |

The following values are tunable and are not stated explicitly in the manuscript text (they appear only in the pseudocode of Fig. 6):

| Parameter | Default | Description |
| --- | --- | --- |
| `w` | `[1 1 1 1 -1 1 1]` | weights of the seven spectral features (crest is negative because it decreases under motion) |
| `stftWin` | 64 | STFT frame length (samples) |
| `stftHop` | 8 | STFT hop (samples) |
| `movWin` | 9 | moving-median/MAD window for robust standardization |

## Requirements

- MATLAB R2019a or later
- Signal Processing Toolbox (`stft`, `emd`, `movmedian`, `movmad`, `smoothdata`, `butter`, `filtfilt`, `findpeaks`, `periodogram`)

## Usage

1. Set `Filename` in `main.m` to the path of a DCA1000 raw `.bin` file.
2. Run `main.m`.

```matlab
>> main
主频: 0.244 Hz, 通带: 0.100-0.494 Hz
呼吸频率: 0.220 Hz (13.18 次/分钟)
```

## Data format

The raw `.bin` files are complex interleaved 16-bit samples recorded by a TI DCA1000 board with the IWR1642BOOST radar: 256 ADC samples per chirp, 2 chirps per frame, 1024 frames, 1 TX / 4 RX, sampled at 10 MS/s. The pipeline retains only RX1 and decimates the slow-time axis by two, yielding a sampling rate of 25 Hz over a 40.96-s window.

## Raw data availability

The raw FMCW radar data that support this study are available from the corresponding author upon reasonable request. The source code is provided in this repository.
