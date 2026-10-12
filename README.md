# Characterizing sleep stages through the complexity–weighted-entropy plane in human intracranial data and in a whole-brain model

Code accompanying the paper:

**Characterizing sleep stages through the complexity–weighted-entropy plane in human intracranial data and in a whole-brain model**
Helena Bordini de Lucas, Leonardo Dalla Porta, Alain Destexhe, Maria V. Sanchez-Vives, Osvaldo A. Rosso, Cláudio R. Mirasso and Fernanda Selingardi Matias.

The paper characterizes wakefulness and the REM, N2 and N3 sleep stages in intracranial EEG (iEEG) from 65 patients, using two symbolic quantifiers computed from ordinal patterns:

- **Weighted permutation entropy (WPE)**, which measures how diverse the ordinal patterns of a signal are, while taking the amplitude of the fluctuations into account.
- **Statistical complexity measure (SCM)**, which measures how structured the pattern distribution is, by comparing it with the equiprobable one.

Plotting the two against each other gives the complexity–entropy plane, where the four states follow the trajectory Wake → REM → N2 → N3. The paper also compares this trajectory with a whole-brain mean-field model in which neuronal adaptation is varied.

## Citation

If you use this code, please cite the paper above.
---

## Repository contents

| File | What it does |
| --- | --- |
| `PE_func.m` | Conventional permutation entropy (Bandt–Pompe), statistical complexity and Fisher information of a time series. |
| `WPE_func.m` | Same as above, but with weighted permutation entropy, where each ordinal vector is weighted by its variance (Eqs. 4, 5 and 7 of the paper). |
| `Cmax_Cmin.m` | Computes and plots the theoretical minimum and maximum complexity curves (black curves in Figs. 5, 6 and 9). |
| `WPE_PE_SCM_separate_channels.m` | Computes entropy, complexity and Fisher information for every channel and every brain state. |
| `WPE_PE_SCM_mean_series.m` | Averages all channels of each brain state into one series and computes the same quantifiers on it. |
| `TVB_Human_SleepStages.ipynb` | Whole-brain simulations in The Virtual Brain (see [Simulations](#simulations)). |
| `README.md` | This file. |

### The quantifier functions

`WPE_func` and `PE_func` share the same interface:

```matlab
[out, info] = WPE_func(x, tauList, D)
[out, info] = PE_func(x, tauList, D)
```

| Argument | Meaning |
| --- | --- |
| `x` | Time series, passed with the zero-amplitude buffers still in place. |
| `tauList` | Vector of embedding delays, in samples. |
| `D` | Embedding dimension. The paper uses `D = 6`. |

`out` has one row per delay and the columns `[H, C, F]`: normalized entropy, statistical complexity and Fisher information. The paper uses H and C only.

Recordings in the database are made of several short segments separated by zero-amplitude buffers. Both functions detect these buffers (runs of at least 10 zeros, or 50 ms at 200 Hz) and use only ordinal vectors that lie entirely inside one continuous segment, so that no artificial transition between separated portions of the recording is created. Segments can also be given explicitly with the `'SegmentId'` or `'ValidMask'` options. `info` reports how many segments were found and how many ordinal vectors were used and discarded for each delay.

### The analysis scripts

Both scripts read `MatlabFile_Filtered.mat` (see [Data](#data)), process the four states in turn, and save one `.mat` file per method and state. Each series is truncated to its first 10,800 valid samples (54 s at 200 Hz), so that all channels and states have the same length.

**`WPE_PE_SCM_separate_channels.m`** treats each of the 1,012 channels as its own series.

- Output: `HC_separate_channels_<method>_<state>.mat`
- Main variable: `results_all_channels`, with dimensions delays × 3 × channels, where the middle dimension is `[H, C, F]`.
- Also saves the number of segments, valid samples and ordinal vectors used and excluded for each channel.

**`WPE_PE_SCM_mean_series.m`** first averages the channels and then computes the quantifiers. The set of channels contributing to the average changes whenever a channel enters or leaves a buffer, so a new segment is started at each change and no ordinal vector spans two different sets of channels.

- Output: `HC_mean_series_<method>_<state>.mat`
- Main variable: `HCF`, with dimensions delays × 3, plus the `info` structure.

In both scripts, `<method>` is `WPE` or `PE` (the conventional version, used as a control analysis), and `<state>` is one of `data_W`, `data_N2`, `data_N3` and `data_R`, for Wake, N2, N3 and REM.

### Theoretical complexity bounds

`Cmax_Cmin.m` builds the curves that bound the statistical complexity for a given embedding dimension: the minimum curve from a family of distributions with one dominant symbol, and the maximum curve from distributions with a block of zero-probability symbols (Martin, Plastino and Rosso, 2006; Lamberti et al., 2004). Data points in the plane always fall between the two curves.

## Data

The recordings come from the public **MNI Open iEEG Atlas** (https://mni-open-ieegatlas.research.mcgill.ca), which contains iEEG from patients with drug-resistant focal epilepsy, restricted to channels with probably normal activity. The paper uses the 1,012 channels (65 patients) that have data in all four states: Wake, REM, N2 and N3.

The analysis scripts expect a preprocessed file, `MatlabFile_Filtered.mat`, that is **not included in this repository** because of its size (about 170 MB). The original atlas data contains up to 1,772 channels from 106 patients, but not every channel was recorded in every state. Only the channels with data in all four states (Wake, REM, N2 and N3) were kept, which leaves 1,012 channels from 65 patients. The file contains four matrices, `data_W`, `data_N2`, `data_N3` and `data_R`, one per state. Each matrix has one row per sample and one column per channel (1,012 columns, in the same channel order in all four matrices).

## Simulations

`TVB_Human_SleepStages.ipynb` runs a 68-region whole-brain model in [The Virtual Brain](https://www.thevirtualbrain.org), using the `ZerlautAdaptationSecondOrder` mean-field model of adaptive exponential integrate-and-fire populations, coupled through a human structural connectome. The level of neuronal adaptation `b` is varied to move the model between asynchronous, wake-like dynamics and slow-wave, sleep-like dynamics. For each value of `b`, the notebook saves a file `DataModel_<b>.pkl` with the simulated time and the firing-rate time series of the 68 regions.

> **Note.** The notebook comes from an earlier, unrevised stage of the work. The adaptation values, the noise settings and some of the model parameters it contains may therefore differ from the final ones reported in the paper. The paper is the reference for the parameters of the published analysis.

To run it you need Python 3.9, NumPy, The Virtual Brain (`tvb-library`) and Jupyter, plus a TVB connectivity file named `Connectivity.zip` with the 68-region connectome. Set `path_c` in the notebook to the location of that file. The simulated series are saved as Python pickles, which MATLAB cannot read directly, so they must be converted before being analyzed with the scripts above.

## Requirements

- **MATLAB** for the analysis scripts and the quantifier functions. They use only base MATLAB, with no extra toolboxes.
- **Python 3.9** with `tvb-library`, NumPy and Jupyter for the simulations.

## How the pieces fit together

```text
MNI Open iEEG Atlas
        |
        |  preprocessing (resampling, buffers, channel selection)
        v
MatlabFile_Filtered.mat
        |
        +--> WPE_PE_SCM_separate_channels.m --> HC_separate_channels_<method>_<state>.mat
        |
        +--> WPE_PE_SCM_mean_series.m -------> HC_mean_series_<method>_<state>.mat
                                                          |
                          Cmax_Cmin.m (bounds) -----------+--> complexity–entropy plane

TVB_Human_SleepStages.ipynb --> DataModel_<b>.pkl --> (conversion) --> same quantifiers
```

## Status

This repository currently contains the quantifier functions, the theoretical bounds, the channel-level and mean-series computations, and the simulation notebook. The scripts that produce the patient-level and regional statistics and the figures of the paper are not included yet.

