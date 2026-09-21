# Binaural Beats and Low-Gamma EEG in Healthy Older Adults

This repository contains the MATLAB code used for the analysis of EEG and
behavioral data reported in:

Effects of Low-Gamma Binaural Beats on Visual-Spatial Working Memory, Selective Attention, and Processing Speed in Healthy Older Adults: An EEG Study

## Study overview

This study investigated the acute and post-stimulation effects of
34–43 Hz binaural-beat stimulation in healthy older adults.

EEG was recorded during three experimental epochs:

- PRE: baseline
- DURING: stimulation
- POST: post-stimulation recovery

## Analyses

The repository contains code for:

1. Whole-scalp 34–43 Hz power analysis
2. Channel-wise 34–43 Hz analysis
3. Exploratory broadband EEG analysis
4. Machine-learning classification using three feature representations
5. Feature-selection stability analysis
6. Welch PSD analysis
7. Whole-record FFT sensitivity analysis

## Feature representations

- Level 1: whole-scalp 34–43 Hz power
- Level 2: 30 channel-resolved 34–43 Hz features
- Level 3: 150 channel × frequency-band features

## Statistical analysis

The analyses include Welch's t-tests, FDR correction,
repeated stratified five-fold cross-validation, and permutation testing.

## Software

All analyses were performed in MATLAB.

MATLAB version: 2024b

Additional toolboxes/packages:
- EEGlab for data preprocessing

## Data availability

Participant-level EEG data are not included in this repository because
of [ethical/privacy/data-sharing restrictions].

## Code availability

The MATLAB scripts required to reproduce the reported analyses are
provided in this repository.
