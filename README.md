# ﻿Super-Heterodyne Receiver Simulation
## Overview
This project simulates the fundamental components of an analog communication system using MATLAB. It focuses on creating an AM (Amplitude Modulation) modulator and a corresponding super-heterodyne receiver, simulating the process from modulation to baseband detection.

### Features
- AM Modulation: Implemented using Double-Sideband Suppressed Carrier (DSB-SC) with a variety of carrier frequencies.
- Super-Heterodyne Receiver: Includes an RF stage, mixer, IF stage, and baseband detector, each simulated using MATLAB.
- Spectrum Analysis: Frequency domain analysis using FFT to estimate signal bandwidth and design filters accordingly.
- Signal Processing: Handles stereo audio signals and implements monophonic reception by combining stereo channels.

### Components
- AM Modulator: Generates modulated signals using custom carrier frequencies.
- RF Stage: Simulated as a tunable band-pass filter.
- Oscillator & Mixer: Creates the intermediate frequency (IF) and mixes it to extract the baseband signal.
- IF Stage: Another band-pass filter centered at the IF.
- Baseband Detection: Final stage of signal recovery using low-pass filtering.

### Instructions
- Audio Input: Stereo signals are provided as .wav files. These are processed in MATLAB using functions like audioread and wavread.
- Simulation Parameters: Ensure carrier frequencies do not violate Nyquist criteria. If necessary, adjust sampling frequency.
- Visualization: Use MATLAB commands such as fft and plot for spectral analysis. Adjust axis scaling to view signal spectra versus frequency.
- Deliverables: The project report should include detailed explanations of each block's design, MATLAB code, and analysis results.

### Usage
Run the MATLAB script provided to simulate the communication system. Follow the comments within the code for detailed steps.


