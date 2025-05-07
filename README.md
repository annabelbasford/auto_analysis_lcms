from pathlib import Path

readme_content = """
# auto_analysis_lcms

Automated analysis pipeline for LC-MS data using UV and MS spectra. This project identifies peaks in UV chromatograms, aligns them with MS scan data, extracts relevant spectra, and visualizes the results.

---

## 📁 Folder Structure

auto_analysis_lcms/
├── *_uv.csv # Input UV chromatograms
├── peak.csv # MS scan data per peak
├── *_uv_peak_rt.csv # Peak summary per sample
├── peak_mz_refined/ # Output: refined MS spectra + plots
└── *.py # Analysis scripts
---

## 🚀 Workflow Overview

1. **UV Peak Detection**
   - Detects peaks in UV data (`*_uv.csv`)
   - Calculates corrected peak area
   - Saves plots and a CSV of peak RTs (`*_uv_peak_rt.csv`)

2. **MS Scan Refinement**
   - For each peak (`*_peak_*.csv`), finds scan closest to `peak_rt_tic`
   - Saves cleaned spectra to `*_refined.csv`
   - Plots m/z vs intensity for each refined scan

---

## 📦 Requirements

- Python 3.7+
- pandas
- numpy
- matplotlib
- scipy
