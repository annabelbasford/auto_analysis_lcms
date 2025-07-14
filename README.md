# auto_analysis_lcms

Automated analysis pipeline for LCMS data using UV and MS spectra. 
This project identifies peaks in UV vs RT trace, aligns them with TIC vs RT trace, then converts the mzML data to extracts relevant mz vs intensity spectra, and visualizes the results.

Following this, each extracted peak mv vs intensity data is screened for possible poc structures and intermediates, assigning by numbe rof each precursor and the number of imine bonds formed.

Currently the project is limited to manual extraction of DAD1 UV trace as a CSV per each Agilent .d file. 
Each Agilent .d file should be converted to an mzML file and kept in /mzML/ folder. This was omitted from github due to large file sizes.

Future work will focus on adding the script used to automate the file conversion between .d and .mzML (code adapted from https://github.com/GreenawayLab/cagey)

---

## 📁 Folder Structure

---

auto_analysis_lcms/
- reaction_planner.json # JSON of each reaction with metadata
- assigned_peaks_data/ # CSVs of the mass_assignment and the processed CSV + annotated spectra plotted of the assigned mz vs intensity data
- DAD_uv/*_uv.csv # Input UV vs RT
- mzML/ *.mzML # mzML files converted from agilent .d
- peaks_mz_csv/ # outputs from mzML to CSV within the rt range for each peak
- peaks_mz_refined/ # filtered CSV for each peak + the plotted data of mz vs intensity for that peak
- pngs/ # Output: annotated plots per peak of UV data
- *_uv_peak_rt.csv # Peak summary per sample for UV and TIC
- *.py # Analysis scripts
"""
---

## 🚀 Workflow Overview

1. **UV Peak Detection**
   - 01_get_peak_rt.py
   - Detects peaks in UV data (`*_uv.csv`)
   - For each peak gives the peak_rt_uv, peak_start_uv, peak_end_uv, peak_area, peak_rt_tic, peak_start_tic, peak_end_tic
   - Saves a CSV of peak RTs (`*_uv_peak_rt.csv`) in the parent folder
   - Plots annotated spectra of RT vs Intensity with the RT and area in /pngs
  
2. **mzMLto CSV**
   - 02_mzML2csv_lcms.py
   - Saves a CSV per peak of 'ScanNumber', 'RetentionTime', 'm/z', 'Intensity' for RT values between the peak_start_tic and peak_end_tic
   - Writes to *_peak_*.csv in /peaks_mz_csv

3. **MS Scan Refinement**
   - 03_refined_cav.py
   - For each peak (`*_peak_*.csv`), finds scan closest to `peak_rt_tic`
   - Saves cleaned spectra to `*_refined.csv` in /peaks_mz_refined
   - Plots m/z vs intensity for each refined scan in /peaks_mz_refined
  
4. **Mass Assignment per Peak**
   - 04_mass_assignment.py
   - loads the reaction data from JSON (reaction_planner.json)
   - For each csv in /peaks_mz_refined it extracts the base_name and peak_number to find a match with the base_name in reactions JSON
   - Calculates the topicity per precursor based of aldehyde and amine SMARTs 
   - Then finds `'num_components'` in JSON and to choose `poc_find_solutions` for 2 components or `poc_find_solutions_ternary` for 3 components
   - Calculates the top 10 isotopes for each solution and screens the refined_csv for any matches, filtering by charge match
   - Outputs in a CSV which is then reordered by the most intense peak
   - Plots the refined_csv of mz vs Intensity with annotations of the matched solutions + rdkit molecules of each precursor

---

## 📦 Requirements

- Python 3.7+
- pandas
- numpy
- matplotlib
- scipy
- pyopenms
- rdkit


## 🧪 Notes

* Peaks with area <1% of the largest are filtered out
* `peak_rt_tic` values are UV RTs adjusted to align with TIC scans
* Raw data folders like `agilent_d/` and `mzML/` are omitted in github due to size

---

## 👩‍🔬 Authors

Annabel Basford
[GitHub: @annabelbasford](https://github.com/annabelbasford)
Benjamin Egleston
[GitHub: @benjybenj](https://github.com//benjybenj)


```

