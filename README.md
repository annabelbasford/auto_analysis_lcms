# auto_analysis_lcms

Automated analysis pipeline for LCMS data using UV and MS spectra. 
This project identifies peaks in UV vs RT trace, aligns them with TIC vs RT trace, then converts the mzML data to extracts relevant mz vs intensity spectra, and visualizes the results.

Currently the project is limited to manual extraction of DAD1 UV trace as a CSV per each Agilent .d file. 
Each Agilent .d file should be converted to an mzML file and kept in /mzML/ folder. This was omitted from github due to large file sizes.

---

## 📁 Folder Structure

auto_analysis_lcms/
─ DAD_uv/*_uv.csv # Input UV vs RT
─ peak.csv # MS scan data per peak
─ *_uv_peak_rt.csv # Peak summary per sample for UV and TIC
─ pngs/ # Output: plots per peak of UV data
─ peak_mz_refined/ # Output: refined CSV spectra + plots per peak
─ *.py # Analysis scripts
---

## 🚀 Workflow Overview

1. **UV Peak Detection**
   - 01_get_peak_rt.py
   - Detects peaks in UV data (`*_uv.csv`)
   - For each peak gives the peak_rt_uv, peak_start_uv, peak_end_uv, peak_area, peak_rt_tic, peak_start_tic, peak_end_tic
   - Saves a CSV of peak RTs (`*_uv_peak_rt.csv`)
  
2. **mzMLto CSV**
   - 02_mzML2csv_lcms.py
   - Saves a CSV per peak of 'ScanNumber', 'RetentionTime', 'm/z', 'Intensity' for RT values between the peak_start_tic and peak_end_tic
   - Writes to *_peak_*.csv

3. **MS Scan Refinement**
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


## 🧪 Notes

* Peaks with area <1% of the largest are filtered out
* `peak_rt_tic` values are UV RTs adjusted to align with TIC scans
* Raw data folders like `agilent_d/` and `mzML/` are ignored due to size

---

## 👩‍🔬 Author

Annabel Basford
[GitHub: @annabelbasford](https://github.com/annabelbasford)

```

---

Let me know if you'd like this written directly to a new `README.md` file or want to add sample plots/screenshots.
```

