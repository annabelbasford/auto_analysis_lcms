import csv
import os
import pandas as pd
from pyteomics import mzml

# Paths
input_folder = '/auto_analysis/mzML/'
output_folder = '/auto_analysis/'
# load all mzML files
filenames = [f for f in os.listdir(input_folder) if f.endswith('.mzML')]


for file in filenames:
    # Load peak ranges from corresponding peak CSV
    base_name = os.path.splitext(file)[0]
    peak_csv_file = os.path.join(output_folder, f'{base_name}_uv_peak_rt.csv')
    peak_df = pd.read_csv(peak_csv_file)
    peak_ranges = list(zip(peak_df['peak_start_tic'], peak_df['peak_end_tic']))

    # Open mzML file
    input_mzml_file = os.path.join(input_folder, f'{base_name}.mzML')
    with mzml.read(input_mzml_file) as mzml_file:
        # Read all spectra once and store them
        spectra = list(mzml_file)

    intensity_threshold = 1000  # Set intensity threshold for filtering

    # For each peak, extract relevant scans and write to a separate CSV
    for i, (start_rt, end_rt) in enumerate(peak_ranges):
        output_csv_file = os.path.join(output_folder, f'{base_name}_peak_{i+1}.csv')
        with open(output_csv_file, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['ScanNumber', 'RetentionTime', 'm/z', 'Intensity'])

            for spectrum in spectra:
                rt = spectrum['scanList']['scan'][0]['scan start time']
                if start_rt <= rt <= end_rt:
                    scan_number = spectrum['id']
                    mzs = spectrum['m/z array']
                    intensities = spectrum['intensity array']
                    filtered = [(mz, intensity) for mz, intensity in zip(mzs, intensities) if intensity >= intensity_threshold]

                    for mz, intensity in filtered:
                        writer.writerow([scan_number, rt, mz, intensity])

        print(f"Saved {base_name}_peak_{i+1}.csv (RT {start_rt:.2f}–{end_rt:.2f} min)")
