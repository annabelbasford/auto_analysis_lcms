import os
import pandas as pd
import csv
import matplotlib.pyplot as plt

# Input/output folders
parent_folder = '/home/abasford/projects/auto_analysis_lcms/'
input_folder = '/home/abasford/projects/auto_analysis_lcms/peaks_mz_csv/'
output_folder = os.path.join(parent_folder, 'peaks_mz_refined')
os.makedirs(output_folder, exist_ok=True)

# Get all per-peak CSV files (e.g., sample_peak_1.csv)
peak_csvs = [f for f in os.listdir(input_folder) if f.endswith('.csv') and '_peak_' in f]

for peak_file in peak_csvs:
    peak_file_path = os.path.join(input_folder, peak_file)
    try:
        # Load peak data
        df = pd.read_csv(peak_file_path)

        # Check required columns exist
        required_columns = {'ScanNumber', 'RetentionTime', 'm/z', 'Intensity'}
        if df.empty or not required_columns.issubset(df.columns):
            print(f"Skipping {peak_file}: missing columns or empty.")
            continue

        # Parse base name and peak index
        base_name = peak_file.split('_peak_')[0]
        peak_idx = int(peak_file.split('_peak_')[1].split('.')[0]) - 1

        # Load peak_rt_tic from UV peak summary
        peak_rt_file = os.path.join(parent_folder, f'{base_name}_uv_peak_rt.csv')
        if not os.path.exists(peak_rt_file):
            print(f"Missing peak_rt file for {base_name}, skipping.")
            continue

        peak_rt_df = pd.read_csv(peak_rt_file)
        if peak_idx >= len(peak_rt_df):
            print(f"Peak index {peak_idx} out of bounds in {peak_rt_file}, skipping.")
            continue

        peak_rt_val = peak_rt_df['peak_rt_tic'].iloc[peak_idx]

        # Group data by scan
        grouped = df.groupby('ScanNumber')

        # Find the scan closest to peak_rt_tic
        rt_diff = [(scan, abs(group['RetentionTime'].iloc[0] - peak_rt_val)) for scan, group in grouped]
        closest_scan = min(rt_diff, key=lambda x: x[1])[0]
        closest_group = grouped.get_group(closest_scan)
        closest_rt = closest_group['RetentionTime'].iloc[0]

        # Save refined CSV with only the closest scan
        refined_csv_path = os.path.join(output_folder, f"{peak_file.replace('.csv', '')}_refined.csv")
        closest_group.to_csv(refined_csv_path, index=False)
        print(f"Saved refined CSV: {refined_csv_path} (Closest RT: {closest_rt:.2f}, Target: {peak_rt_val:.2f})")

        # Plot m/z vs Intensity for this scan
        plt.figure(figsize=(10, 4))
        plt.plot(closest_group['m/z'], closest_group['Intensity'], drawstyle='steps-mid')
        plt.title(f"{base_name} Peak {peak_idx+1} — RT: {closest_rt:.2f} min (Target: {peak_rt_val:.2f})")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.tight_layout()

        plot_path = os.path.join(output_folder, f"{peak_file.replace('.csv', '')}_refined.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"Saved plot: {plot_path}")

    except Exception as e:
        print(f"Error processing {peak_file}: {e}")
