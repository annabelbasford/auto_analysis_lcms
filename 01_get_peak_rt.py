import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import glob
from scipy.signal import find_peaks

# path to the data
data_path = '/auto_analysis_lcms/DAD_uv/'
file_names = [f for f in os.listdir(data_path) if f.endswith('_uv.csv')]

for file in file_names:
    # read the data
    df = pd.read_csv(data_path + file)
    # make data positive
    df['intensity'] = df['intensity'] - df['intensity'].min()
    # drop data points before 2 minutes
    df = df[df['rt'] > 2]

    # find peaks
    peaks, _ = find_peaks(df['intensity'], height=2)

    # collect rt values for peaks
    peak_rt = df['rt'].iloc[peaks].tolist()
    peak_shift = 0.11
    # shift peak to get peak_rt_tic
    peak_rt_tic = [rt + peak_shift for rt in peak_rt]

    # find the start and end of each peak
    peak_start_uv = []
    peak_end_uv = []
    peak_start_tic = []
    peak_end_tic = []
    
    for i in range(len(peaks)):
        start = peaks[i]
        while start > 0 and df['intensity'].iloc[start] > df['intensity'].iloc[start - 1]:
            start -= 1
        peak_start_uv.append(df['rt'].iloc[start])


        end = peaks[i]
        while end < len(df) - 1 and df['intensity'].iloc[end] > df['intensity'].iloc[end + 1]:
            end += 1
        peak_end_uv.append(df['rt'].iloc[end])

    # calculate the corrected peak areas
    peak_area = []
    peak_segments = []  # save peak data segments for plotting
    for i in range(len(peaks)):
        # segment between start and end
        mask = (df['rt'] >= peak_start_uv[i]) & (df['rt'] <= peak_end_uv[i])
        rt_segment = df['rt'][mask]
        intensity_segment = df['intensity'][mask]

        # raw area under the curve
        raw_area = np.trapz(intensity_segment, rt_segment)

        # straight line baseline between start and end
        baseline = np.linspace(intensity_segment.iloc[0], intensity_segment.iloc[-1], len(intensity_segment))
        baseline_area = np.trapz(baseline, rt_segment)

        # corrected area
        corrected_area = raw_area - baseline_area
        peak_area.append(corrected_area)

        # save corrected intensity (subtract baseline)
        corrected_intensity = intensity_segment - baseline
        peak_segments.append((rt_segment, corrected_intensity))

    # Filter out peaks with area < 1% of the largest peak
    max_area = max(peak_area)
    filtered_peak_rt = []
    filtered_peak_start = []
    filtered_peak_end = []
    filtered_peak_area = []
    filtered_peak_segments = []
    filtered_peak_rt_tic = []

    for i in range(len(peak_area)):
        if peak_area[i] >= 0.01 * max_area:
            filtered_peak_rt.append(peak_rt[i])
            filtered_peak_start.append(peak_start_uv[i])
            filtered_peak_end.append(peak_end_uv[i])
            filtered_peak_area.append(peak_area[i])
            filtered_peak_segments.append(peak_segments[i])
            filtered_peak_rt_tic.append(peak_rt_tic[i])

    # Use the filtered lists going forward
    peak_rt_uv     = filtered_peak_rt
    peak_start_uv  = filtered_peak_start
    peak_end_uv    = filtered_peak_end
    peak_area      = filtered_peak_area
    peak_segments  = filtered_peak_segments
    peak_rt_tic    = filtered_peak_rt_tic

    # now shift the peak start and end to get peak_start_tic and peak_end_tic
    for i in range(len(peak_start_uv)):
        start_tic = peak_start_uv[i] + peak_shift
        end_tic = peak_end_uv[i] + peak_shift
        start_tic_rounded = float(f"{start_tic:.4g}")
        end_tic_rounded = float(f"{end_tic:.4g}")
        peak_start_tic.append(start_tic_rounded)
        peak_end_tic.append(end_tic_rounded)

    # Round all values to 4 significant figures
    peak_rt_uv     = [float(f"{x:.4g}") for x in peak_rt_uv]
    peak_start_uv  = [float(f"{x:.4g}") for x in peak_start_uv]
    peak_end_uv    = [float(f"{x:.4g}") for x in peak_end_uv]
    peak_rt_tic    = [float(f"{x:.4g}") for x in peak_rt_tic]
    peak_area      = [float(f"{x:.4g}") for x in peak_area]


    # plot only the corrected peak areas
    plt.figure(figsize=(12, 7))

    colors = ['yellow', 'orange', 'purple', 'pink', 'brown', 'gray']

    for i, (rt_seg, corr_int) in enumerate(peak_segments):
        plt.fill_between(rt_seg, corr_int, color=colors[i % len(colors)], alpha=0.5, label=f'Peak {i+1}')
        
        # Annotate: Peak RT and Area
        peak_rt_val = peak_rt_uv[i]
        peak_area_val = peak_area[i]
        
        # Find max corrected intensity for positioning text
        max_idx = np.argmax(corr_int)
        max_rt = rt_seg.iloc[max_idx]
        max_intensity = corr_int.iloc[max_idx]
        
        # Annotate above the peak
        plt.text(
            max_rt,
            max_intensity + 0.05 * max(corr_int),  # small offset above
            f"RT={peak_rt_val:.2f} min\nArea={peak_area_val:.1f}",
            ha='center',
            va='bottom',
            fontsize=8,
            color='black',
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='gray', facecolor='white', alpha=0.6)
        )

    plt.xlabel('RT (min)')
    plt.ylabel('Corrected Intensity')
    plt.tight_layout()

    # save plot
    png_path = '/home/abasford/projects/BO/auto_analysis_lcms/pngs'
    # save to the png_path
    os.makedirs(png_path, exist_ok=True)
    plt.savefig(os.path.join(png_path, f'{file}.png'), dpi=300)

    # save peak info
    peak_df = pd.DataFrame({'peak_rt_uv': peak_rt_uv, 'peak_start_uv': peak_start_uv, 'peak_end_uv': peak_end_uv, 'peak_area': peak_area, 'peak_rt_tic': peak_rt_tic, 'peak_start_tic': peak_start_tic, 'peak_end_tic': peak_end_tic})
    # split file name to get the base name
    base_name = os.path.splitext(file)[0]
    # save the peak info to a csv file with base name
    peak_df.to_csv(f'{base_name}_peak_rt.csv', index=False)
