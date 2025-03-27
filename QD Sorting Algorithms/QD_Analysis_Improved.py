import pandas as pd
import os
import numpy as np
import glob
import re
from scipy.signal import find_peaks

# Define constants
Cs_wvl_min = 891.6
Cs_wvl_max = 898.2

# Define data folder
folder_name = r"C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\11_03_2025_Test\QD_Text_File_Data_Corrected"
os.chdir(folder_name)

# List to store matching QDs with coordinates
matching_qds = []

# Process each .txt file
for filename in glob.glob('*.txt'):
    print(f"Processing: {filename}")

    # Extract [num num] from the filename
    match = re.search(r'\[(\d+)\s+(\d+)\]', filename)
    if match:
        qd_coords = f"[{match.group(1)} {match.group(2)}]"
    else:
        qd_coords = "UNKNOWN"

    # Read the full data file
    data_frame = pd.read_csv(filename, sep="\t")

    wvlngth_array = np.asarray(data_frame['Wavelen. (nm) '])
    counts_array = np.asarray(data_frame['Count Average'])

    # Find peaks
    spectrum_avg = np.mean(counts_array)
    spectrum_cutoff = 300  # Adjust as needed
    peaks, _ = find_peaks(counts_array, height=spectrum_avg + spectrum_cutoff, distance=75)

    # Check if any peaks fall in Cs transition range
    for peak in peaks:
        if Cs_wvl_min <= wvlngth_array[peak] <= Cs_wvl_max:
            print(f"QD Found in {filename} at {wvlngth_array[peak]:.3f} nm, Coordinates: {qd_coords}")
            matching_qds.append(qd_coords)
            break  # No need to check further in this file

# Output results to a text file
output_filename = "QD_Coordinates.txt"
with open(output_filename, "w") as output_file:
    for qd in matching_qds:
        output_file.write(qd + "\n")

print(f"\nQDs matching Cs D1 transition range saved to {output_filename}.")
print(f"Total matches: {len(matching_qds)}")
