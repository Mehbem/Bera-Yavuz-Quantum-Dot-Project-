% QD Analysis Script 

clc; clear; close all;
warning('off');


% Define constants
Cs_wvl_min = 891.6;
Cs_wvl_max = 898.2;

% Define data folder
folder_name = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\QD_To_Analyze_A";
folder_name_analysis = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\"; 
files = dir(fullfile(folder_name, "*.txt"));

% Initialize array for matching QDs
matching_qds = [];

% Process each .txt file
for i = 1:length(files)
    filename = fullfile(folder_name, files(i).name);
    fprintf("Processing: %s\n", files(i).name);
    
    % Extract [num num] from the filename
    qd_match = regexp(files(i).name, '\[(\d+)\s+(\d+)\]', 'tokens');
    if ~isempty(qd_match)
        qd_coords = sprintf("[%s %s]", qd_match{1}{1}, qd_match{1}{2});
    else
        qd_coords = "UNKNOWN";
    end
    
    % Read data
    data_table = readtable(filename);
    data_array = table2array(data_table); 

    % Convert to arrays
    wvlngth_array = data_array(:,1);
    counts_array = data_array(:,2);

    % Find peaks
    spectrum_avg = mean(counts_array);
    spectrum_cutoff = 300;  % Adjust as needed
    [pks, locs] = findpeaks(counts_array, 'MinPeakHeight', spectrum_avg + spectrum_cutoff, 'MinPeakDistance', 75);

    % Check if any peaks fall in Cs transition range
    if any(wvlngth_array(locs) >= Cs_wvl_min & wvlngth_array(locs) <= Cs_wvl_max)
        fprintf("QD Found in %s at %.3f nm, Coordinates: %s\n", files(i).name, wvlngth_array(locs(1)), qd_coords);
        matching_qds = [matching_qds; qd_coords];  % Append to list
    end
end

% Output results to a text file
output_filename = fullfile(folder_name_analysis, "QD_Coordinates.txt");
fid = fopen(output_filename, 'w');
for j = 1:length(matching_qds)
    fprintf(fid, "%s\n", matching_qds(j));
end
fclose(fid);

fprintf("\nQDs matching Cs D1 transition range saved to %s.\n", output_filename);
fprintf("Total matches: %d\n", length(matching_qds));
