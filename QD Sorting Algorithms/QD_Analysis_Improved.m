% QD Analysis Script 

clc; clear; close all;
warning('off');


% Define constraints
Cs_wvl_min = 894.3;
Cs_wvl_max = 895.3;

% Define data folder (containing all QD to be tested)
%folder_with_QD_tested = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\QD_To_Analyze_A"; (Lab PC)
folder_with_QD_tested = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\Corrected_Wavelengths_Raster"; % (Home pc) 

% Define the base directory
%pathway_with_QD_acquired = "C:\Users\Quantum Dot\Desktop\Bera Yavuz -ANC300 Movement and Images\QD_Data\"; (Lab pc)
pathway_with_QD_acquired = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\QD_Analysis_Results\QD_CS_Range_1000_Count_Up"; %(Home PC)


% Get today's date in YYYY-MM-DD format
today_date = datestr(now, 'mm_dd_yyyy');

% Create the filename with the date
textfile_name = sprintf("QD_Coordinates_%s.txt", today_date);

% Full file path
full_file_path = fullfile(pathway_with_QD_acquired, textfile_name);

files = dir(fullfile(folder_with_QD_tested, "*.txt"));

% Initialize array for matching QDs
matching_qds = [];

% Process each .txt file
for i = 1:length(files)
    filename = fullfile(folder_with_QD_tested, files(i).name);
    %fprintf("Processing: %s\n", files(i).name);
    
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
    %spectrum_avg = 0;
    spectrum_cutoff = 1000;  % Adjust as needed
    [pks, locs] = findpeaks(counts_array, 'MinPeakHeight', spectrum_avg + spectrum_cutoff, 'MinPeakDistance', 75);

    % Check if any peaks fall in Cs transition range
    if any(wvlngth_array(locs) >= Cs_wvl_min & wvlngth_array(locs) <= Cs_wvl_max)
        fprintf("QD Found in %s at %.3f nm, Coordinates: %s\n", files(i).name, wvlngth_array(locs(1)), qd_coords);
        matching_qds = [matching_qds; qd_coords];  % Append to list
    end

    % Progression message
    fprintf("Completed Reading: %d/%d\n",i,length(files))
end

% Output results to a text file
fid = fopen(full_file_path, 'w');
for j = 1:length(matching_qds)
    fprintf(fid, "%s\n", matching_qds(j));
end
fclose(fid);

fprintf("\nQDs matching Cs D1 transition range saved to %s.\n", full_file_path);
disp("*********************************************************************************")
fprintf("Total matches: %d\n", length(matching_qds));
disp("*********************************************************************************")

