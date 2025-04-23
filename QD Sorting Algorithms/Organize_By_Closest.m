% QD Analysis Script 

clc; clear; close all;
warning('off');


% Define Wavelengths of Interest 
Wavelengths_of_Interest = [894.578 894.605]; 

% Variables
spectrum_cutoff = 400; 
Wavelengths_diff_cutoff = 1; 

% Define data folder (containing all QD to be tested)
folder_with_QD_tested = "/Users/bera_yavuz/Desktop/QD Analysis/Corrected_Wavelengths_Raster"; % (Macbook pathway) 

% Define the base directory
main_folder = "/Users/bera_yavuz/Desktop/QD Analysis/QD Closest";
for wvn = 1:numel(Wavelengths_of_Interest) 
    specific_folder = sprintf("QD Ordered Closest to %.3f New Algorithm",Wavelengths_of_Interest(wvn)); 
    pathway_with_QD_acquired = fullfile(main_folder,specific_folder);
    specific_folder_text_file = fullfile(pathway_with_QD_acquired,"text_files"); 

    % Check if folder exists; if not, create it
    if ~exist(specific_folder_text_file, 'dir')
        mkdir(specific_folder_text_file);  
    end

end

% Extracting all QD files to be analyzed 
files = dir(fullfile(folder_with_QD_tested, "*.txt"));

% Initialize array for matching QDs
matching_qds = [];

% Process each .txt file
for i = 1:length(files)
    filename = fullfile(folder_with_QD_tested, files(i).name);
    [~, name, ext] = fileparts(filename);
    specific_file_name = name; 
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

    % Ignore beginning region
    beginnging_region_ignored = floor(numel(wvlngth_array)/4);
    wvlngth_array = wvlngth_array(beginnging_region_ignored:end);
    counts_array = counts_array(beginnging_region_ignored:end);

    % Wavelength range and total number of data points
    wavelength_min = wvlngth_array(1);  % First wavelength value
    wavelength_max = wvlngth_array(end);  % Last wavelength value
    num_points = numel(wvlngth_array);  % Total number of data points
    
    % Calculate the number of data points per nm
    wavelength_range = wavelength_max - wavelength_min;  % Total range in nm
    data_points_per_nm = num_points / wavelength_range;  % Data points per nm

    % Convert 0.3 nm to the number of data points
    MinPeakDistance = round(data_points_per_nm * 0.3);  % Minimum peak distance in data points

    % Find peaks
    spectrum_avg = mean(counts_array);
    % Get the max intensity of the spectrum
    max_intensity = max(counts_array);
    
    % Set dynamic MinPeakProminence (e.g., 5% of peak height)
    min_prom = 0.05 * max_intensity;

    % Skip QD if MinPeakProminence is negative or NaN
    if isnan(min_prom) || min_prom < 0
        fprintf("Skipping %s due to invalid MinPeakProminence.\n", specific_file_name);
        continue;  % Skip to the next QD
    end
    
    % Find peaks with dynamic prominence and narrowness
    % Return peaks with locs in nm
    [pks, locs_nm, widths_nm, proms] = findpeaks(counts_array, wvlngth_array, ...
        'SortStr', 'descend', ...
        'NPeaks', 10, ...
        'MinPeakDistance', MinPeakDistance, ...
        'WidthReference', 'halfheight', ...
        'MinPeakProminence', min_prom,...
        'MinPeakHeight',spectrum_avg+spectrum_cutoff);
    
    % Convert locs (in nm) back to indices in wvlngth_array
    [~, locs_idx] = arrayfun(@(x) min(abs(wvlngth_array - x)), locs_nm);
    
    % Filter sharp peaks
    max_width_nm = 0.6;  % adjust for sharper or looser peaks
    valid_idx = find(widths_nm <= max_width_nm & proms >= min_prom);
    
    % Select final clean peaks
    clean_pks = pks(valid_idx);
    clean_locs = locs_nm(valid_idx);
    clean_locs_idx = locs_idx(valid_idx);  % corresponding indices in wvlngth_array

    [~, sort_idx] = sort(clean_pks, 'descend');
    n_show = min(3, length(sort_idx));
    final_pks = clean_pks(sort_idx(1:n_show));
    final_locs_nm = clean_locs(sort_idx(1:n_show));    % wavelengths in nm
    final_locs_idx = clean_locs_idx(sort_idx(1:n_show));  % corresponding array indices

    for wvn_ID = 1:numel(Wavelengths_of_Interest)
        target_wavelength = Wavelengths_of_Interest(wvn_ID);

        % Find the difference between the target wavelength and detected peaks
        wvnlngth_diff = abs(final_locs_nm - target_wavelength);  % if final_locs are still in nm


        if ~isempty(wvnlngth_diff)
            % Find the closest peak
            [min_diff, idx_min] = min(wvnlngth_diff);
            matched_peak_wavelength = final_locs_nm(idx_min);       % in nm
            matched_peak_index = final_locs_idx(idx_min);           % array index
            matched_peak_count = counts_array(matched_peak_index);  % intensity from counts array


            if min_diff  < Wavelengths_diff_cutoff
            % Store info about this match
            matching_qds(end+1).filename = specific_file_name; 
            matching_qds(end).qd_coords = qd_coords;
            matching_qds(end).peak_wavelength = matched_peak_wavelength;
            matching_qds(end).diff = min_diff;
            matching_qds(end).wvn_ID = wvn_ID;  % to track which wavelength it matched
            matching_qds(end).count = matched_peak_count;
            end
        end
    end

    fprintf("Completed Reading: %d/%d\n", i, length(files));
end

% ---------------------------------------------
% Split results by wavelength and sort by proximity
% ---------------------------------------------

for wvn_ID = 1:numel(Wavelengths_of_Interest)
    % Filter matches for this wavelength
    matches = matching_qds([matching_qds.wvn_ID] == wvn_ID);

    % Sort matches by smallest difference
    [~, sort_idx] = sort([matches.diff]);
    sorted_matches = matches(sort_idx);

    % Generate output file path
    specific_folder = sprintf("QD Ordered Closest to %.3f New Algorithm", Wavelengths_of_Interest(wvn_ID));
    specific_textfile = sprintf("Best_QD_%.3f.txt",Wavelengths_of_Interest(wvn_ID));
    output_txt_path = fullfile(main_folder, specific_folder, specific_textfile);

    % Just QD list
    text_file = sprintf("QD_Coords_%.3f.txt",Wavelengths_of_Interest(wvn_ID));
    QD_list_text = fullfile(main_folder,specific_folder,text_file);

    % Write to text file
    fid = fopen(output_txt_path, 'w');
    for j = 1:length(sorted_matches)
        fprintf(fid, "[%s] | Peak: %.3f nm | Count: %.2f | Δ: %.3f nm | Completion Status:\n", sorted_matches(j).qd_coords, sorted_matches(j).peak_wavelength,sorted_matches(j).count, sorted_matches(j).diff);
    end
    fclose(fid);
    
    % write to QD coord text file 
    fid = fopen(QD_list_text,"w");
    for j = 1:length(sorted_matches)
        fprintf(fid,"%s\n",sorted_matches(j).qd_coords);
    end
    fclose(fid);

    % Move the matching .txt files into the 'text_files.txt' subfolder
    text_file_folder = fullfile(main_folder, specific_folder, "text_files");

    if ~exist(text_file_folder, 'dir')
        mkdir(text_file_folder);
    end

    for j = 1:length(sorted_matches)
        src_file = fullfile(folder_with_QD_tested, sorted_matches(j).filename+ext);
        dest_file = fullfile(text_file_folder, sorted_matches(j).filename+ext);

        % Move file
        if exist(src_file, 'file')
            copyfile(src_file, dest_file);
        end
    end
end

fprintf("\n✨ All matching QDs organized and saved.\n");





