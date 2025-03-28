% Overlap existing plots script
% Bera Yavuz
clear clc; 
% Define the base directory
all_QD_textfile = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\Corrected_Wavelengths_Raster"; %(Home PC)

% Get all PNG files in the directory
files = dir(fullfile(all_QD_textfile, "*.txt"));
QD_plot_filenames = string({files.name}); % Convert to string array for easier comparisons

% Read file with desired QD
Desired_QD_Filename = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\QD_Analysis_Results\QD_CS_Range_1000_Count_UP\QD_Coordinates_03_27_2025.txt";

% Open the file
fid = fopen(Desired_QD_Filename, 'r');

% Read numeric data, ignoring brackets
data = textscan(fid, '[%f %f]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fid);

% QD coordinates
QD_coord = [data{1}, data{2}];
list_length = size(QD_coord, 1);

% Create a figure for stacking plots
figure;
hold on; % Enables multiple plots on the same figure

% Define region of interest
Cs_wvl_min = 894.3;
Cs_wvl_max = 895.3;

% Loop through each QD to find and overlay matching spectra
for QD_ID = 1:list_length
    % Format the expected coordinate pattern
    textfile_desired = sprintf("[%d %d]", QD_coord(QD_ID, :));

    % Find filenames that CONTAIN the desired coordinate pattern
    match_idx = find(contains(QD_plot_filenames, textfile_desired), 1);

    if ~isempty(match_idx)
        % Get the actual filename that contains the desired coordinates
        matching_filename = QD_plot_filenames(match_idx);
    
        % Source path for reading data
        source_file = fullfile(files(match_idx).folder, matching_filename);

        % Read the data (assumes tab or space-separated values)
        data = readmatrix(source_file);

        % Extract columns
        wvlength = data(:, 1); % First column (Wavelength)
        counts = data(:, 2);   % Second column (Counts)

        % Filter counts and wavelength region
        beginning_region_ignored = floor(numel(wvlength) / 4); 
        wvlength = wvlength(beginning_region_ignored:end);
        counts = counts(beginning_region_ignored:end);

        % Plot spectrum on the same figure
        plot(wvlength, counts, '-', 'DisplayName', sprintf("[%d %d]", QD_coord(QD_ID, :)));

        % Print confirmation
        fprintf("Plotted: %s\n", matching_filename);
    end
end

% Add vertical red lines to highlight the region of interest
xline(Cs_wvl_min, 'r--', 'LineWidth', 2);
xline(Cs_wvl_max, 'r--', 'LineWidth', 2);

% Formatting the plot
xlim([min(wvlength), max(wvlength)]);
xlabel('Wavelength [nm]');
ylabel('Arb. Counts');
title('Overlapped QD Spectrum Plots');
legend('show'); % Show legend for each QD plot
grid on;

% Hold off to finalize stacking
hold off;
