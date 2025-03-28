% Wavelength Correction for Multiple Subfolders
% Bera Yavuz 

%Define data folder (containing all QD to be tested)
folder_with_QD_tested = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\Corrected_Wavelengths_Raster"; % (Home pc) 

% Main folder of interest (user needs to update)
main_folder = folder_with_QD_tested; % User inputs full path

% Get all text files from main folder and subfolders
QD_txt_files = dir(fullfile(main_folder, '**', '*.txt'));

% Process each file
for i = 1:numel(QD_txt_files)
    % Get full file path
    file_path = fullfile(QD_txt_files(i).folder, QD_txt_files(i).name);

    % Read the text file
    data = readtable(file_path);
    data = table2array(data); 

    % Extract wavelength and count values
    wavelengths = data(:, 1);
    counts = data(:, 2);

    % ROI for wavelength
    Max_search_wvlngth = 855; 
    idx = find(wavelengths >= Max_search_wvlngth, 1) - 1; % Find first index where condition fails
    wavelength_ROI = wavelengths(1:idx);

    % Determine actual reference wavelength (assuming peak intensity)
    [~, max_idx] = max(counts(1:idx));
    actual_reference = wavelengths(max_idx);

    % Compute the correction factor
    correction_factor = 852.357 - actual_reference;

    % Apply correction
    corrected_wavelengths = wavelengths + correction_factor;

    % Save the corrected data in the same folder as the original file
    corrected_data = [corrected_wavelengths, counts];
    output_file = fullfile(QD_txt_files(i).folder, "Corrected_" + QD_txt_files(i).name);
    writematrix(corrected_data, output_file, 'Delimiter', '\t');

    fprintf("Corrected file saved: %s\n", output_file);
end

disp("Wavelength correction applied to all text files.");
