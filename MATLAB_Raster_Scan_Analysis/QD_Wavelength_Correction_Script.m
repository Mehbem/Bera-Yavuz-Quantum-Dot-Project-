% Wavelength Correction for Multiple Subfolders
% Bera Yavuz
% Main folder of interest containing uncorrected wavelengths (user needs to update)
main_folder = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Mega Folder All QD (UnCorrected)"; % User inputs full path (Home pc)


% Get all text files from main folder and subfolders
QD_txt_files = dir(fullfile(main_folder, '**', '*.txt'));
counter = 0; % keep track of how many textfiles are done so far

% Process each file (change to desired folder) 
corrected_wavelength_output = "C:\Users\yavub\OneDrive\Desktop\Corrected_Wavelengths_Raster"; 

% Defined variables 
max_correction_factor = 0; 
file_of_interest = ''; 


for i = 1:numel(QD_txt_files)
    counter = counter + 1;
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
    correction_factor = abs(852.357 - actual_reference);
    if correction_factor > max_correction_factor && counts(max_idx) > 400 
        max_correction_factor = correction_factor; 
        file_of_interest = QD_txt_files(i).name; 
        reference_of_interest = actual_reference; 
    end
    % Apply correction
    corrected_wavelengths = wavelengths + correction_factor;

    % Define base filename
    base_name = "Corrected_" + QD_txt_files(i).name;
    %mega_folder = "/Users/bera_yavuz/Library/CloudStorage/OneDrive-UniversityofWaterloo/QD Data/Corrected_Wavelengths_Raster"; % macbook
    output_file = fullfile(corrected_wavelength_output, base_name);

    % Check if file exists and create a versioned filename if necessary
    version = 1;
    while isfile(output_file)
        % Generate new filename with version number
        versioned_name = "Corrected_v" + version + "_" + QD_txt_files(i).name;
        output_file = fullfile(corrected_wavelength_output, versioned_name);
        version = version + 1;
    end

    % Write the data
    corrected_data = [corrected_wavelengths(:), counts(:)];
    % Define column headers
    headers = "Wavelength\tCounts";

    % Write headers first
    writelines(headers, output_file);
    writematrix(corrected_data, output_file, 'Delimiter', '\t');
    
    fprintf("Completed: %d/%d\n", counter,numel(QD_txt_files));
end
disp("Wavelength correction applied to all text files.");
disp("The largest amount of shift %.3f nm", max_correction_factor)