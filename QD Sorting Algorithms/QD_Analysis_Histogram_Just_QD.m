% Auto Peak Histogram Generator
% Bera Yavuz

% Main folder containing text files
main_folder = "C:\Users\yavub\OneDrive\Desktop\QD Analysis\Corrected_Wavelengths_Raster";

% Output folder for the histogram
output_folder = "C:\Users\yavub\OneDrive\Desktop\QD Analysis\Peak_Histograms";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get all text files
QD_txt_files = dir(fullfile(main_folder, '**', '*.txt'));

% Initialize array to store all peaks
all_peaks = [];

for QD_ID = 1:numel(QD_txt_files)
    % Read the data
    file_path = fullfile(QD_txt_files(QD_ID).folder, QD_txt_files(QD_ID).name);
    data = readmatrix(file_path);
    
    % Extract columns
    wvlength = data(:, 1);
    counts = data(:, 2);
    
    % Filter data (same as your original)
    beginnging_region_ignored = numel(wvlength)/4; 
    wvlength = wvlength(beginnging_region_ignored:end);
    counts = counts(beginnging_region_ignored:end);
    
    % Find top 3 peaks
    [pks, locs] = findpeaks(counts, wvlength, 'SortStr', 'descend', 'NPeaks', 5, 'MinPeakDistance', 0.3);
    
    % Store the peaks
    all_peaks = [all_peaks; locs(:)];
    
    fprintf("Processed file %d/%d\n", QD_ID, numel(QD_txt_files));
end

% Create histogram of all peaks
figure;
set(gcf, 'Visible', 'off');  % Uncomment to suppress figure display

% Calculate optimal bin width using Freedman-Diaconis rule
if ~isempty(all_peaks)
    h = histogram(all_peaks, 'BinMethod', 'fd', 'FaceColor', 'b', 'EdgeColor', 'k');
    
    % Add labels and title
    xlabel('Wavelength (nm)');
    ylabel('Number of Peaks');
    title('Distribution of Peak Wavelengths Across All QDs');
    grid on;
    
    % Save the histogram
    output_path = fullfile(output_folder, 'QD_peak_wavelength_distribution.png');
    saveas(gcf, output_path);
    fprintf("Histogram saved to: %s\n", output_path);
else
    warning('No peaks found in any files!');
end