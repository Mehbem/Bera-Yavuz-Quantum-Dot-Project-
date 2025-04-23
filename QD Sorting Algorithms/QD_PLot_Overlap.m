% Overlap Plots of Matching QD Spectra from Two Folders
% Author: Bera Yavuz (Modified)
% April 2025

%--- INPUT PATHS ---
large_folder = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Mega Folder All QD (UnCorrected)\150 Gratting Columns 1-14";  % Folder with all QD text files
subset_folder = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Mega Folder All QD (UnCorrected)\150 Gratting Column 1-14 List";  % Folder with a subset of the QD text files
output_folder = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Mega Folder All QD (UnCorrected)\Overlapped Plots";  % Folder to save overlapped plots

% --- CREATE OUTPUT FOLDER IF NOT EXISTS ---
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% --- GET FILE LISTS ---
large_files = dir(fullfile(large_folder, '**', '*.txt'));
subset_files = dir(fullfile(subset_folder, '**', '*.txt'));

% --- PROCESS EACH FILE IN SUBSET ---
counter = 0; 
for i = 1:numel(subset_files)
    counter = counter + 1; 
    % Current file in subset
    subset_name = subset_files(i).name;
    subset_path = fullfile(subset_files(i).folder, subset_name);
    
    % Find match in large dataset
    % Extract prefix (everything before '_Emission')
prefix = erase(subset_name, '_Emission.txt');

% Try to find a matching file in large dataset
match_idx = find(contains({large_files.name}, prefix), 1);
    if isempty(match_idx)
        fprintf("No match found for: %s\n", subset_name);
        continue;
    end
    large_path = fullfile(large_files(match_idx).folder, prefix+".txt");

    % Wavelength Filtering 
    beginnging_region_ignored = 6248/4; 

    % Read data
    subset_data = readtable(subset_path);
    subset_data = table2array(subset_data);

    large_data = readtable(large_path); 
    large_data = table2array(large_data);

    % Extract wavelength and counts
    wavelength_subset = subset_data(:,1);
    wavelength_subset = wavelength_subset(beginnging_region_ignored:end); 
    counts_subset = subset_data(:,2);
    counts_subset = counts_subset(beginnging_region_ignored:end); 

    wavelength_large = large_data(:,1);
    wavelength_large = wavelength_large(beginnging_region_ignored:end);
    counts_large = large_data(:,2);
    counts_large = counts_large(beginnging_region_ignored:end); 

    %% --- PLOT OVERLAPPED SPECTRA ---
    figure('Visible','off'); % Don't show GUI
    plot(wavelength_large, counts_large, 'b-', 'LineWidth', 1.5); hold on;
    plot(wavelength_subset, counts_subset, 'r', 'LineWidth', 1);
    xlim([min([min(wavelength_subset), min(wavelength_large)]), max([max(wavelength_subset), max(wavelength_large)])]);
    legend('Large Set', 'Subset');
    title(['Overlapped Plot: ', subset_name], 'Interpreter', 'none');
    xlabel('Wavelength (nm)');
    ylabel('Counts');
    grid on;

    % Save figure
    [~, name, ~] = fileparts(subset_name);
    save_path = fullfile(output_folder, [name, '_overlap.png']);
    saveas(gcf, save_path);

    % Optional: close figure to save memory
    close(gcf);

    fprintf("Saved overlapped plot for: %d/%d\n", counter,numel(subset_files));
end

disp("All matching plots have been saved.");
%C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Mega Folder All QD (UnCorrected)\150 Gratting Column 1-14 List\QD_[1 4]\Set_1\[1 4]_150mm_gratting_Emission.txt
%C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Mega Folder All QD (UnCorrected)\150 Gratting Columns 1-14\[1 4]_150mm_gratting_Emission.txt