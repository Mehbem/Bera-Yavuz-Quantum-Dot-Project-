% Auto Ploter
% Bera Yavuz 

    % Main folder of interest (user needs to update)
%main_folder = ""; % User inputs full path (LAB pc)
%main_folder = "/Users/bera_yavuz/Library/CloudStorage/OneDrive-UniversityofWaterloo/QD Data/Corrected_Wavelengths_Raster"; % User inputs full path (Macbook) 
main_folder = "C:\Users\yavub\OneDrive\Desktop\QD Analysis\QD_894.600-894.700_nm_Range_1000_Count_Up\text_files"; 

% Out folder with all the plots
%QD_plot_directory = "/Users/bera_yavuz/Library/CloudStorage/OneDrive-UniversityofWaterloo/QD Data/Corrected_Wavelengths_Plots";
QD_plot_directory = "C:\Users\yavub\OneDrive\Desktop\QD Analysis\QD_894.600-894.700_nm_Range_1000_Count_Up\Corrected_Wavelengths_Plots"; 

% Get all text files from main folder and subfolders
QD_txt_files = dir(fullfile(main_folder, '**', '*.txt'));
counter = 0; % keep track of how many textfiles are done so far

% Check if folder exists; if not, create it
if ~exist(QD_plot_directory, 'dir')
    mkdir(QD_plot_directory);  
end

for QD_ID = 1:numel(QD_txt_files)
    counter = counter + 1;
    % Get full file path
    file_path = fullfile(QD_txt_files(QD_ID).folder, QD_txt_files(QD_ID).name);
    [~, name_only, ~] = fileparts(file_path);

    % QD number extraction
    QD_str = regexp(QD_txt_files(QD_ID).name, '\[(\d+)\s+(\d+)\]', 'tokens');
    QD_Coord = [str2num(QD_str{1}{1}) str2num(QD_str{1}{2})]; 
    
    % Read the data (assumes tab or space-separated values)
    data = readmatrix(file_path);
    
    % Extract columns
    wvlength = data(:, 1); % First column
    counts = data(:, 2);       % Second column
    
    % filter counts and wvlength region
    beginnging_region_ignored = numel(wvlength)/4; 
    wvlength = wvlength(beginnging_region_ignored:end);
    counts = counts(beginnging_region_ignored:end);
    
    % Create a new figure
    figure;
    
    % Make the figure invisible
    set(gcf, 'Visible', 'off');

    % plot graph 
    plot(wvlength, counts,'b-');
    
     % Set exact limits for the x-axis
    xlim([min(wvlength), max(wvlength)]);
    
    % Set title and labels
    title_font = sprintf("QD Spectrum Plot: [%d %d]",QD_Coord);
    title(title_font);
    xlabel('Wavelength [nm]');
    ylabel('Arb. Counts');
    
    % finding the main peaks 
    [pks,locs,~,~] = findpeaks(counts,wvlength,'SortStr','descend','NPeaks',3,'MinPeakDistance',0.3);
    
    % Create text strings for top 3 peaks
    peak_text = cell(3,1);
    for n = 1:min(3,length(pks))
        peak_text{n} = sprintf('Peak %d: %.3f nm Abs. Counts: %.2f', n, locs(n),pks(n));
    end
    
    % Add text box in top-right corner
    text(0.95, 0.95, peak_text,...
        'Units', 'normalized',...
        'HorizontalAlignment', 'right',...
        'VerticalAlignment', 'top',...
        'BackgroundColor', [1 1 1 0.7],... % Semi-transparent white
        'EdgeColor', 'k',...
        'FontSize', 10,...
        'Margin', 3);
    
    % Set title and labels
    title_font = sprintf("QD Spectrum Plot: [%d %d]",QD_Coord);
    title(title_font);
    xlabel('Wavelength [nm]');
    ylabel('Arb. Counts');
    
    % plot image
    plot_img = gcf; 
    
    final_plot_full_path = fullfile(QD_plot_directory,name_only); 
    plot_filename = strcat(final_plot_full_path, '_wvl_graph.png'); 
    saveas(gcf, plot_filename);
    fprintf("Completed: %d/%d\n", counter,numel(QD_txt_files));
end