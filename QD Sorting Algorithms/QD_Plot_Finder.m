% QD Plot Finder

clc; clear; close all;
warning('off');

% Define data folder (containing all QD to be tested)
%folder_with_QD_tested = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\QD_To_Analyze_A"; (Lab PC)
folder_with_QD_tested = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\Corrected_Wavelengths_Plots"; % (Home pc) 

% Define the base directory
%pathway_with_QD_acquired = "C:\Users\Quantum Dot\Desktop\Bera Yavuz -ANC300 Movement and Images\QD_Data\"; (Lab pc)
pathway_with_QD_acquired = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\QD_Analysis_Results\QD_CS_Range_1000_Count_Up\Acquired_QD_Plots"; %(Home PC)

% File names
files = dir(fullfile(folder_with_QD_tested, "*.png"));
QD_plot_filenames = {files.name}; 

% read file with desired QD
Desired_QD_Filename = "C:\Users\yavub\OneDrive - University of Waterloo\QD Data\QD_Analysis_Results\QD_CS_Range_1000_Count_UP\QD_Coordinates_03_27_2025.txt";
% Open the file
fid = fopen(Desired_QD_Filename, 'r');

% Read numeric data, ignoring brackets
data = textscan(fid, '[%f %f]');

% QD coordinates
QD_coord = [data{1},data{2}];

list_length = size(QD_coord,1); 
% looping through each qd to find matching points 
for QD_ID = 1:list_length
    % Formatting of file being searched for 
    Plot_desired = sprintf("[%d %d]",QD_coord(QD_ID,:)); 

    % Find filenames that CONTAIN the desired coordinate pattern
    match_idx = find(contains(QD_plot_filenames, Plot_desired), 1); 
    
    if ~isempty(match_idx)
        % Get the actual filename that contains the desired coordinates
        matching_filename = QD_plot_filenames{match_idx};
    
        % Source and destination paths
        source_file = fullfile(files(match_idx).folder, matching_filename);
        destination_file = fullfile(pathway_with_QD_acquired, matching_filename);
        
        % Copy the file to the new location
        copyfile(source_file, destination_file);
        
    end
    
    % Progression message
    fprintf("Completed Reading: %d/%d\n",QD_ID,list_length)
end


