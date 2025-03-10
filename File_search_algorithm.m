% Bera Yavuz 
clear;
clc;


% Define the main (mega) folder to search within
megaFolder = '/Users/bera_yavuz/Desktop/Testing dumb thing folder'; % Change this to your actual folder path

% Define the image file extensions you want to find
imageExtensions = {'*.png', '*.fig','*.txt'}; 

% Initialize a cell array to store found image paths
imageFiles = {};

% Defining desired QD 
desired_QD = [1 1];
search_parameter.square_brackets = sprintf("[%d %d]",desired_QD);
% Defining the folder where everything is put into 
main_folderPath = sprintf('//Users//bera_yavuz//Desktop//All_QD_Organized//Quantum Dot Info_[%d %d]',desired_QD);
folderPath_Dot_position_Imgs = sprintf('%s//Dot_Positions',main_folderPath);
folderPath_ASI_raw_Imgs = sprintf("%s//ASI_IMGs",main_folderPath);
folderPath_Plots_QD = sprintf("%s//Plots_QD",main_folderPath);
folderPath_Rotation_data = sprintf("%s//Rotation_data",main_folderPath);
all_folder_paths = [folderPath_ASI_raw_Imgs,folderPath_Plots_QD,folderPath_Rotation_data];

if ~exist(main_folderPath, 'dir') % Check if folder already exists
    mkdir(main_folderPath);
    mkdir(folderPath_Rotation_data)
    mkdir(folderPath_ASI_raw_Imgs)
    mkdir(folderPath_Plots_QD)
    disp('Associated folders created successfully.');
else
    disp('Folders already exists.');
end

% All folders for different purposes 

% Loop through each image extension and search recursively
for i = 1:length(imageExtensions)
    % Get all image files of the current type recursively
    fileList = dir(fullfile(megaFolder, '**', imageExtensions{i}));
    fileNames = {fileList.name}; 

    % searching for only desired images 
    fileNames_square_brackets = contains(fileNames,search_parameter.square_brackets);
   
    % matched file names for the square bracket format [# #] 
    Matched_filesNames_square_brackets = fileList(fileNames_square_brackets);
    
    

% Copy the desired images to the new folder (square_brackets)
    for j = 1:length(Matched_filesNames_square_brackets)
        srcPath = fullfile(Matched_filesNames_square_brackets(j).folder, Matched_filesNames_square_brackets(j).name);
        destPath = fullfile(all_folder_paths(i), Matched_filesNames_square_brackets(j).name);

        % Handle duplicate file names by adding "_copyN" if needed
        copyCounter = 1;
        [fileDir, fileBase, fileExt] = fileparts(destPath);
        while exist(destPath, 'file')
            newFileName = sprintf('%s_copy%d%s', fileBase, copyCounter, fileExt);
            destPath = fullfile(fileDir, newFileName);
            copyCounter = copyCounter + 1;
        end
   
        % Copy file
        copyfile(srcPath, destPath);
        imageFiles{end+1} = destPath;
    end

end

% Display results
if isempty(imageFiles)
    disp('No image files found.');
else
    disp('Found image files:');
    disp(imageFiles');
end
