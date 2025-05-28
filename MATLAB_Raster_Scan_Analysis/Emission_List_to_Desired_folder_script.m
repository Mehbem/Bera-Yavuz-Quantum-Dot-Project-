% This function searches through a main directory containing quantum dot (QD)
% subfolders organized by columns and grating configurations. Within each QD
% folder (e.g., "QD_[9 11]"), it finds the most recent "Set_*" subfolder and
% copies any ".txt" emission files found there to a specified destination folder.

% The copied files are renamed to include the QD folder name as a prefix to avoid
% filename conflicts.


% Define source and destination directories
baseDir = 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Data_Summary_Reports\NWQD_Raster_Scan_Data\150_ln_mm_Grating_Raster_Scan_Winter_2025\150_ln-mm_Grating_Columns_1-14_List';
destDir = 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Data_Summary_Reports\NWQD_Raster_Scan_Data\150_ln_mm_Grating_Raster_Scan_Winter_2025\150_ln-mm_Grating_Columns_14-50'; % <-- customize this

% Make sure destination exists
if ~exist(destDir, 'dir')
    mkdir(destDir);
end


% Get list of all QD folders
qdFolders = dir(fullfile(baseDir, 'QD_*'));

for i = 1:length(qdFolders)
    qdPath = fullfile(baseDir, qdFolders(i).name);
    
    % Look for Set_* subfolders
    setFolders = dir(fullfile(qdPath, 'Set_*'));
    if isempty(setFolders)
        continue; % Skip if no sets
    end

    % Sort by date to find the most recent one
    [~, idx] = max([setFolders.datenum]);
    latestSet = setFolders(idx);
    latestSetPath = fullfile(qdPath, latestSet.name);

    % Find .txt files in the latest set
    txtFiles = dir(fullfile(latestSetPath, '*.txt'));

    for j = 1:length(txtFiles)
        srcFile = fullfile(latestSetPath, txtFiles(j).name);
        destFile = fullfile(destDir, [qdFolders(i).name '_' txtFiles(j).name]); % add QD info to filename
        copyfile(srcFile, destFile);
    end
end

disp('Text files copied successfully.');
