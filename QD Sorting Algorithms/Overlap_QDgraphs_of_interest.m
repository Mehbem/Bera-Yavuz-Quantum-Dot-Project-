% Overlap known QD points 


% folders to extract points from
QD_folder_1 = "";
QD_folder_2 = "";
QD_folder_3 = ""; % add more if neededed

% matrix of all folders
QD_folders_Interest = [QD_folder_1,QD_folder_2,QD_folder_3];

for folder_ID = 1:numel(QD_folders_Interest) 

    % get all files in current folder of interest
    QD_txt_files = dir(fullfile(QD_folders_Interest{folder_ID}, '**', '*.txt'));
    % Extract filenames from the structure
    QD_txt_filenames = {QD_txt_files.name}; % Convert to a cell array of strings
    QD_txt_filenames = QD_txt_filenames(:); % convert to column array
    QD_txt_filenames = [QD_txt_filenames;QD_txt_filenames]; 
end


