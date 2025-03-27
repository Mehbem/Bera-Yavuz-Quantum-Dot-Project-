% QD Analysis Script 

clc; clear; close all;
warning('off');


% Define data folder
folder_name_list = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\"; 

% Output results to a text file
output_filename = fullfile(folder_name_list, "QD_Search.txt");
fid = fopen(output_filename, 'w');

% Custom Matrix
rows = 1:4:200;
rows = flipud(rows'); 
columns = ones(length(rows),1);
matching_qds = horzcat(rows,columns);

for j = 1:length(matching_qds)
    fprintf(fid, "%d\t%d\n", matching_qds(j,1),matching_qds(j,2));
end
fclose(fid);


fprintf("finished uploading list"); 

file_name = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\QD_Search.txt";
QD_matrix = readmatrix(file_name); 