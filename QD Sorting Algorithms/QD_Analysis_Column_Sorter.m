% Define the text file name
filename = 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\QD_Coordinates.txt';

% Read the file content
fileContent = readlines(filename);

% Initialize an empty matrix to store the coordinates
coordinates = [];

% Loop through each line of the file content
for i = 1:length(fileContent)
    % Use regex to extract the numbers inside the brackets
    match = regexp(fileContent{i}, '\[(\d+)\s+(\d+)\]', 'tokens');
    
    % Check if a valid match is found
    if ~isempty(match)
        % Convert the match to numeric values and store them as a row in the matrix
        coords = str2double(match{1});
        coordinates = [coordinates; coords];
    end
end

% Sort by the second column first (ascending), then by the first column (ascending)
sorted_coordinates = sortrows(coordinates, [2, 1]);

% Find the unique values in the second column (groups)
unique_y_values = unique(sorted_coordinates(:,2));

% Initialize a new matrix to store the modified order
modified_coordinates = [];

% Loop through each unique y-value (column group)
for idx = 1:length(unique_y_values)
    % Extract rows with the same second column value
    group_rows = sorted_coordinates(sorted_coordinates(:,2) == unique_y_values(idx), :);
    
    % If the index is even, flip the order
    if mod(idx, 2) == 0
        group_rows = flipud(group_rows);
    end
    
    % Append to the modified coordinates matrix
    modified_coordinates = [modified_coordinates; group_rows];
end

% Display the modified coordinates
disp('Modified QD Coordinates:');
disp(modified_coordinates);

% Save the modified coordinates back to a file
output_filename = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\QD_Data\Modified_QD_Coordinates.txt";
writematrix(modified_coordinates, output_filename, 'Delimiter', '\t');
