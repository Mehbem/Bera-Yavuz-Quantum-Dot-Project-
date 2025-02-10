% Fill the grid with intensity values based on coordinates
data_table_X = data_table.Positions(:,1);
data_table_Y = data_table.Positions(:,2);

xVals = unique(data_table_X);
yVals = unique(data_table_Y);

data_table_photons = data_table.PeakCounts;

% Determine grid size
numCols = length(xVals);
numRows = length(yVals);

% Create an empty grid for intensity values
Z = nan(numRows, numCols);

for i = 1:size(data_table.PeakCounts,1)
    xIdx = find(xVals == data_table_X(i));
    yIdx = find(yVals == data_table_Y(i)); 
    Z(yIdx, xIdx) = data_table_photons(i);
end

% Flip the Y-axis so (0,0) is at the top-left
Z = flipud(Z);

% Plot the data using imagesc
figure;
imagesc(xVals, yVals, Z);
colormap(jet); % Change colormap if needed
colorbar;
axis equal tight;
ax = gca;
ax.YTickLabel = flip(ax.YTickLabel); % Reverse Y-axis labels
ax.XAxisLocation = 'top'; % Moves X-axis labels to the top  
set(gca, 'YDir', 'Normal'); % Ensure top-left origin
xlabel('X');
ylabel('Y');

title('Raster Scan Photon Count');
