app.StartButton.BackgroundColor = 'g'; 
pause(0.2)
app.StartButton.BackgroundColor = app.DefaultColour; 

if isempty(app.Confirmed_column) | isempty(app.Confirmed_row)
    AppLogger_Specific(app,"Inquired Row/Column Missing",'red')
    return
end

% Parameters
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% X and y factors
Read_XY_factor = app.MyFuncs.XY_Factor_Identifier("","Read","LAB");
x_factor = Read_XY_factor.X_factor; %(Steps/pixels units)
y_factor = Read_XY_factor.Y_factor; %(Steps/pixels units)
Frequency = 20; 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




% Initial Movement Towards Start
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PhotoType = "SteppingPhoto";
readQD_position = app.MyFuncs.QD_tracking_N_Identification("","","Read","LAB","");
QD_counter = readQD_position.lastLine; 
CurrentQD =  QD_counter;
InquiredQD = [app.Confirmed_row app.Confirmed_column]; % user puts what value they want 
[StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking(app.ANC300,PhotoType,QD_counter,app.pyueye_initialization_return);
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% calculating number of rows and column alongside direction needed to travel to inquired QD 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
QD_loca_diff = CurrentQD - InquiredQD; 

if QD_loca_diff(1) > 0 & QD_loca_diff(2) < 0

% direction for each respectiive axis 
direction_row = "topright";
direction_column = "bottomright"; 

% Intervals for each respective axis 
row_interval = [-1 0];
column_interval = [0 1];

elseif QD_loca_diff(1)  > 0 & QD_loca_diff(2) > 0 

% direction for each respectiive axis 
direction_row = "topright";
direction_column = "topleft";

% Intervals for each respective axis 
row_interval = [-1 0];
column_interval = [0 -1];

elseif QD_loca_diff(1)  < 0 & QD_loca_diff(2) < 0 

% direction for each respectiive axis 
direction_row = "bottomleft";
direction_column = "bottomright"; 

% Intervals for each respective axis 
row_interval = [1 0];
column_interval = [0 1];

elseif QD_loca_diff(1)  < 0 & QD_loca_diff(2) > 0 

% direction for each respectiive axis     
direction_row = "bottomleft";
direction_column = "topleft";

% Intervals for each respective axis 
row_interval = [1 0];
column_interval = [0 -1];

elseif QD_loca_diff(1) == 0 | QD_loca_diff(2) == 0
    if QD_loca_diff(1) < 0 
    direction_row = "bottomleft";
    row_interval = [1 0];
    elseif QD_loca_diff(1) > 0 
    direction_row = "topright";
    row_interval = [-1 0];
    elseif QD_loca_diff(2) < 0 
    direction_column = "bottomright"; 
    column_interval = [0 1];
    elseif QD_loca_diff(2) > 0 
    direction_column = "topleft"; 
    column_interval = [0 -1];
    else
        AppLoggerText = sprintf("No movement"); 
        AppLogger_Specific(app,AppLoggerText,'black')
        return
    end

end

% number of rows and columns to move 
num_rows = abs(QD_loca_diff(1));
num_columns = abs(QD_loca_diff(2));
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


%  Movement Towards QD  
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PhotoType = "SteppingPhoto"; 
dots_to_travel_row = QD_loca_diff(1); 
dots_to_travel_column = QD_loca_diff(2);
groups_of_dots_row = 4; 
groups_of_dots_column = 4; 

num_of_breaks_row = floor(dots_to_travel_row / groups_of_dots_row);  % How many times it travels then locks in-between each QD (along the row) 
remaining_QD_row = mod(dots_to_travel_row, groups_of_dots_row);    % What remains that is indiviually travelled (along the row) 

num_of_breaks_column = floor(dots_to_travel_column / groups_of_dots_column);  % How many times it travels then locks in-between each QD (along the column) 
remaining_QD_column = mod(dots_to_travel_column, groups_of_dots_column);    % What remains that is indiviually travelled (along the column) 

Frequency_fast_travel = 100; %whatever frequency we choose for fast travel essentially 

% fast tracking along the QDs 
for column = 1: num_of_breaks_column

    % Moving to the next set of QD along column
    QD_counter = QD_counter + groups_of_dots_column*column_interval; 
    AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
    AppLogger_Specific(app,AppLogger_Specific_text,'black');
    app.MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_column,app.ANC300,Frequency_fast_travel,x_factor,y_factor)

    % Making sure LED is exactly on the dot 
    app.MyFuncs.QD_tracking_N_Identification(QD_counter,"Auto","Write","LAB",direction_column); 
    app.MyFuncs.Precision_Locking(app.ANC300,PhotoType,QD_counter,app.pyueye_initialization_return); % trying to land on the exact dot 
    AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
    AppLogger_Specific(app,AppLogger_Specific_text,app.Darkgreen); 
    
    % Updating all visuals on GUI
    QD_Location_text = sprintf("[%d %d]",QD_counter); 
    app.QD_location_SpecificQD.Text = QD_Location_text; 

    Current_QD_Position_text = sprintf("Nanowire_Photo_[%d %d].jpg",QD_counter); 
    app.Current_QD_Position_Specific.ImageSource = Current_QD_Position_text; 

end

% moving the remaing QD 
QD_counter = QD_counter + remaining_QD_column*column_interval; 
AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
AppLogger_Specific(app,AppLogger_Specific_text,'black');
app.MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_column,app.ANC300,Frequency_fast_travel,x_factor,y_factor)

% Making sure LED is exactly on the dot 
app.MyFuncs.QD_tracking_N_Identification(QD_counter,"Auto","Write","LAB",direction_column); 
[StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking(app.ANC300,PhotoType,QD_counter,app.pyueye_initialization_return); % trying to land on the exact dot 
AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
AppLogger_Specific(app,AppLogger_Specific_text,app.Darkgreen); 

% Updating all visuals on GUI
QD_Location_text = sprintf("[%d %d]",QD_counter); 
app.QD_location_SpecificQD.Text = QD_Location_text; 

Current_QD_Position_text = sprintf("Nanowire_Photo_[%d %d].jpg",QD_counter); 
app.Current_QD_Position_Specific.ImageSource = Current_QD_Position_text; 

for row = 1:num_of_breaks_row

    % Moving to the next set of QD along row
    QD_counter = QD_counter + groups_of_dots_row*row_interval; 
    AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
    AppLogger_Specific(app,AppLogger_Specific_text,'black');
    app.MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_row,app.ANC300,Frequency_fast_travel,x_factor,y_factor)

    % Making sure LED is exactly on the dot 
    app.MyFuncs.QD_tracking_N_Identification(QD_counter,"Auto","Write","LAB",direction_row); 
    app.MyFuncs.Precision_Locking(app.ANC300,PhotoType,QD_counter,app.pyueye_initialization_return); % trying to land on the exact dot 
    AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
    AppLogger_Specific(app,AppLogger_Specific_text,app.Darkgreen); 
    
    % Updating all visuals on GUI
    QD_Location_text = sprintf("[%d %d]",QD_counter); 
    app.QD_location_SpecificQD.Text = QD_Location_text; 

    Current_QD_Position_text = sprintf("Nanowire_Photo_[%d %d].jpg",QD_counter); 
    app.Current_QD_Position_Specific.ImageSource = Current_QD_Position_text; 

end

% moving the remaing QD 
QD_counter = QD_counter + remaining_QD_row*row_interval; 
AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
AppLogger_Specific(app,AppLogger_Specific_text,'black');
app.MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_row,app.ANC300,Frequency_fast_travel,x_factor,y_factor)

% Making sure LED is exactly on the dot 
app.MyFuncs.QD_tracking_N_Identification(QD_counter,"Auto","Write","LAB",direction_row); 
[StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking(app.ANC300,PhotoType,QD_counter,app.pyueye_initialization_return); % trying to land on the exact dot 
AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
AppLogger_Specific(app,AppLogger_Specific_text,app.Darkgreen); 

% Updating all visuals on GUI
QD_Location_text = sprintf("[%d %d]",QD_counter); 
app.QD_location_SpecificQD.Text = QD_Location_text; 

Current_QD_Position_text = sprintf("Nanowire_Photo_[%d %d].jpg",QD_counter); 
app.Current_QD_Position_Specific.ImageSource = Current_QD_Position_text; 

% notifying user the inquired dot has been reached 
AppLogger_Specific_text = sprintf("Arrived to inquired QD: [%d %d]",QD_counter); 
AppLogger_Specific(app,AppLogger_Specific_text,app.Darkgreen); 


if app.PlotDropDown.Value == "No Plot"
    AppLogger_Specific_text = sprintf("No plot was created"); 
    AppLogger_Specific(app,AppLogger_Specific_text,'black'); 

elseif app.PlotDropDown.Value == "Make Plot"
    % Spectroscope Image and Plot
    pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d]',QD_counter); 
    try
        py.asi_func.snap_image(app.asi_initialization_return, pyrun_file_text_Spec)
    catch
        AppLogger_RasterScan(app,"An Error Occured Re-attempting",app.DarkRed)
        terminate(pyenv)
        pyenv("ExecutionMode","OutOfProcess")
        py.asi_func.snap_image(app.asi_initialization_return, pyrun_file_text_Spec)
    end
    AppLogger_Specific_text = sprintf("plot was created and saved"); 
    AppLogger_Specific(app,AppLogger_Specific_text,'black'); 

else
    AppLogger_Specific_text = sprintf("Invalid option for SpectrometerPlot");
    AppLogger_Specific(app,AppLogger_Specific_text,'red'); 
end

% -----------------------------------------------------------------------