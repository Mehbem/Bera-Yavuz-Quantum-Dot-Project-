
% Object with all required functions
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
% Establishing serialconnetion with ANC300 device
ANC300 = serialport("COM7",9600); % Change X to the COM port connected
pyueye_initialization_return = py.pyueye_func.init_camera();
pause(2)

asi_initialization_return = py.asi_func.init_camera();
pause(2)

fprintf(ANC300,"setm 1 stp"); 
fprintf(ANC300,"setm 2 stp"); 
Frequency = 20; 
Read_XY_factor = MyFuncs.XY_Factor_Identifier("","Read","LAB");
x_factor = Read_XY_factor.X_factor;
y_factor = Read_XY_factor.Y_factor;
                


% Initial Movement Towards Start
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
CurrentQD =  [67 1];
QD_counter = CurrentQD; 
InquiredQD = [43 7]; % user puts what value they want
PhotoType = "SteppingPhoto";
fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20");
fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20");
pause(0.3)
[StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = MyFuncs.Precision_Locking(ANC300,PhotoType,QD_counter,pyueye_initialization_return,40);
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
dots_to_travel_row = num_rows; 
dots_to_travel_column = num_columns;
groups_of_dots_row = 4; 
groups_of_dots_column = 4; 


num_of_breaks_row = floor(dots_to_travel_row / groups_of_dots_row);  % How many times it travels then locks in-between each QD (along the row) 
remaining_QD_row = mod(dots_to_travel_row, groups_of_dots_row);    % What remains that is indiviually travelled (along the row) 

num_of_breaks_column = floor(dots_to_travel_column / groups_of_dots_column);  % How many times it travels then locks in-between each QD (along the column) 
remaining_QD_column = mod(dots_to_travel_column, groups_of_dots_column);    % What remains that is indiviually travelled (along the column) 

Frequency_fast_travel = 30; %whatever frequency we choose for fast travel essentially 

% defining the number of steps
Y_Step_num = 90;
X_Step_num = 73;

% additional step number to increase error margin for pause
error_margin_time = 6; % extra time paused = error_margin_time / frequency

% acceptable error margin (pixel distance from dot) for precision locking 
fast_movement_margin = 80; 
locking_on_margin = 40;
final_locking_on_margin = 25; 

% defining the commands depending on the direction being moved 
row_exist  = exist('direction_row', 'var');
column_exist = exist('direction_column', 'var');

if row_exist
if direction_row == "bottomleft" 
    Y_Serial_Comd = sprintf("stepu 2 %d",Y_Step_num);
elseif direction_row == "topright"
    Y_Serial_Comd = sprintf("stepd 2 %d",Y_Step_num); 
end
end

if column_exist
if direction_column == "bottomright" 
    X_Serial_Comd = sprintf("stepu 1 %d",X_Step_num);
elseif direction_column == "topleft"
    X_Serial_Comd = sprintf("stepd 1 %d",X_Step_num); 
end
end
 


% fast tracking along the QDs column
for column = 1: num_of_breaks_column

    % Moving to the next set of QD along column
    QD_counter = QD_counter + groups_of_dots_column*column_interval; 
    % Respective voltage and frequency values (Voltage: 1-20 V Frequency: 1-1000 Hz)
    fprintf(ANC300,"setv 1 30"); 
    fprintf(ANC300,"setf 1 30"); 
    pause(0.3)
    fprintf(ANC300,X_Serial_Comd);
    MyFuncs.StepQueue(X_Step_num+error_margin_time,Frequency_fast_travel); 

    % Making sure LED is exactly on the dot 
    fprintf(ANC300,"setv 1 12"); 
    fprintf(ANC300,"setf 1 20");
    pause(0.3)
    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = MyFuncs.Precision_Locking(ANC300,PhotoType,QD_counter,pyueye_initialization_return,fast_movement_margin); % trying to land on the exact dot 

end

% moving the remaing QD along column
for column_remaining = 1: remaining_QD_column
    % find next point along the column 
    [~,~,~,DistanceBetweenPoints_Rotated] = MyFuncs.OptimizedRasterScan(StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_column);

    % Moving to the next QD along column
    QD_counter = QD_counter + column_interval; 
    MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_column,ANC300,Frequency,x_factor,y_factor)
    
    % Making sure LED is exactly on the dot 
    MyFuncs.QD_tracking_N_IdentificationVer2(QD_counter,"Auto","Write","LAB",direction_column); 
    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = MyFuncs.Precision_Locking(ANC300,PhotoType,QD_counter,pyueye_initialization_return,locking_on_margin); % trying to land on the exact dot 
end

% fast tracking along the QDs row
for row = 1:num_of_breaks_row

    % Moving to the next set of QD along row
    QD_counter = QD_counter + groups_of_dots_row*row_interval; 
    % Respective voltage and frequency values (Voltage: 1-20 V Frequency: 1-1000 Hz)
    fprintf(ANC300,"setv 2 30"); 
    fprintf(ANC300,"setf 2 30"); 
    pause(0.3)
    fprintf(ANC300,Y_Serial_Comd);
    MyFuncs.StepQueue(Y_Step_num+error_margin_time,Frequency_fast_travel); 

    % Making sure LED is exactly on the dot 
    fprintf(ANC300,"setv 2 12"); 
    fprintf(ANC300,"setf 2 20"); 
    pause(0.3)
    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = MyFuncs.Precision_Locking(ANC300,PhotoType,QD_counter,pyueye_initialization_return,fast_movement_margin); % trying to land on the exact dot 
 
end

% moving the remaing QD along row
for row_remaining = 1: remaining_QD_row
    % find next point along the column 
    [~,~,~,DistanceBetweenPoints_Rotated] = MyFuncs.OptimizedRasterScan(StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_row);

    % Moving to the next QD along column
    QD_counter = QD_counter + row_interval; 
    MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_row,ANC300,Frequency,x_factor,y_factor)
    
    % Making sure LED is exactly on the dot 
    MyFuncs.QD_tracking_N_IdentificationVer2(QD_counter,"Auto","Write","LAB",direction_row); 
    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = MyFuncs.Precision_Locking(ANC300,PhotoType,QD_counter,pyueye_initialization_return,locking_on_margin); % trying to land on the exact dot 
end

% final check to make sure exactly on dot 
MyFuncs.OptimizedRasterScan(StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_row,final_locking_on_margin)

pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d]',InquiredQD);
py.asi_func.snap_image(asi_initialization_return, pyrun_file_text_Spec)
fprintf("process complete")
py.pyueye_func.exit_camera(pyueye_initialization_return)
% -----------------------------------------------------------------------