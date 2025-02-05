% Bera Yavuz 

clear;
clc;
% Object with all required functions
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
% Establishing serialconnetion with ANC300 device
ANC300 = serialport("COM7",9600); % Change X to the COM port connected

fprintf(ANC300,"setm 3 gnd"); 
fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp");  
Frequency = 20; 

% Initilize photon detector 
tc = py.ID900_Func.init_ID900();

% defining the smallest movement margin (ANC300 steps)
X_Step = 4;  init_move_mult_X = 8;
Y_Step = 4;  init_move_mult_Y = 8;

% Defining starting position Movement
Init_Y_Pos = init_move_mult_Y*Y_Step;
Y_movement_seria_comd_init = sprintf("stepd 2 %d",Init_Y_Pos);

Init_X_Pos = init_move_mult_X*X_Step;
X_movement_seria_comd_init = sprintf("stepd 1 %d",Init_X_Pos);

% Defining smallest stepping Movemment
Y_movement_seria_comd_back = sprintf("stepd 2 %d",Y_Step);Y_movement_seria_comd_forward = sprintf("stepu 2 %d",Y_Step);
X_movement_seria_comd_back = sprintf("stepd 1 %d",X_Step);X_movement_seria_comd_forward = sprintf("stepu 1 %d",X_Step);

% Defining region to scan
x_rows = init_move_mult_X*2;
y_rows = init_move_mult_Y*2;

% Defining an empty 0 0 matrix for tracking location 
relative_loc = [0 0];

% Going to starting position 
fprintf(ANC300,Y_movement_seria_comd_init)
MyFuncs.StepQueue(Init_Y_Pos,Frequency)
fprintf(ANC300,X_movement_seria_comd_init)
MyFuncs.StepQueue(Init_X_Pos,Frequency)




% get photon count 
init_photon_count = py.ID900_Func.query_photon_counter(tc); 

 

data_table = table(relative_loc,init_photon_count,'VariableNames', {'Positions', 'PeakCounts'}); 



for y_move = 1:y_rows
    if y_move ~= 1
        relative_loc = relative_loc + [0 1]; 
        fprintf(ANC300, Y_movement_seria_comd_forward);
        MyFuncs.StepQueue(Y_Step, Frequency);

        % Read photon count
        photon_count = py.ID900_Func.query_photon_counter(tc); 

        % Update table
        new_row = {relative_loc, photon_count};
        data_table = [data_table; new_row]; 
    end

    for x_move = 1:x_rows
        if mod(y_move, 2) == 0
            relative_loc = relative_loc + [0 -1];

            fprintf(ANC300, X_movement_seria_comd_back);
            MyFuncs.StepQueue(X_Step, Frequency);

        else
            relative_loc = relative_loc + [0 +1];

            fprintf(ANC300, X_movement_seria_comd_forward);
            MyFuncs.StepQueue(X_Step, Frequency);
        end

        % Read photon count
        photon_count = py.ID900_Func.query_photon_counter(tc); 

        % Update table
        new_row = {relative_loc, photon_count};
        data_table = [data_table; new_row]; 
    end
end



[A, i] = max(data_table.PeakCounts);

position_to_return = data_table.Positions(i,:);

