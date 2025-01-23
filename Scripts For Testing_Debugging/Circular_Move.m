% Parameters
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Bera Yavuz
% ANC300 Movement testing 

% % Object with all required functions
% funcs = functionsContainer;
% 
% % Adding needed pathways depending on Device - for lab purposes always use device name LAB 
% funcs.AddPathFunc("LAB") 
% 
% % Establishing serialconnetion with ANC300 device
% ANC300 = serialport("COM7",9600); % Change X to the COM port connected
% 
% % Respective voltage and frequency values (Voltage: 1-20 V Frequency: 1-1000 Hz)
% fprintf(ANC300,"setv 1 10"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
% fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp"); 
% Frequency = 20; 
% 
% % Spectroscope plot saving set environment variables 
% setenv('TCL_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tcl8.6")
% setenv('TK_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tk8.6")
% 
% % Python environment command ('Out-of-order') 
% pyenv("ExecutionMode","OutOfProcess")

% X and y factors
Read_XY_factor = funcs.XY_Factor_Identifier("","Read","LAB");
x_factor = Read_XY_factor.X_factor;
y_factor = Read_XY_factor.Y_factor;


% Both Camera Initilizations 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pyLibraryFolder_Scripts = 'C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Python_Scripts'; 
insert(py.sys.path, int64(0), pyLibraryFolder_Scripts)

% Initialize both ze cameras
pyueye_initialization_return = py.pyueye_func.init_camera();
fprintf('uEye camera has been initialized, yo\n')

pause(2)

asi_initialization_return = py.asi_func.init_camera();
fprintf('ASI camera has been initialized, yo\n')

pause(2)

% Initial Movement Towards Start
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PhotoType = "SteppingPhoto";
readQD_position = funcs.QD_tracking_N_Identification("","","Read","LAB","");
QD_counter = readQD_position.lastLine; 
funcs.Precision_Locking(ANC300,PhotoType,QD_counter,pyueye_initialization_return);


% Parameters for circular movement
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
num_points_circle = 4;  % Number of points around the circle
radius_circle = 50;     % Radius of the circle (in QD units or pixels)

% Defining the matrix used for finding desirable point 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Photon_Count_at_Positions = zeros(1,4);
Positions = zeros(num_points_circle, 2);  % To store positions (x, y) in real units
Current_position = [0 0]; % middle of the dot 

% Circular Movement Loop
theta = linspace(0, 2*pi, num_points_circle + 1);  % Create angles for circle points
theta(end) = [];  % Remove the last point (as it overlaps with the first)

for i = 1:length(theta)
    % Calculate circular coordinates relative to the current QD
    x_move = radius_circle * cos(theta(i));
    y_move = radius_circle * sin(theta(i));
    
    % Convert to actual movement based on ANC300's scaling factors
    x_move_real = round(x_move * x_factor);
    y_move_real = round (y_move * y_factor);
    
    % Save position for tracking
    Positions(i, :) = [x_move_real, y_move_real];
    
    % Calculate relative movement from the current position
    x_displacement = x_move_real - current_position(1);
    y_displacement = y_move_real - current_position(2);
    
    % Update the current position
    current_position = [x_move_real, y_move_real];
    
    % Serial command for movement in X direction
    if x_displacement < 0
        X_Serial_Comd = sprintf("stepu 1 %d", abs(x_displacement));
    else
        X_Serial_Comd = sprintf("stepd 1 %d", abs(x_displacement));
    end
    
    % Serial command for movement in Y direction
    if y_displacement < 0
        Y_Serial_Comd = sprintf("stepu 2 %d", abs(y_displacement));
    else
        Y_Serial_Comd = sprintf("stepe 2 %d", abs(y_displacement));
    end

    % Send the serial commands to move in X and Y directions
    if x_displacement ~= 0
        fprintf(ANC300, X_Serial_Comd);
        funcs.StepQueue(x_displacement, Frequency);
    end

    if y_displacement ~= 0
        fprintf(ANC300, Y_Serial_Comd);
        funcs.StepQueue(y_displacement, Frequency);
    end
    

    % Snap image after each move
    pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d] position#%d', QD_counter, i); 
    py.asi_func.snap_image(asi_initialization_return, pyrun_file_text_Spec);
    
    % Analyze spectrometry and determine desirable value, adding to the matrix
    Photon_Count_at_Positions(i) = photon_count;  % Variable to be assigned
end

% Determine the position with the maximum photon count
[~, best_position_index] = max(Photon_Count_at_Positions);
best_position = Positions(best_position_index, :);

% Trace back to the best position using relative movement
fprintf('Returning to best position: (X = %.2f, Y = %.2f)\n', best_position(1), best_position(2));

% Calculate relative displacement to the best position
x_displacement_to_best = best_position(1) - current_position(1);
y_displacement_to_best = best_position(2) - current_position(2);

% Serial command for returning in X direction
if x_displacement_to_best < 0
    X_Serial_Comd = sprintf("stepu 1 %d", abs(x_displacement_to_best));
else
    X_Serial_Comd = sprintf("stepd 1 %d", abs(x_displacement_to_best));
end

% Serial command for returning in Y direction
if y_displacement_to_best < 0
    Y_Serial_Comd = sprintf("stepu 2 %d", abs(y_displacement_to_best));
else
    Y_Serial_Comd = sprintf("stepe 2 %d", abs(y_displacement_to_best));
end

% Send the commands to return to the best position
if x_displacement_to_best ~= 0
    fprintf(ANC300, X_Serial_Comd);
    funcs.StepQueue(x_displacement_to_best, Frequency);
end


if y_displacement_to_best ~= 0
    fprintf(ANC300, Y_Serial_Comd);
    funcs.StepQueue(y_displacement_to_best, Frequency);
end


fprintf('Reached the best position based on photon count.\n');

