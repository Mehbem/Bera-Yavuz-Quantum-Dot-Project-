% Parameters
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Bera Yavuz
% ANC300 Movement testing 

% Object with all required functions
funcs = functionsContainer;

% Adding needed pathways depending on Device - for lab purposes always use device name LAB 
funcs.AddPathFunc("LAB") 

% Establishing serialconnetion with ANC300 device
ANC300 = serialport("COM7",9600); % Change X to the COM port connected

% Respective voltage and frequency values (Voltage: 1-20 V Frequency: 1-1000 Hz)
fprintf(ANC300,"setv 1 10"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp"); 
Frequency = 20; 

% Spectroscope plot saving set environment variables 
setenv('TCL_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tcl8.6")
setenv('TK_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tk8.6")

% Python environment command ('Out-of-order') 
pyenv("ExecutionMode","OutOfProcess")

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
num_points_circle = 8;  % Number of points around the circle
radius_circle = 50;     % Radius of the circle (in QD units or pixels)



% Circular Movement Loop
theta = linspace(0, 2*pi, num_points_circle + 1);  % Create angles for circle points
theta(end) = [];  % Remove the last point (as it overlaps with the first)

for i = 1:length(theta)
    % Calculate circular coordinates relative to the current QD
    x_move = radius_circle * cos(theta(i));
    y_move = radius_circle * sin(theta(i));
    
    % Convert to actual movement based on ANC300's scaling factors
    x_move_real = x_move * x_factor;
    y_move_real = y_move * y_factor;
    
    % serial command for which direction is needed to be moved towards 
    if x_move_real < 0 % x direction 
    X_Serial_Comd = sprintf("stepu 1 %d",abs(x_move_real));
    else
    X_Serial_Comd = sprintf("stepd 1 %d",abs(x_move_real));
    end 
    
    if  y_move_real < 0 % y direction
    Y_Serial_Comd = sprintf("stepu 2 %d",abs(y_move_real));
    else
    Y_Serial_Comd = sprintf("stepe 2 %d",abs(y_move_real));
    end

    % sending the serial command for movement to occur 
    % X movement 
    if x_move_real ~= 0 
    fprintf(ANC300,X_Serial_Comd);
    end
    StepQueue(obj,x_move_real,Frequency); 
    % Y movement 
     if y_move_real ~= 0 
    fprintf(ANC300,Y_Serial_Comd);
    end
    StepQueue(obj,y_move_real,Frequency); 

   
    
    % Snap image after each move if needed
    pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d] position#%d',QD_counter,i); 
    py.asi_func.snap_image(asi_initialization_return, pyrun_file_text_Spec)
end
