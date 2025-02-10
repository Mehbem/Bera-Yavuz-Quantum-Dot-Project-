% Bera Yavuz 

% Object with all required functions
% --------------------------------------------------------------------------------------------------------
addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 

% Initilizing the necessary devices
% --------------------------------------------------------------------------------------------------------
ANC300 = serialport("COM7",9600); % Establishing serialconnetion with ANC300 device
ell_motor = serialport("COM6",9600); % Establishing a serialconnection with the HWP motor


% initiziling all the needed folders 
py.qd_data_folder_creation.create_qd_data_directories()
            

% defining parameters and initilizing everything 
angle = 0:5:360;
QD_ID = [67 1];


            % Spectroscope plot saving set environment variables 
            terminate(pyenv); % Clear Python 
            setenv('TCL_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tcl8.6")
            setenv('TK_LIBRARY', "C:\Users\Quantum Dot\AppData\Local\Programs\Python\Python311\tcl\tk8.6")

            % Inserting proper pathways needed 
            pyLibraryFolder_Scripts = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools"; 
            insert(py.sys.path, int64(0), pyLibraryFolder_Scripts)
            pause(1)

% taking a snap at angle 0 
asi_initialization_return = py.asi_func.init_camera();
pause(2)
pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d]_0_Degrees',QD_ID); 
py.asi_func.snap_image_FSS(asi_initialization_return, pyrun_file_text_Spec)

for rot_count = 1:length(angle)
% Defining the movement serial code for the rotation 
angle_hxd = dec2hex(floor(mod(angle(rot_count),360)*39822/100), 8);
input_str = "1ma" + angle_hxd; % "2" before ma is to be get from the ELLO software from thorlabs.

% Commiting Command for movement 
fprintf(ell_motor, input_str);


fprintf('Angle: %d \n', angle(rot_count));

pause(1)

% taking a snap at every respective angle 
pyrun_file_text_Spec = sprintf('QD_Spec_Plot_[%d %d]_%d_Degrees',QD_ID,angle(rot_count)); 
py.asi_func.snap_image_FSS(asi_initialization_return, pyrun_file_text_Spec)
end

%[43, 7], [67, 1], 