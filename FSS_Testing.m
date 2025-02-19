% Bera Yavuz 

% Object with all required functions
% --------------------------------------------------------------------------------------------------------
addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 

% Initilizing the necessary devices
% --------------------------------------------------------------------------------------------------------
%ANC300 = serialport("COM7",9600); % Establishing serialconnetion with ANC300 device
ell_motor = serialport("COM10",9600); % Establishing a serialconnection with the HWP motor


            

% defining parameters and initilizing everything 
angle = 0:5:360;
QD_counter = [38 1];
Spectrometer_Gratting = 1800; 


% defining angle of rotation 
angle = 0:5:360;

% for loop for rotating motor 
for rot_count = 1:length(angle)

% Defining the movement serial code for the rotation 
angle_hxd = dec2hex(floor(mod(angle(rot_count),360)*39822/100), 8);
input_str = "1ma" + angle_hxd; % "2" before ma is to be get from the ELLO software from thorlabs.

% Commiting Command for movement 
fprintf(ell_motor, input_str);


Current_angle = sprintf('Angle: %d \n', angle(rot_count));
fprintf(Current_angle)
pause(1)

% taking a snap at every respective angle 
file_name = sprintf('FSS %d',angle(rot_count)); 
file_name_sting = string(file_name);
MyFuncs.ASI_Snap_Img(vid_ASI,src_ASI,"Spectrometer","Yes",Spectrometer_Gratting,QD_counter,file_name_sting);
end


