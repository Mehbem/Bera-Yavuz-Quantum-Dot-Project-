% Bera Yavuz
% ANC300 Settings and Debugging Purpose

clear;
clc;

% Adding all required pathways for functions 
pathway_all_functions = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools\Scripts For Testing_Debugging";
addpath(pathway_all_functions)
% file that contains all the functions 
funcs = functionsContainer; 

% Establishing serialconnetion with ANC300 device
ANC300 = serialport("COM7",9600); % Change X to the COM port connected

% desired ANC300 Settings
% -----------------------------------------------------------------------------------
% X-axis settings
Mode_X = "stp"; % sets the mode of the axis (Possible options: "gnd" or "stp")
Frequency_X = 12; % Possible options: 1 - 10000 (must be integer)
Voltage_X = 20; % Possible options: 1 - 150.0 (can be a float value)

serial_com_Mode_X = sprintf("setm 1 %s",Mode_X); fprintf(ANC300,serial_com_Mode_X)
serial_com_Frequency_X = sprintf("setf 1 %s", num2str(Frequency_X)); fprintf(ANC300,serial_com_Frequency_X)
serial_com_Voltage_X = sprintf("setv 1 %s", num2str(Voltage_X)); fprintf(ANC300,serial_com_Voltage_X)

% Y-axis settings 
Mode_Y = "stp"; % sets the mode of the axis (Possible options: "gnd" or "stp")
Frequency_Y = 30; % Possible options: 1 - 10000 (must be integer)
Voltage_Y = 30; % Possible options: 1 - 150.0 (can be a float value)

serial_com_Mode_Y = sprintf("setm 2 %s",Mode_Y); fprintf(ANC300,serial_com_Mode_Y)
serial_com_Frequency_Y = sprintf("setf 2 %s", num2str(Frequency_Y)); fprintf(ANC300,serial_com_Frequency_Y)
serial_com_Voltage_Y = sprintf("setv 2 %s", num2str(Voltage_Y)); fprintf(ANC300,serial_com_Voltage_Y)
% -----------------------------------------------------------------------------------
