% ANC300 Getting status testing

% Establishing serial connection with ANC300 
ANC300 = serialport("COM8",9600);

% Setting proper voltages, frequencies, and stepping modes
fprintf(ANC300,"setm 3 gnd"); 
fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp");  

% Number of steps to take
step_num = 200; 
serial_comd = sprintf("stepu 1 %d",step_num);

% send serial command to step
fprintf(ANC300,serial_comd)

for i = 1:10
    voltage = query(ANC300, "geto 1");
    fprintf("%s\n",voltage)
    pause(2)
end