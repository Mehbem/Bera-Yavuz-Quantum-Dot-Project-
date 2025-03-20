% ANC300 Getting status testing

% Establishing serial connection with ANC300 
%ANC300 = serialport("COM7",9600);

% Setting proper voltages, frequencies, and stepping modes
%fprintf(ANC300,"setm 3 gnd"); 
%fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
%fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp");  

% Number of steps to take
step_num = 250; 
serial_comd = sprintf("stepd 1 %d",step_num);

% % send serial command to step
fprintf(ANC300,serial_comd);
pause(0.2)
flush(ANC300) % getting rid of buffered text to prevent filling up 


% Waits until movement is safe 
axis = 1;
timeout = 45; 
function step_queue_improved(ANC300,axis,timeout)
% step_queue_improved - Monitors ANC300 piezo movement until voltage reaches zero.
%
% Inputs:
%   ANC300  - Serial object for the ANC300 piezo controller.
%   axis    - Axis to monitor (1 = x, 2 = y).
%   timeout - Maximum wait time in seconds before error.
%
% Description:
%   Continuously queries the ANC300 for the output voltage of the specified axis.  
%   The loop runs until the voltage reaches "0.000000" V, flushing the buffer  
%   after each read to prevent overflow. If movement exceeds the timeout,  
%   an error is triggered. A short delay (0.2s) reduces excessive polling.  
%
% Example:
%   step_queue_improved(ANC300, 1, 10) % Monitor x-axis with a 10s timeout.

    flush(ANC300) % Flush previous undesirable data 
    voltage = ""; % Initalize voltage as an empty string
    tic % Start timer
    while voltage ~= "0.000000"
            serial_comd_get = sprintf("geto %d",axis);
            fprintf(ANC300,serial_comd_get); % Sending command to read current voltage output
            for i = 1:2
                serial_message = fscanf(ANC300); % reading current voltage output
                if i == 2
                    serial_message = string(serial_message);
                    serial_messages = strsplit(serial_message);
                    voltage = serial_messages(3); % extract voltage number from list
                end
            end
            flush(ANC300) % getting rid of buffered text to prevent filling up

            if toc > timeout % checking if stepping is going on longer then expected
            error("movement is taking longer then expected potential issue detected")
            end
            pause(0.2)
    end
    pause(0.5) %extra 0.5 second delay for safe keeping 
end


