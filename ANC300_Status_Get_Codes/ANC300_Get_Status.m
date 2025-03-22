% ANC300 Getting status testing

% Establishing serial connection with ANC300 
%ANC300 = serialport("COM7",9600);

% Setting proper voltages, frequencies, and stepping modes
fprintf(ANC300,"setm 3 gnd"); 
fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp");  

% Number of steps to take
step_num = 200; 
serial_comd = sprintf("stepu 2 %d",step_num);

% % send serial command to step
fprintf(ANC300,serial_comd);
pause(0.2)


% Waits until movement is safe 
axis = 2;
timeout = 45; 
step_queue_improved_v2(ANC300,axis,timeout)


function step_queue_improved(ANC300,axis_ID,timeout)
        % step_queue_improved - Monitors ANC300 piezo movement until voltage reaches zero.
        %
        % Inputs:
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
            pause(0.2)
            voltage = ""; % Initalize voltage as an empty string
            tic % Start timer
            while voltage ~= "0.000000"
                    serial_comd_get = sprintf("geto %d",axis_ID);
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
                     % Display error message
                    error("Movement is taking longer then expected");
                    end
                    pause(0.2)
            end
            pause(0.5) %extra 0.5 second delay for safe keeping 
            disp("voltage hit zero")
        end

        function step_queue_improved_v2(ANC300,axis_ID,timeout)
            flush(ANC300)
            pause(0.1)
            voltage = 2; 
            tic % Start timer
            while voltage ~= 0
                serial_comd_get = sprintf("geto %d",axis_ID);
                fprintf(ANC300,serial_comd_get); % Sending command to read current voltage output
                pause(0.1)
                data = read(ANC300, ANC300.NumBytesAvailable, "uint8"); % Read all available bytes as uint8
                data_text = char(data);
                % Regular expression to extract the number between "voltage =" and "V"
                matches = regexp(data_text, 'voltage\s*=\s*([\d\.\-eE]+)\s*V', 'tokens');
                
                % Extract and convert to numeric format
                if ~isempty(matches)
                    voltage = str2double(matches{1}{1});
                    disp(voltage);
                else
                    disp("Voltage value not found.");
                end

                    flush(ANC300) % getting rid of buffered text to prevent filling up
        
                    if toc > timeout % checking if stepping is going on longer then expected
                     % Display error message
                    error("Movement is taking longer then expected");
                    end
                    pause(0.2)

            end
            pause(0.5) %extra 0.5 second delay for safe keeping 
            disp("voltage hit zero")
        end