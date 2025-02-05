% ID900 Connect Test

% Define the ID900 IP address
ip_address = '169.254.102.137'; % Change this to match your device's IP

% Open a TCP connection
t = tcpclient(ip_address, 5555); % Default SCPI port

% Configure the input channel
writeline(t, 'INPU1:ENAB ON'); % Enable input 1
writeline(t, 'INPU1:THREshold 100mV'); % Set threshold (adjust as needed)
writeline(t, 'INPU1:EDGE RISING'); % Set to trigger on rising edge
writeline(t, 'INPU1:MODE ACCUM'); % Set counter mode to accumulate counts
writeline(t, 'INPU1:RESEt'); % Reset the counter before starting

pause(1); % Give it a second to collect data

% Query the current photon count
writeline(t, 'INPU1:COUNter?'); 
photon_count = str2double(readline(t)); % Read and convert response

% Display the photon count
fprintf('Current photon count: %d\n', photon_count);

% Close the TCP connection
clear t;
