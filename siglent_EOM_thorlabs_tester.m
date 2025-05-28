% Initialize connection to Thorlabs power meter 
meter_list=ThorlabsPowerMeter;                              % Initiate the meter_list
DeviceDescription=meter_list.listdevices;               	% List available device(s)
test_meter=meter_list.connect(DeviceDescription);           % Connect single/the first devices

N_sample = 20;
DC_voltage_vals = linspace(0, 5, N_sample);

% Initialize Siglent waveform 
s_SIGLENT = pyrunfile('siglent_EOM_initialize.py', 's_SIGLENT');
% Set DV voltage value to the EOM
pyrunfile('siglent_EOM_set_vals.py', s_SIGLENT=s_SIGLENT, DC_voltage_val=DC_voltage_vals();


test_meter.disconnect;                                      % Disconnect and release