close all
clear

meter_list=ThorlabsPowerMeter;                              % Initiate the meter_list
DeviceDescription=meter_list.listdevices;               	% List available device(s)
test_meter=meter_list.connect(DeviceDescription);           % Connect single/the first devices

power_meter_readings = zeros(1, 100);

for i=1:1:100   
    test_meter.updateReading(0.5);                          % Update the power reading(with interal period of 0.5s)
    fprintf('%.10f%c\r',test_meter.meterPowerReading,test_meter.meterPowerUnit);
    power_meter_readings(i) = test_meter.meterPowerReading*1e6; % in muW
    truncated_power = fix(power_meter_readings(i) * 100) / 100;
    spectrum_plot_string = 'Spectrometer_power_val_muW_' + truncated_power;
    file_name = sprintf('PowerScan %d', truncated_power); 
    PowerScan_Text = file_name; 
    PowerScan_Text_A = strsplit(PowerScan_Text); 
    Power_Excitation = PowerScan_Text_A{2}; 

end
test_meter.updateReading_V(0.5);                            % To demonstrate that only certain sensors can use this function
                                                            % A warning message is expected here for most of the models
test_meter.disconnect;                                      % Disconnect and release
