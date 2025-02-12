% -------------------- INITIALIZE CAMERA --------------------
clc; clear; close all;
 
% Set up camera
info = imaqhwinfo;
if isempty(info.InstalledAdaptors)
    error('No cameras found. Connect a camera and try again.');
end
adaptor = info.InstalledAdaptors{1};
camInfo = imaqhwinfo(adaptor);
numCameras = numel(camInfo.DeviceIDs); 
% Check if at least two cameras are connected (making sure the UI and ASI
% camera are there) 
if numCameras < 2
    error('Two cameras are required but only %d detected.', numCameras);
end
 
% Identify each camera
for i = 1:numCameras
    deviceID = camInfo.DeviceIDs{i}; % Get device ID
    deviceName = camInfo.DeviceInfo(i).DeviceName; % Get device name
    fprintf('Camera %d: ID = %d, Name = %s\n', i, deviceID, deviceName);
    if contains(deviceName,"ASI") %  checks to see if the detected device is the ASI camera
        ASI_Device_ID = i; 
        fprintf("found ASI (spectrometer camera)\n")
    elseif contains(deviceName,"UI") %  checks to see if the detected device is the ASI camera
        UI_Device_ID = i;
        fprintf("found UI 148x Camera (nanowire camera)\n")
    else % any other camera would have to be an unknown device
        fprintf("unknown device detected\n")
    end
end

% establishing connection and parameters of ASI Device 
vid_ASI = videoinput(adaptor, ASI_Device_ID); % function can take a third input to specify formatting (ASK Sreesh) 
src_ASI = getselectedsource(vid_ASI);
all_props_ASI = propinfo(vid_ASI); 

% establishing connection and parameters of UI Device 
vid_UI = videoinput(adaptor, UI_Device_ID);
src_UI = getselectedsource(vid_UI);
all_props_UI = propinfo(vid_UI); 

% Setting correct Parameters for ASI camera
desired_exposure_time = 2; % unit: seconds 
proper_exposure_setting = round(log2(desired_exposure_time)); %2^-15 seconds to 2^11 seconds 
set(src_ASI,"ExposureMode","manual", "Exposure", proper_exposure_setting); %exposure is calculated as 2^n seconds by the camera (n = [-15 11])
set(src_ASI,"GainMode", "manual" ,"Gain", 300); % setting gain 
set(src_ASI, "GammaMode","manual","Gamma",50); % setting gamma
set(src_ASI,"Brightness",50); % setting brightness



% adding an itial delay to make sure everything is working as expected 
pause(1)

preview(vid)
 
% Capture initial image
img = getsnapshot(vid);
imshow(img)

closepreview(vid)
stop(vid)
imaqreset
delete(vid)
clear vid
