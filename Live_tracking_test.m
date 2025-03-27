% Bera Yavuz Live tracking testing 

% % Inserting required pathways
addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools\Scripts For Testing_Debugging\Testing Images")
my_obj = functionsContainer; 

 % finding current imaging devices connected 
info = imaqhwinfo;
if isempty(info.InstalledAdaptors)
    error('No cameras found. Connect a camera and try again.');
end
adaptor = info.InstalledAdaptors{1};
camInfo = imaqhwinfo(adaptor);
numCameras = numel(camInfo.DeviceIDs); 

% Identify each camera
for i = 1:numCameras
    deviceID = camInfo.DeviceIDs{i}; % Get device ID
    deviceName = camInfo.DeviceInfo(i).DeviceName; % Get device name
    fprintf('Camera %d: ID = %d, Name = %s\n', i, deviceID, deviceName);
    if contains(deviceName,"ASI") %  checks to see if the detected device is the ASI camera
        ASI_Device_ID = i;
        ASI_name = deviceName; 
        camInfo.DeviceInfo(i).DefaultFormat = 'RGB8_6248x4176';%'RGB8_1280x960'; 
        fprintf("found ASI (spectrometer camera)\n")
    elseif contains(deviceName,"UI") %  checks to see if the detected device is the ASI camera
        UI_Device_ID = i;
        UI_name = deviceName; 
        fprintf("found UI 148x Camera (nanowire camera)\n")
    else % any other camera would have to be an unknown device
        fprintf("unknown device detected\n")
    end
end

% Initialize camera (adjust as needed)
vid_UI = videoinput(adaptor, UI_name);
src_UI = getselectedsource(vid_UI);

 % Vertically Flipping the UI camera to the proper orientation   
src_UI.VerticalFlip = 'on'; 
src_UI.ExposureMode = "manual"; 
src_UI.GainMode = "manual";
src_UI.ContrastMode = "manual";
vid_UI.ReturnedColorSpace = 'grayscale'; % rgb, grayscale, bayer


% Assuming vid_UI is correctly initialized
screenSize = get(0, 'ScreenSize'); % gets pixel resolution of screen
dimensions_UI = vid_UI.VideoResolution; % Get the resolution of the video input
width = dimensions_UI(1); % Video width
height = dimensions_UI(2); % Video height

% Shrink the figure to a fraction of the original resolution
figWidth = round(width / 3.5);
figHeight = round(height / 3.5);

% Position the window on the screen (shifted right and centered vertically)
left = (screenSize(3) - figWidth) * 0.75; 
bottom = (screenSize(4) - figHeight) / 2;

figure('Units', 'pixels', 'Position', [left, bottom, figWidth, figHeight],"NumberTitle","off");
% Create axes within the figure (slightly smaller for padding)
axes('Units', 'pixels', 'Position', [10, 10, figWidth - 20, figHeight - 20]);

% Get the video resolution and number of bands (for RGB or grayscale)
vidRes = get(vid_UI, 'VideoResolution');
nBands = get(vid_UI, 'NumberOfBands');

% Initialize the image handle in the specified axes, with the same resolution as the video
hImage = image(zeros(vidRes(2), vidRes(1), nBands));

% Start the preview on the video input object, passing the axes handle for display
preview(vid_UI, hImage); 
pause(2); % Allow camera to adjust

 
% Capture the first frame
firstFrame = getsnapshot(vid_UI);
 
% Filtering Settings
scaling = 0.5; 
skyBlue = [0.53, 0.81, 0.92]; 
radii_big_circle = 500; 
center_big_circle = [650,500]; 
sigma_flatfield = 35; 
Salt_pepper_pixel_factor = [20 20]; 
Min_circle_area = 200;
radiusQD = 25;
num_lines = 9; % Num lines on each side of the original line
sep_red = 345; % Separation distance red lines
sep_blue = 370; % Seperation distance blue lines

% Extracting QD
[grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = my_obj.Finalized_Analyzed_QDBinaryImg_With_Dots(scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,firstFrame);
[allNextPts,allPerpPts] = my_obj.Finalized_MainAxes(Img,CopyCentroid,radiusQD);
[x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = my_obj.Finalized_GridLines(allNextPts,allPerpPts); 
[VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = my_obj.Finalized_VirtualQD(num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled); 

figure;
imshow(firstFrame); hold on;
plot(FullQDList_sorted(:,1), FullQDList_sorted(:,2), 'go', 'MarkerSize', 10, 'LineWidth', 1.5);
title("First Instance");

% Initialize Point Tracker
tracker = vision.PointTracker('MaxBidirectionalError', 2);
initialize(tracker, FullQDList_sorted, firstFrame);
 
figure;
for i = 1:100 % Track for 100 frames (adjust as needed)
    frame = getsnapshot(vid_UI);
    % Track dots
    [trackedPoints, isValid] = step(tracker, frame);
    % Display results
    imshow(frame); hold on;
    plot(trackedPoints(isValid,1), trackedPoints(isValid,2), 'go', 'MarkerSize', 10, 'LineWidth', 1.5);
    title(['Frame ', num2str(i)]);
    hold off;
 
    pause(0.5); % Adjust based on frame rate
end
 
clear vid_UI; % Release camera