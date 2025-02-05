      addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 

ASI_Settings.desired_exposure_time = 2;
ASI_Settings.Gain = 300;
ASI_Settings.Gamma = 50; 
ASI_Settings.Brightness = 50; 
UI_Settings = ""; 


[vid_ASI,src_ASI,vid_UI,src_UI,CamInfo] = MyFuncs.ASI_UI_CameraInit(ASI_Settings,UI_Settings);



src_UI.VerticalFlip = 'on'; 

exposure = propinfo(src_UI,"Exposure");
src_UI.ExposureMode = "manual"; 
src_UI.GainMode = "manual";
src_UI.ContrastMode = "manual";

src_UI.BacklightCompensation = "on"; 
vid_UI.ReturnedColorSpace = 'grayscale'; 


% Assuming vid_UI is correctly initialized
screenSize = get(0, 'ScreenSize'); % gets pixel resolution of screen
dimensions_UI = vid_UI.VideoResolution; % Get the resolution of the video input
width = dimensions_UI(1); % Video width
height = dimensions_UI(2); % Video height

% Shrink the figure to a fraction of the original resolution
figWidth = width / 3.5;
figHeight = height / 3.5;

% Position the window on the screen (shifted right and centered vertically)
left = (screenSize(3) - figWidth) * 0.75; 
bottom = (screenSize(4) - figHeight) / 2;

% Create a figure with specified dimensions and position
UI_fig = figure('Units', 'pixels', 'Position', [left, bottom, figWidth, figHeight],"NumberTitle","off");

% Create axes within the figure (slightly smaller for padding)
axes('Units', 'pixels', 'Position', [10, 10, figWidth - 20, figHeight - 20]);

% Get the video resolution and number of bands (for RGB or grayscale)
vidRes = get(vid_UI, 'VideoResolution');
nBands = get(vid_UI, 'NumberOfBands');

% Initialize the image handle in the specified axes, with the same resolution as the video
hImage = image(zeros(vidRes(2), vidRes(1), nBands));

% Start the preview on the video input object, passing the axes handle for display
preview(vid_UI, hImage); 

imaqreset