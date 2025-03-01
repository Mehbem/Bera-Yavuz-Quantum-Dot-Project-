% Bera Yavuz 

% Object with all required functions
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
% Establishing serialconnetion with ANC300 device
%ANC300 = serialport("COM7",9600); % Change X to the COM port connected
% 
% fprintf(ANC300,"setm 3 gnd"); 
% fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
% fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp");  
%Frequency = 20;



  % Assigning the default settings to ASI
            ASI_Settings.desired_exposure_time = 2;
            ASI_Settings.Gain = 300;
            ASI_Settings.Gamma = 50; 
            ASI_Settings.Brightness = 50; 

            % Assigning the default settings to UI
            UI_Settings = ""; % no settings as of yet due to only require basic snapping function 
            [vid_ASI,src_ASI,vid_UI,src_UI] = MyFuncs.ASI_UI_CameraInit(ASI_Settings,UI_Settings); 

                      % Vertically Flipping the UI camera to the proper orientation   
            src_UI.VerticalFlip = 'on'; 
            src_UI.ExposureMode = "manual"; 
            src_UI.GainMode = "manual";
            src_UI.ContrastMode = "manual";

            % UI camera capturing parameters
            vid_UI.FramesPerTrigger = 1;
            vid_UI.TriggerRepeat = Inf;
            vid_UI.Timeout = 20; 
            triggerconfig(vid_UI,'manual')
            start(vid_UI)

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
                
                % Create a figure with specified dimensions and position
                app.UI_Fig = figure('Units', 'pixels', 'Position', [left, bottom, figWidth, figHeight],"NumberTitle","off");
                
                % Create axes within the figure (slightly smaller for padding)
                axes('Units', 'pixels', 'Position', [10, 10, figWidth - 20, figHeight - 20]);
                
                % Get the video resolution and number of bands (for RGB or grayscale)
                vidRes = get(vid_UI, 'VideoResolution');
                nBands = get(vid_UI, 'NumberOfBands');
                
                % Initialize the image handle in the specified axes, with the same resolution as the video
                hImage = image(zeros(vidRes(2), vidRes(1), nBands));
                
                % Start the preview on the video input object, passing the axes handle for display
                preview(vid_UI, hImage); 
