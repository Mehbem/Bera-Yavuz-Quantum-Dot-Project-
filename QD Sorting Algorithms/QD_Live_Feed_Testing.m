% Bera Yavuz Live Spectrum Generator from Live Feed

% --------------------------- INIT SETUP ---------------------------

% Add paths and initialize ASI camera
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 

ASI_Settings.desired_exposure_time = 0.2;
ASI_Settings.Gain = 300;
ASI_Settings.Gamma = 50; 
ASI_Settings.Brightness = 50; 

UI_Settings = ""; % no UI settings required
[vid_ASI, ~, ~, ~] = MyFuncs.ASI_UI_CameraInit(ASI_Settings, UI_Settings); 

% --------------------------- LIVE FEED + CALLBACK ---------------------------

image = preview(vid_ASI); 

