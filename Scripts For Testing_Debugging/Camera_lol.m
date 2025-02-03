% Background Initilization
% Assigning the default settings to ASI
ASI_Settings.desired_exposure_time = 2;
ASI_Settings.Gain = 300;
ASI_Settings.Gamma = 50; 
ASI_Settings.Brightness = 50; 
UI_Settings = ""; 
addpath("C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools")
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 

[vid_ASI,src_ASI,vid_UI,src_UI,CamInfo] = MyFuncs.ASI_UI_CameraInit(ASI_Settings,UI_Settings);
fig = figure();
camsize = vid_UI.VideoResolution;
im = image(zeros(camsize));

preview(vid_UI);     %pass the image object to preview!
