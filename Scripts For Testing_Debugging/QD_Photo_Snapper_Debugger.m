% Bera Yavuz 
% Photo Debug Testing 
clear;
clc; 
tic
directoryPath_Scripts = "C:\Users\Quantum Dot\Desktop\Bera Yavuz - ANC300 Movement and Images\Scripts_&_Debugging_Tools"; 
addpath(directoryPath_Scripts)

% Create Folders for the day 
py.qd_data_folder_creation.create_qd_data_directories()

% Adding all required pathways for functions 
pathway_all_functions = "C:\Users\Quantum Dot\Desktop\Bera_Yavuz_GitHub\AttoCube-Project-Stuff";
addpath(pathway_all_functions)


% file that contains all the functions 
funcs = functionsContainer; 
funcs.AddPathFunc("LAB")


% Assigning the default settings to ASI
ASI_Settings.desired_exposure_time = 2;
ASI_Settings.Gain = 300;
ASI_Settings.Gamma = 50; 
ASI_Settings.Brightness = 50; 
UI_Settings = ""; 


[vid_ASI,src_ASI,vid_UI,src_UI,CamInfo] = funcs.ASI_UI_CameraInit(ASI_Settings,UI_Settings);


[UI_Position_Img] = funcs.UI_Snap_Img(vid_UI,src_UI,"Yes",[1 1]); 
UI_Position_Img = flipud(UI_Position_Img); 


% Filtering Settings
scaling = 0.5; 
skyBlue = [0.53, 0.81, 0.92]; 
radii_big_circle = 500; 
center_big_circle = [650,500]; 
sigma_flatfield = 35; 
Salt_pepper_pixel_factor = [20 20]; 
Min_circle_area = 200;
radiusQD = 25;

% LED 
LEDSpotCentroid = funcs.LED_Cooordinate_Identifer("","Read","LAB"); 
LEDSpotCentroid = round(LEDSpotCentroid); 
LEDSpotCentroidX = LEDSpotCentroid(1);
LEDSpotCentroidY = LEDSpotCentroid(2);

% Define the number of parallel lines and their separation
num_lines = 9; % Num lines on each side of the original line
sep_red = 345; % Separation distance red lines
sep_blue = 370; % Seperation distance blue lines

angle = 45; %43.696474776653320; % use the angle found by the algorithm in the XYANC300Axes script 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


% Processing raw photo for QD location data
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = funcs.Finalized_Analyzed_QDBinaryImg_With_Dots(scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,UI_Position_Img);
[allNextPts,allPerpPts] = funcs.Finalized_MainAxes(Img,CopyCentroid,radiusQD);
[x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = funcs.Finalized_GridLines(allNextPts,allPerpPts); 
[VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = funcs.Finalized_VirtualQD(num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled); 
[VirtualQDList_rotated,RealQD_rotated,LEDSpotCentroid_rotated,Table_FullQDList_sorted,Rotated_Table_FullQDList_sorted,rotated_image] = funcs.Rotated_Img_n_Pts(angle,MaskedImage,FullQDList_sorted,VirtualQDList,RealQDCentroids,LEDSpotCentroid);
toc
% PLotting Original Photo
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_title = sprintf("LED spot highlighted On Original Image: %s",UI_Position_Img); 
figure("Name",figure_title,"NumberTitle","off","Color",skyBlue)
imshow(MaskedImage)
hold on; 

% Plot Virtual QD
plot(VirtualQDList(:,1), VirtualQDList(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r',"LineStyle","none");

% Plot Real QD 
plot(centroidx, centroidy, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none");

% Plot LED 
plot(LEDSpotCentroidX,LEDSpotCentroidY,"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o")
hold off; 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Plotting rotated image 
figure_title = sprintf("Virtual and Real QD rotated: %s",UI_Position_Img); 
figure("Name",figure_title,"NumberTitle","off","Color",skyBlue)
imshow(rotated_image);
hold on;

% Plot Virtual Roated QD 
plot(VirtualQDList_rotated(:,1), VirtualQDList_rotated(:,2), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r',"LineStyle","none");

% Plot Real Rotated QD 
plot(RealQD_rotated(:,1), RealQD_rotated(:,2), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none");

% Plot Rotated LED 
plot(LEDSpotCentroid_rotated(1),LEDSpotCentroid_rotated(2),"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o")

title("All Real and Virtual QD Rotated")
hold off; 

