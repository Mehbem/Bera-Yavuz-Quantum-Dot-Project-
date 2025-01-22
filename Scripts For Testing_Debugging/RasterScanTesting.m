clear;
clc;
workspace;

% Parameters
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Object with all required functions
my_obj = functionsContainer; 

% Adding needed pathways depending on Device - for lab purposes always use device name LAB 
my_obj.AddPathFunc("LAB")


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
LEDSpotCentroid = [1411	886];
LEDSpotCentroidX = LEDSpotCentroid(1);
LEDSpotCentroidY = LEDSpotCentroid(2);

% Define the number of parallel lines and their separation
num_lines = 9; % Num lines on each side of the original line
sep_red = 345; % Separation distance red lines
sep_blue = 370; % Seperation distance blue lines

angle = 45; %43.696474776653320; % use the angle found by the algorithm in the XYANC300Axes script 

I = "Nanowire_Photo_[26 6].jpg";  
[grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = my_obj.Finalized_Analyzed_QDBinaryImg_With_Dots(scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,I);
[allNextPts,allPerpPts] = my_obj.Finalized_MainAxes(Img,CopyCentroid,radiusQD);
[x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = my_obj.Finalized_GridLines(allNextPts,allPerpPts); 
[VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = my_obj.Finalized_VirtualQD(num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled); 
[VirtualQDList_rotated,RealQD_rotated,LEDSpotCentroid_rotated,Table_FullQDList_sorted,Rotated_Table_FullQDList_sorted,rotated_image] = my_obj.Rotated_Img_n_Pts(angle,MaskedImage,FullQDList_sorted,VirtualQDList,RealQDCentroids,LEDSpotCentroid);
[StartingQD_rotated,ShortestDistance,ID] = my_obj.FindClosestPt(Rotated_Table_FullQDList_sorted.("QD Coordinates"),LEDSpotCentroid_rotated); 
StartingQD = Table_FullQDList_sorted.("QD Coordinates")(ID,:); 


% PLotting Original Photo 
figure_title = sprintf("LED spot highlighted On Original Image: %s",I); 
figure("Name",figure_title,"NumberTitle","off","Color",skyBlue)
imshow(MaskedImage)
hold on; 

% Plot Virtual QD
plot(VirtualQDList(:,1), VirtualQDList(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r',"LineStyle","none");

% Plot Real QD 
plot(centroidx, centroidy, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none");

% Plot LED 
plot(LEDSpotCentroidX,LEDSpotCentroidY,"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o","LineStyle","none")
direction_QD = "bottomleft"; % manually edit this for debugging accordingly 
repeatingNum = 1; 
%StartingQD = [1217.92 1032.41]; 
%[DistanceBetweenPoints,DistanceBetweenPoints_Rotated,NextPt_Coords,Rotated_NextPt,AllPointsOnDiagonal_Original,AllPointsOnDiagonal_Rotated] = my_obj.Finalized_RasterScanPatt(radiusQD,StartingQD,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted,direction_QD,repeatingNum);
[NextPoint,Rotated_NextPoint,XY_Difference,Rotated_XY_Difference] = my_obj.OptimizedRasterScan(StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_QD);

% next point 
plot(NextPoint(1),NextPoint(2),"MarkerSize",8,"MarkerEdgeColor",[0 0 1],"MarkerFaceColor",[0 0 1],"Marker","o","LineStyle","none")

%point to start from 
plot(StartingQD(1),StartingQD(2),"MarkerSize",8,"MarkerEdgeColor",[1 0 1],"MarkerFaceColor",[1 0 1],"Marker","o","LineStyle","none")
legend("VirtualQD","RealQD","LED Point","Next Point", "Starting Point")