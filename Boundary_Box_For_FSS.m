% Boundary Box Around QD 
% Parameters
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Object with all required functions
my_obj = functionsContainer; 

% Adding needed pathways depending on Device - for lab purposes always use device name LAB 
my_obj.AddPathFunc("LAB")
Photon_count_init = py.ID900_Func.init_ID900();

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
LEDSpotCentroid = my_obj.LED_Cooordinate_Identifer("","Read","LAB"); 
LEDSpotCentroidX = LEDSpotCentroid(1);
LEDSpotCentroidY = LEDSpotCentroid(2);

% x y factors for small movement 
Read_XY_factor = my_obj.XY_Factor_Identifier("","Read","LAB");
x_factor = Read_XY_factor.X_factor; %(Steps/pixels units)
y_factor = Read_XY_factor.Y_factor; %(Steps/pixels units)

% Defining Frequency
Frequency = 20; 

% Define the number of parallel lines and their separation
num_lines = 9; % Num lines on each side of the original line
sep_red = 345; % Separation distance red lines
sep_blue = 370; % Seperation distance blue lines

angle = 45; %43.696474776653320; % use the angle found by the algorithm in the XYANC300Axes script 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
QD_counter = [67 1];
    
[StartingQD,StartingQD_Rotated,~,~] = my_obj.Precision_Locking_Matlab(ANC300,QD_counter,vid_UI,src_UI,20);
chosen_QD_Coords = StartingQD; 
I = getsnapshot(vid_UI); 
[grayImage,centroid,CopyCentroid,mask,Contrast_adjusted_Img,Img,mask_scaled,MaskedImage] = my_obj.Finalized_Analyzed_QDBinaryImg_With_Dots(scaling,radii_big_circle,center_big_circle,sigma_flatfield,Salt_pepper_pixel_factor,Min_circle_area,I);
[allNextPts,allPerpPts] = my_obj.Finalized_MainAxes(Img,CopyCentroid,radiusQD);
[x_main,y,x_perp,y_perp,b,perp_b,m,perp_m] = my_obj.Finalized_GridLines(allNextPts,allPerpPts); 
[VirtualQDList,FullQDList_sorted,AllPossibleQDList,RealQDCentroids,centroidx,centroidy] = my_obj.Finalized_VirtualQD(num_lines,sep_red,sep_blue,m,perp_m,b,perp_b,CopyCentroid,Img,mask_scaled); 
[VirtualQDList_rotated,RealQD_rotated,LEDSpotCentroid_rotated,Table_FullQDList_sorted,Rotated_Table_FullQDList_sorted,rotated_image] = my_obj.Rotated_Img_n_Pts(angle,MaskedImage,FullQDList_sorted,VirtualQDList,RealQDCentroids,LEDSpotCentroid);


% Coordinates of upper and lower boundary lines 
y_bottom = chosen_QD_Coords(2) + 2.5*radiusQD;
y_top = chosen_QD_Coords(2) - 2.5*radiusQD;

% Coordinates of left and right boundary lines 
x_left = chosen_QD_Coords(1) - 2.5*radiusQD; 
x_right = chosen_QD_Coords(1) + 2.5*radiusQD;


num_points_horz = length(x_left:x_right);
num_points_vert = length(y_top:y_bottom);

y_top_pts = ones(1,num_points_horz) * y_top;
y_btm_pts = ones(1,num_points_horz) * y_bottom;

x_left_pts = ones(1,num_points_vert) * x_left; 
x_right_pts = ones(1,num_points_vert) * x_right; 

horz_line = x_left:x_right;
vert_line = y_top:y_bottom;

% Creating mesh grid for each spot to check and go to 
% Define step size (adjust as needed)
step_size = radiusQD/1.4; % Distance between points

% Generate grid points using meshgrid
[x_grid, y_grid] = meshgrid(x_left:step_size:x_right, y_top:step_size:y_bottom);
[grid_size,~] = size(x_grid); 

% Flatten the grid matrices into vectors for plotting
x_points = x_grid(:);
y_points = y_grid(:);

% Define empty Matrix with coordinates 
raster_order_points = zeros(length(x_points),2);

% Get Unique X and Y vals
Unique_x_vals = unique(x_points);
Unique_y_vals = unique(y_points);

ID = 1; 
for columns = 1:grid_size
    for rows = 1:grid_size
        if mod(columns,2) == 0
            raster_order_points(ID,1) = Unique_x_vals(end+1-rows);
        else
            raster_order_points(ID,1) = Unique_x_vals(rows);
        end
        ID = ID + 1; 
    end
    raster_order_points(columns*grid_size-(grid_size-1):columns*grid_size,2) = Unique_y_vals(columns);
end

% PLotting Original Photo 
figure_title = sprintf("LED spot highlighted On Original Image: %s",I); 
figure("Name",figure_title,"NumberTitle","off","Color",skyBlue)
imshow(MaskedImage)
[height,width] = size(MaskedImage); 
hold on; 

% Plot Virtual QD
plot(VirtualQDList(:,1), VirtualQDList(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r',"LineStyle","none");

% Plot Real QD 
plot(centroidx, centroidy, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none");

% Plot LED 
plot(LEDSpotCentroidX,LEDSpotCentroidY,"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o")

% Plot Chosen QD

plot(chosen_QD_Coords(1),chosen_QD_Coords(2),"MarkerSize",5,"MarkerEdgeColor",[0.2 0 0.7],"MarkerFaceColor",[0.2 0 0.7],"Marker","o")


% plotting boundary lines
plot(horz_line,y_top_pts,"MarkerSize",1,"MarkerEdgeColor",[0.2 0 0.7],"MarkerFaceColor",[0.2 0 0.7],"Marker","o") % top boundary
plot(horz_line,y_btm_pts,"MarkerSize",1,"MarkerEdgeColor",[0.2 0 0.7],"MarkerFaceColor",[0.2 0 0.7],"Marker","o") % bottom boundary
plot(x_left_pts,vert_line,"MarkerSize",1,"MarkerEdgeColor",[0.2 0 0.7],"MarkerFaceColor",[0.2 0 0.7],"Marker","o") % left boundary
plot(x_right_pts,vert_line,"MarkerSize",1,"MarkerEdgeColor",[0.2 0 0.7],"MarkerFaceColor",[0.2 0 0.7],"Marker","o") % right boundary

% Plot Grid Points
plot(x_points, y_points, 'k.', 'MarkerSize', 5); % Black dots for grid points
hold off 

% PLotting Original Photo 
figure_title = sprintf("Raster Scan Patttern: %s",I); 
figure("Name",figure_title,"NumberTitle","off","Color",skyBlue)
imshow(MaskedImage)
[height,width] = size(MaskedImage); 
hold on; 

% Plot the raster scan
plot(raster_order_points(:,1), raster_order_points(:,2), '-o','MarkerSize', 5,"MarkerFaceColor",'green');
% Plot Virtual QD
plot(VirtualQDList(:,1), VirtualQDList(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r',"LineStyle","none");
% Plot Real QD 
plot(centroidx, centroidy, 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none");
% Plot LED 
plot(LEDSpotCentroidX,LEDSpotCentroidY,"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o")
% Plot Chosen QD
plot(chosen_QD_Coords(1),chosen_QD_Coords(2),"MarkerSize",5,"MarkerEdgeColor",[0.2 0 0.7],"MarkerFaceColor",[0.2 0 0.7],"Marker","o")

hold off; 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Rotating points
theta = -angle * pi / 180; % Convert angle to radians
center_pt = [size(MaskedImage, 2)/2, size(MaskedImage, 1)/2]; % Image center

% Virtual points rotated
raster_order_points_X_rotated = center_pt(1) + (raster_order_points(:,1) - center_pt(1)) * cos(theta) - (raster_order_points(:,2) - center_pt(2)) * sin(theta);
raster_order_points_Y_rotated = center_pt(2) + (raster_order_points(:,1) - center_pt(1)) * sin(theta) + (raster_order_points(:,2) - center_pt(2)) * cos(theta);
raster_order_points_rotated = horzcat(raster_order_points_X_rotated,raster_order_points_Y_rotated); 


% Plotting rotated image 
figure_title = sprintf("Virtual and Real QD rotated: %s",I); 
figure("Name",figure_title,"NumberTitle","off","Color",skyBlue)
imshow(rotated_image);
hold on;
plot(VirtualQDList_rotated(:,1), VirtualQDList_rotated(:,2), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r',"LineStyle","none");
plot(RealQD_rotated(:,1), RealQD_rotated(:,2), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g',"LineStyle","none");
plot(LEDSpotCentroid_rotated(1),LEDSpotCentroid_rotated(2),"MarkerSize",5,"MarkerEdgeColor",[1 0.5 0],"MarkerFaceColor",[1 0.5 0],"Marker","o")
plot(raster_order_points_X_rotated,raster_order_points_Y_rotated,'k.', 'MarkerSize', 5)



title("All Real and Virtual QD Rotated")



hold off; 
