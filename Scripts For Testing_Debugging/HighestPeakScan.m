% Bera Yavuz 

% Object with all required functions
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
% Establishing serialconnetion with ANC300 device
ANC300 = serialport("COM7",9600); % Change X to the COM port connected
ell_motor = serialport("COM10",9600); % Establishing a serialconnection with the HWP motor

fprintf(ANC300,"setm 3 gnd"); 
fprintf(ANC300,"setv 1 12"); fprintf(ANC300,"setf 1 20"); fprintf(ANC300,"setm 1 stp"); 
fprintf(ANC300,"setv 2 12"); fprintf(ANC300,"setf 2 20"); fprintf(ANC300,"setm 2 stp");  
Frequency = 20; 

% Initilize photon detector 
tc = py.ID900_Func.init_ID900();

% Defining different movement options 
X_Step_num_big = 35; X_Step_num_small = 5;
Y_Step_num_big = 30; Y_Step_num_small = 5;
loop_interval_X = X_Step_num_big/X_Step_num_small;
loop_interval_Y = Y_Step_num_big/Y_Step_num_small; 



% Defining smallest stepping Movemment
Y_movement_seria_comd_forward = sprintf("stepu 2 %d",Y_Step_num_small);
X_movement_seria_comd_forward = sprintf("stepu 1 %d",X_Step_num_small);

Y_movement_seria_comd_back = sprintf("stepd 2 %d",Y_Step_num_small);
X_movement_seria_comd_back = sprintf("stepd 1 %d",X_Step_num_small+1);


% Defining an empty 0 0 matrix for tracking location 
relative_loc = [0 0];
max_photon_count = 0; 
percentage_complete = 0;
QD_counter = [1 1];
Spectrometer_Gratting = 1800; 
total_amount = loop_interval_Y*4*loop_interval_X+loop_interval_Y*2; 
break_all = false; 

MyFuncs.Precision_Locking_Matlab(ANC300,QD_counter,vid_UI,src_UI,20)

% get initial photon count 
init_photon_count = py.ID900_Func.query_photon_counter(tc); 

 
% initial data table 
data_table = table(relative_loc,init_photon_count,'VariableNames', {'Positions', 'PeakCounts'}); 


for num_scans = 1:2 
    %Defining starting position Movement
    Y_move_forward = sprintf("stepd 2 %d",Y_Step_num_small);
    for i = 1:loop_interval_Y
        fprintf(ANC300,Y_move_forward);
        MyFuncs.StepQueue(Y_Step_num_small,Frequency)
        pause(0.1)
    end
    
    X_move_forward = sprintf("stepd 1 %d",X_Step_num_small);
    for i = 1:loop_interval_X
        fprintf(ANC300,X_move_forward);
        MyFuncs.StepQueue(X_Step_num_small,Frequency);
        pause(0.1)
    end

    for y_move = 1:loop_interval_Y*2
        if y_move ~= 1
            relative_loc = relative_loc + [0 1]; 
            fprintf(ANC300, Y_movement_seria_comd_forward);
            MyFuncs.StepQueue(X_Step_num_small, Frequency);
            pause(0.1)
    
            % Read photon count
            photon_count = py.ID900_Func.query_photon_counter(tc); 
    
            % Update table
            new_row = {relative_loc, photon_count};
            data_table = [data_table; new_row]; 
            percentage_complete = percentage_complete + 1;
    
            % Update Message
            fprintf("Completed %d/%d",percentage_complete,total_amount)
            if photon_count > 0.8*max_photon_count & max_photon_count ~= 0 
                break_all = true; 
                break
            end
        end
    
        for x_move = 1:loop_interval_X*2
            if mod(y_move, 2) == 0
                relative_loc = relative_loc + [-1 0];
    
                fprintf(ANC300, X_movement_seria_comd_back);
                MyFuncs.StepQueue(X_Step_num_small, Frequency);
                pause(0.1)
            else
                relative_loc = relative_loc + [+1 0];
    
                fprintf(ANC300, X_movement_seria_comd_forward);
                MyFuncs.StepQueue(X_Step_num_small, Frequency);
                pause(0.1)
            end
    
            % Read photon count
            photon_count = py.ID900_Func.query_photon_counter(tc); 
    
            % Update table
            new_row = {relative_loc, photon_count};
            data_table = [data_table; new_row]; 
            percentage_complete = percentage_complete + 1; 
    
             % Update Message
            fprintf("Completed %d/%d",percentage_complete,total_amount)
            if photon_count > 0.8*max_photon_count & max_photon_count ~= 0 
                break_all = true; 
                break
            end
        end
        if break_all == true
            break
        end
    end
    if break_all == true
        break
    end
    % find maximum photon count from first scan 
    max_photon_count = max(data_table.PeakCounts); 
    original_data_table = data_table; 
    if num_scans < 2

        MyFuncs.Precision_Locking_Matlab(ANC300,QD_counter,vid_UI,src_UI,20); % trying to land on the exact dot
        pause(0.2)
    end
end


            % defining angle of rotation 
            angle = 0:5:360;

            % for loop for rotating motor 
            for rot_count = 1:length(angle)

            % Defining the movement serial code for the rotation 
            angle_hxd = dec2hex(floor(mod(angle(rot_count),360)*39822/100), 8);
            input_str = "1ma" + angle_hxd; % "2" before ma is to be get from the ELLO software from thorlabs.
            
            % Commiting Command for movement 
            fprintf(ell_motor, input_str);
            
            
            Current_angle = sprintf('Angle: %d \n', angle(rot_count));
            fprintf(Current_angle)
            pause(1)
            
            % taking a snap at every respective angle 
            file_name = sprintf('FSS %d',angle(rot_count)); 
            MyFuncs.ASI_Snap_Img(vid_ASI,src_ASI,"Spectrometer","No",Spectrometer_Gratting,QD_counter,file_name)
            end





% Extract data from table
X_Coords = original_data_table{:,'Positions'}(:,1);
Y_Coords = original_data_table{:,'Positions'}(:,2);
Photon_count_vals = original_data_table{:,'PeakCounts'};


% Extract unique X and Y values (assuming evenly spaced grid)
xVals = unique(X_Coords);
yVals = unique(Y_Coords); 

% Determine grid size
numCols = length(xVals);
numRows = length(yVals);

% Create an empty grid for intensity values
Z = nan(numRows, numCols);

% Fill the grid with intensity values based on coordinates
for i = 1:size(X_Coords,1)
    xIdx = find(xVals == X_Coords(i));
    yIdx = find(yVals == Y_Coords(i));
    Z(yIdx, xIdx) = Photon_count_vals(i);
end

% Flip the Y-axis so (0,0) is at the top-left
Z = flipud(Z);

% Plot the data using imagesc
figure;
imagesc(xVals, yVals, Z);
colormap(jet); % Change colormap if needed
colorbar;
axis equal tight;
ax = gca;
ax.YTickLabel = flip(ax.YTickLabel); % Reverse Y-axis labels
ax.XAxisLocation = 'top'; % Moves X-axis labels to the top  
set(gca, 'YDir', 'Normal'); % Ensure top-left origin
xlabel('X');
ylabel('Y');

title('Raster Scan Photon Count');
            file_name = sprintf("C:\\Users\\Quantum Dot\\Desktop\\Bera Yavuz - ANC300 Movement and Images\\FSS_Photon_Count_Scan\\Photon_Count_Scan__[%d %d]",QD_counter); 
            saveas(gca,file_name)

