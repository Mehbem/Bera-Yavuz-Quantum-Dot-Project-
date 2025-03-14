function Specific_Movement(app,text_type)
% Parameters
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                Frequency = 20; 
                QD_Progress = 0; %variable for tracking current number of QD done
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                
                % reinitilizing spectrometer camera after a given number of
                % QD
                QD_Progress = QD_Progress + 1; 
                if mod(QD_Progress,50) == 0
                    reinitialize_ASI_Camera(app)
                end
                
                
                % Initial Movement Towards Start
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                readQD_position = app.MyFuncs.QD_tracking_N_IdentificationVer2("","","Read","LAB","");
                QD_counter = readQD_position.lastLine; 
                CurrentQD =  QD_counter;
                InquiredQD = app.QD_matrix(QD_Scanned,:); 
                readQD_position = app.MyFuncs.QD_tracking_N_IdentificationVer2("","","Read","LAB","");
                [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking_Matlab(app.ANC300,QD_counter,app.vid_UI,app.src_UI,25);
          
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                
                % calculating number of rows and column alongside direction needed to travel to inquired QD 
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                QD_loca_diff = CurrentQD - InquiredQD; 
                app.MyFuncs.FSS_Folder_Creation(InquiredQD)
                
                if QD_loca_diff(1) > 0 & QD_loca_diff(2) < 0
                
                % direction for each respectiive axis 
                direction_row = "topright";
                direction_column = "bottomright"; 
                
                % Intervals for each respective axis 
                
                row_interval = [-1 0];
                column_interval = [0 1];
                
                elseif QD_loca_diff(1)  > 0 & QD_loca_diff(2) > 0 
                
                % direction for each respectiive axis 
                direction_row = "topright";
                direction_column = "topleft";
                
                % Intervals for each respective axis 
                row_interval = [-1 0];
                column_interval = [0 -1];
                
                elseif QD_loca_diff(1)  < 0 & QD_loca_diff(2) < 0 
                
                % direction for each respectiive axis 
                direction_row = "bottomleft";
                direction_column = "bottomright"; 
                
                % Intervals for each respective axis 
                row_interval = [1 0];
                column_interval = [0 1];
                
                elseif QD_loca_diff(1)  < 0 & QD_loca_diff(2) > 0 
                
                % direction for each respectiive axis     
                direction_row = "bottomleft";
                direction_column = "topleft";
                
                % Intervals for each respective axis 
                row_interval = [1 0];
                column_interval = [0 -1];
                
                elseif QD_loca_diff(1) == 0 | QD_loca_diff(2) == 0
                    if QD_loca_diff(1) < 0 
                    direction_row = "bottomleft";
                    row_interval = [1 0];
                    elseif QD_loca_diff(1) > 0 
                    direction_row = "topright";
                    row_interval = [-1 0];
                    elseif QD_loca_diff(2) < 0 
                    direction_column = "bottomright"; 
                    column_interval = [0 1];
                    elseif QD_loca_diff(2) > 0 
                    direction_column = "topleft"; 
                    column_interval = [0 -1];
                    else
                        AppLoggerText = sprintf("No movement"); 
                        AppLogger_FSS(app,AppLoggerText,'black')
                    end
                
                end
                
                % number of rows and columns to move 
                num_rows = abs(QD_loca_diff(1));
                num_columns = abs(QD_loca_diff(2));
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                
                
                %  Movement Towards QD  
                % --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                dots_to_travel_row = num_rows; 
                dots_to_travel_column = num_columns;
                groups_of_dots_row = 4; 
                groups_of_dots_column = 4; 
    
                num_of_breaks_row = floor(dots_to_travel_row / groups_of_dots_row);  % How many times it travels then locks in-between each QD (along the row) 
                remaining_QD_row = mod(dots_to_travel_row, groups_of_dots_row);    % What remains that is indiviually travelled (along the row) 
                
                num_of_breaks_column = floor(dots_to_travel_column / groups_of_dots_column);  % How many times it travels then locks in-between each QD (along the column) 
                remaining_QD_column = mod(dots_to_travel_column, groups_of_dots_column);    % What remains that is indiviually travelled (along the column) 
                
                Frequency_fast_travel = 30; %whatever frequency we choose for fast travel essentially 
                
                % defining the number of steps
                Y_Step_num = app.data.Y_Fast_Movement_Step; 
                X_Step_num = app.data.X_Fast_Movement_Step;
    
                % additional step number to increase error margin for pause
                error_margin_time = 15; % extra time paused = error_margin_time / frequency
    
                
                % acceptable error margin (pixel distance from dot) for precision locking 
                fast_movement_margin = 80; 
                locking_on_margin = 40;
                final_locking_on_margin = 20; 
    
                % defining the commands depending on the direction being moved 
                row_exist  = exist('direction_row', 'var');
                column_exist = exist('direction_column', 'var');
                
                if row_exist
                if direction_row == "bottomleft" 
                    Y_Serial_Comd = sprintf("stepu 2 %d",Y_Step_num);
                elseif direction_row == "topright"
                    Y_Serial_Comd = sprintf("stepd 2 %d",Y_Step_num); 
                end
                end
                
                if column_exist
                if direction_column == "bottomright" 
                    X_Serial_Comd = sprintf("stepu 1 %d",X_Step_num);
                elseif direction_column == "topleft"
                    X_Serial_Comd = sprintf("stepd 1 %d",X_Step_num); 
                end
                end

                AppLogger_FSS(app,"Starting Movement to ",'black'); 
                 
                % fast tracking along the QDs column
                for column = 1: num_of_breaks_column
                
                    % Moving to the next set of QD along column
                    QD_counter = QD_counter + groups_of_dots_column*column_interval; 
                    % Respective voltage and frequency values (Voltage: 1-20 V Frequency: 1-1000 Hz)
                    fprintf(appANC300,"setv 1 30"); 
                    fprintf(app.ANC300,"setf 1 30"); 
                    pause(0.3)
                    AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,'black');
                    fprintf(app.ANC300,X_Serial_Comd);
                    app.MyFuncs.StepQueue(X_Step_num+error_margin_time,Frequency_fast_travel); 
                
                    % Making sure LED is exactly on the dot 
                    fprintf(app.ANC300,"setv 1 12"); 
                    fprintf(app.ANC300,"setf 1 20");
                    pause(0.3)
                    app.MyFuncs.QD_tracking_N_IdentificationVer2(QD_counter,"Auto","Write","LAB",direction_column); 
                    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking_Matlab(app.ANC300,QD_counter,app.vid_UI,app.src_UI,fast_movement_margin);
                    AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,app.Darkgreen); 
    
                    % Updating all visuals on GUI
                    QD_Location_text = sprintf("[%d %d]",QD_counter); 
                    app.QD_location_SpecificQD.Text = QD_Location_text; 
    
                end
    
                % moving the remaing QD along column
                for column_remaining = 1: remaining_QD_column
                    % find next point along the column 
                    [~,~,~,DistanceBetweenPoints_Rotated] = app.MyFuncs.OptimizedRasterScan(StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_column);
                
                    % Moving to the next QD along column
                    QD_counter = QD_counter + column_interval; 
                    AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,'black');
                    app.MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_column,app.ANC300,Frequency,app.x_factor*3/4,app.y_factor*3/4,app.X_Fact_Back*3/4,app.Y_Fact_Back*3/4)

                    
                    % Making sure LED is exactly on the dot 
                    app.MyFuncs.QD_tracking_N_IdentificationVer2(QD_counter,"Auto","Write","LAB",direction_column); 
                    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking_Matlab(app.ANC300,QD_counter,app.vid_UI,app.src_UI,locking_on_margin); % trying to land on the exact dot 
                    AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,app.Darkgreen); 
    
                    % Updating all visuals on GUI
                    QD_Location_text = sprintf("[%d %d]",QD_counter); 
                    app.QD_location_SpecificQD.Text = QD_Location_text; 
                end
    
    
    
                
                % fast tracking along the QDs row
                for row = 1:num_of_breaks_row
                
                    % Moving to the next set of QD along row
                    QD_counter = QD_counter + groups_of_dots_row*row_interval; 
                    % Respective voltage and frequency values (Voltage: 1-20 V Frequency: 1-1000 Hz)
                    fprintf(app.ANC300,"setv 2 30"); 
                    fprintf(app.ANC300,"setf 2 30"); 
                    pause(0.3)
                    AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,'black');
                    fprintf(app.ANC300,Y_Serial_Comd);
                    app.MyFuncs.StepQueue(Y_Step_num+error_margin_time,Frequency_fast_travel); 
                
                    % Making sure LED is exactly on the dot 
                    fprintf(app.ANC300,"setv 2 12"); 
                    fprintf(app.ANC300,"setf 2 20"); 
                    pause(0.3)
                    app.MyFuncs.QD_tracking_N_IdentificationVer2(QD_counter,"Auto","Write","LAB",direction_row); 
                    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking_Matlab(app.ANC300,QD_counter,app.vid_UI,app.src_UI,fast_movement_margin);
                    
                    % trying to land on the exact dot 
                    AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,app.Darkgreen); 
    
                    % Updating all visuals on GUI
                    QD_Location_text = sprintf("[%d %d]",QD_counter); 
                    app.QD_location_SpecificQD.Text = QD_Location_text; 
                end
    
                % moving the remaing QD along row
                for row_remaining = 1: remaining_QD_row
                    % find next point along the column 
                    [~,~,~,DistanceBetweenPoints_Rotated] = app.MyFuncs.OptimizedRasterScan(StartingQD, Table_FullQDList_sorted, Rotated_Table_FullQDList_sorted,direction_row);
                
                    % Moving to the next QD along column
                    QD_counter = QD_counter + row_interval; 
                    AppLogger_Specific_text = sprintf("Moving to QD [%d %d]",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,'black');
                    app.MyFuncs.Dual_ANC300_Movement(DistanceBetweenPoints_Rotated(1),DistanceBetweenPoints_Rotated(2),direction_row,app.ANC300,Frequency,app.x_factor*3/4,app.y_factor*3/4,app.X_Fact_Back*3/4,app.Y_Fact_Back*3/4)

                    
                    % Making sure LED is exactly on the dot 
                    app.MyFuncs.QD_tracking_N_IdentificationVer2(QD_counter,"Auto","Write","LAB",direction_row); 
                    [StartingQD,~,Rotated_Table_FullQDList_sorted,Table_FullQDList_sorted] = app.MyFuncs.Precision_Locking_Matlab(app.ANC300,QD_counter,app.vid_UI,app.src_UI,locking_on_margin); % trying to land on the exact dot 
                    AppLogger_Specific_text = sprintf("Currently on [%d %d] QD",QD_counter);
                    AppLogger_FSS(app,AppLogger_Specific_text,app.Darkgreen); 
    
                    % Updating all visuals on GUI
                    QD_Location_text = sprintf("[%d %d]",QD_counter); 
                    app.QD_location_SpecificQD.Text = QD_Location_text; 
                end
    

                AppLogger_Specific_text = sprintf("Arrived to inquired QD: [%d %d]",QD_counter); 
                AppLogger_FSS(app,AppLogger_Specific_text,app.Darkgreen); 
end