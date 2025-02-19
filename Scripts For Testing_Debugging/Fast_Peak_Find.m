% Bera Yavuz
% Object with all required functions
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
% Establishing serialconnetion with ANC300 device
%ANC300 = serialport("COM7",9600); % Change X to the COM port connected

check_x_move_step = 4;
check_y_move_step = 4; 
Frequency = 20; 

            % Init
            Photon_count_init = py.ID900_Func.init_ID900();


            % Get the original photon count
            original_count = py.ID900_Func.query_photon_counter(Photon_count_init); 
            
            % Movement commands (forward and backward) - using fprintf with anc400
            move_x_back_cmd = sprintf("stepd 1 %d", check_x_move_step); move_x_forward_cmd = sprintf("stepu 1 %d", check_x_move_step);
            move_y_back_cmd = sprintf("stepd 2 %d", check_y_move_step); move_y_forward_cmd = sprintf("stepu 2 %d", check_y_move_step);

            move_x = {@() fprintf(ANC300, move_x_back_cmd), ...
                      @() fprintf(ANC300, move_x_forward_cmd)}; % X+ (right), X- (left)
            move_y = {@() fprintf(ANC300, move_y_back_cmd), ...
                      @() fprintf(ANC300, move_y_forward_cmd)}; % Y+ (up), Y- (down)
            
            % Directions: {Move forward, Move backward}
            directions = {move_x{1}, move_x{2}, move_y{1}, move_y{2}};
            reverse_directions = {move_x{2}, move_x{1}, move_y{2}, move_y{1}}; % Reverse each direction
            step_types = [check_x_move_step check_y_move_step]; 
            
            % main while loop that iterates each 4 direction movement
            no_movement = 0; 
            while no_movement ~= 4
            
            % counter to see how many times a movement didn't happen for each direct
            no_movement = 0; 
            
                % for loop used for going in each direction  
                for check = 1:4
                
                    while true
                        % Move in the specified direction
                        directions{check}();
                        if check == 1 || check == 2
                            step_num = step_types(1);
                        else
                            step_num = step_types(2);
                        end
                        MyFuncs.StepQueue(step_num,Frequency);
                
                        
                        % Get the new photon count
                        new_count = py.ID900_Func.query_photon_counter(Photon_count_init); 
                        pause(1)
                        % If photon count increases, update original count and continue moving
                        if new_count > original_count
                            original_count = new_count; 
                        else
                            % Move back to the previous position and break out of loop
                            reverse_directions{check}();
                            no_movement = no_movement+1; 
                            break;
                        end
                    end
                end
            end