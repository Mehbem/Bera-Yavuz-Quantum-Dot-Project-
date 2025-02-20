 function Find_High_Photon_fast_algorithm(obj,check_x_move_step,check_y_move_step,Photon_count_init,ANC300,Frequency)
            
            % Get the original photon count
            original_count = py.ID900_Func.query_photon_counter(Photon_count_init); 
            
            % Movement commands (forward and backward) - using fprintf with anc400
            move_fwd_x = sprintf("stepu 1 %d", check_x_move_step);move_back_x = sprintf("stepd 1 %d", check_x_move_step); 
            move_fwd_y = sprintf("stepu 2 %d", check_y_move_step);move_back_y = sprintf("stepd 2 %d", check_y_move_step); 
            move_x = {@() fprintf(ANC300, move_fwd_x), ...
                      @() fprintf(ANC300, move_back_x)}; % X+ (right), X- (left)
            move_y = {@() fprintf(ANC300, move_fwd_y), ...
                      @() fprintf(ANC300, move_back_y)}; % Y+ (up), Y- (down)
            
            % Directions: {Move forward, Move backward}
            directions = {move_x{1}, move_x{2}, move_y{1}, move_y{2}};
            reverse_directions = {move_x{2}, move_x{1}, move_y{2}, move_y{1}}; % Reverse each direction
            step_types = [check_x_move_step check_y_move_step]; 
            
            no_movement = 0; 
            % main while loop that iterates each 4 direction movement
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
                        StepQueue(obj,step_num,Frequency);
                
                        
                        % Get the new photon count
                        pause(2)
                        new_count = py.ID900_Func.query_photon_counter(Photon_count_init); 
                        
                        % If photon count increases, update original count and continue moving
                        if new_count > original_count
                            original_count = new_count; 
                        else
                            % Move back to the previous position and break out of loop
                            reverse_directions{check}();
                            StepQueue(obj,step_num,Frequency);
                            no_movement = no_movement+1; 
                            pause(2)
                            original_count = py.ID900_Func.query_photon_counter(Photon_count_init);
                            break;
                        end
                    end
                end
            end
        end