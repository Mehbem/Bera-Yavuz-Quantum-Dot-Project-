function Find_Highest_Photon_Spot_UpHill(obj, step_size_x, step_size_y, Photon_count_init, ANC300, Frequency)
            % Hill-climbing algorithm 
            
            % Get the original photon count
            best_count = py.ID900_Func.query_photon_counter(Photon_count_init);
            
            % Movement commands (forward and backward) - using fprintf with anc400
            move_fwd_x = sprintf("stepu 1 %d", step_size_x);move_back_x = sprintf("stepd 1 %d", step_size_x); 
            move_fwd_y = sprintf("stepu 2 %d", step_size_y);move_back_y = sprintf("stepd 2 %d", step_size_y); 
            move_x = {@() fprintf(ANC300, move_fwd_x), ...
                      @() fprintf(ANC300, move_back_x)}; % X+ (right), X- (left)
            move_y = {@() fprintf(ANC300, move_fwd_y), ...
                      @() fprintf(ANC300, move_back_y)}; % Y+ (up), Y- (down)

            reverse_x = {move_x{2}, move_x{1}}; % Reverse for X
            reverse_y = {move_y{2}, move_y{1}}; % Reverse for Y
            
            % Direction order: [X+, X-, Y+, Y-]
            directions = {move_x{1}, move_x{2}, move_y{1}, move_y{2}};
            reverse_directions = {reverse_x{1}, reverse_x{2}, reverse_y{1}, reverse_y{2}};
        
            % Position tracking
            visited_positions = []; % Store past positions and their photon counts
            step_types = [step_size_x, step_size_y]; % Step sizes
        
            while true
                best_direction = -1; % Index of the best direction
                max_count = best_count; % Best photon count found
                
                % Try moving in all 4 directions
                for check = 1:4
                    % Move in the specified direction
                    directions{check}();
                    if check <= 2
                        step_num = step_types(1); % X-movement
                    else
                        step_num = step_types(2); % Y-movement
                    end
                    StepQueue(obj, step_num, Frequency);
        
                    % Get the new photon count
                    new_count = py.ID900_Func.query_photon_counter(Photon_count_init);
                    
                    % Check if it's the highest count so far
                    if new_count > max_count
                        max_count = new_count;
                        best_direction = check;
                    end
                    
                    % Store position data
                    visited_positions = [visited_positions; check, new_count]; 
                    
                    % Move back to the original position
                    reverse_directions{check}();
                end
                
                % Stop if no improvement is found
                if best_direction == -1
                    break;
                end
                
                % Move to the new best position
                directions{best_direction}();
                StepQueue(obj, step_types(mod(best_direction,2)+1), Frequency);
                
                % Update the best count
                best_count = max_count;
            end
        end

function Find_Highest_Photon_Spot_Stoichastic(obj, step_size_x, step_size_y, Photon_count_init, ANC300, Frequency)
            % Get the initial photon count
            best_count = py.ID900_Func.query_photon_counter(Photon_count_init);
            
            % Define step sizes
            step_sizes = [step_size_x, step_size_y];
            
            % Search directions (X+, X-, Y+, Y-)
            move_fwd_x = sprintf("stepu 1 %d", step_size_x);move_back_x = sprintf("stepd 1 %d", step_size_x); 
            move_fwd_y = sprintf("stepu 2 %d", step_size_y);move_back_y = sprintf("stepd 2 %d", step_size_y); 
            move_x = {@() fprintf(ANC300, move_fwd_x), ...
                      @() fprintf(ANC300, move_back_x)}; % X+ (right), X- (left)
            move_y = {@() fprintf(ANC300, move_fwd_y), ...
                      @() fprintf(ANC300, move_back_y)}; % Y+ (up), Y- (down)
            directions = {move_x{1}, move_x{2}, move_y{1}, move_y{2}};
            
            % Grid search initialization
            num_samples = 8; % Number of random movements before refining
            best_position = [0, 0]; % Tracking relative movements
            
            % **Step 1: Broad Random Sampling**
            for i = 1:num_samples
                dir_idx = randi(4); % Pick a random direction
                directions{dir_idx}();
                StepQueue(obj,step_size_x, Frequency);
                
                % Get photon count
                new_count = py.ID900_Func.query_photon_counter(Photon_count_init);
                
                % If the photon count improves, update the best position
                if new_count > best_count
                    best_count = new_count;
                    best_position = [dir_idx, new_count]; % Store best direction & count
                end
            end
            
            % **Step 2: Move to Best Found Position**
            if best_position(1) ~= 0
                directions{best_position(1)}(); 
                StepQueue(obj,step_size_x, Frequency);
            end
            
            % **Step 3: Refine with Smaller Steps**
            step_sizes = step_sizes / 2; % Reduce step size for precision
            improvement = true;
            
            while improvement
                improvement = false;
                for check = 1:4
                    directions{check}();
                    StepQueue(obj, step_sizes(mod(check,2)+1), Frequency);
                    
                    new_count = py.ID900_Func.query_photon_counter(Photon_count_init);
                    
                    if new_count > best_count
                        best_count = new_count;
                        improvement = true;
                    else
                        % No improvement, ignore bad movement
                    end
                end
            end
        end

function Find_Highest_Photon_RSRA(obj, step_size_x, step_size_y, Photon_count_init, ANC300, Frequency, max_attempts)
            % Get the initial photon count
            best_count = py.ID900_Func.query_photon_counter(Photon_count_init);
            
            % Define movement functions
            move_fwd_x = sprintf("stepu 1 %d", step_size_x);move_back_x = sprintf("stepd 1 %d", step_size_x); 
            move_fwd_y = sprintf("stepu 2 %d", step_size_y);move_back_y = sprintf("stepd 2 %d", step_size_y); 
            move_x = {@() fprintf(ANC300, move_fwd_x), ...
                      @() fprintf(ANC300, move_back_x)}; % X+ (right), X- (left)
            move_y = {@() fprintf(ANC300, move_fwd_y), ...
                      @() fprintf(ANC300, move_back_y)}; % Y+ (up), Y- (down)
            directions = {move_x{1}, move_x{2}, move_y{1}, move_y{2}}; % List of movements
        
            % Search parameters
            best_position = [0, 0]; % Assume initial position as reference
            step_reduction_factor = 0.7; % Reduce step size in later iterations
            attempts = 0;
            
            while attempts < max_attempts
                % Randomly shuffle movement order
                move_order = randperm(4);  
                found_better = false;
                
                for i = 1:4
                    dir_idx = move_order(i);
                    directions{dir_idx}();
                    StepQueue(obj,step_size_x, Frequency); % Queue step to sync piezo
                    
                    % Get new photon count
                    new_count = py.ID900_Func.query_photon_counter(Photon_count_init);
                    
                    % Check if this is the best count found
                    if new_count > best_count
                        best_count = new_count;
                        best_position = [dir_idx, new_count]; % Store best position info
                        found_better = true;
                    end
                end
                
                % Reduce step size if we are improving
                if found_better
                    step_size_x = round(step_size_x * step_reduction_factor);
                    step_size_y = round(step_size_y * step_reduction_factor);
                else
                    % If no better spot is found, stop searching
                    break;
                end
            end
        end
