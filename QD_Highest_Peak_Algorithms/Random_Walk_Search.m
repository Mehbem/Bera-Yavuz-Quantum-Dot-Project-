
function [x_pos_new, y_pos_new] = new_stochastic_movement(obj,x_pos, y_pos, dir_x, dir_y, ANC300, Frequency)
            
            % Calculate the random perturbation (random step size)
            step_size_x = randi([3 6],1); % Random step size between 3 and 6
            step_size_y = randi([3 6],1); % Random step size between 3 and 6
        
        
            % step direction command 
            if dir_x == 1
            stochastic_dir_x = sprintf("stepu 1 %d",step_size_x);
            x_factor = 1; 
            else
            stochastic_dir_x = sprintf("stepd 1 %d",step_size_x);
            x_factor = -1; 
            end
            
            if dir_y == 1
            stochastic_dir_y = sprintf("stepu 2 %d",step_size_y);
            y_factor = 1;
            else
            stochastic_dir_y = sprintf("stepd 2 %d",step_size_y);
            y_factor = -1; 
            end

            % Update position
            x_pos_new = x_pos + x_factor*step_size_x;
            y_pos_new = y_pos + y_factor*step_size_y;
        
            % applying actual movement 
            fprintf(ANC300,stochastic_dir_x)
            StepQueue(obj,step_size_x,Frequency);
        
            fprintf(ANC300,stochastic_dir_y)
            StepQueue(obj,step_size_y,Frequency);
            
    
end

       

% Random walker parameters
num_steps = 200; % Number of steps for the random walk
start_position = [0, 0]; % Start from an arbitrary point away from the center
positions = start_position; % Initialize positions array
Photon_count_init = py.ID900_Func.init_ID900();

% Simulate the biased random walk
for step = 1:num_steps
    fprintf("Step %d/%d\n",step,num_steps)
    % Get current position
    current_pos = positions(end, :);
    x_pos = current_pos(1);
    y_pos = current_pos(2);

    % Randomly choose direction for x and y
    dir_x = 2*randi([0,1])-1; 
    dir_y = 2*randi([0,1])-1; 

    % Calculate new position
    [x_pos_new, y_pos_new] = new_stochastic_movement(x_pos, y_pos, dir_x, dir_y, ANC300,20);

    % Get the change in photon count 
    photon_counts = py.ID900_Func.query_photon_counter(Photon_count_init);

    % If intensity decreases, retrace step
    if photon_counts < 0
        [x_pos_new, y_pos_new] = new_stochastic_movement(x_pos_new, y_pos_new, -dir_x, -dir_x,ANC300,20);
    end

    % If intensity increases or stays the same, accept the new position
    if photon_counts >= 0
        positions = [positions; x_pos_new, y_pos_new]; % Append new position
    end
end
