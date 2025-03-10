% random walker algorithm 

% Random walker parameters
num_steps = 200; % Number of steps for the random walk
start_position = [0, 0]; % Start from an arbitrary point away from the center
positions = start_position; % Initialize positions array
Photon_count_init = py.ID900_Func.init_ID900();
ANC300 = serialport("COM7",9600); 
MyFuncs = functionsContainer;
old_photon_count = 0; 

% intilizing the highest peak counter 
highest_count = 0; 

% checking basis
check_after = 50; 
breaking_cond = false; 

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
    [x_pos_new, y_pos_new] = MyFuncs.new_stochastic_movement(x_pos, y_pos, dir_x, dir_y, ANC300,20);

    % Get the change in photon count 
    photon_counts = py.ID900_Func.query_photon_counter(Photon_count_init);
    
    % stop stepping since highest has been achieved 
    if photon_counts > highest_count*0.98 & breaking_cond
        break
    end

    % If intensity decreases, retrace step
    if photon_counts < old_photon_count
        [x_pos_new, y_pos_new] = MyFuncs.new_stochastic_movement(x_pos_new, y_pos_new, -dir_x, -dir_y,ANC300,20);
    end

    % If intensity increases or stays the same, accept the new position
    if photon_counts > highest_count
        highest_count = photon_counts; % Append new position
        prev_highest_count = highest_count; 
    end

    if mod(step,check_after) == 0
        if highest_count < 1.1 * prev_highest_count
            breaking_cond = true; 
        end 
    end

    old_photon_count = photon_counts; 
end

