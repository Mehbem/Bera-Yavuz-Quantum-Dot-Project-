% Plotting two different axes across each other to find angle 


% First Image 
points_1 = [1618.9 712.62
1427.49 887.507
1239.75 1056.92
1056.55 1231.38];

points_2 = [1372.08 856.075
1184.74 1026.32
1000.54 1193.32
821.378 1367.92];


% Normalize the data (scale between 0 and 1)
all_points = [points_1; points_2];
min_vals = min(all_points);
max_vals = max(all_points);
range_vals = max_vals - min_vals;

points_1_norm = (points_1 - min_vals) ./ range_vals;
points_2_norm = (points_2 - min_vals) ./ range_vals;

% Extract normalized x and y values
x1 = points_1_norm(:,1); y1 = points_1_norm(:,2);
x2 = points_2_norm(:,1); y2 = points_2_norm(:,2);

% Get coefficients of best-fit lines
coeff_1 = polyfit(x1, y1, 1);
coeff_2 = polyfit(x2, y2, 1);

% Generate fitted lines
xFit = linspace(0, 1, 1000);
yFit1 = polyval(coeff_1, xFit);
yFit2 = polyval(coeff_2, xFit);

% Calculate angle between the two lines
slope1 = coeff_1(1);
slope2 = coeff_2(1);
theta_rad = atan(abs((slope2 - slope1) / (1 + slope1 * slope2))); % Angle in radians
theta_deg = rad2deg(theta_rad); % Convert to degrees

% Plot normalized data and best-fit lines
figure;
hold on;
plot(x1, y1, 'bo', 'MarkerSize', 8, 'DisplayName', 'Dataset 1');
plot(xFit, yFit1, 'b-', 'LineWidth', 2, 'DisplayName', 'Best Fit 1');
plot(x2, y2, 'ro', 'MarkerSize', 8, 'DisplayName', 'Dataset 2');
plot(xFit, yFit2, 'r-', 'LineWidth', 2, 'DisplayName', 'Best Fit 2');

% Display angle
text(0.05, 0.9, sprintf('Angle: %.2f°', theta_deg), 'FontSize', 12, 'FontWeight', 'bold');

grid on;
xlabel('Normalized X');
ylabel('Normalized Y');
legend('show');
title('Comparison of Best Fit Lines with Angle Between Them');
hold off;
