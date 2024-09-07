close all
clear
clc

% Define the corners and boundaries of the geometry
x_left = [0, 0]; % X-coordinates for the left side
y_left = [0, 300]; % Y-coordinates for the left side

x_right = [300, 300]; % X-coordinates for the right side
y_right = [0, 200]; % Y-coordinates for the right side

x_bottom = [0, 300]; % X-coordinates for the bottom side
y_bottom = [0, 0]; % Y-coordinates for the bottom side

% Define the top side using the sine function
x_top = linspace(0, 300, 100); % X-coordinates from 0 to 300
y_top = 50 * cos(pi * x_top / 100) + 250; % Y-coordinates based on the function

% Combine all boundary points in the correct order
boundary_x = [x_left(2), x_top, x_right(2), x_right(1), x_left(1)]; % X-coordinates of the entire boundary
boundary_y = [y_left(2), y_top, y_right(2), y_right(1), y_left(1)]; % Y-coordinates of the entire boundary

%%
% Load the colormap from the MAT file
load('brocO.mat');  % Ensure 'brocO' is the correct variable name in the MAT file

positive_colors =flip( [
    brocO(76,:);  
    brocO(106,:)  
]);

negative_colors = flip([
    brocO(176,:); 
    brocO(186,:); 
    brocO(196,:); 
    brocO(206,:); 
    brocO(216,:); 
    brocO(226,:)  
]);

% Combine the selected colors into one colormap
selected_colors = [negative_colors; positive_colors];
% Load the CSV file while preserving the original column headers
data = readtable('ZeroStress.csv', 'VariableNamingRule', 'preserve');

% Extract the necessary columns
x = data.('X'); % X column for x-coordinates
y = data.('Y'); % Y column for y-coordinates
values = data.('S-S11'); % S-S11 column for the values to contour

% Calculate the min and max for x and y
min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);

% Normalize x and y to the range [0, 300]
x = 300 * (x - min_x) / (max_x - min_x);
y = 300 * (y - min_y) / (max_y - min_y);

% Create a grid of points where you want to evaluate the contour
[xq, yq] = meshgrid(linspace(0, 300, 300), linspace(0,300, 300));

% Interpolate the values on this grid
vq = griddata(x, y, values, xq, yq, 'cubic');

% Create a logical mask for the area inside the boundary
mask = inpolygon(xq, yq, boundary_x, boundary_y);

% Apply the mask to the contour data
vq(~mask) = NaN; % Set the values outside the boundary to NaN


% Plot the masked contour map
figure;

% Plot the original contour
contourf(xq, yq, vq, [-15 -10 -5 0 1 2 3 4 5], 'LineColor', 'none'); 
hold on;

% Mirror the contour data along the x=0 axis and plot it next to the original
contourf(-xq, yq, vq, [-15 -10 -5 0 1 2 3 4 5], 'LineColor', 'none');

% Apply the custom colormap
colormap((selected_colors));  
colorbar; % Add a color bar to indicate the values
clim([-15 5]); % Set the color axis limits

% Plot the original geometry boundary on top of the contour
plot(boundary_x, boundary_y, 'k-', 'LineWidth', 1); 

% Mirror the boundary data and plot it
plot(-boundary_x, boundary_y, 'k-', 'LineWidth', 1); 

hold off;

% Set the axes limits and remove grid
axis([-300 300 0 300]);
axis equal;
%set(gca, 'XColor', 'none', 'YColor', 'none'); % Remove X and Y axis lines
grid off; % Remove the grid
% Save the figure as a vector-based EPS file
print('zero0_S11', '-depsc', '-r300','-vector');
%%
% Load the CSV file while preserving the original column headers
data = readtable('Compre0.csv', 'VariableNamingRule', 'preserve');

% Extract the necessary columns
x = data.('X')*5; % X column for x-coordinates
y = data.('Y')*5; % Y column for y-coordinates
values = data.('S-S11'); % S-S11 column for the values to contour

% Calculate the min and max for x and y
min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);

% Normalize x and y to the range [0, 300]
x = 300 * (x - min_x) / (max_x - min_x);
y = 300 * (y - min_y) / (max_y - min_y);

% Create a grid of points where you want to evaluate the contour
[xq, yq] = meshgrid(linspace(0, 300, 300), linspace(0,300, 300));

% Interpolate the values on this grid
vq = griddata(x, y, values, xq, yq, 'cubic');

% Create a logical mask for the area inside the boundary
mask = inpolygon(xq, yq, boundary_x, boundary_y);

% Apply the mask to the contour data
vq(~mask) = NaN; % Set the values outside the boundary to NaN


% Plot the masked contour map
figure;

% Plot the original contour
contourf(xq, yq, vq, [-15 -10 -5 0 1 2 3 4 5], 'LineColor', 'none'); 
hold on;

% Mirror the contour data along the x=0 axis and plot it next to the original
contourf(-xq, yq, vq, [-15 -10 -5 0 1 2 3 4 5], 'LineColor', 'none');

% Apply the custom colormap
colormap((selected_colors));  
colorbar; % Add a color bar to indicate the values
clim([-15 5]); % Set the color axis limits

% Plot the original geometry boundary on top of the contour
plot(boundary_x, boundary_y, 'k-', 'LineWidth', 1); 

% Mirror the boundary data and plot it
plot(-boundary_x, boundary_y, 'k-', 'LineWidth', 1); 

hold off;

% Set the axes limits and remove grid
axis([-300 300 0 300]);
axis equal;
%set(gca, 'XColor', 'none', 'YColor', 'none'); % Remove X and Y axis lines
grid off; % Remove the grid
legend('off');
% Save the figure as a vector-based EPS file
print('compre0_S11', '-depsc', '-r300','-vector');
%%
% Load the CSV file while preserving the original column headers
data = readtable('tensile0.csv', 'VariableNamingRule', 'preserve');

% Extract the necessary columns
x = data.('X')*5; % X column for x-coordinates
y = data.('Y')*5; % Y column for y-coordinates
values = data.('S-S11'); % S-S11 column for the values to contour

% Calculate the min and max for x and y
min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);

% Normalize x and y to the range [0, 300]
x = 300 * (x - min_x) / (max_x - min_x);
y = 300 * (y - min_y) / (max_y - min_y);

% Create a grid of points where you want to evaluate the contour
[xq, yq] = meshgrid(linspace(0, 300, 300), linspace(0,300, 300));

% Interpolate the values on this grid
vq = griddata(x, y, values, xq, yq, 'cubic');

% Create a logical mask for the area inside the boundary
mask = inpolygon(xq, yq, boundary_x, boundary_y);

% Apply the mask to the contour data
vq(~mask) = NaN; % Set the values outside the boundary to NaN


% Plot the masked contour map
figure;

% Plot the original contour
contourf(xq, yq, vq, [-15 -10 -5 0 1 2 3 4 5], 'LineColor', 'none'); 
hold on;

% Mirror the contour data along the x=0 axis and plot it next to the original
contourf(-xq, yq, vq, [-15 -10 -5 0 1 2 3 4 5], 'LineColor', 'none');

% Apply the custom colormap
colormap((selected_colors));  
colorbar; % Add a color bar to indicate the values
clim([-15 5]); % Set the color axis limits

% Plot the original geometry boundary on top of the contour
plot(boundary_x, boundary_y, 'k-', 'LineWidth', 1); 

% Mirror the boundary data and plot it
plot(-boundary_x, boundary_y, 'k-', 'LineWidth', 1); 

hold off;

% Set the axes limits and remove grid
axis([-300 300 0 300]);
axis equal;
%set(gca, 'XColor', 'none', 'YColor', 'none'); % Remove X and Y axis lines
grid off; % Remove the grid
legend('off');
% Save the figure as a vector-based EPS file
print('tensile0_S11', '-depsc', '-r300','-vector');
