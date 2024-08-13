% Load the data
load('All_Pier_Data_30_Min_Bins.mat');

% Apply data quality control
valid_idx = Cz_p <= 0.01; 

valid_idx =  r2_uw_p >= .8 & valid_idx; 

% Filter the data for wind coming 90 degrees from onshore
onshore90_idx = onshore90_p == 1 & valid_idx;
filtered_Wdir_measured_p = Wdir_measured_p(onshore90_idx);
filtered_Sdir_measured = Sdir_measured_p(onshore90_idx);

% Shift wind direction data
shifted_Wdir_measured_p = filtered_Wdir_measured_p;
shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) = shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) + 360;

% Find the mean of the shifted wind direction data
mean_shifted_Wdir = mean(shifted_Wdir_measured_p);

% Determine the shift needed to center data between -90 and 90 degrees
shift_needed = 90 - max(shifted_Wdir_measured_p);
shifted_Wdir_measured_p = shifted_Wdir_measured_p + shift_needed;

% Define the wind direction ranges
alongshore_idx = (shifted_Wdir_measured_p >= -90 & shifted_Wdir_measured_p < -45) | ...
                 (shifted_Wdir_measured_p > 45 & shifted_Wdir_measured_p <= 90);
onshore_idx = shifted_Wdir_measured_p >= -45 & shifted_Wdir_measured_p <= 45;

% Extract stress direction data for each category
alongshore_Wdir = shifted_Wdir_measured_p(alongshore_idx);
alongshore_Sdir = filtered_Sdir_measured(alongshore_idx);
onshore_Wdir = shifted_Wdir_measured_p(onshore_idx);
onshore_Sdir = filtered_Sdir_measured(onshore_idx);

% Bin the wind direction data
bin_size = 5; 
bin_edges = -90:bin_size:90;
bin_centers = bin_edges(1:end-1) + bin_size / 2;

% Calculate the mean stress direction for each bin and count the data points
mean_alongshore_Sdir = zeros(size(bin_centers));
mean_onshore_Sdir = zeros(size(bin_centers));
count_alongshore = zeros(size(bin_centers));
count_onshore = zeros(size(bin_centers));

for i = 1:length(bin_centers)
    % Alongshore
    bin_idx_alongshore = alongshore_Wdir >= bin_edges(i) & alongshore_Wdir < bin_edges(i+1);
    count_alongshore(i) = sum(bin_idx_alongshore);
    if count_alongshore(i) > 0
        mean_alongshore_Sdir(i) = mean(alongshore_Sdir(bin_idx_alongshore));
    else
        mean_alongshore_Sdir(i) = NaN;
    end
    
    % Onshore
    bin_idx_onshore = onshore_Wdir >= bin_edges(i) & onshore_Wdir < bin_edges(i+1);
    count_onshore(i) = sum(bin_idx_onshore);
    if count_onshore(i) > 0
        mean_onshore_Sdir(i) = mean(onshore_Sdir(bin_idx_onshore));
    else
        mean_onshore_Sdir(i) = NaN;
    end
end

% Create the plot
figure;
hold on;

% Plot the mean stress direction for alongshore and onshore
plot(bin_centers, mean_alongshore_Sdir, 'r-o', 'LineWidth', 2, 'DisplayName', 'Alongshore Mean Stress Direction');
plot(bin_centers, mean_onshore_Sdir, 'b-o', 'LineWidth', 2, 'DisplayName', 'Onshore Mean Stress Direction');

% Add grid and labels
title('Mean Stress Direction vs Wind Direction');
xlabel('Wind Direction (degrees)');
ylabel('Mean Stress Direction');
legend show;
grid on;

% Add reference lines and text
y_limits = ylim;
plot([-45 -45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([45 45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([0 0], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

% Add North and South text annotations
x_limits = xlim;
text(x_limits(2), y_limits(1), 'North', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
text(x_limits(1), y_limits(1), 'South', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

hold off;
