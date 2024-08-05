% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Convert date column to MATLAB datenum
date_nums = datenum(data.Date);

% Helper function to convert angles to unit vectors
angles_to_vectors = @(angles) [cosd(angles), sind(angles)];

% Helper function to convert unit vectors to angles
vectors_to_angles = @(vectors) atan2d(vectors(:,2), vectors(:,1));

% Convert wave and current direction data for awc6 to unit vectors
waveMeanDirection_awc6 = str2double(data.('waveMeanDirection_awc6'));
currentDirection_awc6 = str2double(data.('currentDirection_awc6'));

wave_vectors_awc6 = angles_to_vectors(waveMeanDirection_awc6);
current_vectors_awc6 = angles_to_vectors(currentDirection_awc6);

% Interpolate missing data for awc6 using unit vectors
wave_time_idx_awc6 = find(~isnan(waveMeanDirection_awc6) & waveMeanDirection_awc6 ~= 9999);
current_time_idx_awc6 = find(~isnan(currentDirection_awc6) & currentDirection_awc6 ~= 9999);

wave_vectors_awc6_interp = interp1(wave_time_idx_awc6, wave_vectors_awc6(wave_time_idx_awc6,:), 1:length(waveMeanDirection_awc6), 'pchip', 'extrap');
current_vectors_awc6_interp = interp1(current_time_idx_awc6, current_vectors_awc6(current_time_idx_awc6,:), 1:length(currentDirection_awc6), 'pchip', 'extrap');

% Convert interpolated unit vectors back to angles
waveMeanDirection_awc6_interp = vectors_to_angles(wave_vectors_awc6_interp);
currentDirection_awc6_interp = vectors_to_angles(current_vectors_awc6_interp);

% Repeat the process for awc4.5
waveMeanDirection_awc45 = str2double(data.('waveMeanDirection_awc4.5'));
currentDirection_awc45 = str2double(data.('currentDirection_awc4.5'));

wave_vectors_awc45 = angles_to_vectors(waveMeanDirection_awc45);
current_vectors_awc45 = angles_to_vectors(currentDirection_awc45);

wave_time_idx_awc45 = find(~isnan(waveMeanDirection_awc45) & waveMeanDirection_awc45 ~= 9999);
current_time_idx_awc45 = find(~isnan(currentDirection_awc45) & currentDirection_awc45 ~= 9999);

wave_vectors_awc45_interp = interp1(wave_time_idx_awc45, wave_vectors_awc45(wave_time_idx_awc45,:), 1:length(waveMeanDirection_awc45), 'pchip', 'extrap');
current_vectors_awc45_interp = interp1(current_time_idx_awc45, current_vectors_awc45(current_time_idx_awc45,:), 1:length(currentDirection_awc45), 'pchip', 'extrap');

waveMeanDirection_awc45_interp = vectors_to_angles(wave_vectors_awc45_interp);
currentDirection_awc45_interp = vectors_to_angles(current_vectors_awc45_interp);

% Initialize combined data arrays
waveMeanDirection_combined = waveMeanDirection_awc6_interp;
currentDirection_combined = currentDirection_awc6_interp;

% Replace missing values in awc6 with values from awc4.5 for specified dates
dates_awc6_missing = date_nums >= datenum(2022, 4, 1) & date_nums <= datenum(2022, 7, 31);
waveMeanDirection_combined(dates_awc6_missing) = waveMeanDirection_awc45_interp(dates_awc6_missing);
currentDirection_combined(dates_awc6_missing) = currentDirection_awc45_interp(dates_awc6_missing);

% Apply data quality control
valid_idx = ~strcmp(data.('Wdir.range'), 'Offshore') & data.('Wdir.deg') ~= 9999;
valid_idx = valid_idx & data.Cdz <= 0.01;
valid_idx = valid_idx & ~isnan(data.Cdz) & ~isnan(data.('Wdir.deg')) & ~isnan(data.('Sdir.measured')) & ~isnan(data.Uz);
valid_idx = valid_idx & data.r2_uw >= 0.9;

% Filtered data
filtered_waveMeanDirection = waveMeanDirection_combined(valid_idx);
filtered_currentDirection = currentDirection_combined(valid_idx);
filtered_windDirection = data.('Wdir.deg')(valid_idx);
filtered_stressDirection = data.('Sdir.measured')(valid_idx);
filtered_dates = date_nums(valid_idx);
filtered_wind_speed = data.Uz(valid_idx);

% Interpolate wave height and current speed
wave_height_awc6 = str2double(data.waveHs_awc6);
current_speed_awc6 = str2double(data.currentSpeed_awc6);

wave_height_awc6_interp = interp1(wave_time_idx_awc6, wave_height_awc6(wave_time_idx_awc6), 1:length(wave_height_awc6), 'pchip', 'extrap');
current_speed_awc6_interp = interp1(current_time_idx_awc6, current_speed_awc6(current_time_idx_awc6), 1:length(current_speed_awc6), 'pchip', 'extrap');

filtered_wave_height = wave_height_awc6_interp(valid_idx);
filtered_current_speed = current_speed_awc6_interp(valid_idx);

% Calculate wave speed
g = 9.81; % gravity in m/s^2
wave_period = 30 * 60; % 30 minutes in seconds
filtered_wave_speed = sqrt(g * filtered_wave_height / (2 * pi / wave_period));

% Adjust wave and current directions to the shoreward reference frame
coastline_direction_deg = 72; % Coastline direction in degrees

adjust_direction = @(direction) arrayfun(@(x) ...
    (x >= 0 && x <= 72) * (x - coastline_direction_deg) + ...
    (x > 72 && x <= 360) * (x - coastline_direction_deg) + ...
    (x < 0) * (360 + x - coastline_direction_deg), direction);

filtered_waveMeanDirection = adjust_direction(filtered_waveMeanDirection);
filtered_currentDirection = adjust_direction(filtered_currentDirection);

% Ensure directions are within the range [-180, 180)
normalize_direction = @(x) mod(x + 180, 360) - 180;
filtered_waveMeanDirection = normalize_direction(filtered_waveMeanDirection);
filtered_currentDirection = normalize_direction(filtered_currentDirection);

% Flip current and wave direction so that north is negative and south is positive
filtered_currentDirection = -filtered_currentDirection;
filtered_waveMeanDirection = -filtered_waveMeanDirection;

% Ensure directions are within the range [-180, 180) again after flipping
filtered_waveMeanDirection = normalize_direction(filtered_waveMeanDirection);
filtered_currentDirection = normalize_direction(filtered_currentDirection);

% Transpose the row vectors to column vectors
filtered_wave_height = filtered_wave_height';
filtered_current_speed = filtered_current_speed';
filtered_wave_speed = filtered_wave_speed';

% Combine data into one matrix and remove rows with NaN values
combined_data = [filtered_dates, filtered_windDirection, filtered_waveMeanDirection, filtered_currentDirection, filtered_stressDirection, filtered_wind_speed, filtered_wave_height, filtered_current_speed, filtered_wave_speed];
combined_data = combined_data(~any(isnan(combined_data), 2), :);

% Save the combined data to a new CSV file
output_table = array2table(combined_data, 'VariableNames', {'Date', 'WindDirection', 'WaveMeanDirection', 'CurrentDirection', 'StressDirection', 'WindSpeed', 'WaveHeight', 'CurrentSpeed', 'WaveSpeed'});
writetable(output_table, 'Corrected_Combined_Pier_Data_BE.csv');

% Determine individual radius limits for each plot
max_wind_speed = max(filtered_wind_speed);
max_wave_speed = max(filtered_wave_speed);
max_current_speed = max(filtered_current_speed);

%% THREE GRAPHS
% Extract necessary columns
filtered_Wdir_measured = combined_data(:, 2);
filtered_Sdir_measured = combined_data(:, 5);
adjusted_wave_direction = combined_data(:, 3);
adjusted_current_direction = combined_data(:, 4);

% Define bins for wind direction
bin_size = 15;
bin_edges = -90:bin_size:90;
bin_centers = bin_edges(1:end-1) + bin_size / 2;
binned_95th_percentile = zeros(size(bin_centers));
binned_5th_percentile = zeros(size(bin_centers));
binned_mean = zeros(size(bin_centers));
binned_sem = zeros(size(bin_centers)); % Standard error of the mean

for i = 1:length(bin_centers)
    bin_idx = filtered_Wdir_measured >= bin_edges(i) & filtered_Wdir_measured < bin_edges(i+1);
    bin_data = filtered_Sdir_measured(bin_idx);

    if ~isempty(bin_data)
        binned_95th_percentile(i) = prctile(bin_data, 90);
        binned_5th_percentile(i) = prctile(bin_data, 10);
        binned_mean(i) = mean(bin_data);
        binned_sem(i) = std(bin_data) / sqrt(length(bin_data)); % Calculate the standard error of the mean
    end
end

% Create figure with 3 subplots
figure;

% Subplot 1: Filled plot showing 90% data range and mean
subplot(3, 1, 1);
hold on;
fill([bin_centers, fliplr(bin_centers)], ...
     [binned_5th_percentile, fliplr(binned_95th_percentile)], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', '90% Data Range');
errorbar(bin_centers, binned_mean, binned_sem, 'k-o', 'MarkerFaceColor', 'k', 'DisplayName', 'Mean'); % Plot with error bars
title('Stress Direction vs Wind Direction');
xlabel('Wind Direction (degrees)');
ylabel('Stress Direction');
ylim([-25 25]); 
legend show;
grid on;
y_limits = ylim;
plot([-45 -45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([45 45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([0 0], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
text(-67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(2, y_limits(1), 'Onshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
x_limits = xlim;
text(x_limits(2), y_limits(1), 'South', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
text(x_limits(1), y_limits(1), 'North', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
hold off;

% Subplot 2: Histogram of wave direction
subplot(3, 1, 2);
histogram(adjusted_wave_direction, 'BinEdges', -180:10:180, 'Normalization', 'probability');
title('Histogram of Wave Direction');
xlabel('Direction (degrees)');
ylabel('Frequency');
xlim([-120 120]);

% Subplot 3: Histogram of current direction
subplot(3, 1, 3);
histogram(adjusted_current_direction, 'BinEdges', -180:10:180, 'Normalization', 'probability');
title('Histogram of Current Direction');
xlabel('Direction (degrees)');
ylabel('Frequency');
xlim([-120 120]);

% Link x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');
%% CIRCLE AND TIME SERIES
% Plotting the data
figure;

% Polar plot for Wind Direction and Speed
subplot(1, 3, 1);
polarplot(deg2rad(filtered_windDirection), filtered_wind_speed, 'r.');
title('Wind Direction and Speed');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
rlim([0 max_wind_speed]);

% Polar plot for Wave Mean Direction and Speed
subplot(1, 3, 2);
polarplot(deg2rad(filtered_waveMeanDirection), filtered_wave_speed, 'b.');
title('Wave Mean Direction and Speed');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
rlim([0 max_wave_speed]);

% Polar plot for Current Direction and Speed
subplot(1, 3, 3);
polarplot(deg2rad(filtered_currentDirection), filtered_current_speed, 'g.');
title('Current Direction and Speed');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
rlim([0 max_current_speed]);

% Add annotations for Upcoast, Downcoast, North, and South
for i = 1:3
    subplot(1, 3, i);
    ax = gca;
    r_max = ax.RLim(2);
    text(deg2rad(-90), r_max, 'Upcoast', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(deg2rad(90), r_max, 'Downcoast', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(deg2rad(-72), r_max, 'N', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(deg2rad(108), r_max, 'S', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Save the figure
saveas(gcf, 'Polar_Plots.png');

% Plotting time series
figure;

% Time series for Wind Direction
subplot(3, 1, 1);
plot(filtered_dates, filtered_windDirection, 'r');
datetick('x', 'mmm dd, yyyy');
xlabel('Date');
ylabel('Wind Direction (deg)');
title('Time Series of Wind Direction');

% Time series for Wave Mean Direction
subplot(3, 1, 2);
plot(filtered_dates, filtered_waveMeanDirection, 'b');
datetick('x', 'mmm dd, yyyy');
xlabel('Date');
ylabel('Wave Mean Direction (deg)');
title('Time Series of Wave Mean Direction');

% Time series for Current Direction
subplot(3, 1, 3);
plot(filtered_dates, filtered_currentDirection, 'g');
datetick('x', 'mmm dd, yyyy');
xlabel('Date');
ylabel('Current Direction (deg)');
title('Time Series of Current Direction');

%% DIFF CURRENTS

% Extract necessary columns
filtered_Wdir_measured = combined_data(:, 2);
filtered_Sdir_measured = combined_data(:, 5);
adjusted_wave_direction = combined_data(:, 3);
adjusted_current_direction = combined_data(:, 4);

% Filter data where adjusted_current_direction > 0
number = 0;
idx_current_gt_0 = adjusted_current_direction > number;
Wdir_gt_0 = filtered_Wdir_measured(idx_current_gt_0);
Sdir_gt_0 = filtered_Sdir_measured(idx_current_gt_0);

% Plot wind direction versus stress direction for current direction > 0
createFilledPlot(Wdir_gt_0, Sdir_gt_0, 'Wind Direction vs. Stress Direction (Current Direction = South)');

% Filter data where adjusted_current_direction < 0
idx_current_lt_0 = adjusted_current_direction < number;
Wdir_lt_0 = filtered_Wdir_measured(idx_current_lt_0);
Sdir_lt_0 = filtered_Sdir_measured(idx_current_lt_0);

% Plot wind direction versus stress direction for current direction < 0
createFilledPlot(Wdir_lt_0, Sdir_lt_0, 'Wind Direction vs. Stress Direction (Current Direction = North)');

%% MULTIVARIATE

% Load the combined data
data = readtable('Corrected_Combined_Pier_Data_BE.csv');

% Extract relevant variables
windDirection = data.WindDirection;
waveMeanDirection = data.WaveMeanDirection;
currentDirection = data.CurrentDirection;
stressDirection = data.StressDirection;

% Calculate the true stress direction (wind direction + stress direction)
true_stressDirection = windDirection + stressDirection;

% Prepare the design matrix for regression
X = [windDirection, waveMeanDirection, currentDirection];
X = [ones(size(X,1),1) X]; % Add a column of ones for the intercept term

% Perform multiple linear regression
[b,~,~,~,stats] = regress(true_stressDirection, X);

% Display the regression coefficients and statistics
disp('Regression Coefficients:');
disp(b);
disp('R-squared:');
disp(stats(1));
disp('F-statistic and p-value:');
disp([stats(2), stats(3)]);
disp('Standard error of the estimate:');
disp(stats(4));

% Check for multicollinearity using VIF
vif = @(X) diag(inv(corrcoef(X))).^(1/2);
vif_values = vif(X(:,2:end)); % Exclude the intercept for VIF calculation
disp('Variance Inflation Factors (VIF):');
disp(vif_values);

% Calculate residuals
residuals = true_stressDirection - X * b;

% Create residual plots
figure;
subplot(2,1,1);
scatter(X(:,2), residuals, 'filled');
xlabel('Wind Direction');
ylabel('Residuals');
title('Residuals vs Wind Direction');
grid on;

subplot(2,1,2);
scatter(X(:,3), residuals, 'filled');
xlabel('Wave Mean Direction');
ylabel('Residuals');
title('Residuals vs Wave Mean Direction');
grid on;

% Save the figure
saveas(gcf, 'Residuals_Plot.png');

% Plot actual vs predicted values
predicted_stressDirection = X * b;

figure;
scatter(true_stressDirection, predicted_stressDirection, 'filled');
hold on;
plot([min(true_stressDirection), max(true_stressDirection)], [min(true_stressDirection), max(true_stressDirection)], 'r--');
xlabel('True Stress Direction');
ylabel('Predicted Stress Direction');
title('True Stress Direction vs Predicted Stress Direction');
legend('Data Points', 'Ideal Fit Line', 'Location', 'Best');
grid on;

%% LINE GRAPH
% Load the new combined data
combined_data = readtable('Corrected_Combined_Pier_Data_BE.csv');

% Extract relevant columns from the combined data
filtered_dates = combined_data.Date;
filtered_windDirection = combined_data.WindDirection;
filtered_waveMeanDirection = combined_data.WaveMeanDirection;
filtered_currentDirection = combined_data.CurrentDirection;
filtered_stressDirection = combined_data.StressDirection;
filtered_wind_speed = combined_data.WindSpeed;
filtered_wave_height = combined_data.WaveHeight;
filtered_current_speed = combined_data.CurrentSpeed;
filtered_wave_speed = combined_data.WaveSpeed;

% Apply data quality control based on the criteria used in the preset
valid_idx = filtered_wind_speed >= 0 & ...
            filtered_wave_speed >= 0 & ...
            filtered_current_speed >= 0 & ...
            ~isnan(filtered_windDirection) & ...
            ~isnan(filtered_waveMeanDirection) & ...
            ~isnan(filtered_currentDirection) & ...
            ~isnan(filtered_stressDirection);

% Filter data
filtered_dates = filtered_dates(valid_idx);
filtered_windDirection = filtered_windDirection(valid_idx);
filtered_waveMeanDirection = filtered_waveMeanDirection(valid_idx);
filtered_currentDirection = filtered_currentDirection(valid_idx);
filtered_stressDirection = filtered_stressDirection(valid_idx);
filtered_wind_speed = filtered_wind_speed(valid_idx);
filtered_wave_height = filtered_wave_height(valid_idx);
filtered_current_speed = filtered_current_speed(valid_idx);
filtered_wave_speed = filtered_wave_speed(valid_idx);

% Normalize direction to be between -180 and 180 degrees
normalize_direction = @(x) mod(x + 180, 360) - 180;
filtered_windDirection = normalize_direction(filtered_windDirection);

% Define the wind direction ranges
alongshore_idx = (filtered_windDirection >= -90 & filtered_windDirection < -45) | ...
                 (filtered_windDirection > 45 & filtered_windDirection <= 90);
onshore_idx = filtered_windDirection >= -45 & filtered_windDirection <= 45;

% Extract stress direction data for each category
alongshore_windDirection = filtered_windDirection(alongshore_idx);
alongshore_stressDirection = filtered_stressDirection(alongshore_idx);
onshore_windDirection = filtered_windDirection(onshore_idx);
onshore_stressDirection = filtered_stressDirection(onshore_idx);

% Bin the wind direction data
bin_size = 5; 
bin_edges = -90:bin_size:90;
bin_centers = bin_edges(1:end-1) + bin_size / 2;

% Calculate the mean stress direction for each bin and count the data points
mean_alongshore_stressDirection = zeros(size(bin_centers));
mean_onshore_stressDirection = zeros(size(bin_centers));
count_alongshore = zeros(size(bin_centers));
count_onshore = zeros(size(bin_centers));

for i = 1:length(bin_centers)
    % Alongshore
    bin_idx_alongshore = alongshore_windDirection >= bin_edges(i) & alongshore_windDirection < bin_edges(i+1);
    count_alongshore(i) = sum(bin_idx_alongshore);
    if count_alongshore(i) > 0
        mean_alongshore_stressDirection(i) = mean(alongshore_stressDirection(bin_idx_alongshore));
    else
        mean_alongshore_stressDirection(i) = NaN;
    end
    
    % Onshore
    bin_idx_onshore = onshore_windDirection >= bin_edges(i) & onshore_windDirection < bin_edges(i+1);
    count_onshore(i) = sum(bin_idx_onshore);
    if count_onshore(i) > 0
        mean_onshore_stressDirection(i) = mean(onshore_stressDirection(bin_idx_onshore));
    else
        mean_onshore_stressDirection(i) = NaN;
    end
end

% Create the plot
figure;
hold on;

% Plot the mean stress direction for alongshore and onshore
plot(bin_centers, mean_alongshore_stressDirection, 'r-o', 'LineWidth', 2, 'DisplayName', 'Alongshore Mean Stress Direction');
plot(bin_centers, mean_onshore_stressDirection, 'b-o', 'LineWidth', 2, 'DisplayName', 'Onshore Mean Stress Direction');

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

% Add North and South text annotations, switched as per the requirement
x_limits = xlim;
text(90, y_limits(1), 'South', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
text(-90, y_limits(1), 'North', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

% Set x-axis limits to -90 to 90
xlim([-90 90]);

hold off;