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
coastline_direction_deg = 72; % Convert to shoreward reference
filtered_waveMeanDirection = filtered_waveMeanDirection - coastline_direction_deg;
filtered_currentDirection = filtered_currentDirection - coastline_direction_deg;

% Ensure directions are within the range [-180, 180)
normalize_direction = @(x) mod(x + 180, 360) - 180;
filtered_waveMeanDirection = normalize_direction(filtered_waveMeanDirection);
filtered_currentDirection = normalize_direction(filtered_currentDirection);

% Correctly center wave and current directions around 0 with -90 as upcoast and 90 as downcoast
filtered_waveMeanDirection = normalize_direction(filtered_waveMeanDirection - 90);
filtered_currentDirection = normalize_direction(filtered_currentDirection - 90);

% Apply additional +90 degree rotation to correctly align directions
filtered_waveMeanDirection = normalize_direction(filtered_waveMeanDirection + 90);
filtered_currentDirection = normalize_direction(filtered_currentDirection + 90);

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

% Extract necessary columns
filtered_Wdir_measured = combined_data(:, 2);
filtered_Sdir_measured = combined_data(:, 5);
adjusted_wave_direction = combined_data(:, 3);
adjusted_current_direction = combined_data(:, 4);

% Define bins for wind direction
bin_size = 9;
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
