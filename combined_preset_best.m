% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Convert date column to MATLAB datenum
date_nums = datenum(data.Date);

% Convert wave and current direction data for awc6
waveMeanDirection_awc6 = str2double(data.('waveMeanDirection_awc6'));
currentDirection_awc6 = str2double(data.('currentDirection_awc6'));

% Step 1: Rotate the directions to align with the desired orientation
% This is done by adjusting the coordinate system such that north is -72°,
% south is 108°, upcoast is -90°, and downcoast is 90°
waveMeanDirection_awc6 = wrapTo360(waveMeanDirection_awc6 - 72);
currentDirection_awc6 = wrapTo360(currentDirection_awc6 - 72);

% Step 2: Normalize directions to the range [-180, 180)
waveMeanDirection_awc6 = mod(waveMeanDirection_awc6 + 180, 360) - 180;
currentDirection_awc6 = mod(currentDirection_awc6 + 180, 360) - 180;

% Filter out invalid directions (9999 or 'NA')
valid_wave_direction_idx_awc6 = ~isnan(waveMeanDirection_awc6) & waveMeanDirection_awc6 ~= 9999;
valid_current_direction_idx_awc6 = ~isnan(currentDirection_awc6) & currentDirection_awc6 ~= 9999;

% Interpolate missing data for awc6 using unit vectors
wave_vectors_awc6 = [cosd(waveMeanDirection_awc6(valid_wave_direction_idx_awc6)), sind(waveMeanDirection_awc6(valid_wave_direction_idx_awc6))];
current_vectors_awc6 = [cosd(currentDirection_awc6(valid_current_direction_idx_awc6)), sind(currentDirection_awc6(valid_current_direction_idx_awc6))];

wave_time_idx_awc6 = find(valid_wave_direction_idx_awc6);
current_time_idx_awc6 = find(valid_current_direction_idx_awc6);

waveMeanDirection_awc6_interp = nan(size(waveMeanDirection_awc6));
currentDirection_awc6_interp = nan(size(currentDirection_awc6));

if ~isempty(wave_time_idx_awc6)
    wave_vectors_awc6_interp = interp1(wave_time_idx_awc6, wave_vectors_awc6, 1:length(waveMeanDirection_awc6), 'pchip', 'extrap');
    waveMeanDirection_awc6_interp = atan2d(wave_vectors_awc6_interp(:,2), wave_vectors_awc6_interp(:,1));
end

if ~isempty(current_time_idx_awc6)
    current_vectors_awc6_interp = interp1(current_time_idx_awc6, current_vectors_awc6, 1:length(currentDirection_awc6), 'pchip', 'extrap');
    currentDirection_awc6_interp = atan2d(current_vectors_awc6_interp(:,2), current_vectors_awc6_interp(:,1));
end

% Repeat the process for awc4.5
waveMeanDirection_awc45 = str2double(data.('waveMeanDirection_awc4.5'));
currentDirection_awc45 = str2double(data.('currentDirection_awc4.5'));

% Rotate the directions to align with the desired orientation
waveMeanDirection_awc45 = wrapTo360(waveMeanDirection_awc45 - 72);
currentDirection_awc45 = wrapTo360(currentDirection_awc45 - 72);

% Normalize directions to the range [-180, 180)
waveMeanDirection_awc45 = mod(waveMeanDirection_awc45 + 180, 360) - 180;
currentDirection_awc45 = mod(currentDirection_awc45 + 180, 360) - 180;

valid_wave_direction_idx_awc45 = ~isnan(waveMeanDirection_awc45) & waveMeanDirection_awc45 ~= 9999;
valid_current_direction_idx_awc45 = ~isnan(currentDirection_awc45) & currentDirection_awc45 ~= 9999;

% Interpolate missing data for awc45 using unit vectors
wave_vectors_awc45 = [cosd(waveMeanDirection_awc45(valid_wave_direction_idx_awc45)), sind(waveMeanDirection_awc45(valid_wave_direction_idx_awc45))];
current_vectors_awc45 = [cosd(currentDirection_awc45(valid_current_direction_idx_awc45)), sind(currentDirection_awc45(valid_current_direction_idx_awc45))];

wave_time_idx_awc45 = find(valid_wave_direction_idx_awc45);
current_time_idx_awc45 = find(valid_current_direction_idx_awc45);

waveMeanDirection_awc45_interp = nan(size(waveMeanDirection_awc45));
currentDirection_awc45_interp = nan(size(currentDirection_awc45));

if ~isempty(wave_time_idx_awc45)
    wave_vectors_awc45_interp = interp1(wave_time_idx_awc45, wave_vectors_awc45, 1:length(waveMeanDirection_awc45), 'pchip', 'extrap');
    waveMeanDirection_awc45_interp = atan2d(wave_vectors_awc45_interp(:,2), wave_vectors_awc45_interp(:,1));
end

if ~isempty(current_time_idx_awc45)
    current_vectors_awc45_interp = interp1(current_time_idx_awc45, current_vectors_awc45, 1:length(currentDirection_awc45), 'pchip', 'extrap');
    currentDirection_awc45_interp = atan2d(current_vectors_awc45_interp(:,2), current_vectors_awc45_interp(:,1));
end

% Combine awc6 and awc45 data
waveMeanDirection_combined = waveMeanDirection_awc6_interp;
currentDirection_combined = currentDirection_awc6_interp;

% Replace missing values in awc6 with values from awc45 for specified dates
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

wave_height_awc6_interp = interp1(find(~isnan(wave_height_awc6)), wave_height_awc6(~isnan(wave_height_awc6)), 1:length(wave_height_awc6), 'pchip', 'extrap');
current_speed_awc6_interp = interp1(find(~isnan(current_speed_awc6)), current_speed_awc6(~isnan(current_speed_awc6)), 1:length(current_speed_awc6), 'pchip', 'extrap');

filtered_wave_height = wave_height_awc6_interp(valid_idx);
filtered_current_speed = current_speed_awc6_interp(valid_idx);

% Calculate wave speed
g = 9.81; % gravity in m/s^2
wave_period = 30 * 60; % 30 minutes in seconds
filtered_wave_speed = sqrt(g * filtered_wave_height / (2 * pi / wave_period));

% Ensure all filtered variables have the same number of rows
nRows = length(filtered_dates);
filtered_waveMeanDirection = filtered_waveMeanDirection(:);
filtered_currentDirection = filtered_currentDirection(:);
filtered_windDirection = filtered_windDirection(:);
filtered_stressDirection = filtered_stressDirection(:);
filtered_wind_speed = filtered_wind_speed(:);
filtered_wave_height = filtered_wave_height(:);
filtered_current_speed = filtered_current_speed(:);
filtered_wave_speed = filtered_wave_speed(:);

% Combine data into one matrix and remove rows with NaN values
combined_data = [filtered_dates, filtered_windDirection, filtered_waveMeanDirection, ...
                 filtered_currentDirection, filtered_stressDirection, ...
                 filtered_wind_speed, filtered_wave_height, filtered_current_speed, filtered_wave_speed];
combined_data = combined_data(~any(isnan(combined_data), 2), :);

% Save the combined data to a new CSV file
output_table = array2table(combined_data, 'VariableNames', ...
    {'Date', 'WindDirection', 'WaveMeanDirection', 'CurrentDirection', ...
     'StressDirection', 'WindSpeed', 'WaveHeight', 'CurrentSpeed', 'WaveSpeed'});
writetable(output_table, 'Corrected_Combined_Pier_Data_BE.csv');

% Determine individual radius limits for each plot
max_wind_speed = max(filtered_wind_speed);
max_wave_speed = max(filtered_wave_speed);
max_current_speed = max(filtered_current_speed);

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
    text(deg2rad(-90), r_max, 'Upcoast', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(deg2rad(90), r_max, 'Downcoast', 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(deg2rad(-72), r_max, 'N', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(deg2rad(108), r_max, 'S', 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Final adjustments to figure
sgtitle('Polar Plots of Wind, Wave, and Current Directions with Speeds');

% Define bins for wind direction
bin_size = 15;
bin_edges = -90:bin_size:90;
bin_centers = bin_edges(1:end-1) + bin_size / 2;
binned_95th_percentile = zeros(size(bin_centers));
binned_5th_percentile = zeros(size(bin_centers));
binned_mean = zeros(size(bin_centers));
binned_sem = zeros(size(bin_centers)); % Standard error of the mean

for i = 1:length(bin_centers)
    bin_idx = filtered_windDirection >= bin_edges(i) & filtered_windDirection < bin_edges(i+1);
    bin_data = filtered_stressDirection(bin_idx);

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
histogram(filtered_waveMeanDirection, 'BinEdges', -180:5:180, 'Normalization', 'probability');
title('Histogram of Wave Direction');
xlabel('Direction (degrees)');
ylabel('Frequency');
xlim([-120 120]);

% Subplot 3: Histogram of current direction
subplot(3, 1, 3);
histogram(filtered_currentDirection, 'BinEdges', -180:5:180, 'Normalization', 'probability');
title('Histogram of Current Direction');
xlabel('Direction (degrees)');
ylabel('Frequency');
xlim([-120 120]);

% Link x-axes of all subplots
linkaxes(findall(gcf, 'Type', 'axes'), 'x');