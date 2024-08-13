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

% Save the combined data to a new CSV file
output_table = array2table(combined_data, 'VariableNames', {'Date', 'WindDirection', 'WaveMeanDirection', 'CurrentDirection', 'StressDirection', 'WindSpeed', 'WaveHeight', 'CurrentSpeed', 'WaveSpeed'});
writetable(output_table, 'Corrected_Combined_Pier_Data_BE.csv');

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

% Create a correlation matrix
correlation_matrix = corrcoef([windDirection, waveMeanDirection, currentDirection, stressDirection]);
disp('Correlation Matrix:');
disp(correlation_matrix);

% Visualize the correlation matrix
figure;
imagesc(correlation_matrix);
colorbar;
xticks(1:4);
yticks(1:4);
xticklabels({'WindDir', 'WaveDir', 'CurrentDir', 'StressDir'});
yticklabels({'WindDir', 'WaveDir', 'CurrentDir', 'StressDir'});
title('Correlation Matrix');
grid on;

% Save the figure
saveas(gcf, 'Correlation_Matrix.png');

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