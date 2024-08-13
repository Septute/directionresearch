% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Convert date column to MATLAB datenum
% Assuming your date column is named 'Date'
date_nums = datenum(data.Date);

% Interpolate wave and current direction data for awc6
waveMeanDirection_awc6 = str2double(data.('waveMeanDirection_awc6'));
currentDirection_awc6 = str2double(data.('currentDirection_awc6'));
wave_time_idx_awc6 = find(~isnan(waveMeanDirection_awc6) & waveMeanDirection_awc6 ~= 9999);
current_time_idx_awc6 = find(~isnan(currentDirection_awc6) & currentDirection_awc6 ~= 9999);

% Interpolate missing data using pchip method
waveMeanDirection_awc6_interp = interp1(wave_time_idx_awc6, waveMeanDirection_awc6(wave_time_idx_awc6), 1:length(waveMeanDirection_awc6), 'pchip', 'extrap');
currentDirection_awc6_interp = interp1(current_time_idx_awc6, currentDirection_awc6(current_time_idx_awc6), 1:length(currentDirection_awc6), 'pchip', 'extrap');

% Interpolate wave and current direction data for awc4.5
waveMeanDirection_awc45 = str2double(data.('waveMeanDirection_awc4.5'));
currentDirection_awc45 = str2double(data.('currentDirection_awc4.5'));
wave_time_idx_awc45 = find(~isnan(waveMeanDirection_awc45) & waveMeanDirection_awc45 ~= 9999);
current_time_idx_awc45 = find(~isnan(currentDirection_awc45) & currentDirection_awc45 ~= 9999);

% Interpolate missing data using pchip method
waveMeanDirection_awc45_interp = interp1(wave_time_idx_awc45, waveMeanDirection_awc45(wave_time_idx_awc45), 1:length(waveMeanDirection_awc45), 'pchip', 'extrap');
currentDirection_awc45_interp = interp1(current_time_idx_awc45, currentDirection_awc45(current_time_idx_awc45), 1:length(currentDirection_awc45), 'pchip', 'extrap');

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

% Adjust directions to align with the coastline (shoreward)
coastline_direction_deg = 72;
filtered_waveMeanDirection = filtered_waveMeanDirection - coastline_direction_deg;
filtered_currentDirection = filtered_currentDirection - coastline_direction_deg;

% Flip current and wave direction so that north is negative and south is positive
filtered_currentDirection = -filtered_currentDirection;
filtered_waveMeanDirection = -filtered_waveMeanDirection;

% Ensure all vectors are column vectors
filtered_windDirection = filtered_windDirection(:);
filtered_waveMeanDirection = filtered_waveMeanDirection(:);
filtered_currentDirection = filtered_currentDirection(:);
filtered_stressDirection = filtered_stressDirection(:);

% Combine data into one matrix and remove rows with NaN values
combined_data = [filtered_windDirection, filtered_waveMeanDirection, filtered_currentDirection, filtered_stressDirection];
combined_data = combined_data(~any(isnan(combined_data), 2), :);

% Separate predictors and response after removing NaN values
X = combined_data(:, 1:3);
y = combined_data(:, 4);

% Perform statistical tests
% Fit multiple linear regression model
mdl = fitlm(X, y);

% Display the model
disp(mdl);

% Calculate and display the R-squared value
R2 = mdl.Rsquared.Ordinary;
disp(['R-squared: ', num2str(R2)]);

% Perform correlation analysis
[R, P] = corrcoef(combined_data);

% Display the correlation matrix and p-values
disp('Correlation matrix:');
disp(R);
disp('P-values:');
disp(P);