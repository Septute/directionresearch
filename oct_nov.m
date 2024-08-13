% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Remove rows where Wdir.range is 'Offshore' or Wdir.deg is 9999
valid_idx = ~strcmp(data.('Wdir.range'), 'Offshore') & data.('Wdir.deg') ~= 9999;

% Apply data quality control
valid_idx = valid_idx & data.Cdz <= 0.01;
valid_idx = valid_idx & ~isnan(data.Cdz) & ~isnan(data.('Wdir.deg')) & ~isnan(data.('Sdir.measured')) & ~isnan(data.Uz);
valid_idx = valid_idx & data.r2_uw >= 0.9;

% Remove data from April 2022 to July 2022
date_nums = datenum(data.year, data.month, data.day);
remove_dates = date_nums >= datenum(2022, 4, 1) & date_nums <= datenum(2022, 7, 31);
valid_idx = valid_idx & ~remove_dates;

% Define specific date ranges (3 days around the specified dates)
october_dates = (date_nums >= datenum(2022, 10, 5)) & (date_nums <= datenum(2022, 10, 11));
november_dates = (date_nums >= datenum(2022, 10, 30)) & (date_nums <= datenum(2022, 11, 4));

% Apply specific date filtering to the valid data
october_valid_idx = valid_idx & october_dates;
november_valid_idx = valid_idx & november_dates;

october_data = data(october_valid_idx, :);
november_data = data(november_valid_idx, :);

% Extract relevant columns and convert to numeric where needed
time_october = datetime(october_data.year, october_data.month, october_data.day, october_data.hour, october_data.minute, 0);
time_november = datetime(november_data.year, november_data.month, november_data.day, november_data.hour, november_data.minute, 0);

wind_direction_october = october_data.('Wdir.deg');
wind_direction_november = november_data.('Wdir.deg');

stress_direction_october = october_data.('Sdir.measured');
stress_direction_november = november_data.('Sdir.measured');

wave_direction_october = str2double(october_data.waveMeanDirection_awc6);
wave_direction_november = str2double(november_data.waveMeanDirection_awc6);

current_direction_october = str2double(october_data.currentDirection_awc6);
current_direction_november = str2double(november_data.currentDirection_awc6);

% Adjust wave and current directions
wave_direction_october = mod(wave_direction_october - 72 + 360, 360);
wave_direction_october(wave_direction_october > 180) = wave_direction_october(wave_direction_october > 180) - 360;

wave_direction_november = mod(wave_direction_november - 72 + 360, 360);
wave_direction_november(wave_direction_november > 180) = wave_direction_november(wave_direction_november > 180) - 360;

current_direction_october = mod(current_direction_october - 72 + 360, 360);
current_direction_october(current_direction_october > 180) = current_direction_october(current_direction_october > 180) - 360;

current_direction_november = mod(current_direction_november - 72 + 360, 360);
current_direction_november(current_direction_november > 180) = current_direction_november(current_direction_november > 180) - 360;

% Combine wind and stress directions
combined_stress_wind_direction_october = mod(wind_direction_october + stress_direction_october, 360);
combined_stress_wind_direction_november = mod(wind_direction_november + stress_direction_november, 360);

% Filter out NaN values
valid_wave_idx_october = ~isnan(wave_direction_october);
valid_wave_idx_november = ~isnan(wave_direction_november);

valid_current_idx_october = ~isnan(current_direction_october);
valid_current_idx_november = ~isnan(current_direction_november);

% Create plots
figure;

% Plot for October
subplot(2, 1, 1);
hold on;
plot(time_october, combined_stress_wind_direction_october, 'r-o', 'DisplayName', 'Wind+Stress Direction');
plot(time_october(valid_wave_idx_october), wave_direction_october(valid_wave_idx_october), 'b-s', 'DisplayName', 'Wave Direction');
plot(time_october(valid_current_idx_october), current_direction_october(valid_current_idx_october), 'g-^', 'DisplayName', 'Current Direction');
xlabel('Time');
ylabel('Direction (degrees)');
title('Directions over Time (October)');
legend('show');
hold off;

% Plot for November
subplot(2, 1, 2);
hold on;
plot(time_november, combined_stress_wind_direction_november, 'r-o', 'DisplayName', 'Wind+Stress Direction');
plot(time_november(valid_wave_idx_november), wave_direction_november(valid_wave_idx_november), 'b-s', 'DisplayName', 'Wave Direction');
plot(time_november(valid_current_idx_november), current_direction_november(valid_current_idx_november), 'g-^', 'DisplayName', 'Current Direction');
xlabel('Time');
ylabel('Direction (degrees)');
title('Directions over Time (November)');
legend('show');
hold off;