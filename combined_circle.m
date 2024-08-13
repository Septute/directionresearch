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
    
    % Ensure directions are within the range [-180, 180)
    normalize_direction = @(x) mod(x + 180, 360) - 180;
    filtered_waveMeanDirection = normalize_direction(filtered_waveMeanDirection);
    filtered_currentDirection = normalize_direction(filtered_currentDirection);
    
    % Adjust wave and current directions to the shoreward reference frame
    coastline_direction_deg = 72; % Coastline direction in degrees
    
    adjust_direction = @(direction) arrayfun(@(x) ...
        (x >= 0 && x <= 72) * (x - coastline_direction_deg) + ...
        (x > 72 && x <= 360) * (x - coastline_direction_deg) + ...
        (x < 0) * (360 + x - coastline_direction_deg), direction);
    
    filtered_waveMeanDirection = adjust_direction(filtered_waveMeanDirection);
    filtered_currentDirection = adjust_direction(filtered_currentDirection);
    
    
    % Flip current and wave direction so that north is negative and south is positive
    filtered_currentDirection = -filtered_currentDirection;
    filtered_waveMeanDirection = -filtered_waveMeanDirection;
    
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