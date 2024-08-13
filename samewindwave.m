% Load the corrected combined data
data = readtable('Corrected_Combined_Pier_Data_BE.csv');
date_nums = datenum(data.Date);

% Extract relevant columns
waveDirection = data.WaveMeanDirection;
windDirection = data.WindDirection;
currentDirection = data.CurrentDirection;
stressDirection = data.StressDirection;

% Compute the true stress direction (stress + wind direction)
true_stress_direction = stressDirection + windDirection;

% Normalize directions to the range [0, 360)
normalize_direction = @(x) mod(x, 360);
waveDirection = normalize_direction(waveDirection);
windDirection = normalize_direction(windDirection);
currentDirection = normalize_direction(currentDirection);
true_stress_direction = normalize_direction(true_stress_direction);

% Identify instances where wave and wind directions are the same (or close to each other within a tolerance)
tolerance = 1; % Define a tolerance for comparison (e.g., 1 degree)
same_wave_wind_idx = abs(waveDirection - windDirection) <= tolerance;

% Initialize a logical array to store results
is_between_wind_and_current = false(size(true_stress_direction));

% Check if true stress direction is between wind and current direction
for i = 1:length(true_stress_direction)
    if same_wave_wind_idx(i)
        wind_dir = windDirection(i);
        current_dir = currentDirection(i);
        true_stress_dir = true_stress_direction(i);
        
        % Adjust for circular range
        if current_dir < wind_dir
            current_dir = current_dir + 360;
        end
        
        if true_stress_dir < wind_dir
            true_stress_dir = true_stress_dir + 360;
        end
        
        % Check if true stress direction is between wind and current direction
        if true_stress_dir >= wind_dir && true_stress_dir <= current_dir
            is_between_wind_and_current(i) = true;
        end
    end
end

% Calculate the percentage of instances
percentage_between = sum(is_between_wind_and_current) / sum(same_wave_wind_idx) * 100;

% Output the result
fprintf('Percentage of time true stress direction is between wind and current direction when wave and wind directions are the same: %.2f%%\n', percentage_between);