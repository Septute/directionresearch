% Load the data
data = readtable('Pier_Data_BE.csv', 'VariableNamingRule', 'preserve');

% Remove rows where Wdir.range is 'Offshore' or Wdir.deg is 9999
valid_idx = ~strcmp(data.('Wdir.range'), 'Offshore') & data.('Wdir.deg') ~= 9999;

% Apply data quality control
valid_idx = valid_idx & data.Cdz <= 0.01;
valid_idx = valid_idx & ~isnan(data.Cdz) & ~isnan(data.('Wdir.deg')) & ~isnan(data.('Sdir.measured')) & ~isnan(data.Uz);
valid_idx = valid_idx & data.r2_uw >= 0.8;

% Remove data from April 2022 to July 2022
date_nums = datenum(data.year, data.month, data.day);
remove_dates = date_nums >= datenum(2022, 4, 1) & date_nums <= datenum(2022, 7, 31);
valid_idx = valid_idx & ~remove_dates;

% Define specific date ranges (3 days around the specified dates)
specific_date_range1 = (date_nums >= datenum(2022, 10, 5)) & (date_nums <= datenum(2022, 10, 11));
specific_date_range2 = (date_nums >= datenum(2022, 10, 30)) & (date_nums <= datenum(2022, 11, 4));
specific_dates = specific_date_range1 | specific_date_range2;

% Apply specific date filtering to the valid data
filtered_valid_idx = valid_idx & specific_dates;
filtered_data = data(filtered_valid_idx, :);

% Extract relevant columns and convert to numeric where needed
wind_direction = filtered_data.('Wdir.deg');
stress_direction = filtered_data.('Sdir.measured');
wave_direction = str2double(filtered_data.waveMeanDirection_awc6);
current_direction = str2double(filtered_data.currentDirection_awc6);

% Adjust wave and current directions
wave_direction = mod(wave_direction - 72 + 360, 360);
wave_direction(wave_direction > 180) = wave_direction(wave_direction > 180) - 360;

current_direction = mod(current_direction - 72 + 360, 360);
current_direction(current_direction > 180) = current_direction(current_direction > 180) - 360;

% Filter out NaN values
valid_wave_idx = ~isnan(wave_direction);
valid_current_idx = ~isnan(current_direction);

% Create a single plot
figure;
hold on;

% Plot Wind vs. Stress Direction
scatter(wind_direction, stress_direction, 50, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Stress Direction');

% Plot Wind vs. Wave Direction
scatter(wind_direction(valid_wave_idx), wave_direction(valid_wave_idx), 50, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Wave Direction');

% Plot Wind vs. Current Direction
scatter(wind_direction(valid_current_idx), current_direction(valid_current_idx), 50, 'filled', 'MarkerFaceColor', 'g', 'DisplayName', 'Current Direction');

% Add labels and title
xlabel('Wind Direction (degrees)');
ylabel('Direction (degrees)');
title('Interaction of Stress, Wave, and Current Directions with Wind Direction');
legend('show');

% Add reference lines and text
y_limits = ylim;
plot([-45 -45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([45 45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([0 0], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
text(-67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(2, y_limits(1), 'Onshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
x_limits = xlim;
text(x_limits(2), y_limits(1), 'North', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
text(x_limits(1), y_limits(1), 'South', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

hold off;