% Load the data
load('All_Pier_Data_30_Min_Bins.mat');

% Apply data quality control
valid_idx = Cz_p <= 0.01; 

mvalid_idx = r2_uw_p >= .9 & valid_idx 

% Filter the data for wind coming 90 degrees from onshore
onshore90_idx = onshore90_p == 1 & mvalid_idx;
filtered_Wdir_measured_p = Wdir_measured_p(onshore90_idx);
filtered_Sdir_measured = Sdir_measured_p(onshore90_idx);
filtered_Uz_p = Uz_p(onshore90_idx);


% Shift wind direction data
shifted_Wdir_measured_p = filtered_Wdir_measured_p;
shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) = shifted_Wdir_measured_p(shifted_Wdir_measured_p < 60) + 360;

% Find the mean of the shifted wind direction data
mean_shifted_Wdir = mean(shifted_Wdir_measured_p);

% Determine the shift needed to center data between -90 and 90 degrees
shift_needed = 90 - max(shifted_Wdir_measured_p);
shifted_Wdir_measured_p = shifted_Wdir_measured_p + shift_needed;

% Bin the data
bin_size = 15; % Change this value to adjust bin size
bin_edges = -90:bin_size:90;
bin_centers = bin_edges(1:end-1) + bin_size / 2;
binned_95th_percentile = zeros(size(bin_centers));
binned_5th_percentile = zeros(size(bin_centers));
binned_mean = zeros(size(bin_centers));

for i = 1:length(bin_centers)
    bin_idx = shifted_Wdir_measured_p >= bin_edges(i) & shifted_Wdir_measured_p < bin_edges(i+1);
    bin_data = filtered_Sdir_measured(bin_idx);
    
    if ~isempty(bin_data)
        binned_95th_percentile(i) = prctile(bin_data, 90);
        binned_5th_percentile(i) = prctile(bin_data, 10);
        binned_mean(i) = mean(bin_data);
    end
end

% Create the plot
figure;
hold on;

% Plot the 90% data range (5th to 95th percentiles)
fill([bin_centers, fliplr(bin_centers)], ...
     [binned_5th_percentile, fliplr(binned_95th_percentile)], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', '90% Data Range');

% Plot the mean line
plot(bin_centers, binned_mean, 'k-o', 'MarkerFaceColor', 'k', 'DisplayName', 'Mean');

% Add grid and labels
title('Stress Direction vs Wind Direction');
xlabel('Wind Direction (degrees)');
ylabel('Stress Direction');
legend show;
grid on;

% Add reference lines and text
y_limits = ylim;
plot([-45 -45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([45 45], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
plot([0 0], y_limits, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
text(-67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(2, y_limits(1), 'Onshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
% Add North and South text annotations
x_limits = xlim;
text(x_limits(2), y_limits(1), 'North', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
text(x_limits(1), y_limits(1), 'South', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

hold off;

