% Load the data
load('All_Pier_Data_30_Min_Bins.mat');

% Apply data quality control
valid_idx = Cz_p <= 0.01; 


mvalid_idx =  r2_uw_p >= .9 & valid_idx; 

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

% Create the plot
figure;
hold on;

% Plot boxplots for each bin
for i = 1:length(bin_centers)
    bin_idx = shifted_Wdir_measured_p >= bin_edges(i) & shifted_Wdir_measured_p < bin_edges(i+1);
    bin_data = filtered_Sdir_measured(bin_idx);
    
    if ~isempty(bin_data)
        % Calculate the position for the boxplot
        positions = repmat(bin_centers(i), size(bin_data));
        
        % Plot boxplot
        boxplot(bin_data, 'Positions', bin_centers(i), 'Widths', bin_size * 0.8, 'Colors', 'k', 'Symbol', 'r+');
    end
end

% Adjust x-axis tick labels
set(gca, 'XTick', bin_centers);
set(gca, 'XTickLabel', arrayfun(@num2str, bin_centers, 'UniformOutput', false));

% Add grid and labels
title('Stress Direction vs Wind Direction');
xlabel('Wind Direction (degrees)');
ylabel('Stress Direction');
grid on;

% Add reference lines and text
y_limits = ylim;
plot([-45 -45], y_limits, '--k', 'LineWidth', 1);
plot([45 45], y_limits, '--k', 'LineWidth', 1);
plot([0 0], y_limits, '--k', 'LineWidth', 1);
text(-67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(2, y_limits(1), 'Onshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
text(67.5, y_limits(1), 'Alongshore', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'k', 'FontName', 'Arial');
