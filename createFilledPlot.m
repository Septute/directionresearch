function createFilledPlot(windDir, stressDir, plotTitle)
    bin_size = 15;
    bin_edges = -90:bin_size:90;
    bin_centers = bin_edges(1:end-1) + bin_size / 2;
    binned_95th_percentile = zeros(size(bin_centers));
    binned_5th_percentile = zeros(size(bin_centers));
    binned_mean = zeros(size(bin_centers));
    binned_sem = zeros(size(bin_centers)); % Standard error of the mean

    for i = 1:length(bin_centers)
        bin_idx = windDir >= bin_edges(i) & windDir < bin_edges(i+1);
        bin_data = stressDir(bin_idx);

        if ~isempty(bin_data)
            binned_95th_percentile(i) = prctile(bin_data, 90);
            binned_5th_percentile(i) = prctile(bin_data, 10);
            binned_mean(i) = mean(bin_data);
            binned_sem(i) = std(bin_data) / sqrt(length(bin_data)); % Calculate the standard error of the mean
        end
    end

    figure;
    hold on;

    fill([bin_centers, fliplr(bin_centers)], ...
         [binned_5th_percentile, fliplr(binned_95th_percentile)], ...
         [0.8 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', '90% Data Range');

    errorbar(bin_centers, binned_mean, binned_sem, 'k-o', 'MarkerFaceColor', 'k', 'DisplayName', 'Mean'); % Plot with error bars

    title(plotTitle);
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
end