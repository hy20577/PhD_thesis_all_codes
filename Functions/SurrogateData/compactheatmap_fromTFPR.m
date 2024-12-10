function compactheatmap_fromTFPR(truePositiveRate,falsePositiveRate, label)

% Calculate the size scaling factor for the circles
maxCircleSize = 1e3; % You can adjust this based on your preference
circleSize = maxCircleSize * falsePositiveRate;

% Create the heatmap of true positive rate
h = imagesc(truePositiveRate);
% h.GridVisible = "off"; 
% colormap hot; % Apply the colormap here
xlabel('$pc_{2}$', 'Interpreter', 'latex');
ylabel('$pc_{1}$', 'Interpreter', 'latex');
% Set the font size for x-axis and y-axis labels

c = colorbar('Location', 'eastoutside');
c.Label.String = 'TPR'; % Set the colorbar label
c.Label.FontSize = 30;
c.Label.Interpreter = 'latex';
% set(cb.Label, 'Location', 'manual');
% Adjust the Position property to position the label just above the colormap
% set(get(c, 'Title'), 'String', 'TPR');
% title(cb, 'TPR', 'Fontsize', 20, 'interpreter', 'latex');
clim([0, 1]);
% Get the size of the matrices
[numRows, numCols] = size(truePositiveRate);

hold on;

% Loop through each cell in the matrix and plot circles
for row = 1:numRows
    for col = 1:numCols
        centerX = col; % + 0.5; % X-coordinate of the cell center
        centerY = row; %+ 0.5; % Y-coordinate of the cell center
        radius = circleSize(row, col); % Circle radius
        
        if radius > 0
            scatter(centerX, centerY, radius, 'r', 'filled');
        end
    end
end

xticks(1:numCols);
xticklabels(num2str(linspace(0, 1, numCols)', '%.2f'));

yticks(1:numRows);
yticklabels(num2str(linspace(1, 0, numRows)', '%.2f'));


% Set square aspect ratio
axis equal;

for row = 1:numRows
    for col = 1:numCols
        if all([truePositiveRate(row, col) == 1, falsePositiveRate(row,col) == 0])
            centerX = col; % + 0.5; % X-coordinate of the cell center
            centerY = row; %+ 0.5; % Y-coordinate of the cell center
            radius = 1e2; % Circle radius
                if radius > 0
                    scatter(centerX, centerY, radius, 'g', 'filled');
                end
        end
    end
end

% ax = gca;
% ax.XAxis.Label.FontSize = 50;  % Set the font size for the x-axis label
% ax.YAxis.Label.FontSize = 50;

hold off;

% Adjust axis limits and labels
axis([0.5 numCols+0.5 0.5 numRows+0.5]);

lg = legend('FPR', 'Location', 'northoutside');
lg.Interpreter = 'latex';
set(lg, 'Fontsize', 30);

% set(findobj(obj, 'type', 'line'), 'MarkerSize', 150)
% icons_m = findobj(icons, 'type', 'line');
% set(icons_m,'MarkerSize',20);
legend boxoff 
% fontsize(lg, 20, 'points')

% saveas(gcf, sprintf('%s.fig', label))
% saveas(gcf, sprintf('%s', label), 'epsc')
end