
%% plot learning curves
load('logfiles/Part2/PychReview_Part2_Sim3_InductiveBias_6P3F_v2_h100_r20');

plotCurve = 4;

plotSEM = 1;

fig1 = figure(1);
set(fig1, 'Position', [100 100 750 180]);

maxIterationsPlot = 30;

mean_0shared = mean(batch_log.MSE_log_0shared(:,1:maxIterationsPlot));
mean_1shared1 = mean(batch_log.MSE_log_1shared1(:,1:maxIterationsPlot));
mean_1shared2 = mean(batch_log.MSE_log_1shared2(:,1:maxIterationsPlot));
mean_1shared3 = mean(batch_log.MSE_log_1shared3(:,1:maxIterationsPlot));
% mean_2shared = mean(batch_log.MSE_log_2shared(:,1:maxIterationsPlot));
% mean_3shared = mean(batch_log.MSE_log_3shared(:,1:maxIterationsPlot));

sem_0shared = std(batch_log.MSE_log_0shared(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared1 = std(batch_log.MSE_log_1shared1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared2 = std(batch_log.MSE_log_1shared2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared3 = std(batch_log.MSE_log_1shared3(:,1:maxIterationsPlot)) / sqrt(replications);
% sem_2shared = std(batch_log.MSE_log_2shared(:,1:maxIterationsPlot)) / sqrt(replications);
% sem_3shared = std(batch_log.MSE_log_3shared(:,1:maxIterationsPlot)) / sqrt(replications);

plotSettings;

colors(5,:) = 1-repmat(0.75,1,3);
colors(1,:) = 1-repmat(0.5,1,3);
colors(7,:) = 1-repmat(0.25,1,3);
colors(8,:) = 1-repmat(0,1,3);

currentColor = [100 181 203]/255;
switch plotCurve
    
    case 1
        colors(5,:) = currentColor;
    case 2
        colors(1,:) = currentColor;
    case 3
        colors(7,:) = currentColor;
    case 4
        colors(8,:) = currentColor;
end

if(~plotSEM)
    plot(1:maxIterationsPlot, mean_0shared, 'LineWidth', 3, 'color', colors(5,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared2, 'LineWidth', 3, 'color', colors(7,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared3, 'LineWidth', 3, 'color', colors(8,:)); hold on;
%     plot(1:maxIterationsPlot, mean_2shared, 'LineWidth', 3, 'color', colors(2,:)); hold on;
%     plot(1:maxIterationsPlot, mean_3shared, 'LineWidth', 3, 'color', colors(3,:)); hold off;
else
    if(plotCurve > 0)
        errorbar(1:maxIterationsPlot, mean_0shared, sem_0shared, 'LineWidth', 3, 'color', colors(5,:)); hold on;
    end
    if(plotCurve > 1)
        errorbar(1:maxIterationsPlot, mean_1shared1, sem_1shared1,'LineWidth', 3, 'color', colors(1,:)); hold on;
    end
    if(plotCurve > 2)
        errorbar(1:maxIterationsPlot, mean_1shared2, sem_1shared2, 'LineWidth', 3, 'color', colors(7,:)); hold on;
    end
    if(plotCurve > 3)
        errorbar(1:maxIterationsPlot, mean_1shared3, sem_1shared3, 'LineWidth', 3, 'color', colors(8,:)); hold on;
    end
%     errorbar(1:maxIterationsPlot, mean_2shared, sem_2shared, 'LineWidth', 3, 'color', colors(2,:)); hold on;
%     errorbar(1:maxIterationsPlot, mean_3shared, sem_3shared, 'LineWidth', 3, 'color', colors(3,:)); hold off;
end

xlabel('Training Iterations','FontSize', fontSize_xlabel-1);
ylabel({'Mean Squared Error on', 'Target Single Tasks'}, 'FontSize', fontSize_ylabel-1, 'FontName', fontName);
leg = legend('No Pre-Training', 'Pre-Trained on 1  Task', 'Pre-Trained on 2 Tasks', 'Pre-Trained on 3 Tasks','Location', 'northeast');
set(leg, 'FontName', fontName, 'fontSize', fontSize_legend+1);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca-1);

set(leg, 'TextColor', 'w');
set(leg, 'Color', 'none');
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

%% MDS PLOTS
plotSettings;

width = 200;
height = 180;
fontSize_title = 11;
fontSize_xlabel = 16;
fontSize_ylabel = fontSize_xlabel;

repetitionIdx  = 4; % 14 17
scalar_0 = 5;
scalar_1 = 5;
scalar_2 = scalar_1;
scalar_3 = scalar_1;

scalar_rest = 5;

markerSize = 100;
markerLineWidth = 2;
markerLineWidth_pretrained = 1;

fig = figure(1);
set(fig, 'Position', [100 100 width height]);

% plot single task training;
xlimit = [-1 1] * scalar_0;
ylimit = [-1 1] * scalar_0;

x = batch_log.MSE_hidden_0shared{repetitionIdx}(:,1);
y = batch_log.MSE_hidden_0shared{repetitionIdx}(:,2);

color = nan(numel(x,1),3);
color(1,:) = plotSettings.colors(cContrast1,:);
color(2,:) = plotSettings.colors(cContrast2,:);
color(3,:) = plotSettings.colors(cContrast3,:);

scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth);
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

% plot 1 shared
xlimit = [-1 1] * scalar_1;
ylimit = [-1 1] * scalar_1;

fig = figure(2);
set(fig, 'Position', [100 300 width height]);

x = batch_log.MSE_hidden_1shared1{repetitionIdx}(2:end,1);
y = batch_log.MSE_hidden_1shared1{repetitionIdx}(2:end,2);
x_pretrained = batch_log.MSE_hidden_1shared1{repetitionIdx}(1:1,1);
y_pretrained = batch_log.MSE_hidden_1shared1{repetitionIdx}(1:1,2);

color = zeros(numel(x,1),3);
color(1,:) = plotSettings.colors(cContrast1,:);
color(2,:) = plotSettings.colors(cContrast2,:);
color(3,:) = plotSettings.colors(cContrast3,:);
scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth); hold on;

color = zeros(numel(x_pretrained,1),3);
color(1,:) = plotSettings.colors(cContrast1,:);
scatter(x_pretrained, y_pretrained, markerSize, color, 'LineWidth', markerLineWidth_pretrained); hold off;
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

% plot 2 shared
xlimit = [-1 1] * scalar_2;
ylimit = [-1 1] * scalar_2;

fig = figure(3);
set(fig, 'Position', [100 500 width height]);

x = batch_log.MSE_hidden_1shared2{repetitionIdx}(3:end,1);
y = batch_log.MSE_hidden_1shared2{repetitionIdx}(3:end,2);
x_pretrained = batch_log.MSE_hidden_1shared2{repetitionIdx}(1:2,1);
y_pretrained = batch_log.MSE_hidden_1shared2{repetitionIdx}(1:2,2);

color = zeros(numel(x,1),3);
color(1,:) = plotSettings.colors(cContrast1,:);
color(2,:) = plotSettings.colors(cContrast2,:);
color(3,:) = plotSettings.colors(cContrast3,:);
scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth); hold on;

color = zeros(numel(x_pretrained,1),3);
color(1,:) = plotSettings.colors(cContrast1,:);
color(2,:) = plotSettings.colors(cContrast2,:);
scatter(x_pretrained, y_pretrained, markerSize, color, 'LineWidth', markerLineWidth_pretrained); hold off;
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

% plot 3 shared
xlimit = [-1 1] * scalar_3;
ylimit = [-1 1] * scalar_3;

fig = figure(4);
set(fig, 'Position', [100 700 width height]);

x = batch_log.MSE_hidden_1shared3{repetitionIdx}(4:end,1);
y = batch_log.MSE_hidden_1shared3{repetitionIdx}(4:end,2);
x_pretrained = batch_log.MSE_hidden_1shared3{repetitionIdx}(1:3,1);
y_pretrained = batch_log.MSE_hidden_1shared3{repetitionIdx}(1:3,2);

color = zeros(numel(x,1),3);
color(1,:) = plotSettings.colors(cContrast1,:);
color(2,:) = plotSettings.colors(cContrast2,:);
color(3,:) = plotSettings.colors(cContrast3,:);

scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth); hold on;
scatter(x_pretrained, y_pretrained, markerSize, color, 'LineWidth', markerLineWidth_pretrained); hold off;
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

