%% data settings
clear all;
clc;
close all;

datafilePrefix = 'PsychReview_Part1_Sim4_PRP_LS2000_2P8F_4tasks_22_h100_r';

Nhidden = 100;
correlation_threshold = 0.5;
plotBlack = 0;

% load data
folder = 'logfiles/Part1/';
files = dir(folder);

% get valid file names
validFileNames = {};
for fileIndex =1:length(files)
    % check if this is a desired data file
    if(~isempty(strfind(files(fileIndex).name, datafilePrefix)))
        validFileNames{end+1} = files(fileIndex).name;
    end
end
if(~isempty(validFileNames))
    numRepetitionsNet = length(validFileNames);
    % load every valid file
    for fileIndex = 1:length(validFileNames);
        disp(['loading ' validFileNames{fileIndex} '...']);
        load(strcat(folder, validFileNames{fileIndex}));
        % initial setup of batch_log_tmp
        if(fileIndex == 1)
            numHiddenUnitsTested = size(batch_log,1);
            numTausTested = length(tau_range);
            batch_log_tmp = repmat(batch_log(1,1), numHiddenUnitsTested, numRepetitionsNet, numTausTested);
        end
        batch_log_tmp(:, fileIndex, :) = batch_log(:, rep, :);
    end
else
    error('No valid file names found');
end
batch_log = batch_log_tmp;

numHiddenIdx = find(numHiddenUnits == Nhidden,1);
if(isempty(numHiddenIdx))
    warning('Requested number of hidden units does not exist in data set');
end

corrThreshIdx = find(correlation_thresholds == correlation_threshold,1);
if(isempty(corrThreshIdx))
    warning('Requested correlation threshold for MIS extraction does not exist in data set');
end

numCorrelationThresholds = length(correlation_thresholds);

% load plotSettings
loadPlotSettings;

disp('Done loading.');

% PREP DATA

% prepare MIS data

for tauIdx = 1:length(tau_range)

     for rep =1:numRepetitionsNet
         
         if(isfield(batch_log, 'funcDependence_inc') & isfield(batch_log, 'funcDependence_con'))
             
             funcDependence_inc.task1RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task1RT_log;
             funcDependence_inc.task2RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task2RT_log;
         
             funcDependence_con.task1RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task1RT_log;
             funcDependence_con.task2RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task2RT_log;
         
         end
     end
 
end
 
%% PLOTS - CONGRUENT VS. INCONGRUENT (POWER POINT)

plotBlack = 0;
plotTask1 = 1;
plottedTaus = [0.1];
plottedTausAll = [0.1];

plotSettings;
if(plotBlack)
    colors = [254 209 103; ...
                   254 137 45; ...
                   47 205 58; ...
                   252 13 27]/255;
else
    colors = [255 140 41; ... % orange
              31 104 172; ... % blue
              44 155 37; ... % green
              255 0 0; ... % red
              0     0   0  ; ... % black
            142 142 142; ... % grey
            253 120 21; ... % orange
            233 100 1; ... % orange
            255 255 255] / 255; % white 
end

fig = figure(1);
if(plotBlack)
    set(fig, 'Position', [100 100 400 200]);
else
    set(fig, 'Position', [100 100 450 220]);
end
x = SOA_range * 0.1;
ylimit = [0 1]; % [0 0.5];
lineWidth = 3;
plotSEM = 1;


if(plotBlack)
    cContrast1 = 4;
    cContrast2 = 3;
    color_gradient1 = nan(length(plottedTausAll), 3);
    color_gradient2 = nan(length(plottedTausAll), 3);

    offsetWhite = 1;
    whiteContrast1 = min(1, colors(cContrast1, :) + offsetWhite);
    whiteContrast2 = min(1, colors(cContrast2, :) + offsetWhite);
    color_gradient1(:,1) = linspace(colors(cContrast1,1), whiteContrast1(1), size(color_gradient1,1));
    color_gradient1(:,2) = linspace(colors(cContrast1,2), whiteContrast1(2), size(color_gradient1,1));
    color_gradient1(:,3) = linspace(colors(cContrast1,3), whiteContrast1(3), size(color_gradient1,1));
    color_gradient2(:,1) = linspace(colors(cContrast2,1), whiteContrast2(1), size(color_gradient2,1));
    color_gradient2(:,2) = linspace(colors(cContrast2,2), whiteContrast2(2), size(color_gradient2,1));
    color_gradient2(:,3) = linspace(colors(cContrast2,3), whiteContrast2(3), size(color_gradient2,1));
else
    cContrast1 = 1;
    cContrast2 = 2;
    color_gradient1 = nan(length(plottedTausAll), 3);
    color_gradient2 = nan(length(plottedTausAll), 3);

    offsetWhite = -0.5;
    whiteContrast1 = max(0, colors(cContrast1, :) + offsetWhite);
    whiteContrast2 = max(0, colors(cContrast2, :) + offsetWhite);
    color_gradient1(:,1) = linspace(colors(cContrast1,1), whiteContrast1(1), size(color_gradient1,1));
    color_gradient1(:,2) = linspace(colors(cContrast1,2), whiteContrast1(2), size(color_gradient1,1));
    color_gradient1(:,3) = linspace(colors(cContrast1,3), whiteContrast1(3), size(color_gradient1,1));
    color_gradient2(:,1) = linspace(colors(cContrast2,1), whiteContrast2(1), size(color_gradient2,1));
    color_gradient2(:,2) = linspace(colors(cContrast2,2), whiteContrast2(2), size(color_gradient2,1));
    color_gradient2(:,3) = linspace(colors(cContrast2,3), whiteContrast2(3), size(color_gradient2,1));
end

legLabels = {};
% functional dependence
for plotTauIdx = 1:length(plottedTaus)
    tauIdx = find(tau_range == plottedTaus(plotTauIdx));
    tau = tau_range(tauIdx);
    
    if(~plotTask1)
    y_mean = squeeze(mean(funcDependence_inc.task2RT_log(:, tauIdx, :) ));
    y_sem = squeeze(std(funcDependence_inc.task2RT_log(:, tauIdx, :) )) / sqrt(numRepetitionsNet);
    else
    y_mean = squeeze(mean(funcDependence_inc.task1RT_log(:, tauIdx, :) ));
    y_sem = squeeze(std(funcDependence_inc.task1RT_log(:, tauIdx, :) )) / sqrt(numRepetitionsNet);    
    end
    
    colorIdx = find(plottedTausAll == tau);
    if(plotSEM)
        errorbar(x, y_mean, y_sem, '-', 'Color', color_gradient1(colorIdx,:), 'LineWidth', lineWidth); hold on;
    else
        plot(x, y_mean, y_sem, '-', 'Color', color_gradient1(colorIdx,:), 'LineWidth', lineWidth); hold on;
    end
    
    legLabels{end+1} = ['Incongruent, p = ' num2str(1-tau)];
end

% independence
for plotTauIdx = 1:length(plottedTaus)
    tauIdx = find(tau_range == plottedTaus(plotTauIdx));
    tau = tau_range(tauIdx);
    
    if(~plotTask1)
    y_mean = squeeze(mean(funcDependence_con.task2RT_log(:, tauIdx, :) ));
    y_sem = squeeze(std(funcDependence_con.task2RT_log(:, tauIdx, :) )) / sqrt(numRepetitionsNet);
    else
    y_mean = squeeze(mean(funcDependence_con.task1RT_log(:, tauIdx, :) ));
    y_sem = squeeze(std(funcDependence_con.task1RT_log(:, tauIdx, :) )) / sqrt(numRepetitionsNet);    
    end
    
    if(tau == 1 & plotBlack)
        lineType = '.';
    else
        lineType = '-';
    end
    
    colorIdx = find(plottedTausAll == tau);
    if(plotSEM)
        errorbar(x, y_mean, y_sem, lineType, 'Color', color_gradient2(colorIdx,:), 'LineWidth', lineWidth); hold on;
    else
        plot(x, y_mean, y_sem, lineType, 'Color', color_gradient2(colorIdx,:), 'LineWidth', lineWidth); hold on;
    end
    
    legLabels{end+1} = ['Congruent, p = ' num2str(1-tau)];
end

hold off;

ylim(ylimit);
if(~plotTask1)
    ylabel({'Reaction Time' 'of Task 1 (s)'}, 'fontSize', fontSize_ylabel-2);
else
    ylabel({'Reaction Time' 'of Task 2 (s)'}, 'fontSize', fontSize_ylabel-2);
end
xlabel('SOA (s)', 'fontSize', fontSize_ylabel-2);
leg = legend(legLabels, 'Location', 'eastoutside');
set(leg, 'FontSize', fontSize_legend);
set(gca, 'FontSize', fontSize_gca);

if(plotBlack)
set(leg, 'TextColor', 'w');
set(leg, 'Color', 'none');
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');
end



%% plot learned representations

plotSettings;

M = nan(numRepetitionsNet, 1, length(tasksToPerform), length(tasksToPerform));
for rep = 1:size(batch_log,2)
    for strengthIdx = 1:size(batch_log,1)
        M(rep, strengthIdx,:,:) = batch_log(strengthIdx, rep).hiddenCorr;
    end
end

M = nanmean(M);
if(size(M,2) > 1)
    M = squeeze(mean(M));
else
    M = squeeze(M);
end

fig2 = figure(2);
set(fig2, 'Position', [100 100 230 260]);
imagesc(M);

map = colormap(copper);
map(1:end, 1) = linspace(0, 1, size(map,1));
map(1:end, 2) = linspace(0, 1, size(map,1));
map(1:end, 3) = linspace(0, 1, size(map,1));
map = map(fliplr(1:size(map,1)),:);
colormap(map);
colorbar('northoutside');
caxis([0 1]);
set(gca, 'FontSize', 12);
set(gca, 'XTick', 1:size(M,1));
set(gca, 'YTick', 1:size(M,2));
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D', 'E'});
set(gca, 'YTickLabel', {'A', 'B', 'C', 'D', 'E'});
ylabel('Tasks','FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Tasks','FontSize', fontSize_xlabel, 'FontName', fontName);
% title({'Correlation of Task Representations' ' at Associative Layer'},'FontSize', fontSize_xlabel-2, 'FontName', fontName);
