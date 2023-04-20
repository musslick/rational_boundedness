%% data settings
clear all;
clc;
close all;

% NOTE: This data file needs to be generated first using Manuscript_Simulation_3_PRP_BCE_tau_fun.m
datafilePrefix = 'PsychReview_Part1_Sim4_PRP_BCE_tau_3P2F_4tasks_3_h100_r';

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
         
         funcDependence.task1RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task1RT_log;
         task1RT_distribution = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task1RT_distribution_log{1};
         funcDependence.task1RT_distribution(rep, tauIdx, :) = task1RT_distribution(:);
         funcDependence.task1Accuracy_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task1Accuracy_log;
         funcDependence.task1Iterations_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task1Iterations_log;
         funcDependence.task2RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task2RT_log;
         task2RT_distribution = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task2RT_distribution_log{1};
         funcDependence.task2RT_distribution(rep, tauIdx, :) = task2RT_distribution(:);
         funcDependence.task2Accuracy_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task2Accuracy_log;
         funcDependence.task2Iterations_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_inc(tauIdx).task2Iterations_log;
         
         independence.task1RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task1RT_log;
         task1RT_distribution = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task1RT_distribution_log{1};
         independence.task1RT_distribution(rep, tauIdx, :) = task1RT_distribution(:);
         independence.task1Accuracy_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task1Accuracy_log;
         independence.task1Iterations_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task1Iterations_log;
         independence.task2RT_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task2RT_log;
         task2RT_distribution = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task2RT_distribution_log{1};
         independence.task2RT_distribution(rep, tauIdx, :) = task2RT_distribution(:);
         independence.task2Accuracy_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task2Accuracy_log;
         independence.task2Iterations_log(rep, tauIdx, :) = batch_log(numHiddenIdx, rep).funcDependence_con(tauIdx).task2Iterations_log;
         
     end
 
end

%% PLOTS - SUMMARY

plotSettings;

fig = figure(1);
set(fig, 'Position', [100 100 240 350]);
x = SOA_range * 0.1;
ylimit = [0 0.7];
plottedTaus = [0.1 0.2 0.5 1];
tauIdx = 1;
lineWidth = 3;
plotSEM = 1;
plotTask1 = 0;
plotDependence = 1;
plotIndependence = 1;

color_gradient1 = nan(length(plottedTaus), 3);
color_gradient2 = nan(length(plottedTaus), 3);

colors = [255 140 41; ... % orange
              31 104 172; ... % blue
              44 155 37; ... % green
              255 0 0; ... % red
              0     0   0  ; ... % black
            142 142 142; ... % grey
            253 120 21; ... % orange
            233 100 1; ... % orange
            255 255 255] / 255; % white 
% colors(cContrast1, :) = max(0,colors(cContrast1, :)-0.2);
whiteContrast1 = min(1, colors(cContrast1, :) + 0.4);
whiteContrast2 = min(1, colors(cContrast2, :) + 0.4);
color_gradient1(:,1) = linspace(colors(cContrast1,1), whiteContrast1(1), size(color_gradient1,1));
color_gradient1(:,2) = linspace(colors(cContrast1,2), whiteContrast1(2), size(color_gradient1,1));
color_gradient1(:,3) = linspace(colors(cContrast1,3), whiteContrast1(3), size(color_gradient1,1));
color_gradient2(:,1) = linspace(colors(cContrast2,1), whiteContrast2(1), size(color_gradient2,1));
color_gradient2(:,2) = linspace(colors(cContrast2,2), whiteContrast2(2), size(color_gradient2,1));
color_gradient2(:,3) = linspace(colors(cContrast2,3), whiteContrast2(3), size(color_gradient2,1));

% for each simulation, compute RT quantiles
task1_con_RT= mean(independence.task1RT_log(:,tauIdx,1));
task1_inc_RT= mean(funcDependence.task1RT_log(:,tauIdx,1));
task2_con_RT= mean(independence.task2RT_log(:,tauIdx,1));
task2_inc_RT= mean(funcDependence.task2RT_log(:,tauIdx,1));

task1_con_RT_sem= std(independence.task1RT_log(:,tauIdx,1))/sqrt(numRepetitionsNet);
task1_inc_RT_sem= std(funcDependence.task1RT_log(:,tauIdx,1))/sqrt(numRepetitionsNet);
task2_con_RT_sem= std(independence.task2RT_log(:,tauIdx,1))/sqrt(numRepetitionsNet);
task2_inc_RT_sem= std(funcDependence.task2RT_log(:,tauIdx,1))/sqrt(numRepetitionsNet);

x = [1 2];
y_mean = [task1_con_RT task1_inc_RT];
y_sem = [task1_con_RT_sem task1_inc_RT_sem];
bar(x(1), y_mean(1), 'FaceColor', colors(cContrast3, :)); hold on;
bar(x(2), y_mean(2), 'FaceColor', colors(cContrast1, :)); hold on;
errorbar(x, y_mean, y_sem, '.k');

x = [3 4];
y_mean = [task2_con_RT task2_inc_RT];
y_sem = [task2_con_RT_sem task2_inc_RT_sem];
bar(x(1), y_mean(1), 'FaceColor', colors(cContrast3, :)); hold on;
bar(x(2), y_mean(2), 'FaceColor', colors(cContrast1, :)); hold on;
errorbar(x, y_mean, y_sem, '.k');

hold off;

%% PLOTS - HOMMEL 1998

plotSettings;

fig = figure(1);
set(fig, 'Position', [100 100 240 350]);
x = SOA_range * 0.1;
ylimit = [0 0.7];
plottedTaus = [0.1 0.2 0.5 1];
tauIdx = 1;
lineWidth = 3;
plotSEM = 1;
plotTask1 = 0;
plotDependence = 1;
plotIndependence = 1;

color_gradient1 = nan(length(plottedTaus), 3);
color_gradient2 = nan(length(plottedTaus), 3);

colors = [255 140 41; ... % orange
              31 104 172; ... % blue
              44 155 37; ... % green
              255 0 0; ... % red
              0     0   0  ; ... % black
            142 142 142; ... % grey
            253 120 21; ... % orange
            233 100 1; ... % orange
            255 255 255] / 255; % white 
% colors(cContrast1, :) = max(0,colors(cContrast1, :)-0.2);
whiteContrast1 = min(1, colors(cContrast1, :) + 0.4);
whiteContrast2 = min(1, colors(cContrast2, :) + 0.4);
color_gradient1(:,1) = linspace(colors(cContrast1,1), whiteContrast1(1), size(color_gradient1,1));
color_gradient1(:,2) = linspace(colors(cContrast1,2), whiteContrast1(2), size(color_gradient1,1));
color_gradient1(:,3) = linspace(colors(cContrast1,3), whiteContrast1(3), size(color_gradient1,1));
color_gradient2(:,1) = linspace(colors(cContrast2,1), whiteContrast2(1), size(color_gradient2,1));
color_gradient2(:,2) = linspace(colors(cContrast2,2), whiteContrast2(2), size(color_gradient2,1));
color_gradient2(:,3) = linspace(colors(cContrast2,3), whiteContrast2(3), size(color_gradient2,1));

% for each simulation, compute RT quantiles
task1_con_RT= nan(numRepetitionsNet, 5);
task1_inc_RT= nan(numRepetitionsNet, 5);
task2_con_RT= nan(numRepetitionsNet, 5);
task2_inc_RT= nan(numRepetitionsNet, 5);

p = [0.2:0.2:1];
for rep = 1:numRepetitionsNet
    task1_con_RT(rep,:) = quantile(squeeze(independence.task1RT_distribution(rep, tauIdx, :)), p);
    task1_inc_RT(rep,:) = quantile(squeeze(funcDependence.task1RT_distribution(rep, tauIdx, :)), p);
    task2_con_RT(rep,:) = quantile(squeeze(independence.task2RT_distribution(rep, tauIdx, :)), p);
    task2_inc_RT(rep,:) = quantile(squeeze(funcDependence.task2RT_distribution(rep, tauIdx, :)), p);
end
  

% TASK 1 RT CON
x = 1:length(p);
y_mean = mean(task1_con_RT);
y_sem = std(task1_con_RT,[],1)/sqrt(numRepetitionsNet);
errorbar(x, y_mean, y_sem, '-', 'Color', colors(cContrast3,:), 'LineWidth', lineWidth); hold on;

% TASK 1 RT INC
y_mean = mean(task1_inc_RT);
y_sem = std(task1_inc_RT,[],1)/sqrt(numRepetitionsNet);
errorbar(x, y_mean, y_sem, '--', 'Color', colors(cContrast3,:), 'LineWidth', lineWidth); hold on

% TASK 2 RT CON
x = 1:length(p);
y_mean = mean(task2_con_RT);
y_sem = std(task2_con_RT,[],1)/sqrt(numRepetitionsNet);
errorbar(x, y_mean, y_sem, '-', 'Color', colors(cContrast2,:), 'LineWidth', lineWidth); hold on;

% TASK 2 RT INC
y_mean = mean(task2_inc_RT);
y_sem = std(task2_inc_RT,[],1)/sqrt(numRepetitionsNet);
errorbar(x, y_mean, y_sem, '--', 'Color', colors(cContrast2,:), 'LineWidth', lineWidth); hold on

hold off


 
