%% data settings
clear all;
clc;
close all;

datafilePrefix = 'PsychReview_Part1_Sim3_TaskSwitching_3P3F_5tasks_11_h100_r';

Nhidden = 100;
correlation_threshold = 0.5;
plotBlack = 0;

% load data
folder = 'logfiles/Part1/';
files = dir(folder);

% get valid file names
validFileNames = {};
for i =1:length(files)
    % check if this is a desired data file
    if(~isempty(strfind(files(i).name, datafilePrefix)))
        validFileNames{end+1} = files(i).name;
    end
end
if(~isempty(validFileNames))
    numRepetitions = length(validFileNames);
    % load every valid file
    for i = 1:length(validFileNames);
        disp(['loading ' validFileNames{i} '...']);
        load(strcat(folder, validFileNames{i}));
        % initial setup of batch_log_tmp
        if(i == 1)
            numHiddenUnitsTested = size(batch_log,1);
            batch_log_tmp = repmat(batch_log(1,1), numHiddenUnitsTested, numRepetitions);
        end
        batch_log_tmp(:, i) = batch_log(:, rep);
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

 for rep =1:numRepetitions
        
     structDependence.ER_repCon(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.ER_repCon);
     structDependence.ER_switchCon(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.ER_switchCon);
     structDependence.ER_repInc(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.ER_repInc);
     structDependence.ER_switchInc(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.ER_switchInc);
     
     funcDependence.ER_repCon(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.ER_repCon);
     funcDependence.ER_switchCon(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.ER_switchCon);
     funcDependence.ER_repInc(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.ER_repInc);
     funcDependence.ER_switchInc(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.ER_switchInc);
   
     independence.ER_repCon(rep) = mean(batch_log(numHiddenIdx, rep).independence.ER_repCon);
     independence.ER_switchCon(rep) = mean(batch_log(numHiddenIdx, rep).independence.ER_switchCon);
     independence.ER_repInc(rep) = mean(batch_log(numHiddenIdx, rep).independence.ER_repInc);
     independence.ER_switchInc(rep) = mean(batch_log(numHiddenIdx, rep).independence.ER_switchInc);
     
     structDependence.RT_repCon(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.RT_repCon);
     structDependence.RT_switchCon(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.RT_switchCon);
     structDependence.RT_repInc(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.RT_repInc);
     structDependence.RT_switchInc(rep) = mean(batch_log(numHiddenIdx, rep).structDependence.RT_switchInc);
     
     funcDependence.RT_repCon(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.RT_repCon);
     funcDependence.RT_switchCon(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.RT_switchCon);
     funcDependence.RT_repInc(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.RT_repInc);
     funcDependence.RT_switchInc(rep) = mean(batch_log(numHiddenIdx, rep).funcDependence.RT_switchInc);
   
     independence.RT_repCon(rep) = mean(batch_log(numHiddenIdx, rep).independence.RT_repCon);
     independence.RT_switchCon(rep) = mean(batch_log(numHiddenIdx, rep).independence.RT_switchCon);
     independence.RT_repInc(rep) = mean(batch_log(numHiddenIdx, rep).independence.RT_repInc);
     independence.RT_switchInc(rep) = mean(batch_log(numHiddenIdx, rep).independence.RT_switchInc);
 end
 
%% PLOTS - Interaction Transition x Congruency
singleLine = 1;

x = [1, 2];
ylimit = [0 0.5];
% xlimit = [0.5 2.5];
lineWidth = 3;
showTitle = 0;
ERInPercent = 1;
color_struct = [254 209 103]/255;
color_func = [252 9 27]/255;
color_ind = [47 205 58]/255;


if(ERInPercent)
   ylimit = ylimit * 100;
   scale = 100;
else
   scale = 1; 
end

% ERROR RATES

% structural dependence
fig3 = figure(3);
set(fig3, 'Position', [100 400 300 200]);
y_switch_mean = [mean(structDependence.ER_switchCon) mean(structDependence.ER_switchInc)];
y_switch_sem = [std(structDependence.ER_switchCon) std(structDependence.ER_switchInc)] / sqrt(numRepetitions);
y_rep_mean = [mean(structDependence.ER_repCon) mean(structDependence.ER_repInc)];
y_rep_sem = [std(structDependence.ER_repCon) std(structDependence.ER_repInc)] / sqrt(numRepetitions);

if(~singleLine)
errorbar(x, y_switch_mean*scale, y_switch_sem*scale, '-',  'LineWidth', lineWidth, 'Color', color_struct); hold on;
end
errorbar(x, y_rep_mean*scale, y_rep_sem*scale, '--', 'LineWidth', lineWidth, 'Color', color_struct); hold off;
ylim(ylimit);
xlim([0.8 2.2]);

if(ERInPercent)
    ylabel('Error Rate (%)', 'fontSize', fontSize_ylabel-2);
else
    ylabel('Error Rate', 'fontSize', fontSize_ylabel-2);
end
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Congruent', 'Incongruent'});
if(showTitle)
    title({'Structural Dependence'}, 'FontSize', fontSize_title);
end
if(~singleLine)
    leg =  legend({'Switch', 'Repeat'}, 'Location', 'north');
else
    leg =  legend({'Repeat'}, 'Location', 'north');
end
set(gca, 'FontSize', fontSize_gca);

set(leg, 'TextColor', 'w');
set(leg, 'Color', 'none');
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

% functional dependence
fig1 = figure(1);
set(fig1, 'Position', [100 400 300 200]);
y_switch_mean = [mean(funcDependence.ER_switchCon) mean(funcDependence.ER_switchInc)];
y_switch_sem = [std(funcDependence.ER_switchCon) std(funcDependence.ER_switchInc)] / sqrt(numRepetitions);
y_rep_mean = [mean(funcDependence.ER_repCon) mean(funcDependence.ER_repInc)];
y_rep_sem = [std(funcDependence.ER_repCon) std(funcDependence.ER_repInc)] / sqrt(numRepetitions);

errorbar(x, y_switch_mean*scale, y_switch_sem*scale, '-', 'LineWidth', lineWidth, 'Color', color_func); hold on;
errorbar(x, y_rep_mean*scale, y_rep_sem*scale, '--', 'LineWidth', lineWidth, 'Color', color_func); hold off;
ylim(ylimit);
xlim([0.8 2.2]);

if(ERInPercent)
    ylabel('Error Rate (%)', 'fontSize', fontSize_ylabel-2);
else
    ylabel('Error Rate', 'fontSize', fontSize_ylabel-2);
end
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Congruent', 'Incongruent'});
if(showTitle)
    title({'Functional Dependence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Switch', 'Repeat'}, 'Location', 'north');
set(gca, 'FontSize', fontSize_gca);

set(leg, 'TextColor', 'w');
set(leg, 'Color', 'none');
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

% independence dependence
fig2 = figure(2);
set(fig2, 'Position', [100 400 300 200]);
y_switch_mean = [mean(independence.ER_switchCon) mean(independence.ER_switchInc)];
y_switch_sem = [std(independence.ER_switchCon) std(independence.ER_switchInc)] / sqrt(numRepetitions);
y_rep_mean = [mean(independence.ER_repCon) mean(independence.ER_repInc)];
y_rep_sem = [std(independence.ER_repCon) std(independence.ER_repInc)] / sqrt(numRepetitions);

errorbar(x, y_switch_mean*scale, y_switch_sem*scale, '-', 'LineWidth', lineWidth, 'Color', color_ind); hold on;
errorbar(x, y_rep_mean*scale, y_rep_sem*scale, '--', 'LineWidth', lineWidth, 'Color', color_ind); hold off;
ylim(ylimit);
xlim([0.8 2.2]);

if(ERInPercent)
    ylabel('Error Rate (%)', 'fontSize', fontSize_ylabel-2);
else
    ylabel('Error Rate', 'fontSize', fontSize_ylabel-2);
end
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Congruent', 'Incongruent'});
if(showTitle)
    title({'Independence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Switch', 'Repeat'}, 'Location', 'north');
set(gca, 'FontSize', fontSize_gca);

set(leg, 'TextColor', 'w');
set(leg, 'Color', 'none');
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'xColor', 'w');
set(gca, 'yColor', 'w');
set(gca, 'zColor', 'w');

%% Reaction Time (s)S

fig = figure(2);
set(fig, 'Position', [100 100 900 200]);
x = [1, 2];
ylimit = [0 0.5];
% xlimit = [0.5 2.5];
lineWidth = 3;

% structural dependence
subplot(1, 3, 1);
y_switch_mean = [mean(structDependence.RT_switchCon) mean(structDependence.RT_switchInc)];
y_switch_sem = [std(structDependence.RT_switchCon) std(structDependence.RT_switchInc)] / sqrt(numRepetitions);
y_rep_mean = [mean(structDependence.RT_repCon) mean(structDependence.RT_repInc)];
y_rep_sem = [std(structDependence.RT_repCon) std(structDependence.RT_repInc)] / sqrt(numRepetitions);

errorbar(x, y_switch_mean, y_switch_sem, '-k', 'LineWidth', lineWidth); hold on;
errorbar(x, y_rep_mean, y_rep_sem, '--k', 'LineWidth', lineWidth); hold off;
ylim(ylimit);
xlim([0.8 2.2]);

ylabel('Reaction Time (s)', 'fontSize', fontSize_ylabel-2);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Congruent', 'Incongruent'});
if(showTitle)
    title({'Structural Dependence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Switch', 'Repeat'}, 'Location', 'north');
set(gca, 'FontSize', fontSize_gca);

% functional dependence
subplot(1, 3, 2);
y_switch_mean = [mean(funcDependence.RT_switchCon) mean(funcDependence.RT_switchInc)];
y_switch_sem = [std(funcDependence.RT_switchCon) std(funcDependence.RT_switchInc)] / sqrt(numRepetitions);
y_rep_mean = [mean(funcDependence.RT_repCon) mean(funcDependence.RT_repInc)];
y_rep_sem = [std(funcDependence.RT_repCon) std(funcDependence.RT_repInc)] / sqrt(numRepetitions);

errorbar(x, y_switch_mean, y_switch_sem, '-k', 'LineWidth', lineWidth); hold on;
errorbar(x, y_rep_mean, y_rep_sem, '--k', 'LineWidth', lineWidth); hold off;
ylim(ylimit);
xlim([0.8 2.2]);

ylabel('Reaction Time (s)', 'fontSize', fontSize_ylabel-2);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Congruent', 'Incongruent'});
if(showTitle)
    title({'Functional Dependence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Switch', 'Repeat'}, 'Location', 'north');
set(gca, 'FontSize', fontSize_gca);

% independence dependence
subplot(1, 3, 3);
y_switch_mean = [mean(independence.RT_switchCon) mean(independence.RT_switchInc)];
y_switch_sem = [std(independence.RT_switchCon) std(independence.RT_switchInc)] / sqrt(numRepetitions);
y_rep_mean = [mean(independence.RT_repCon) mean(independence.RT_repInc)];
y_rep_sem = [std(independence.RT_repCon) std(independence.RT_repInc)] / sqrt(numRepetitions);

errorbar(x, y_switch_mean, y_switch_sem, '-k', 'LineWidth', lineWidth); hold on;
errorbar(x, y_rep_mean, y_rep_sem, '--k', 'LineWidth', lineWidth); hold off;
ylim(ylimit);
xlim([0.8 2.2]);

ylabel('Reaction Time (s)', 'fontSize', fontSize_ylabel-2);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Congruent', 'Incongruent'});
if(showTitle)
    title({'Independence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Switch', 'Repeat'}, 'Location', 'north');
set(gca, 'FontSize', fontSize_gca);



%% stats

% use variance from task switching study to fit model
