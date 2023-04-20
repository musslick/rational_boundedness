%% data settings
clear all;
clc;
close all;

datafilePrefix = 'PsychReview_Part1_Sim3_TaskSwitchingUnivalent_3P3F_5tasks_15_h100_r';

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
            numTausTested = size(batch_log,3);
            batch_log_tmp = repmat(batch_log(1,1), numHiddenUnitsTested, numRepetitions, numTausTested);
        end
        batch_log_tmp(:, i, :) = batch_log(:, rep, :);
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

for tauIdx = 1:length(tau_range)

     for rep =1:numRepetitions
         
         % structural dependence
         structDependence.ER_univalent_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_univalent_rep);
         structDependence.ER_univalent_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_univalent_switch);
         structDependence.ER_bivalent_con_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_bivalent_con_rep);
         structDependence.ER_bivalent_con_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_bivalent_con_switch);
         structDependence.ER_bivalent_inc_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_bivalent_inc_rep);
         structDependence.ER_bivalent_inc_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_bivalent_inc_switch);
         structDependence.ER_bivalent_con_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_bivalent_con_switchCost);
         structDependence.ER_bivalent_inc_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).ER_bivalent_inc_switchCost);
         
         structDependence.RT_univalent_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_univalent_rep);
         structDependence.RT_univalent_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_univalent_switch);
         structDependence.RT_bivalent_con_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_bivalent_con_rep);
         structDependence.RT_bivalent_con_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_bivalent_con_switch);
         structDependence.RT_bivalent_inc_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_bivalent_inc_rep);
         structDependence.RT_bivalent_inc_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_bivalent_inc_switch);
         structDependence.RT_bivalent_con_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_bivalent_con_switchCost);
         structDependence.RT_bivalent_inc_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).structDependence(tauIdx).RT_bivalent_inc_switchCost);
         
         % functional dependence

         funcDependence.ER_univalent_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_univalent_rep);
         funcDependence.ER_univalent_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_univalent_switch);
         funcDependence.ER_bivalent_con_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_bivalent_con_rep);
         funcDependence.ER_bivalent_con_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_bivalent_con_switch);
         funcDependence.ER_bivalent_inc_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_bivalent_inc_rep);
         funcDependence.ER_bivalent_inc_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_bivalent_inc_switch);
         funcDependence.ER_bivalent_con_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_bivalent_con_switchCost);
         funcDependence.ER_bivalent_inc_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).ER_bivalent_inc_switchCost);
         
         funcDependence.RT_univalent_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_univalent_rep);
         funcDependence.RT_univalent_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_univalent_switch);
         funcDependence.RT_bivalent_con_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_bivalent_con_rep);
         funcDependence.RT_bivalent_con_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_bivalent_con_switch);
         funcDependence.RT_bivalent_inc_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_bivalent_inc_rep);
         funcDependence.RT_bivalent_inc_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_bivalent_inc_switch);
         funcDependence.RT_bivalent_con_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_bivalent_con_switchCost);
         funcDependence.RT_bivalent_inc_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).funcDependence(tauIdx).RT_bivalent_inc_switchCost);
         
         % independence

         independence.ER_univalent_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_univalent_rep);
         independence.ER_univalent_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_univalent_switch);
         independence.ER_bivalent_con_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_bivalent_con_rep);
         independence.ER_bivalent_con_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_bivalent_con_switch);
         independence.ER_bivalent_inc_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_bivalent_inc_rep);
         independence.ER_bivalent_inc_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_bivalent_inc_switch);
         independence.ER_bivalent_con_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_bivalent_con_switchCost);
         independence.ER_bivalent_inc_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).ER_bivalent_inc_switchCost);
         
         independence.RT_univalent_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_univalent_rep);
         independence.RT_univalent_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_univalent_switch);
         independence.RT_bivalent_con_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_bivalent_con_rep);
         independence.RT_bivalent_con_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_bivalent_con_switch);
         independence.RT_bivalent_inc_rep(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_bivalent_inc_rep);
         independence.RT_bivalent_inc_switch(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_bivalent_inc_switch);
         independence.RT_bivalent_con_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_bivalent_con_switchCost);
         independence.RT_bivalent_inc_switchCost(rep, tauIdx) = mean(batch_log(numHiddenIdx, rep).independence(tauIdx).RT_bivalent_inc_switchCost);
         
     end
 
end

 
%% PLOTS SUMMARY FOR SPECIFIC TAU

% parameters
plotBlack = 0;
showTitle = 1;
plotRT = 1;
ERInPercent = 1;
plotSEM = 0;
selectedTau = 0.1;
tauIdx = find(tau_range == selectedTau);

singleLine = 0;

x = [1, 2, 3];
if(~plotRT)
    ylimit = [0 0.04];
else
    ylimit = [0 0.5];
end
lineWidth = 1;
markerSize = 8;
if(plotBlack)
    color_struct = [254 209 103]/255;
    color_func = [252 9 27]/255;
    color_ind = [47 205 58]/255;
else
    color_struct = colors(cContrast3,:);
    color_func = colors(cContrast1,:);
    color_ind = colors(cContrast2,:);
end


if(~plotRT)
    if(ERInPercent)
       ylimit = ylimit * 100;
       scale = 100;
    else
       scale = 1; 
    end
else
    scale = 1;
end

% ERROR RATES

% structural dependence
fig3 = figure(3);
set(fig3, 'Position', [100 400 300 200]);
if(~plotRT)
    y_switch_mean = [mean(structDependence.ER_univalent_switch(:, tauIdx)) mean(structDependence.ER_bivalent_con_switch(:, tauIdx)) mean(structDependence.ER_bivalent_inc_switch(:, tauIdx))];
    y_switch_sem = [std(structDependence.ER_univalent_switch(:, tauIdx)) std(structDependence.ER_bivalent_con_switch(:, tauIdx)) std(structDependence.ER_bivalent_inc_switch(:, tauIdx))] / sqrt(numRepetitions);
    y_rep_mean = [mean(structDependence.ER_univalent_rep(:, tauIdx)) mean(structDependence.ER_bivalent_con_rep(:, tauIdx)) mean(structDependence.ER_bivalent_inc_rep(:, tauIdx))];
    y_rep_sem = [std(structDependence.ER_univalent_rep(:, tauIdx)) std(structDependence.ER_bivalent_con_rep(:, tauIdx)) std(structDependence.ER_bivalent_inc_rep(:, tauIdx))] / sqrt(numRepetitions);
else
    y_switch_mean = [mean(structDependence.RT_univalent_switch(:, tauIdx)) mean(structDependence.RT_bivalent_con_switch(:, tauIdx)) mean(structDependence.RT_bivalent_inc_switch(:, tauIdx))];
    y_switch_sem = [std(structDependence.RT_univalent_switch(:, tauIdx)) std(structDependence.RT_bivalent_con_switch(:, tauIdx)) std(structDependence.RT_bivalent_inc_switch(:, tauIdx))] / sqrt(numRepetitions);
    y_rep_mean = [mean(structDependence.RT_univalent_rep(:, tauIdx)) mean(structDependence.RT_bivalent_con_rep(:, tauIdx)) mean(structDependence.RT_bivalent_inc_rep(:, tauIdx))];
    y_rep_sem = [std(structDependence.RT_univalent_rep(:, tauIdx)) std(structDependence.RT_bivalent_con_rep(:, tauIdx)) std(structDependence.RT_bivalent_inc_rep(:, tauIdx))] / sqrt(numRepetitions);
end


if(plotSEM)
    errorbar(x, y_switch_mean*scale, y_switch_sem*scale, 's',  'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold on;
    errorbar(x, y_rep_mean*scale, y_rep_sem*scale, '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold off;
else
    plot(x, y_switch_mean*scale, 's',  'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold on;
    plot(x, y_rep_mean*scale, '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold off;    
end
ylim(ylimit);
xlim([0.8 3.2]);

if(plotRT)
    ylabel('Reaction Time (s)', 'fontSize', fontSize_ylabel-2);
else
    if(ERInPercent)
        ylabel('Error Rate (%)', 'fontSize', fontSize_ylabel-2);
    else
        ylabel('Error Rate', 'fontSize', fontSize_ylabel-2);
    end
end
set(gca, 'XTick', x);
set(gca, 'XTickLabel', {'Univalent', 'Congruent', 'Incongruent'});
if(showTitle)
    title({'Structural Dependence'}, 'FontSize', fontSize_title);
end
if(~singleLine)
    leg =  legend({'Task Switch', 'Task Repetition'}, 'Location', 'south');
    set(leg, 'FontSize', fontSize_legend);
else
    leg =  legend({'Repeat'}, 'Location', 'south');
end
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

% functional dependence
fig1 = figure(1);
set(fig1, 'Position', [100 400 300 200]);
if(~plotRT)
    y_switch_mean = [mean(funcDependence.ER_univalent_switch(:, tauIdx)), mean(funcDependence.ER_bivalent_con_switch(:, tauIdx)) mean(funcDependence.ER_bivalent_inc_switch(:, tauIdx))];
    y_switch_sem = [std(funcDependence.ER_univalent_switch(:, tauIdx)) std(funcDependence.ER_bivalent_con_switch(:, tauIdx)) std(funcDependence.ER_bivalent_inc_switch(:, tauIdx))] / sqrt(numRepetitions);
    y_rep_mean = [mean(funcDependence.ER_univalent_rep(:, tauIdx)) mean(funcDependence.ER_bivalent_con_rep(:, tauIdx)) mean(funcDependence.ER_bivalent_inc_rep(:, tauIdx))];
    y_rep_sem = [std(funcDependence.ER_univalent_rep(:, tauIdx)) std(funcDependence.ER_bivalent_con_rep(:, tauIdx)) std(funcDependence.ER_bivalent_inc_rep(:, tauIdx))] / sqrt(numRepetitions);
else
    y_switch_mean = [mean(funcDependence.RT_univalent_switch(:, tauIdx)) mean(funcDependence.RT_bivalent_con_switch(:, tauIdx)) mean(funcDependence.RT_bivalent_inc_switch(:, tauIdx))];
    y_switch_sem = [std(funcDependence.RT_univalent_switch(:, tauIdx)) std(funcDependence.RT_bivalent_con_switch(:, tauIdx)) std(funcDependence.RT_bivalent_inc_switch(:, tauIdx))] / sqrt(numRepetitions);
    y_rep_mean = [mean(funcDependence.RT_univalent_rep(:, tauIdx)) mean(funcDependence.RT_bivalent_con_rep(:, tauIdx)) mean(funcDependence.RT_bivalent_inc_rep(:, tauIdx))];
    y_rep_sem = [std(funcDependence.RT_univalent_rep(:, tauIdx)) std(funcDependence.RT_bivalent_con_rep(:, tauIdx)) std(funcDependence.RT_bivalent_inc_rep(:, tauIdx))] / sqrt(numRepetitions);
end

if(plotSEM)
    errorbar(x, y_switch_mean*scale, y_switch_sem*scale, 's',  'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold on;
    errorbar(x, y_rep_mean*scale, y_rep_sem*scale, '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold off;
else
    plot(x, y_switch_mean*scale, 's',  'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_func); hold on;
    plot(x, y_rep_mean*scale, '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_func); hold off;    
end

ylim(ylimit);
xlim([0.8 3.2]);

if(plotRT)
    ylabel('Reaction Time (s)', 'fontSize', fontSize_ylabel-2);
else
    if(ERInPercent)
        ylabel('Error Rate (%)', 'fontSize', fontSize_ylabel-2);
    else
        ylabel('Error Rate', 'fontSize', fontSize_ylabel-2);
    end
end
set(gca, 'XTick', x);
set(gca, 'XTickLabel', {'Univalent', 'Congruent', 'Incongruent'});
if(showTitle)
    title({'Functional Dependence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Task Switch', 'Task Repetition'}, 'Location', 'south');
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

% independence dependence
fig2 = figure(2);
set(fig2, 'Position', [100 400 300 200]);
if(~plotRT)
    y_switch_mean = [mean(independence.ER_univalent_switch(:, tauIdx)) mean(independence.ER_bivalent_con_switch(:, tauIdx)) mean(independence.ER_bivalent_inc_switch(:, tauIdx))];
    y_switch_sem = [std(independence.ER_univalent_switch(:, tauIdx)) std(independence.ER_bivalent_con_switch(:, tauIdx)) std(independence.ER_bivalent_inc_switch(:, tauIdx))] / sqrt(numRepetitions);
    y_rep_mean = [mean(independence.ER_univalent_rep(:, tauIdx)) mean(independence.ER_bivalent_con_rep(:, tauIdx)) mean(independence.ER_bivalent_inc_rep(:, tauIdx))];
    y_rep_sem = [std(independence.ER_univalent_rep(:, tauIdx)) std(independence.ER_bivalent_con_rep(:, tauIdx)) std(independence.ER_bivalent_inc_rep(:, tauIdx))] / sqrt(numRepetitions);
else
    y_switch_mean = [mean(independence.RT_univalent_switch(:, tauIdx)) mean(independence.RT_bivalent_con_switch(:, tauIdx)) mean(independence.RT_bivalent_inc_switch(:, tauIdx))];
    y_switch_sem = [std(independence.RT_univalent_switch(:, tauIdx)) std(independence.RT_bivalent_con_switch(:, tauIdx)) std(independence.RT_bivalent_inc_switch(:, tauIdx))] / sqrt(numRepetitions);
    y_rep_mean = [mean(independence.RT_univalent_rep(:, tauIdx)) mean(independence.RT_bivalent_con_rep(:, tauIdx)) mean(independence.RT_bivalent_inc_rep(:, tauIdx))];
    y_rep_sem = [std(independence.RT_univalent_rep(:, tauIdx)) std(independence.RT_bivalent_con_rep(:, tauIdx)) std(independence.RT_bivalent_inc_rep(:, tauIdx))] / sqrt(numRepetitions);
end

if(plotSEM)
    errorbar(x, y_switch_mean*scale, y_switch_sem*scale, 's',  'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold on;
    errorbar(x, y_rep_mean*scale, y_rep_sem*scale, '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_struct); hold off;
else
    plot(x, y_switch_mean*scale, 's',  'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_ind); hold on;
    plot(x, y_rep_mean*scale, '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', color_ind); hold off;    
end
ylim(ylimit);
xlim([0.8 3.2]);

if(plotRT)
    ylabel('Reaction Time (s)', 'fontSize', fontSize_ylabel-2);
else
    if(ERInPercent)
        ylabel('Error Rate (%)', 'fontSize', fontSize_ylabel-2);
    else
        ylabel('Error Rate', 'fontSize', fontSize_ylabel-2);
    end
end
set(gca, 'XTick', x);
set(gca, 'XTickLabel', {'Univalent', 'Congruent', 'Incongruent'});
if(showTitle)
    title({'Independence'}, 'FontSize', fontSize_title);
end
leg =  legend({'Task Switch', 'Task Repetition'}, 'Location', 'south');
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

 
%% PLOTS - SWITCH COSTS

plotBlack = 0;
showTitle = 1;
plotRT = 1;

fig = figure(1);
set(fig, 'Position', [100 100 350 200]);
x = 1-tau_range;
ylimit = [0 0.12];
xlimit = [0 1];
lineWidth = 3;

if(~plotRT)
    % ERROR RATES

    % structural dependence
    y_switch_struct_mean = mean(structDependence.ER_bivalent_inc_switchCost);
    y_switch_struct_sem = std(structDependence.ER_bivalent_inc_switchCost) / sqrt(numRepetitions);

    % functional dependence
    y_switch_func_mean = mean(funcDependence.ER_bivalent_inc_switchCost);
    y_switch_func_sem = std(funcDependence.ER_bivalent_inc_switchCost) / sqrt(numRepetitions);

    % independence
    y_switch_indep_mean = mean(independence.ER_bivalent_inc_switchCost);
    y_switch_indep_sem = std(independence.ER_bivalent_inc_switchCost) / sqrt(numRepetitions);
else
    % RTs

    % structural dependence
    y_switch_struct_mean = mean(structDependence.RT_bivalent_inc_switchCost);
    y_switch_struct_sem = std(structDependence.RT_bivalent_inc_switchCost) / sqrt(numRepetitions);

    % functional dependence
    y_switch_func_mean = mean(funcDependence.RT_bivalent_inc_switchCost);
    y_switch_func_sem = std(funcDependence.RT_bivalent_inc_switchCost) / sqrt(numRepetitions);

    % independence
    y_switch_indep_mean = mean(independence.RT_bivalent_inc_switchCost);
    y_switch_indep_sem = std(independence.RT_bivalent_inc_switchCost) / sqrt(numRepetitions);
end

errorbar(x, y_switch_struct_mean, y_switch_struct_sem, '-', 'Color', colors(cContrast3,:), 'LineWidth', lineWidth); hold on;
errorbar(x, y_switch_func_mean, y_switch_func_sem, '-', 'Color', colors(cContrast1,:), 'LineWidth', lineWidth); 
errorbar(x, y_switch_indep_mean, y_switch_indep_sem, '-', 'Color', colors(cContrast2,:), 'LineWidth', lineWidth); hold off;
hold on;
ylim(ylimit);

if(~plotRT)
    ylabel({'Switch Costs' 'in Error Rate (%)'}, 'fontSize', fontSize_ylabel-2);
else
    ylabel({'Switch Costs' 'in Reaction Time (s)'}, 'fontSize', fontSize_ylabel-2);
end

xlabel('Persistence p', 'fontSize', fontSize_ylabel-2);
% set(gca, 'XTick', tau_range);
% set(gca, 'XTickLabel', tau_range);
leg = legend('Structural Dependence', 'Functional Dependence', 'Independence', 'Location', 'northwest');
set(leg, 'FontSize', fontSize_legend);
set(gca, 'FontSize', fontSize_gca);
hold off;


%% stats

