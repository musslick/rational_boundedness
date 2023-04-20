
% PAPER RELEVANT

clear all;
clc;
load('logfiles/Part2/PsychReview_Part2_Exp1_modifiedStroop_3P4F_4tasks_13_h100_ortho0.mat');

% PREPARE DATA

numTaskStrengths = length(interferenceTaskStrength);

LCA_funcDependence = nan(numRepetitions, numTaskStrengths);

% to account for excluded subjects
numRepetitions = 21;

maxSize = 0;
for rep = 1:numRepetitions
    for taskStrengthIdx = 1:length(interferenceTaskStrength)
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).Independence.A_B)]);
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).Independence.A_B_1)]);
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).Independence.min_A_B)]);

        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).funcDependence.A_B)]);
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).funcDependence.A_B_1)]);
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).funcDependence.min_A_B)]);
    end
end

A_B_matrix_type0_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
min_A_B_matrix_type0_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
A_B_1_matrix_type0_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
C_CNLP{taskStrengthIdx}.data = nan(numRepetitions, maxSize);

A_B_matrix_type4_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
min_A_B_matrix_type4_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
A_B_1_matrix_type4_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
C_CNWP{taskStrengthIdx}.data = nan(numRepetitions, maxSize);

for rep = 1:numRepetitions
    
   for taskStrengthIdx = 1:length(interferenceTaskStrength)
       LCA_funcDependence(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy;
       LCA_funcDependence_inc(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy_inc;
       LCA_funcDependence_con(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy_con;
       
       LCA_Independence(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy;
       LCA_Independence_inc(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy_inc;
       LCA_Independence_con(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy_con;
       
       if(isfield(batch_log, 'taskA'))
           LCA_taskA(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskA.LCA_Accuracy;
           LCA_taskA_inc(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskA.LCA_Accuracy_inc;
           LCA_taskA_con(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskA.LCA_Accuracy_con;
       end
       if(isfield(batch_log, 'taskB'))
           LCA_taskB(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskB.LCA_Accuracy;
           LCA_taskB_inc(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskB.LCA_Accuracy_inc;
           LCA_taskB_con(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskB.LCA_Accuracy_con;
       end
       if(isfield(batch_log, 'taskC'))
           LCA_taskC(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy;
           LCA_taskC_inc(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy_inc;
           LCA_taskC_con(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy_con;
       end

        % add townsend data (independence)
        numElements = length(batch_log(taskStrengthIdx, rep).Independence.A_B);
        A_B_matrix_type0_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.A_B;

        numElements = length(batch_log(taskStrengthIdx, rep).Independence.A_B_1);
        A_B_1_matrix_type0_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.A_B_1;

        numElements = length(batch_log(taskStrengthIdx, rep).Independence.min_A_B);
        min_A_B_matrix_type0_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.min_A_B;
        
        numElements = length(batch_log(taskStrengthIdx, rep).Independence.C);
        C_CNLP{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.C;

        % add townsend data (functional dependence)
        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.A_B);
        A_B_matrix_type4_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.A_B;

        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.A_B_1);
        A_B_1_matrix_type4_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.A_B_1;

        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.min_A_B);
        min_A_B_matrix_type4_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.min_A_B;
        
        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.C);
        C_CNWP{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.C;

        % extract MDS data
        MDS_data(taskStrengthIdx, rep, :, :) = batch_log(taskStrengthIdx, rep).MDS_data;
        
   end
end

%% plot learned representations

plotSettings;

M = nan(numRepetitions, length(interferenceTaskStrength), length(tasksToPerform), length(tasksToPerform));
for rep = 1:size(batch_log,2)
    for strengthIdx = 1:size(batch_log,1)
        M(rep, strengthIdx,:,:) = batch_log(strengthIdx, rep).R_hidden;
    end
end

M = mean(M);
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
set(gca, 'XTickLabel', {'CN', 'WR', 'WP', 'LM'});
set(gca, 'YTickLabel', {'CN', 'WR', 'WP', 'LM'});
ylabel('Tasks','FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Tasks','FontSize', fontSize_xlabel, 'FontName', fontName);
% title({'Correlation of Task Representations' ' at Associative Layer'},'FontSize', fontSize_xlabel-2, 'FontName', fontName);


%% FINAL ACCURACY PLOT (SEBASTIAN)

plotAccuracy = 1;

nSubj = numRepetitions;

fig1 = figure(20);
set(fig1, 'Position', [100 100 900 250]);

% colors
plotSettings;
fontSize_ylabel =17;
barColor = colors(cWeak, :);
scatterColor = colors(cContrast1, :);

% labels
xLabels = {'CN', 'LM', 'WM', 'CN+LM', 'CN+WM'};

% plot
for plotID = 1:2
    
    subplot(1,2,plotID);

    x = 1:5;

    scatterdata_x = repmat(x, nSubj, 1);
    scatterdata_y = nan(size(scatterdata_x));
    
    if(plotID == 1)
        bardata_mean = [nanmean(LCA_taskA_con(:, 1)), ...
                         nanmean(LCA_taskC_con(:, 1)), ...
                         nanmean(LCA_taskB_con(:, 1)), ...
                         nanmean(LCA_Independence_con(:, 1)), ...
                         nanmean(LCA_funcDependence_con(:, 1))] * 100; 

        bardata_sem =   [nanstd(LCA_taskA_con(:, 1)), ...
                         nanstd(LCA_taskC_con(:, 1)), ...
                         nanstd(LCA_taskB_con(:, 1)), ...
                         nanstd(LCA_Independence_con(:, 1)), ...
                         nanstd(LCA_funcDependence_con(:, 1))] * 100 / sqrt(nSubj);  
                     
        scatterdata_y(:,1) = LCA_taskA_con(:, 1) * 100;
        scatterdata_y(:,2) = LCA_taskC_con(:, 1) * 100;
        scatterdata_y(:,3) = LCA_taskB_con(:, 1) * 100;
        scatterdata_y(:,4) = LCA_Independence_con(:, 1) * 100;
        scatterdata_y(:,5) = LCA_funcDependence_con(:, 1) * 100;
    
    elseif(plotID == 2)            
          
        bardata_mean = [nanmean(LCA_taskA_inc(:, 1)), ...
                         nanmean(LCA_taskC_inc(:, 1)), ...
                         nanmean(LCA_taskB_inc(:, 1)), ...
                         nanmean(LCA_Independence_inc(:, 1)), ...
                         nanmean(LCA_funcDependence_inc(:, 1))] * 100; 

        bardata_sem =   [nanstd(LCA_taskA_inc(:, 1)), ...
                         nanstd(LCA_taskC_inc(:, 1)), ...
                         nanstd(LCA_taskB_inc(:, 1)), ...
                         nanstd(LCA_Independence_inc(:, 1)), ...
                         nanstd(LCA_funcDependence_inc(:, 1))] * 100 / sqrt(nSubj);  
                     
        scatterdata_y(:,1) = LCA_taskA_inc(:, 1) * 100;
        scatterdata_y(:,2) = LCA_taskC_inc(:, 1) * 100;
        scatterdata_y(:,3) = LCA_taskB_inc(:, 1) * 100;
        scatterdata_y(:,4) = LCA_Independence_inc(:, 1) * 100;
        scatterdata_y(:,5) = LCA_funcDependence_inc(:, 1) * 100;
        
    end

    hold on;
    bar(x, bardata_mean, 'FaceColor', barColor);
    errorbar(x, bardata_mean, bardata_sem, '.k');
    scatter(scatterdata_x(:), scatterdata_y(:), markerSize, scatterColor, 'LineWidth', 1);
    plot(scatterdata_x', scatterdata_y', 'Color', [scatterColor 0.2]);
    hold off;

    ylim([0 100]);
    xlim([min(x)-0.5 max(x)+0.5]);
    set(gca, 'XTick', x);
    set(gca, 'XTickLabels', xLabels);
    set(gca, 'FontSize',  fontSize_ylabel);
    ylabel('Accuracy (%)','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
    if(plotID == 1)
        title('Congruent Stroop Stimuli');
    elseif(plotID == 2)
        title('Incongruent Stroop Stimuli');
    end

end

%% TOWNSEND & WENGER ANALYSIS
plotSEM = 1;

% PEFFORM ANALYSIS

taskStrengthIdx = 1;
% maxRT = 1.5;
% RT_x = 0.00:0.01:maxRT;
RT_x = CDFs{1}.x;
nSubj = numRepetitions;

CNLP_C_all = nan(nSubj, length(RT_x));
CNWP_C_all = nan(nSubj, length(RT_x));

for subj_idx = 1:nSubj
    
    % CNLP
    CNLP_C_all(subj_idx, :) = C_CNLP{taskStrengthIdx}.data(subj_idx,:);
    CNWP_C_all(subj_idx, :) = C_CNWP{taskStrengthIdx}.data(subj_idx,:);
    
end

% remove invalid values
CNLP_C_all(CNLP_C_all == -Inf) = nan;
CNWP_C_all(CNWP_C_all == -Inf) = nan;
CNLP_C_all(CNLP_C_all == Inf) = nan;
CNWP_C_all(CNWP_C_all == Inf) = nan;

% plot
close all;
fig2 = figure(2);
set(fig2, 'Position', [100 100 400 200]);
plotSettings;

CNLP_mean = nanmean(CNLP_C_all);
CNLP_sem = nanstd(CNLP_C_all)./sqrt(sum(~isnan(CNLP_C_all)));
CNLP_sem(CNLP_sem == 0) = nan;

CNWP_mean = nanmean(CNWP_C_all);
CNWP_sem = nanstd(CNWP_C_all)./sqrt(sum(~isnan(CNWP_C_all)));
CNWP_sem(CNWP_sem == 0) = nan;

hold on;
if(plotSEM)
%     errorbar(RT_x, CNLP_mean, CNLP_sem, '-', 'LineWidth', 1, 'Color', [0 0 0 0.1]);
%     errorbar(RT_x, CNWP_mean, CNWP_sem, '-', 'LineWidth', 1, 'Color', [0 0 0 0.1]);
        s = shadedErrorBar(RT_x,CNLP_C_all,{@nanmean,@nansem},'lineprops', {'-', 'Color',  colors(cContrast2,:)});
        s.patch.FaceColor = colors(cContrast2,:);
        set(s.edge,'LineWidth',1,'LineStyle','-')
        
        s = shadedErrorBar(RT_x,CNWP_C_all,{@nanmean,@nansem},'lineprops', {'-', 'Color',  colors(cContrast3,:)});
        s.patch.FaceColor = colors(cContrast3,:);
        set(s.patch, 'EdgeColor', colors(cContrast3,:));
        set(s.edge,'LineWidth',1,'LineStyle','-')
end
p1 = plot(RT_x, CNLP_mean, '-k', 'LineWidth', 3, 'Color', colors(cContrast2,:));
p2 = plot(RT_x, CNWP_mean, '-k', 'LineWidth', 3, 'Color', colors(cContrast3,:));
hold off;

xlim([0, max(RT_x)]);
ylabel('Capacity C','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
xlabel('Time t in Seconds','FontSize', fontSize_ylabel, 'FontName', fontName, 'FontSize', fontSize_ylabel);
leg = legend([p1 p2], 'CN+LP', 'CN+WM', 'location', 'northwest');
set(leg, 'FontSize', fontSize_ylabel); 
set(gca, 'FontSize', fontSize_gca)
xlim([0 0.4]);

%% MDS PLOT

width = 200; % 200
height = 400; % 200
fontSize_title = 11;
fontSize_xlabel = 12;
fontSize_ylabel = fontSize_xlabel;

repetitionIdx  = 2;
xlimit = [-10 10];
ylimit = [-10 10];

markerSize = 38;
markerLineWidth = 1;

fig = figure(1);
set(fig, 'Position', [100 100 width height]);

x = MultitaskConTraining_MDS_hidden_Phase1(repetitionIdx,1,:,1);
y = MultitaskConTraining_MDS_hidden_Phase1(repetitionIdx,1,:,2);

color = nan(numel(x,1),3);
for i =1:length(x)
    color(i, :) = colors(ceil(i/NPathways),:);
% color(i, :) = colors(cSingle,:);
end

scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth);
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
title({'100% Single Task Training'},'FontSize', fontSize_title, 'FontName', fontName, 'Color', colors(cSingle,:));
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

%% MDS PLOT OLD
plotSettings;

plotLegend = 1;

width = 200;
height = 200;
fontSize_title = 11;
fontSize_xlabel = 14;
fontSize_ylabel = fontSize_xlabel;

repetitionIdx  = 1; 
taskStrengthIdx = 1;
scalar_0 = 10;

scalar_rest = 5;

markerSize = 70;
markerLineWidth = 1;

fig = figure(2);
set(fig, 'Position', [100 100 width height]);

% plot projection of single task representations
xlimit = [-1 1] * scalar_0;
ylimit = [-1 1] * scalar_0;

x = squeeze(MDS_data(taskStrengthIdx, repetitionIdx, :,1));
y = squeeze(MDS_data(taskStrengthIdx, repetitionIdx, :,2));

color = nan(numel(x,1),3);
color(1,:) = colors(cContrast1,:);
color(2,:) = [112 48 160]/255;
color(3,:) = colors(cContrast2,:);
color(4,:) = colors(cContrast3,:);

hold on;
scatter(x(1), y(1), markerSize, color(1,:), 'LineWidth', markerLineWidth);
scatter(x(2), y(2), markerSize, color(2,:), 'LineWidth', markerLineWidth);
scatter(x(3), y(3), markerSize, color(3,:), 'LineWidth', markerLineWidth);
scatter(x(4), y(4), markerSize, color(4,:), 'LineWidth', markerLineWidth);
hold off;

if(plotLegend)
    leg = legend('CN', 'WR', 'WP', 'LM', 'Location', 'north');
    set(leg, 'fontSize', fontSize_xlabel);
end
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
set(gca, 'FontSize', fontSize_xlabel);
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

%% STATS

% set up data matrix
numDataPoints = 2 * 2 * numRepetitions;
X = nan(numDataPoints, 4);

taskStrengthIdx = 1;

% collect data
for i = [1:4:numDataPoints]
    
    rep = ceil(i/4);
    
    % define multitasking condition
    X(i:(i+1),2) = 1; % CNWP
    X((i+2):(i+3),2) = 0; % CNLP
    
    % define congruency condition
    X([i (i + 2)],3) = 1; % incongruent
    X([(i+1) (i+3)],3) = 0; % congruent
    
    % define grouping variable (network)
    X(i:(i+3),4) = rep;
    
    % fill in accuracy
    X(i,1) = LCA_funcDependence_inc(rep, taskStrengthIdx); % CNWP inc
    X(i+1,1) = LCA_funcDependence_con(rep, taskStrengthIdx); % CNWP con
    X(i+2,1) = LCA_Independence_inc(rep, taskStrengthIdx); % CNLP inc
    X(i+3,1) = LCA_Independence_con(rep, taskStrengthIdx); % CNLP con
    
end

% define table
tbl = table(X(:,1), X(:,2), X(:,3), X(:,4), 'VariableNames', {'accuracy', 'multicondition', 'congruency', 'network'});

lme = fitlme(tbl,'accuracy~multicondition*congruency+(1|network)');

pValue = lme.Coefficients.pValue;
estimate = lme.Coefficients.Estimate;
SE = lme.Coefficients.SE;

disp('++++');
idx = 2;
disp(['multicondition: $(\beta = ' num2str(estimate(idx)) ',~SEM = ' num2str(SE(idx)) ',~p = ' num2str(pValue(idx)) ')$']);
idx = 3;
disp(['congruency: $(\beta = ' num2str(estimate(idx)) ',~SEM = ' num2str(SE(idx)) ',~p = ' num2str(pValue(idx)) ')$']);
idx = 4;
disp(['multicondition:congruency: $(\beta = ' num2str(estimate(idx)) ',~SEM = ' num2str(SE(idx)) ',~p = ' num2str(pValue(idx)) ')$']);
