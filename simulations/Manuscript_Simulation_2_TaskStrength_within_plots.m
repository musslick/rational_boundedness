
clear all;
clc;

% results with fixed task-hidden weights (set to 1)
load('logfiles/Part1/PsychReview_Part1_Sim2_V5_TS_within_inc__3P3F_5tasks_11_h100_ortho0_r20.mat');

% PREPARE DATA

numTaskStrengths = length(interferenceTaskStrength);

LCA_funcDependence = nan(numRepetitions, numTaskStrengths);

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

A_B_matrix_type4_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
min_A_B_matrix_type4_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);
A_B_1_matrix_type4_reps{taskStrengthIdx}.data = nan(numRepetitions, maxSize);

for rep = 1:numRepetitions
    
   for taskStrengthIdx = 1:length(interferenceTaskStrength)
       LCA_funcDependence(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy;
       LCA_Independence(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy;
       if(isfield(batch_log, 'taskD'))
           LCA_taskD(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskD.LCA_Accuracy;
       end
       if(isfield(batch_log, 'taskE'))
           LCA_taskE(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskE.LCA_Accuracy;
       end
       if(isfield(batch_log, 'taskC'))
           LCA_taskC(rep, taskStrengthIdx)  = batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy;
       end

        % add townsend data (independence)
        numElements = length(batch_log(taskStrengthIdx, rep).Independence.A_B);
        A_B_matrix_type0_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.A_B;

        numElements = length(batch_log(taskStrengthIdx, rep).Independence.A_B_1);
        A_B_1_matrix_type0_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.A_B_1;

        numElements = length(batch_log(taskStrengthIdx, rep).Independence.min_A_B);
        min_A_B_matrix_type0_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).Independence.min_A_B;

        % add townsend data (functional dependence)
        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.A_B);
        A_B_matrix_type4_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.A_B;

        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.A_B_1);
        A_B_1_matrix_type4_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.A_B_1;

        numElements = length(batch_log(taskStrengthIdx, rep).funcDependence.min_A_B);
        min_A_B_matrix_type4_reps{taskStrengthIdx}.data(rep,1:numElements) = batch_log(taskStrengthIdx, rep).funcDependence.min_A_B;

   end
end

%% plot learned representations (in Figure)

plotSettings;

sharedness_A_D = nan(numRepetitions, length(interferenceTaskStrength));
sharedness_B_E = nan(numRepetitions, length(interferenceTaskStrength));

M = nan(numRepetitions, length(interferenceTaskStrength), length(tasksToPerform), length(tasksToPerform));
for rep = 1:size(batch_log,2)
    for strengthIdx = 1:size(batch_log,1)
        M(rep, strengthIdx,:,:) = batch_log(strengthIdx, rep).hiddenCorr;
        sharedness_A_D(rep, strengthIdx) = M(rep, strengthIdx,1,4);
        sharedness_B_E(rep, strengthIdx) = M(rep, strengthIdx,2,5);
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
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D', 'E'});
set(gca, 'YTickLabel', {'A', 'B', 'C', 'D', 'E'});
ylabel('Tasks','FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Tasks','FontSize', fontSize_xlabel, 'FontName', fontName);
% title({'Correlation of Task Representations' ' at Associative Layer'},'FontSize', fontSize_xlabel-2, 'FontName', fontName);

%% plot learned representations for particular task strength (not in paper)

plotSettings;

sharedness_A_D = nan(numRepetitions, length(interferenceTaskStrength));
sharedness_B_E = nan(numRepetitions, length(interferenceTaskStrength));

M = nan(numRepetitions, length(interferenceTaskStrength), length(tasksToPerform), length(tasksToPerform));
for rep = 1:size(batch_log,2)
    for strengthIdx = 1:size(batch_log,1)
        M(rep, strengthIdx,:,:) = batch_log(strengthIdx, rep).hiddenCorr;
        sharedness_A_D(rep, strengthIdx) = M(rep, strengthIdx,1,4);
        sharedness_B_E(rep, strengthIdx) = M(rep, strengthIdx,2,5);
    end
end

strength = 16;
M = mean(M(:, strength, :, :));
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



%% PLOT MULTITASKING ACCURACY

plotSettings;

if(isfield(batch_log, 'taskD'))

plotSettings;

plotBars = 0;

% set up figure
fig = figure(1);

if(plotBars)
    set(fig, 'Position', [100 100 250 250]);
else
    set(fig, 'Position', [100 100 300 300]);
end

% plot error bar data

grey = [0.8 0.8 0.8];

scale = 100;

% gather data
xLabels = {};


if(plotBars)
    % plot bars                      
    bar(bardata_mean, 'FaceColor', grey, 'EdgeColor', colors(cSingle,:)); hold on;
    % plot error bars
    errorbar(bardata_mean, bardata_sem, '.', 'LineWidth', lineWidth, 'Color', colors(cSingle,  :));
    hold off;
else
    x = interferenceTaskStrength*100;

    % TASK A  
    bardata_mean = nan(1, 0 + numTaskStrengths);
    bardata_sem = nan(1, 0 + numTaskStrengths);
    for i = 1:numTaskStrengths
        bardata_mean(0+i) = mean(LCA_taskD(:, i));
        bardata_sem(0+i) = std(LCA_taskD(:, i)) / sqrt(numRepetitions);
        xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
    end
    bardata_mean = bardata_mean * 100;
    bardata_sem = bardata_sem * 100;
    errorbar(x, bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(3,  :)); hold on;

    % TASK B
    bardata_mean = nan(1, 0 + numTaskStrengths);
    bardata_sem = nan(1, 0 + numTaskStrengths);
    for i = 1:numTaskStrengths
        bardata_mean(0+i) = mean(LCA_taskE(:, i));
        bardata_sem(0+i) = std(LCA_taskE(:, i)) / sqrt(numRepetitions);
        xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
    end
    bardata_mean = bardata_mean * 100;
    bardata_sem = bardata_sem * 100;
    errorbar(x, bardata_mean, bardata_sem, '--', 'LineWidth', lineWidth, 'Color', colors(cSingle,  :)); hold on;

%     % TASK C
%     bardata_mean = nan(1, 0 + numTaskStrengths);
%     bardata_sem = nan(1, 0 + numTaskStrengths);
%     for i = 1:numTaskStrengths
%         bardata_mean(0+i) = mean(LCA_taskC(:, i));
%         bardata_sem(0+i) = std(LCA_taskC(:, i)) / sqrt(numRepetitions);
%         xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
%     end
%     bardata_mean = bardata_mean * 100;
%     bardata_sem = bardata_sem * 100;
%     errorbar(x, bardata_mean, bardata_sem, '--', 'LineWidth', lineWidth, 'Color', colors(2,  :)); hold on;

    % FUNCTIONAL DEPENDENCE    
    bardata_mean = nan(1, 0 + numTaskStrengths);
    bardata_sem = nan(1, 0 + numTaskStrengths);
    for i = 1:numTaskStrengths
        bardata_mean(0+i) = mean(LCA_funcDependence(:, i));
        bardata_sem(0+i) = std(LCA_funcDependence(:, i)) / sqrt(numRepetitions);
        xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
    end
    bardata_mean = bardata_mean * 100;
    bardata_sem = bardata_sem * 100;
    errorbar(x, bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(1,  :)); hold on;

    % INDEPENDENCE
    bardata_mean = nan(1, 0 + numTaskStrengths);
    bardata_sem = nan(1, 0 + numTaskStrengths);
    for i = 1:numTaskStrengths
        bardata_mean(0+i) = mean(LCA_Independence(:, i));
        bardata_sem(0+i) = std(LCA_Independence(:, i)) / sqrt(numRepetitions);
        xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
    end
    bardata_mean = bardata_mean * 100;
    bardata_sem = bardata_sem * 100;
    errorbar(x, bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(2,  :)); hold off;
end

% set axes
ylim([0 1 * scale]);
if(plotBars)
    xlim([0.5 numTaskStrengths+0.5]);
end
if(plotBars)
    set(gca,'XTickLabel',{' '},'FontSize', fontSize_xlabel, 'FontName', fontName);
else
%     set(gca, 'XTick', 1:length(interferenceTaskStrength));
%     set(gca,'XTickLabel',xLabels,'FontSize', fontSize_xlabel-1, 'FontName', fontName);
end
set(gca, 'FontSize', fontSize_ylabel);
ylabel('Accuracy (%)','FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel({'% Training on Tasks D and E' ' Compared to Tasks A, B and C'},'FontSize', fontSize_ylabel, 'FontName', fontName);
leg = legend('Peforming Task D Alone','Performing Task E Alone','Multitasking Tasks A and B', 'Multitasking Tasks A and C', 'Location', 'northoutside');
set(leg, 'FontSize', fontSize_ylabel); 


end



%% Fig 5: Response Time Inequality Series and Task Dependence. 

plotSettings;

plotSEM = 0;
if(~trainMultitasking)
    xlimit = [0.25 0.8];
else
    xlimit = [0.25 0.45];
end
ylimit = [-1 1];
fontSize_title = fontSize_gca;

corrThreshIdx_tmp = 1;

% x_data = 1:size(A_B_matrix_type0_reps{corrThreshIdx_tmp}.data,2);
x_data = CDFs{1}.x;

% type 0 (independence), low task strength
taskStrength = 1.0; % 0.1
taskStrengthIdx = find(interferenceTaskStrength == taskStrength);

y_type0_A_B_mean = nanmean(A_B_matrix_type0_reps{taskStrengthIdx}.data);
y_type0_A_B_1_mean = nanmean(A_B_1_matrix_type0_reps{taskStrengthIdx}.data);
y_type0_min_A_B_mean = nanmean(min_A_B_matrix_type0_reps{taskStrengthIdx}.data);

y_type0_A_B_sem = nanstd(A_B_matrix_type0_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_matrix_type0_reps{taskStrengthIdx}.data));
y_type0_A_B_1_sem = nanstd(A_B_1_matrix_type0_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_1_matrix_type0_reps{taskStrengthIdx}.data));
y_type0_min_A_B_sem = nanstd(min_A_B_matrix_type0_reps{taskStrengthIdx}.data) ./ sqrt(nansum(min_A_B_matrix_type0_reps{taskStrengthIdx}.data));

% type 4 (symmetric dependence)
y_type4_A_B_mean = nanmean(A_B_matrix_type4_reps{taskStrengthIdx}.data);
y_type4_A_B_1_mean = nanmean(A_B_1_matrix_type4_reps{taskStrengthIdx}.data);
y_type4_min_A_B_mean = nanmean(min_A_B_matrix_type4_reps{taskStrengthIdx}.data);

y_type4_A_B_sem = nanstd(A_B_matrix_type4_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_matrix_type4_reps{taskStrengthIdx}.data));
y_type4_A_B_1_sem = nanstd(A_B_1_matrix_type4_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_1_matrix_type4_reps{taskStrengthIdx}.data));
y_type4_min_A_B_sem = nanstd(min_A_B_matrix_type4_reps{taskStrengthIdx}.data) ./ sqrt(nansum(min_A_B_matrix_type4_reps{taskStrengthIdx}.data));

% plot

fig1 = figure(1);
set(fig1, 'Position', [100 100 350 400]);

% INDEPENDENCE

subplot(2, 1, 2);

if(length(x_data) > length(y_type0_A_B_mean))
    x_data = x_data(1:length(y_type0_A_B_mean));
elseif(length(x_data) < length(y_type0_A_B_mean))
    numAdditionalElements = length(y_type0_A_B_mean) - length(x_data);
    x_data(end+(1:numAdditionalElements)) = x_data(end) + (x_data(end)-x_data(end-1)) * (1:numAdditionalElements);
end

if(plotSEM)
    errorbar(x_data, y_type0_A_B_mean, y_type0_A_B_sem, '-k', 'LineWidth',3); hold on;
    errorbar(x_data, y_type0_A_B_1_mean, y_type0_A_B_1_sem, '-k', 'LineWidth',1);
    errorbar(x_data, y_type0_min_A_B_mean, y_type0_min_A_B_sem, '--k', 'LineWidth',1); 
    hold off;
else
    plot(x_data, y_type0_A_B_mean, '-k', 'LineWidth',3, 'Color', colors(2,:)); hold on;
    plot(x_data, y_type0_A_B_1_mean, '-k', 'LineWidth',1, 'Color', colors(2,:));
    plot(x_data, y_type0_min_A_B_mean, '--k', 'LineWidth',1, 'Color', colors(2,:)); 
    hold off;    
end

ylabel({'Probability of',  'Response Before t'}, 'fontSize', fontSize_ylabel-2);
xlabel({'Time t in Seconds'}, 'fontSize', fontSize_xlabel-2);
if(trainMultitasking)
    title({'Tasks A and C'}, 'FontSize', fontSize_title);
end
% leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
leg =  legend({'A AND C', 'A + C - 1', 'min(A, C)'}, 'Location', 'southeast');
ylim(ylimit);
xlim(xlimit);
set(gca, 'FontSize', fontSize_gca);

% FUNCTIONAL DEPENDENCE

subplot(2, 1, 1);

if(plotSEM)
    errorbar(x_data, y_type4_A_B_mean, y_type4_A_B_sem, '-k', 'LineWidth',3); hold on;
    errorbar(x_data, y_type4_A_B_1_mean, y_type4_A_B_1_sem, '-k', 'LineWidth',1);
    errorbar(x_data, y_type4_min_A_B_mean, y_type4_min_A_B_sem, '--k', 'LineWidth',1); 
    hold off;
else
    plot(x_data, y_type4_A_B_mean, '-k', 'LineWidth',3, 'Color', colors(1,:)); hold on;
    plot(x_data, y_type4_A_B_1_mean, '-k', 'LineWidth',1, 'Color', colors(1,:));
    plot(x_data, y_type4_min_A_B_mean, '--k', 'LineWidth',1, 'Color', colors(1,:)); 
    hold off;    
end

ylabel({'Probability of',  'Response Before t'}, 'fontSize', fontSize_ylabel-2);
xlabel({'Time t in Seconds'}, 'fontSize', fontSize_xlabel-2);
if(trainMultitasking)
    title({'Tasks A and B'}, 'FontSize', fontSize_title);
end
% leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
leg =  legend({'A AND B', 'A + B - 1', 'min(A, B)'}, 'Location', 'southeast');
ylim(ylimit);
xlim(xlimit);
set(gca, 'FontSize', fontSize_gca);

%% STATS
disp('Fig. 5: Townsend & Wenger analysis.');

if(splitDependence)
    
    smallestCommonWindow = ~isnan(y_type0_A_B_mean) & ~isnan(y_type0_A_B_1_mean) & ~isnan(y_type0_min_A_B_mean) ...
                                             & ~isnan(y_type3_A_B_mean) & ~isnan(y_type3_A_B_1_mean) & ~isnan(y_type3_min_A_B_mean) ...
                                             & ~isnan(y_type4_A_B_mean) & ~isnan(y_type4_A_B_1_mean) & ~isnan(y_type4_min_A_B_mean);


    type0_aboveCapacity = sum(poslin(A_B_matrix_type0_reps{taskStrengthIdx}.data - min_A_B_matrix_type0_reps{taskStrengthIdx}.data), 2);
    type3_aboveCapacity = sum(poslin(A_B_matrix_type3_reps{taskStrengthIdx}.data - min_A_B_matrix_type3_reps{taskStrengthIdx}.data), 2);
    type4_aboveCapacity = sum(poslin(A_B_matrix_type3_reps{taskStrengthIdx}.data - min_A_B_matrix_type3_reps{taskStrengthIdx}.data), 2);

    type0_belowCapacity = sum(poslin(A_B_1_matrix_type0_reps{taskStrengthIdx}.data - A_B_matrix_type0_reps{taskStrengthIdx}.data), 2);
    type3_belowCapacity = sum(poslin(A_B_1_matrix_type3_reps{taskStrengthIdx}.data - A_B_matrix_type3_reps{taskStrengthIdx}.data), 2);
    type4_belowCapacity = sum(poslin(A_B_1_matrix_type4_reps{taskStrengthIdx}.data - A_B_matrix_type4_reps{taskStrengthIdx}.data), 2);

    [H,P,CI, STATS] = ttest(type0_aboveCapacity, type3_aboveCapacity, 'tail', 'right');
    disp(['t-test type 0 is more above capacity than. type 3, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
    [H,P,CI, STATS] = ttest(type0_aboveCapacity, type4_aboveCapacity, 'tail', 'right');
    disp(['t-test type 0 is more above capacity than. type 4, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
    [H,P,CI, STATS] = ttest(type3_aboveCapacity, type4_aboveCapacity, 'tail', 'right');
    disp(['t-test type 3 is more above capacity than. type 4, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);

    [H,P,CI, STATS] = ttest(type0_belowCapacity, type3_belowCapacity, 'tail', 'left');
    disp(['t-test type 0 is less below capacity than. type 3, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
    [H,P,CI, STATS] = ttest(type0_belowCapacity, type4_belowCapacity, 'tail', 'left');
    disp(['t-test type 0 is less below capacity than. type 4, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
    [H,P,CI, STATS] = ttest(type3_belowCapacity, type4_belowCapacity, 'tail', 'left');
    disp(['t-test type 3 is less below capacity than. type 4, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);

else
    
    smallestCommonWindow = ~isnan(y_type0_A_B_mean) & ~isnan(y_type0_A_B_1_mean) & ~isnan(y_type0_min_A_B_mean) ...
                                             & ~isnan(y_type34_A_B_mean) & ~isnan(y_type34_A_B_1_mean) & ~isnan(y_type34_min_A_B_mean);


    type0_aboveCapacity = sum(poslin(A_B_matrix_type0_reps{taskStrengthIdx}.data - min_A_B_matrix_type0_reps{taskStrengthIdx}.data), 2);
    type34_aboveCapacity = sum(poslin(A_B_matrix_type34_reps{taskStrengthIdx}.data - min_A_B_matrix_type34_reps{taskStrengthIdx}.data), 2);

    type0_belowCapacity = sum(poslin(A_B_1_matrix_type0_reps{taskStrengthIdx}.data - A_B_matrix_type0_reps{taskStrengthIdx}.data), 2);
    type34_belowCapacity = sum(poslin(A_B_1_matrix_type34_reps{taskStrengthIdx}.data - A_B_matrix_type34_reps{taskStrengthIdx}.data), 2);

    [H,P,CI, STATS] = ttest(type0_aboveCapacity, type34_aboveCapacity, 'tail', 'right');
    disp(['t-test type 0 is more above capacity than. type 34, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);

    [H,P,CI, STATS] = ttest(type0_belowCapacity, type34_belowCapacity, 'tail', 'left');
    disp(['t-test type 0 is less below capacity than. type 34, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
    
end


%% PLOT MULTITASKING ACCURACY (POWER POINT)

plotBlack = 0;

plotSettings;

if(isfield(batch_log, 'taskD'))

plotSettings;

plotBars = 0;

% set up figure
close all;

if(plotBlack)
colors = [254 209 103; ...
               254 137 45; ...
               47 205 58; ...
               252 13 27]/255;
end


for plotIdx = 0:4
    
    fig = figure();

    if(plotBars)
        set(fig, 'Position', [100 100 250 250]);
    else
        if(plotIdx == 4)
            set(fig, 'Position', [100 100 300 300]);
        else
            set(fig, 'Position', [100 100 300 200]);
        end
    end

    % plot error bar data

    grey = [0.8 0.8 0.8];

    scale = 100;

    % gather data
    xLabels = {};


    if(plotBars)
        % plot bars                      
        bar(bardata_mean, 'FaceColor', grey, 'EdgeColor', colors(cSingle,:)); hold on;
        % plot error bars
        errorbar(bardata_mean, bardata_sem, '.', 'LineWidth', lineWidth, 'Color', colors(cSingle,  :));
        hold off;
    else
        x = interferenceTaskStrength*100;

        if(plotIdx >= 2 || plotIdx == 0)
            % TASK D  
            bardata_mean = nan(1, 0 + numTaskStrengths);
            bardata_sem = nan(1, 0 + numTaskStrengths);
            for i = 1:numTaskStrengths
                bardata_mean(0+i) = mean(LCA_taskD(:, i));
                bardata_sem(0+i) = std(LCA_taskD(:, i)) / sqrt(numRepetitions);
                xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
            end
            bardata_mean = bardata_mean * 100;
            bardata_sem = bardata_sem * 100;
            if(plotBlack)
                colorIdx = 1;
            else
                colorIdx = 6;
            end
            errorbar(x, bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(colorIdx,  :)); hold on;

        end
        if(plotIdx >= 2 || plotIdx == 0)
            % TASK E
            bardata_mean = nan(1, 0 + numTaskStrengths);
            bardata_sem = nan(1, 0 + numTaskStrengths);
            for i = 1:numTaskStrengths
                bardata_mean(0+i) = mean(LCA_taskE(:, i));
                bardata_sem(0+i) = std(LCA_taskE(:, i)) / sqrt(numRepetitions);
                xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
            end
            bardata_mean = bardata_mean * 100;
            bardata_sem = bardata_sem * 100;
            if(plotBlack)
                colorIdx = 2;
            else
                colorIdx = cSingle;
            end
            errorbar(x, bardata_mean, bardata_sem, '--', 'LineWidth', lineWidth, 'Color', colors(colorIdx,  :)); hold on;

        %     % TASK C
        %     bardata_mean = nan(1, 0 + numTaskStrengths);
        %     bardata_sem = nan(1, 0 + numTaskStrengths);
        %     for i = 1:numTaskStrengths
        %         bardata_mean(0+i) = mean(LCA_taskC(:, i));
        %         bardata_sem(0+i) = std(LCA_taskC(:, i)) / sqrt(numRepetitions);
        %         xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
        %     end
        %     bardata_mean = bardata_mean * 100;
        %     bardata_sem = bardata_sem * 100;
        %     errorbar(x, bardata_mean, bardata_sem, '--', 'LineWidth', lineWidth, 'Color', colors(2,  :)); hold on;

        end
        
        if(plotIdx >= 3 || plotIdx == 0)
            % INDEPENDENCE
            bardata_mean = nan(1, 0 + numTaskStrengths);
            bardata_sem = nan(1, 0 + numTaskStrengths);
            for i = 1:numTaskStrengths
                bardata_mean(0+i) = mean(LCA_Independence(:, i));
                bardata_sem(0+i) = std(LCA_Independence(:, i)) / sqrt(numRepetitions);
                xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
            end
            bardata_mean = bardata_mean * 100;
            bardata_sem = bardata_sem * 100;
            if(plotBlack)
                colorIdx = 3;
            else
                colorIdx = 2;
            end
            errorbar(x, bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(colorIdx,  :)); hold on;
        end
        
        if(plotIdx >= 4 || plotIdx == 0)
            % FUNCTIONAL DEPENDENCE    
            bardata_mean = nan(1, 0 + numTaskStrengths);
            bardata_sem = nan(1, 0 + numTaskStrengths);
            for i = 1:numTaskStrengths
                bardata_mean(0+i) = mean(LCA_funcDependence(:, i));
                bardata_sem(0+i) = std(LCA_funcDependence(:, i)) / sqrt(numRepetitions);
                xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
            end
            bardata_mean = bardata_mean * 100;
            bardata_sem = bardata_sem * 100;
            if(plotBlack)
                colorIdx = 4;
            else
                colorIdx = 1;
            end
            errorbar(x, bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(colorIdx,  :)); hold off;

        end
        
            
    end

    % set axes
    ylim([0 1 * scale]);
    if(plotBars)
        xlim([0.5 numTaskStrengths+0.5]);
    end
    if(plotBars)
        set(gca,'XTickLabel',{' '},'FontSize', fontSize_xlabel, 'FontName', fontName);
    else
    %     set(gca, 'XTick', 1:length(interferenceTaskStrength));
    %     set(gca,'XTickLabel',xLabels,'FontSize', fontSize_xlabel-1, 'FontName', fontName);
    end
    set(gca, 'FontSize', fontSize_ylabel);
    ylabel('Accuracy (%)','FontSize', fontSize_ylabel, 'FontName', fontName);
    xlabel({'% Training on Tasks D and E' ' Compared to Tasks A, B and C'},'FontSize', fontSize_ylabel, 'FontName', fontName);
    if(plotIdx == 4)
        leg = legend('Peforming Task D Alone','Performing Task E Alone','Multitasking Tasks A and C', 'Multitasking Tasks A and B', 'Location', 'northoutside');
        set(leg, 'FontSize', fontSize_ylabel); 
        if(plotBlack)
            set(leg, 'TextColor', 'w');
            set(leg, 'Color', 'none');
        end
    end
    if(plotBlack)
        set(gcf, 'Color', 'k');
        set(gca, 'Color', 'k');
        set(gca, 'xColor', 'w');
        set(gca, 'yColor', 'w');
        set(gca, 'zColor', 'w');
    end
    
    
    hold off;
    end

end


%% Fig 5: Response Time Inequality Series and Task Dependence. (POWERPOINT)

close all;
dontShow = 0;
plotBlack = 0;
plotSettings;

if(plotBlack)
    colors = [254 209 103; ...
                   254 137 45; ...
                   47 205 58; ...
                   252 13 27]/255;
end

plotSEM = 0;
if(~trainMultitasking)
    xlimit = [0.25 0.8];
else
    xlimit = [0.25 0.45];
end
ylimit = [-1 1];
fontSize_title = fontSize_gca;

x_data = CDFs{1}.x;

% type 0 (independence), low task strength
taskStrength = 1;
taskStrengthIdx = find(interferenceTaskStrength == taskStrength);

y_type0_A_B_mean = nanmean(A_B_matrix_type0_reps{taskStrengthIdx}.data);
y_type0_A_B_1_mean = nanmean(A_B_1_matrix_type0_reps{taskStrengthIdx}.data);
y_type0_min_A_B_mean = nanmean(min_A_B_matrix_type0_reps{taskStrengthIdx}.data);

y_type0_A_B_sem = nanstd(A_B_matrix_type0_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_matrix_type0_reps{taskStrengthIdx}.data));
y_type0_A_B_1_sem = nanstd(A_B_1_matrix_type0_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_1_matrix_type0_reps{taskStrengthIdx}.data));
y_type0_min_A_B_sem = nanstd(min_A_B_matrix_type0_reps{taskStrengthIdx}.data) ./ sqrt(nansum(min_A_B_matrix_type0_reps{taskStrengthIdx}.data));

% type 4 (symmetric dependence)
y_type4_A_B_mean = nanmean(A_B_matrix_type4_reps{taskStrengthIdx}.data);
y_type4_A_B_1_mean = nanmean(A_B_1_matrix_type4_reps{taskStrengthIdx}.data);
y_type4_min_A_B_mean = nanmean(min_A_B_matrix_type4_reps{taskStrengthIdx}.data);

y_type4_A_B_sem = nanstd(A_B_matrix_type4_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_matrix_type4_reps{taskStrengthIdx}.data));
y_type4_A_B_1_sem = nanstd(A_B_1_matrix_type4_reps{taskStrengthIdx}.data) ./ sqrt(nansum(A_B_1_matrix_type4_reps{taskStrengthIdx}.data));
y_type4_min_A_B_sem = nanstd(min_A_B_matrix_type4_reps{taskStrengthIdx}.data) ./ sqrt(nansum(min_A_B_matrix_type4_reps{taskStrengthIdx}.data));

% plot

fig1 = figure(1);
if(plotBlack)
set(fig1, 'Position', [100 100 350 200]);
else
set(fig1, 'Position', [100 100 350 160]);    
end

% INDEPENDENCE


if(length(x_data) > length(y_type0_A_B_mean))
    x_data = x_data(1:length(y_type0_A_B_mean));
elseif(length(x_data) < length(y_type0_A_B_mean))
    numAdditionalElements = length(y_type0_A_B_mean) - length(x_data);
    x_data(end+(1:numAdditionalElements)) = x_data(end) + (x_data(end)-x_data(end-1)) * (1:numAdditionalElements);
end

if(plotSEM)
    if(~dontShow)
    errorbar(x_data, y_type0_A_B_mean, y_type0_A_B_sem, '-k', 'LineWidth',3); 
    end
    hold on;
    if(plotBlack)
    errorbar(x_data, y_type0_A_B_1_mean, y_type0_A_B_1_sem, '-w', 'LineWidth',1);
    errorbar(x_data, y_type0_min_A_B_mean, y_type0_min_A_B_sem, '--w', 'LineWidth',1); 
    else
    errorbar(x_data, y_type0_A_B_1_mean, y_type0_A_B_1_sem, '-k', 'LineWidth',1);
    errorbar(x_data, y_type0_min_A_B_mean, y_type0_min_A_B_sem, '--k', 'LineWidth',1); 
    end
    hold off;
else
	if(~dontShow)
        if(plotBlack)
            plot(x_data, y_type0_A_B_mean, '-','Color',colors(3,:), 'LineWidth',3); 
        else
            plot(x_data, y_type0_A_B_mean, '-','Color',colors(2,:), 'LineWidth',3); 
        end
    end
    hold on;
    if(plotBlack)
    plot(x_data, y_type0_A_B_1_mean, '-w', 'LineWidth',1);
    plot(x_data, y_type0_min_A_B_mean, '--w', 'LineWidth',1); 
    else
    plot(x_data, y_type0_A_B_1_mean, '-k', 'LineWidth',1);
    plot(x_data, y_type0_min_A_B_mean, '--k', 'LineWidth',1);     
    end
    hold off;    
end

ylabel({'Probability of',  'Response Before t'}, 'fontSize', fontSize_ylabel-2);
xlabel({'Time t in Seconds'}, 'fontSize', fontSize_xlabel-2);
if(trainMultitasking & plotBlack)
    title({'Tasks A and C'}, 'FontSize', fontSize_title);
end
% leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
if(~dontShow)
leg =  legend({'A AND C', 'A + C - 1', 'min(A, C)'}, 'Location', 'southeast');
else
leg =  legend({'A + C - 1', 'min(A, C)'}, 'Location', 'southeast');    
end
set(leg, 'FontSize', fontSize_ylabel); 
ylim(ylimit);
xlim(xlimit);
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

fig2 = figure(2);
if(plotBlack)
set(fig2, 'Position', [100 100 350 200]);
else
set(fig2, 'Position', [100 100 350 160]);    
end

% FUNCTIONAL DEPENDENCE

if(plotSEM)
    if(~dontShow)
        if(plotBlack)
            errorbar(x_data, y_type4_A_B_mean, y_type4_A_B_sem, '-w', 'LineWidth',3); 
        else
            errorbar(x_data, y_type4_A_B_mean, y_type4_A_B_sem, '-k', 'LineWidth',3); 
        end
    end
    hold on;
    errorbar(x_data, y_type4_A_B_1_mean, y_type4_A_B_1_sem, '-k', 'LineWidth',1);
    errorbar(x_data, y_type4_min_A_B_mean, y_type4_min_A_B_sem, '--k', 'LineWidth',1); 
    hold off;
else
    if(~dontShow)
        if(plotBlack)
            plot(x_data, y_type4_A_B_mean, '-','Color',colors(4,:), 'LineWidth',3); 
        else
            plot(x_data, y_type4_A_B_mean, '-','Color',colors(1,:), 'LineWidth',3); 
        end
    end
    hold on;
    if(plotBlack)
    plot(x_data, y_type4_A_B_1_mean, '-w', 'LineWidth',1);
    plot(x_data, y_type4_min_A_B_mean, '--w', 'LineWidth',1); 
    else
    plot(x_data, y_type4_A_B_1_mean, '-k', 'LineWidth',1);
    plot(x_data, y_type4_min_A_B_mean, '--k', 'LineWidth',1);     
    end
    hold off;    
end

ylabel({'Probability of',  'Response Before t'}, 'fontSize', fontSize_ylabel-2);
xlabel({'Time t in Seconds'}, 'fontSize', fontSize_xlabel-2);
if(trainMultitasking)
    title({'Tasks A and B'}, 'FontSize', fontSize_title);
end
% leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
if(~dontShow)
    leg =  legend({'A AND B', 'A + B - 1', 'min(A, B)'}, 'Location', 'southeast');
else
    leg =  legend({ 'A + B - 1', 'min(A, B)'}, 'Location', 'southeast');
end
set(leg, 'FontSize', fontSize_ylabel); 
ylim(ylimit);
xlim(xlimit);
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