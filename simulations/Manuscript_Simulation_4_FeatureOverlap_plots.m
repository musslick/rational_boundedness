%% data settings
clear all;
clc;

datafilePrefix ='PychReview_Part2_Sim1_FeatureOverlap';       

correlation_threshold = 0.5;
performanceMeasure = 'LCA';   % alternatives: 'Accuracy'
optimalPerformanceThreshold = 0.9;
useWeightSimilarity = 0;

repetitionIdx = 1;

% load data
% logFolder = 'logfiles/param_set 1 complete/';
logFolder = 'logfiles/Part2/';
files = dir(logFolder);

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
    for fileIdx = 1:length(validFileNames);
        disp(['loading ' validFileNames{fileIdx} '...']);
        load(strcat(logFolder, validFileNames{fileIdx}));
        % initial setup of batch_log_tmp
        if(fileIdx == 1)
            batch_log_tmp = repmat(batch_log(1,1), numRepetitions, 1);
        end
        batch_log_tmp(fileIdx) = batch_log(1,1);
    end
else
    error('No valid file names found');
end
batch_log = batch_log_tmp;

corrThreshIdx = find(corr_thresholds == correlation_threshold,1);
if(isempty(corrThreshIdx))
    warning('Requested correlation threshold for MIS extraction does not exist in data set');
end

optPerformanceThreshIdx = find(round(goodPerformanceThresholds*100)/100 == optimalPerformanceThreshold,1);
if(isempty(optPerformanceThreshIdx))
    warning('Requested optimal performance threshold does not exist in data set');
end

numCorrThresholds = length(corr_thresholds);
numGoodPerformanceThresholds = length(goodPerformanceThresholds);

% plot settings

fontSize_title = 14;
fontSize_gca = 14;
fontSize_xlabel = 14;
fontSize_ylabel = 14;
fontSize_legend = 14;

fontName = 'Helvetica';
markerSize = 50;
sem_width = 2;
sem_marker = '-';

lineWidth = 3;

colors = [253 120 21; ... % orange
              31 104 172; ... % blue
              44 155 37; ... % green
              0     0   0  ; ... % black
            142 142 142; ... % grey 
            255 255 255] / 255; % white 
        
cContrast1 = 1;
cContrast2 = 2;
cContrast3 = 3;
cSingle = 4;
cWeak = 5;
cWhite = 6;

% plot graph analysis results

close all;

%% prepare data

genPerformance = nan(numRepetitions, length(taskSimilarities));
multiPerformance  = nan(numRepetitions, length(taskSimilarities), NPathways);
learningPerformance = nan(numRepetitions, length(taskSimilarities));
MIS = nan(numRepetitions, length(taskSimilarities), numCorrThresholds);
taskCardinality = nan(numRepetitions, length(taskSimilarities), numGoodPerformanceThresholds);
weightCorr_final = nan(numRepetitions, length(taskSimilarities));

basis_template = eye(nTasks,nTasks);
for row = 1:size(basis_template,1)
    basis_template(row, (ceil(row/NPathways)-1)*NPathways+(1:NPathways)) = 1;
    basis_template(row, row:end) = 0;
end

for i = 1:length(taskSimilarities)
    
    for rep = 1:numRepetitions
        
        learnIdx = find(batch_log(rep).MSE_log(i,:) > 0);
        learnIdx = learnIdx(end);
        learningPerformance(rep, i) = learnIdx;
        
        if(isfield(batch_log(rep), 'hiddenTaskRepWeights') && useWeightSimilarity)
            hiddenTaskRepWeights = squeeze(batch_log(rep).hiddenTaskRepWeights(i, learnIdx, :,:));
            hiddenTaskRep = corr(hiddenTaskRepWeights);
        else
            hiddenTaskRep = squeeze(batch_log(rep).hiddenTaskRep(i, :,:));
        end
        weightCorr_final(rep, i) = mean(hiddenTaskRep(basis_template == 1));
        
        % merge data
        switch performanceMeasure
            case 'MSE'
                % generalization performance
                genPerformance(rep, i) = batch_log(rep).test_MSE(i, learnIdx(end));
                % multitasking performance
                for cap = 2:NPathways
                    multiPerformance(rep, i, cap) = batch_log(rep).MSE_multi{i, learnIdx, cap};
                end
                % optimal performance
                if(isfield(batch_log, 'taskCardinalityLog_meanAcc'))
                    for perf_thresh_idx = 1:numGoodPerformanceThresholds
                        taskCardinality(rep, i, perf_thresh_idx) = batch_log(rep).taskCardinalityLog_meanAcc(i, learnIdx, perf_thresh_idx);
                    end
                else
                    warning('No task cardinality data available.');
                    taskCardinality = [];
                end

            case 'Accuracy'
                % generalization performance
                genPerformance(rep, i) = batch_log(rep).test_Accuracy(i, learnIdx(end));
                % multitasking performance
                for cap = 2:NPathways
                    multiPerformance(rep, i, cap) = batch_log(rep).Accuracy_multi_inc{i, learnIdx, cap};
                end
                % optimal performance
                if(isfield(batch_log, 'taskCardinalityLog_meanAcc'))
                    for perf_thresh_idx = 1:numGoodPerformanceThresholds
                        taskCardinality(rep, i, perf_thresh_idx) = batch_log(rep).taskCardinalityLog_meanAcc(i, learnIdx, perf_thresh_idx);
                    end
                else
                    warning('No task cardinality data available.');
                    taskCardinality = [];
                end
                
            case 'PCorrect'
                % generalization performance
                genPerformance(rep, i) = batch_log(rep).test_PCorrect(i, learnIdx(end));
                % multitasking performance
                for cap = 2:NPathways
                    multiPerformance(rep, i, cap) = batch_log(rep).PCorrect_multi_inc{i, learnIdx, cap};
                end
                % optimal performance
                if(isfield(batch_log, 'taskCardinalityLog_meanAcc'))
                    for perf_thresh_idx = 1:numGoodPerformanceThresholds
                        taskCardinality(rep, i, perf_thresh_idx) = batch_log(rep).taskCardinalityLog_respProb(i, learnIdx, perf_thresh_idx);
                    end
                else
                    warning('No task cardinality data available.');
                    taskCardinality = [];
                end
                
            case 'LCA'

                % multitasking performance
                for cap = 2:NPathways
                    multiPerformance(rep, i, cap) = batch_log(rep).LCA_multi{i, cap};
                end
   
        end
       
    end
    
end

%% Fig 1: Learned Weight Similarity vs. Feature Overlap
plotColors = 0;

% create scatter plot
x = taskSimilarities;
y = nanmean(weightCorr_final,1);
y_sem = nanstd(weightCorr_final,[],1)/sqrt(size(weightCorr_final,1));

cmap = colormap;
colors = nan(length(taskSimilarities), 3);
for i = 1:size(colors,1)
    cIndex = round((size(cmap,1)-1)*taskSimilarities(i)/max(taskSimilarities));
    colors(i,:) = cmap(cIndex+1,:);
end
if(plotColors)
    scatter(x, y, markerSize, colors, 'filled', 'LineWidth',1, 'MarkerEdgeColor',[0 0 0]); hold on;
else
    errorbar(x, y, y_sem,'-k','LineWidth',2);
end
%caxis([]);
set(gca, 'FontSize', fontSize_gca-2, 'FontName', fontName);
ylabels = (get(gca,'YTickLabel'));

xlabel('Feature Overlap Between Tasks','FontSize', fontSize_gca, 'FontName', fontName, 'Color', 'k');
ylabel({'Learned Task Similarity', 'in Hidden Layer'},'FontSize', fontSize_gca, 'FontName', fontName, 'Color', 'k');

set(gcf, 'Color', 'k');
set(gcf, 'Position', [500 500 220 200]);

if(isfield(batch_log(rep), 'hiddenTaskRepWeights'))
    xlim([0 1]);
end

% plot SEM
plotSEM = 1;
sem_marker = '-k';
if(plotSEM & plotColors)
    for i = 1:length(taskSimilarities)
        plot([x(i) x(i)], [y(i)  y(i)+y_sem(i)], sem_marker, 'LineWidth', sem_width);
        plot([x(i) x(i)], [y(i)  y(i)-y_sem(i)], sem_marker, 'LineWidth', sem_width);
    end
end
hold off;

% added for black figures

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'xColor', 'k');
set(gca, 'yColor', 'k');

% preserve background color when saving figure
fig = gcf;
fig.InvertHardcopy = 'off';

%% Fig 2: Learning Speed vs. Feature Overlap

colorMultitasking = 1;
cap = 2;

% create scatter plot
x = taskSimilarities;

y = nanmean(learningPerformance(:,:),1);
y_sem = nanstd(learningPerformance, [], 1)/sqrt(size(learningPerformance,1));

if(~colorMultitasking)
    cmap = colormap;
    colors = nan(length(taskSimilarities), 3);
    for i = 1:size(colors,1)
        cIndex = round((size(cmap,1)-1)*taskSimilarities(i)/max(taskSimilarities));
        colors(i,:) = cmap(cIndex+1,:);
    end
else
    multitaskingAccuracy = nanmean(multiPerformance(:,:,cap),1);
    multitaskingAccuracy = multitaskingAccuracy - min(multitaskingAccuracy);
    multitaskingAccuracy = multitaskingAccuracy / max(multitaskingAccuracy);
    cmap = colormap;
    colors = nan(length(multitaskingAccuracy), 3);
    for i = 1:size(colors,1)
        cIndex = round((size(cmap,1)-1)*multitaskingAccuracy(i));
        colors(i,:) = cmap(cIndex+1,:);
    end
end

scatter(x, y, markerSize, colors, 'filled', 'LineWidth',1, 'MarkerEdgeColor',[0 0 0]); hold on;
%caxis([]);
set(gca, 'FontSize', fontSize_gca-2, 'FontName', fontName);
ylabels = (get(gca,'YTickLabel'));

xlabel('Feature Overlap Between Tasks','FontSize', fontSize_gca, 'FontName', fontName, 'Color', 'k');
ylabel('Iterations Required to Train','FontSize', fontSize_gca, 'FontName', fontName, 'Color', 'k');

set(gcf, 'Color', 'k');
set(gcf, 'Position', [500 500 250 200]);

if(isfield(batch_log(rep), 'hiddenTaskRepWeights'))
    xlim([0 1]);
end

% plot SEM
plotSEM = 1;
sem_marker = '-k';
if(plotSEM & ~colorMultitasking)
    for i = 1:length(taskSimilarities)
        plot([x(i) x(i)], [y(i)  y(i)+y_sem(i)], sem_marker, 'LineWidth', sem_width);
        plot([x(i) x(i)], [y(i)  y(i)-y_sem(i)], sem_marker, 'LineWidth', sem_width);
    end
end

if(colorMultitasking)
    multitaskingAccuracy = nanmean(multiPerformance(:,:,cap),1) * 100;
    cbh = colorbar;
    set(cbh,'XTickLabel',round(linspace(min(multitaskingAccuracy), max(multitaskingAccuracy), 6)*10)/10)
    errorbar(x, y, y_sem,'-k','LineWidth',1); hold on;
end

hold off;

% added for black figures

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'xColor', 'k');
set(gca, 'yColor', 'k');

% preserve background color when saving figure
fig = gcf;
fig.InvertHardcopy = 'off';

%% PRINT MULTI PERFORMANCE
cap = 2;

y = nanmean(multiPerformance(:,:,cap),1);
y_sem = nanstd(multiPerformance(:,:,cap), [], 1);
disp(round(y*10000)/100);
disp(y_sem);

%% Fig 3:  Multitasking Performance vs. Learning Speed vs. Feature Overlap

cap = 2;
plotSEM = 1;

% plot data
figure(1)

x = nanmean(learningPerformance(:,:),1);
x_sem = nanstd(learningPerformance, [], 1)/sqrt(size(learningPerformance,1));
y = nanmean(multiPerformance(:,:,cap),1);
y_sem = nanstd(multiPerformance(:,:,cap), [], 1)/sqrt(size(multiPerformance(:,:,cap),1));
z = taskSimilarities;
disp(y);

% create scatter plot
if(plotSEM)
    e_v = errorbar(x, y, y_sem, '.', 'Color', 'k'); hold on;
    e_v.LineStyle = 'none';
    e_v.LineWidth = lineWidth-2;
    e_h = herrorbar(x, y, x_sem, x_sem, '.'); hold on;
    e_h(1).LineWidth = lineWidth-2;
    e_h(1).Color = [0 0 0];
end

% compute axis parameters
switch performanceMeasure
    case 'MSE'
        measurementLabel = 'MSE';
        ygranularity = 1000;
        ymargin = 0.01;
    case 'Accuracy'
        measurementLabel = 'Accuracy';
        ygranularity = 10000;
        ymargin = 0.001;
    case 'PCorrect'   
        %measurementLabel = 'P(Correct)';
        y = y*100;
        y_sem = y_sem*100;
        measurementLabel = 'Accuracy (%)';
        ygranularity = 1000;
        ymargin = 0.04;
end
% ylimit = [floor(min(y-y_sem)*ygranularity)/ygranularity ceil(max(y+y_sem)*ygranularity)/ygranularity] .* [1-ymargin 1+ymargin];
% xgranularity = 10;
% xmargin = 0.02;
% xlimit = [floor(min(x)*xgranularity)/xgranularity ceil(max(x)*xgranularity)/xgranularity] .* [1-xmargin 1+xmargin];

scatter(x, y, markerSize, z, 'filled', 'MarkerEdgeColor',[0 0 0]);
hold off;


set(gca, 'FontSize', fontSize_gca-2, 'FontName', fontName);
% xlim(xlimit);
% ylim(ylimit);

ylabel({'Multitasking Accuracy (%)'},'FontSize', fontSize_gca, 'FontName', fontName, 'Color', 'k');
xlabel(['Iterations Required To Train'],'FontSize', fontSize_gca, 'Color', 'k');

hcb = colorbar;
caxis([0 1]);
set(hcb,'Color','k');
ylabel(hcb, 'Feature Overlap','FontSize', fontSize_gca-2, 'Color', 'k');
set(gcf, 'Color', 'k');
set(gcf, 'Position', [500 500 300 200]);

hold off;

% added for black figures

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'xColor', 'k');
set(gca, 'yColor', 'k');

% preserve background color when saving figure
fig = gcf;
fig.InvertHardcopy = 'off';

