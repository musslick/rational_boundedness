%% data settings

%  Manuscript FIG_DEPENDENCY_RESULTS

clear all;
clc;
close all;

datafilePrefix = 'PsychReview_Part1_Sim1_sharedRepresentation_5P3F_10tasks_20_h100_r';

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
    corrThreshIdx = 1;
    warning(['Requested correlation threshold for MIS extraction does not exist in data set, using first correlation threshold with value ' num2str(correlation_thresholds(corrThreshIdx))]);
end

numCorrelationThresholds = length(correlation_thresholds);

% load plotSettings
loadPlotSettings;

%% PREP DATA

% prepare MIS data

if(isfield(batch_log(numHiddenIdx,1), 'maxCarryingCapacity'))

    data_extractedMIS = nan(numRepetitions, length(correlation_thresholds));

    for rep =1:numRepetitions
        data_extractedMIS(rep,:) = str2num(batch_log(numHiddenIdx,rep).maxCarryingCapacity);
    end
    
    numCorrThresholds = size(data_extractedMIS,2);
else
    numCorrThresholds = length(correlation_thresholds);
end

% prepare performance curve data

maxN = NPathways;
maximumPerformanceCurveLog.LCA = nan(numRepetitions, 1, maxN);

for rep = 1:numRepetitions
   
    maxCardinaility = size(batch_log(numHiddenIdx, rep).maximumPerformanceCurveLog.LCA,2);
    maximumPerformanceCurveLog.LCA(rep,1,1:maxCardinaility) = batch_log(numHiddenIdx, rep).maximumPerformanceCurveLog.LCA;

end

% prepare dependency data

depencencyLogLCA = nan([size(batch_log(1, 1).depencencyLogLCA_mean{1}) numRepetitions]);
depencencyLogAvgLCA = nan([size(batch_log(1, 1).depencencyLogAvgLCA_mean{1}) numRepetitions]);

for rep = 1:numRepetitions
    for corrIdx = 1:numCorrThresholds
        depencencyLogLCA(:,:, rep) = batch_log(numHiddenIdx, rep).depencencyLogLCA_mean{corrThreshIdx};
        depencencyLogAvgLCA(:,:, rep) = batch_log(numHiddenIdx, rep).depencencyLogAvgLCA_mean{corrThreshIdx};
    end
end

depencencyLogLCA_summary = nan(size(depencencyLogLCA,1)-1, 2, size(depencencyLogLCA,3));
for setSize = 1:size(depencencyLogLCA_summary,1)
    depencencyLogLCA_summary(setSize,1,:) = depencencyLogLCA(setSize+1,1,:);
    depencencyLogLCA_summary(setSize,2,:) = nanmean(depencencyLogLCA(setSize+1,:,:),2);
end

depencencyLogLCA_summary_mean = nanmean(depencencyLogLCA_summary, 3);
depencencyLogLCA_summary_sem= nanstd(depencencyLogLCA_summary, [], 3) ./ sum(~isnan(depencencyLogLCA_summary),3);

depencencyLogLCA_mean = nanmean(depencencyLogLCA, 3);
depencencyLogLCA_sem= nanstd(depencencyLogLCA, [], 3) ./ sum(~isnan(depencencyLogLCA),3);

depencencyLogAvgLCA_mean = nanmean(depencencyLogAvgLCA, 3);
depencencyLogAvgLCA_sem= nanstd(depencencyLogAvgLCA, [], 3) ./ sum(~isnan(depencencyLogAvgLCA),3);

% get dependency data for activation plots

for corrThreshIdx_tmp = 1:numCorrThresholds
    for rep =1:numRepetitions

        if(isfield(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}, 'taskDependencies'))
        
            % get indices for different interference types
            type0Idx =  find(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies(3,:) == 0);
            type3Idx =  find(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies(3,:) == 3);
            type4Idx =  find(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies(3,:) == 4);
            type34Idx =  find(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies(3,:) == 3 ...
                                        | batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies(3,:) == 4);

            % congruent
            dependencyData{corrThreshIdx_tmp}.activation_type0_con(rep, :) = mean(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.responseActivation_con(type0Idx,:),1);
            dependencyData{corrThreshIdx_tmp}.activation_type3_con(rep, :) = mean(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.responseActivation_con(type3Idx,:),1);
            dependencyData{corrThreshIdx_tmp}.activation_type4_con(rep, :) = mean(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.responseActivation_con(type4Idx,:),1);

            % incongruent
            dependencyData{corrThreshIdx_tmp}.activation_type0_inc(rep, :) = mean(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.responseActivation_inc(type0Idx,:),1);
            dependencyData{corrThreshIdx_tmp}.activation_type3_inc(rep, :) = mean(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.responseActivation_inc(type3Idx,:),1);
            dependencyData{corrThreshIdx_tmp}.activation_type4_inc(rep, :) = mean(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.responseActivation_inc(type4Idx,:),1);

            if(isfield(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}, 'A_B'))

                 % determine maximum size of CDFs
                maxSize = 0;
                for combIdx = 1:(length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B))
                    maxSize = max([maxSize length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{combIdx})]);
                    maxSize = max([maxSize length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B_1{combIdx})]);
                    maxSize = max([maxSize length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.min_A_B{combIdx})]);
                end

                A_B_matrix_type0 = nan(length(type0Idx), maxSize);
                A_B_matrix_type3 = nan(length(type3Idx), maxSize);
                A_B_matrix_type4 = nan(length(type4Idx), maxSize);
                A_B_matrix_type34 = nan(length(type34Idx), maxSize);
                
                A_B_1_matrix_type0 = nan(length(type0Idx), maxSize);
                A_B_1_matrix_type3 = nan(length(type3Idx), maxSize);
                A_B_1_matrix_type4 = nan(length(type4Idx), maxSize);
                A_B_1_matrix_type34 = nan(length(type34Idx), maxSize);
                
                min_A_B_matrix_type0 = nan(length(type0Idx), maxSize);
                min_A_B_matrix_type3 = nan(length(type3Idx), maxSize);
                min_A_B_matrix_type4 = nan(length(type4Idx), maxSize);
                min_A_B_matrix_type34 = nan(length(type34Idx), maxSize);

                % collect CDFs for given network and correlation threshold (type 0)
                for combIdx = 1:length(type0Idx)
                    currentSize = length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type0Idx(combIdx)});
                    A_B_matrix_type0(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type0Idx(combIdx)};
                    A_B_1_matrix_type0(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B_1{type0Idx(combIdx)};
                    min_A_B_matrix_type0(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.min_A_B{type0Idx(combIdx)};
                end
                
                % collect CDFs for given network and correlation threshold (type 3)
                for combIdx = 1:length(type3Idx)
                    currentSize = length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type3Idx(combIdx)});
                    A_B_matrix_type3(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type3Idx(combIdx)};
                    A_B_1_matrix_type3(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B_1{type3Idx(combIdx)};
                    min_A_B_matrix_type3(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.min_A_B{type3Idx(combIdx)};
                end
                
                % collect CDFs for given network and correlation threshold (type 4)
                for combIdx = 1:length(type4Idx)
                    currentSize = length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type4Idx(combIdx)});
                    A_B_matrix_type4(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type4Idx(combIdx)};
                    A_B_1_matrix_type4(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B_1{type4Idx(combIdx)};
                    min_A_B_matrix_type4(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.min_A_B{type4Idx(combIdx)};
                end
                
                % collect CDFs for given network and correlation threshold (type 34)
                for combIdx = 1:length(type34Idx)
                    currentSize = length(batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type34Idx(combIdx)});
                    A_B_matrix_type34(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B{type34Idx(combIdx)};
                    A_B_1_matrix_type34(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.A_B_1{type34Idx(combIdx)};
                    min_A_B_matrix_type34(combIdx, 1:currentSize) = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.min_A_B{type34Idx(combIdx)};
                end

            end
        
        end
        
        A_B_type0_reps{rep} = nanmean(A_B_matrix_type0);
        A_B_type3_reps{rep} = nanmean(A_B_matrix_type3);
        A_B_type4_reps{rep} = nanmean(A_B_matrix_type4);
        A_B_type34_reps{rep} = nanmean(A_B_matrix_type34);

        A_B_1_type0_reps{rep} = nanmean(A_B_1_matrix_type0);
        A_B_1_type3_reps{rep} = nanmean(A_B_1_matrix_type3);
        A_B_1_type4_reps{rep} = nanmean(A_B_1_matrix_type4);
        A_B_1_type34_reps{rep} = nanmean(A_B_1_matrix_type34);

        min_A_B_type0_reps{rep} = nanmean(min_A_B_matrix_type0);
        min_A_B_type3_reps{rep} = nanmean(min_A_B_matrix_type3);
        min_A_B_type4_reps{rep} = nanmean(min_A_B_matrix_type4);
        min_A_B_type34_reps{rep} = nanmean(min_A_B_matrix_type34);

    end
    
    % determine maximum size of CDFs
    maxSize = 0;
    for rep = 1:numRepetitions
        maxSize = max([maxSize length(A_B_type0_reps{rep})]);
        maxSize = max([maxSize length(A_B_type3_reps{rep})]);
        maxSize = max([maxSize length(A_B_type4_reps{rep})]);
        maxSize = max([maxSize length(A_B_type34_reps{rep})]);
        maxSize = max([maxSize length(A_B_1_type0_reps{rep})]);
        maxSize = max([maxSize length(A_B_1_type3_reps{rep})]);
        maxSize = max([maxSize length(A_B_1_type4_reps{rep})]);
        maxSize = max([maxSize length(A_B_1_type34_reps{rep})]);
        maxSize = max([maxSize length(min_A_B_type0_reps{rep})]);
        maxSize = max([maxSize length(min_A_B_type3_reps{rep})]);
        maxSize = max([maxSize length(min_A_B_type4_reps{rep})]);
        maxSize = max([maxSize length(min_A_B_type34_reps{rep})]);
    end

    A_B_matrix_type0_reps{corrThreshIdx_tmp}.data = nan(numRepetitions, maxSize);
    A_B_matrix_type3_reps{corrThreshIdx_tmp}.data = nan(numRepetitions, maxSize);
    A_B_matrix_type4_reps{corrThreshIdx_tmp}.data = nan(numRepetitions, maxSize);
    A_B_matrix_type34_reps{corrThreshIdx_tmp}.data = nan(numRepetitions, maxSize);

    A_B_1_matrix_type0_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    A_B_1_matrix_type3_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    A_B_1_matrix_type4_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    A_B_1_matrix_type34_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);

    min_A_B_matrix_type0_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    min_A_B_matrix_type3_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    min_A_B_matrix_type4_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    min_A_B_matrix_type34_reps{corrThreshIdx_tmp}.data  = nan(numRepetitions, maxSize);
    
    for rep = 1:numRepetitions
        
        % collect CDFs for given network and correlation threshold (type 0)
        currentSize = length(A_B_type0_reps{rep});
        A_B_matrix_type0_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_type0_reps{rep};
        A_B_1_matrix_type0_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_1_type0_reps{rep};
        min_A_B_matrix_type0_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = min_A_B_type0_reps{rep};
        
        % collect CDFs for given network and correlation threshold (type 3)
        currentSize = length(A_B_type3_reps{rep});
        A_B_matrix_type3_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_type3_reps{rep};
        A_B_1_matrix_type3_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_1_type3_reps{rep};
        min_A_B_matrix_type3_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = min_A_B_type3_reps{rep};
        
         % collect CDFs for given network and correlation threshold (type 4)
        currentSize = length(A_B_type4_reps{rep});
        A_B_matrix_type4_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_type4_reps{rep};
        A_B_1_matrix_type4_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_1_type4_reps{rep};
        min_A_B_matrix_type4_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = min_A_B_type4_reps{rep};
        
        % collect CDFs for given network and correlation threshold (type 4)
        currentSize = length(A_B_type34_reps{rep});
        A_B_matrix_type34_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_type34_reps{rep};
        A_B_1_matrix_type34_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = A_B_1_type34_reps{rep};
        min_A_B_matrix_type34_reps{corrThreshIdx_tmp}.data(rep,1:currentSize) = min_A_B_type34_reps{rep};
        
    end
    
end

% fix task dependencies
for corrThreshIdx_tmp = 1:numCorrThresholds
    for rep = 1:numRepetitions
        dualTaskCombs = batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies;
        A_tasksIdx = batch_log(numHiddenIdx, rep).A_tasksIdx{1};
        R_hidden = batch_log(numHiddenIdx, rep).R_hidden_avg;
        tasksToPerform = batch_log(numHiddenIdx, rep).tasksToPerform;
        [dualTaskCombs, similarity] = recomputeTaskDependencies(A_tasksIdx, dualTaskCombs, R_hidden, tasksToPerform);
         
        batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.taskDependencies = dualTaskCombs;
        batch_log(numHiddenIdx, rep).taskCombData{corrThreshIdx_tmp}.similarity = similarity;
    end
end

%% FIG 1: Parallel Processing Capacity -- Performance as a function of predicted MIS (for specific correlation threshold)

% parameters 
optimalPerformanceThreshold = 0.90; % threshold on performance to decide whether the multitasking attempt was successful or not (only necessary for effective PPC metric)
cardinality_window = -2:2;  % plotted cardinality relative to predicted MIS
performanceMetric = 'LCA';         
performanceType = 'max';                               
showSingleGraph = 0; % 200
maxN = 5;
maxMIS = 4;

% select corresponding performance
switch performanceMetric
    case 'LCA'
        switch performanceType
            case 'max'
                performanceData = maximumPerformanceCurveLog.LCA;
        end
end

% first transform data such that performance curve is aligned with predicted N
accuracyCurve = nan(numRepetitions, length(cardinality_window));

performanceData = squeeze(nanmean(performanceData, 2)); % take mean performance over all repetitions per graph
MISData = squeeze(mode(data_extractedMIS, 2));  % take mode of MIS across all repetitions per graph
MISData_specific = squeeze(MISData(:, corrThreshIdx)); % pick MIS corresponding to selected correlation threshold
middleIndex = find(cardinality_window == 0);

for graphIdx = 1:size(performanceData, 1)
    
    minIndexOrg = max(1, MISData_specific(graphIdx) + cardinality_window(1));
    maxIndexOrg = min(maxN, MISData_specific(graphIdx) + cardinality_window(end));
    
    minIndexNormalized = middleIndex - min(MISData_specific(graphIdx), middleIndex) + 1;
    maxIndexNormalized = middleIndex + min(maxN - MISData_specific(graphIdx),  middleIndex-1);
    
    accuracyCurve(graphIdx, minIndexNormalized:maxIndexNormalized) = performanceData(graphIdx,minIndexOrg:maxIndexOrg);
    
end

% plot all performance curves
fig1 = figure(1);
set(fig1, 'Position', [500 500 550 250]);
plotColors = colormap(copper(maxN+1));

colorOffset = 0.2;
for row = 1:size(plotColors, 1)
    for col = 1:size(plotColors, 2)
        plotColors(row, col) = min(plotColors(row, col) + colorOffset, 1);
    end
end
plotColors = plotColors(fliplr(1:size(plotColors,1)),:);

if(strcmp(performanceMetric, 'MSE'))
    scalar = 1;
else
    scalar = 100;
end

if(any(MISData_specific == 0))
    warning('Found MIS of 0. Fixing to 1.');
    MISData_specific(MISData_specific == 0) = 1;
end

if(showSingleGraph == 0)
    for graphIdx = 1:size(accuracyCurve, 1)
        plot(accuracyCurve(graphIdx,:) * scalar, '-', 'LineWidth', lineWidth-1, 'color', plotColors(MISData_specific(graphIdx),:)); hold on;
    end
else
    graphIdx = showSingleGraph;
    plot(accuracyCurve(graphIdx,:) * scalar, '-', 'LineWidth', lineWidth-1, 'color', 'k'); hold on;
end

lineColor = '--k';
plot([middleIndex middleIndex], [0 100], lineColor, 'LineWidth', lineWidth);

% legend
if(showSingleGraph == 0 | showSingleGraph ~= 0)
    h = zeros(maxMIS, 1);
    legendLabels = {};
    for i = 1:maxMIS
        h(i) = plot(NaN, NaN, '-', 'LineWidth', lineWidth, 'color', plotColors(i,:)); 
        legendLabels{i} = ['MIS =' num2str(i)];
    end
    leg = legend(h, legendLabels,'Location','eastoutside');
    set(leg, 'FontSize', fontSize_gca);
end

% axes
xTickLabels = {};
if(showSingleGraph == 0)
    for i = 1:length(cardinality_window)
        if(sign(cardinality_window(i)) > 0)
            signLabel = ['+' num2str(cardinality_window(i))];
        elseif(sign(cardinality_window(i)) < 0)
            signLabel = num2str(cardinality_window(i));
        else
            signLabel = '';
        end
        xTickLabels{i} = ['MIS' signLabel];
    end
else
    for i = 1:length(cardinality_window)
        xTickLabels{i} = num2str(MISData_specific(graphIdx) + min(cardinality_window) + i -1);
    end
end


set(gca, 'XTick', 1:length(cardinality_window));
set(gca, 'XTickLabel', xTickLabels);
if(strcmp(performanceMetric, 'LCA'))
    ylim([0 100]);
    ylabel('Multitasking Accuracy (%)', 'FontSize', fontSize_ylabel);
elseif(strcmp(performanceMetric, 'MSE'))
    ylim([0 0.08]);
    ylabel('MSE', 'FontSize', fontSize_ylabel);
else
    ylabel('Performance (%)', 'FontSize', fontSize_ylabel);
end
xlim([0.5, length(cardinality_window) + 0.5]);
% title([performanceMetric ' ' performanceType], 'FontSize', fontSize_title);
xlabel('Task Set Size', 'FontSize', fontSize_ylabel);
set(gca, 'FontSize', fontSize_gca);

% plot line to indicate MIS
if(plotBlack)
    lineColor = '--w';
else
    lineColor = '--k';
end

hold off;

if(plotBlack)
    if(showSingleGraph == 0)
        set(leg, 'TextColor', 'w');
        set(leg, 'Color', 'k');
    end
    set(gcf, 'Color', 'k');
    set(gca, 'Color', 'k');
    set(gca, 'xColor', 'w');
    set(gca, 'yColor', 'w');
    set(gca, 'zColor', 'w');
end


% PERFORM STATISTICAL TEST

% select only graphs with 1 data point before and 1 data point after MIS
testCurve = accuracyCurve;
removeIdx = [];
middle = round(size(accuracyCurve,2)/2);
for graphIdx = 1:size(accuracyCurve, 1)
    if(any(isnan(accuracyCurve(graphIdx, [middle-1 middle middle+1]))))
        removeIdx = [removeIdx graphIdx];
    end
end
testCurve(removeIdx,:) = [];

% compute bias of fitted sigmoid for each curve
testBiases = nan(1, size(testCurve,1));
x_full = (1:size(testCurve,2))-middle;
for graphIdx = 1:size(testCurve, 1)
    
    nanIdx = isnan(testCurve(graphIdx,:));
    x = x_full;
    y = testCurve(graphIdx,:);
    x(nanIdx) = [];
    y(nanIdx) = [];
    
    performanceMax = max(y);
    performanceMin = min(y);
    
    [param] = sigm_fit(x,y, [performanceMin, performanceMax , NaN , NaN], [performanceMin performanceMax 0 -2]);
    testBiases(graphIdx) = param(3);
    
end

% perform two one-sided t-tests
disp('Fig1: Testing sigmoidal fit to performance curve:');
[H,P,CI, STATS] = ttest(testBiases, 0, 'tail', 'right');
disp(['t-test bias above 0, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
[H,P,CI, STATS] = ttest(testBiases, 1, 'tail', 'left');
disp(['t-test bias below 1, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p = ' num2str(P)]);
    

%% Fig 2: Parallel Processing Performance as a Function of Number of Task Dependencies
fig1 = figure;
set(fig1, 'Position', [470 400 400 250]);
performanceMetric = 'LCA';         

bardata_mean = depencencyLogLCA_mean;       % first row: good MSE, 2nc column: bad MSE, rows indicate capacities
bardata_sem = depencencyLogLCA_sem;

bardata_mean(end,:) = [];
bardata_sem(end,:) = [];
legText = {};

plotColors = colormap(copper(size(bardata_mean,2)));
colorOffset = 0.2;
for row = 1:size(plotColors, 1)
    for col = 1:size(plotColors, 2)
        plotColors(row, col) = min(plotColors(row, col) + colorOffset, 1);
    end
end
plotColors = plotColors(fliplr(1:size(plotColors,1)),:);


if(~isempty(bardata_mean))
    xlim([0.5 size(bardata_mean,1)+0.5]);
    bar_handle = errorbar_groups(bardata_mean', bardata_sem','FigID', fig1, 'bar_colors', plotColors,'errorbar_width',0.5,'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',1.5});
    set(gca,'XTickLabel',{' ', ' ',' ',' '},'FontSize', fontSize_gca, 'FontName', fontName);
    performanceMetricLabel = performanceMetric;
    performanceMetricLabel(1) = upper(performanceMetric(1));
    if(strcmp(performanceMetric, 'LCA')) 
        performanceMetricLabel = {'Accuracy (%)'};
        ylimit = [0 1];
    end
    ylabel(performanceMetricLabel,'FontSize', fontSize_ylabel, 'FontName', fontName);
    xlabel('Number of Tasks Performed','FontSize', fontSize_xlabel);
    
    bardata_sum = nansum(bardata_mean);
    lastNonNan = find(bardata_sum > 0, 1, 'last');
    for legIdx = 1:lastNonNan
        legText{legIdx} = [num2str(legIdx-1) ' Dependencies Among Tasks'];
    end
    ylim(ylimit);
    set(gca, 'XTickLabel', 1:size(bardata_mean,1));
    set(gcf, 'Color', [1 1 1])
    
    % draw in chance performance
    hold on;
    for set_size = 1:size(bardata_mean, 1)
        chance = (1/NFeatures)^set_size;
        x = [0.5 6.5] + (set_size-1)*6;
        y = [chance, chance];
        plot(x, y, '--k', 'LineWidth', 2);
    end
    hold off;
    legText{end+1} = 'Chance Performance';
    [l, objh] = legend(legText, 'Location', 'northoutside');
    lineh = findobj(objh,'type','line')';
    set(lineh,'color','k')
    set(l, 'FontName', fontName);
end

% regression
disp('Fig2: Regression across all task set sizes:');

X_dependencies = repmat([1: size(depencencyLogAvgLCA,2)], 1, 1, numRepetitions);

X = X_dependencies(:);
Y = depencencyLogAvgLCA(:);
toRemove = isnan(Y);
X(toRemove) = [];
Y(toRemove) = [];
lm = fitlm(X,Y,'linear');
an = anova(lm,'summary');
disp(['regression: performance  ~ dependencies, b = ' num2str(lm.Coefficients.Estimate(2))...
                                                                                                    ', t(' num2str(an{1, 2}) ...
                                                                                                    ') = ' num2str(lm.Coefficients.tStat(2)) ...
                                                                                                    ', p = ' num2str(lm.Coefficients.pValue(2))]);

                                                                                                

