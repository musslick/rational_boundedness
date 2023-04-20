%% data settings

clear all;
clc;

datafilePrefix = 'PychReview_Part2_Sim2_BasisVsTensorTraining';   % 100 % feature overlap, weight initialization

performanceMeasure = 'MSE';

iterations_MultiTrain = 1000;
multiTrainingProportion = [0 0.2 0.4 0.6 0.8];
thresh = 0.01;
NPathways = 3;

% load data
% logFolder = 'logfiles/';
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
            batch_log_tmp = repmat(batch_log, numRepetitions, 1);
        end
        batch_log_tmp(fileIdx) = batch_log;
    end
else
    error('No valid file names found');
end
batch_log = batch_log_tmp;

numMultiConditions = length(multiTrainingProportion);
selectedMultiConditions = multiTrainingProportion;

% plot settings

plotSettings;

colorGradientCon = zeros(length(selectedMultiConditions), 3);
colorGradientCon(:,1) = linspace(0, colors(cContrast1,1), length(selectedMultiConditions));
colorGradientCon(:,2) = linspace(0, colors(cContrast1,2), length(selectedMultiConditions));
colorGradientCon(:,3) = linspace(0, colors(cContrast1,3), length(selectedMultiConditions));

colorGradientInc = zeros(length(selectedMultiConditions), 3);
colorGradientInc(:,1) = linspace(0, colors(cContrast2,1), length(selectedMultiConditions));
colorGradientInc(:,2) = linspace(0, colors(cContrast2,2), length(selectedMultiConditions));
colorGradientInc(:,3) = linspace(0, colors(cContrast2,3), length(selectedMultiConditions));


% plot graph analysis results

close all;


%% prepare data

Performance_MultitaskConTraining_Single_timeCourse = nan(numRepetitions, numMultiConditions, iterations_MultiTrain);
Performance_MultitaskIncTraining_Single_timeCourse = nan(numRepetitions, numMultiConditions, iterations_MultiTrain);
Performance_MultitaskConTraining_Multi_timeCourse = nan(numRepetitions, numMultiConditions, iterations_MultiTrain);
Performance_MultitaskIncTraining_Multi_timeCourse = nan(numRepetitions, numMultiConditions, iterations_MultiTrain);

TaskCorr_MultitaskConTraining = nan(numRepetitions, numMultiConditions, iterations_MultiTrain);
TaskCorr_MultitaskIncTraining = nan(numRepetitions, numMultiConditions, iterations_MultiTrain);

MultitaskConTraining_LCA_Multi = nan(numRepetitions, numMultiConditions);
MultitaskIncTraining_LCA_Multi = nan(numRepetitions, numMultiConditions);

singleTaskLearning_MultiCon = nan(numRepetitions, numMultiConditions);
singleTaskLearning_MultiInc = nan(numRepetitions, numMultiConditions);

for rep = 1:numRepetitions
    
    for multiCond = 1:numMultiConditions

        MultitaskConTraining_similarTasksCorr_hidden(rep, multiCond, :) = batch_log(rep).MultitaskConTraining_similarTasksCorr_hidden(multiCond,:);
        MultitaskIncTraining_similarTasksCorr_hidden(rep, multiCond, :) = batch_log(rep).MultitaskIncTraining_similarTasksCorr_hidden(multiCond,:);

        MultitaskConTraining_MDS_hidden_Phase1(rep, multiCond, :, :) = batch_log(rep).MultitaskConTraining_MDS_hidden_Phase1{multiCond};
        MultitaskIncTraining_MDS_hidden_Phase1(rep, multiCond, :, :) = batch_log(rep).MultitaskIncTraining_MDS_hidden_Phase1{multiCond};

        MultitaskConTraining_LCA_Multi(rep, multiCond) = batch_log(rep).MultitaskConTraining_LCA_Multi(multiCond);
        MultitaskIncTraining_LCA_Multi(rep, multiCond) = batch_log(rep).MultitaskIncTraining_LCA_Multi(multiCond);
        
    switch performanceMeasure   
        
        case 'MSE'
  
            Performance_MultitaskConTraining_Single_timeCourse(rep, multiCond, :) = batch_log(rep).MultitaskConTraining_MSE_Single_timeCourse(multiCond,:);
            Performance_MultitaskIncTraining_Single_timeCourse(rep, multiCond, :) = batch_log(rep).MultitaskIncTraining_MSE_Single_timeCourse(multiCond,:);
            Performance_MultitaskConTraining_Multi_timeCourse(rep, multiCond, :) = batch_log(rep).MultitaskConTraining_MSE_Multi_timeCourse(multiCond,:);
            Performance_MultitaskIncTraining_Multi_timeCourse(rep, multiCond, :) = batch_log(rep).MultitaskIncTraining_MSE_Multi_timeCourse(multiCond,:);
            
            singleTaskLearning_MultiCon(rep, multiCond) = find(batch_log(rep).MultitaskConTraining_MSE_Single_timeCourse(multiCond,:) <= thresh, 1, 'first');
            singleTaskLearning_MultiInc(rep, multiCond) = find(batch_log(rep).MultitaskIncTraining_MSE_Single_timeCourse(multiCond,:) <= thresh, 1, 'first');

    end
    
    end
    
end

%% Fig 1:  similarity between critical tasks over the course of multitask training 

cla;
plotSEM = 0;

fig = figure(1);
set(fig, 'Position', [100 100 900 250]);
x = 1:iterations_MultiTrain;

subplot(1,2,1);

legendText = {};

for selectedMultiCond = 1:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    y_MultiCon = squeeze(nanmean(MultitaskConTraining_similarTasksCorr_hidden(:, multiCond,:),1));
    y_MultiCon_sem = squeeze(nanstd(MultitaskConTraining_similarTasksCorr_hidden(:, multiCond,:),[],1)) / sqrt(numRepetitions);

    if(plotSEM)
        errorbar(x, y_MultiCon(x), y_MultiCon_sem(x), '-', 'LineWidth', lineWidth-2, 'Color', colorGradientCon(selectedMultiCond,:));
    else
        plot(x, y_MultiCon(x), '-', 'LineWidth', lineWidth+1, 'Color', colorGradientCon(selectedMultiCond,:));
    end
    hold on;
    
    legendText{selectedMultiCond} = [num2str(sglProp) '% Single Task / ' num2str(multiProp) '% Multitask Training'];

end
hold off;

ylim([0 1]);

ylabel('Task Representation Similarity', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Training Iterations','FontSize', fontSize_xlabel, 'FontName', fontName);
title('Mixed Congruent Multitask Training','FontSize', fontSize_title, 'FontName', fontName);
                
l = legend(legendText,'Location', 'southwest');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend-3);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca);

% INCONGRUENT

subplot(1,2,2);

legendText = {};

for selectedMultiCond = 1:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    y_MultiCon = squeeze(nanmean(MultitaskIncTraining_similarTasksCorr_hidden(:, multiCond,:),1));
    y_MultiCon_sem = squeeze(nanstd(MultitaskIncTraining_similarTasksCorr_hidden(:, multiCond,:),[],1)) / sqrt(numRepetitions);

    if(plotSEM)
        errorbar(x, y_MultiCon(x), y_MultiCon_sem(x), '-', 'LineWidth', lineWidth-2, 'Color', colorGradientInc(selectedMultiCond,:));
    else
        plot(x, y_MultiCon(x), '-', 'LineWidth', lineWidth+1, 'Color', colorGradientInc(selectedMultiCond,:));
    end
    hold on;
    
    legendText{selectedMultiCond} = [num2str(sglProp) '% Single Task / ' num2str(multiProp) '% Multitask Training'];

end
hold off;

ylim([0 1]);

ylabel('Task Representation Similarity', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Training Iterations','FontSize', fontSize_xlabel, 'FontName', fontName);
title('Mixed Incongruent Multitask Training','FontSize', fontSize_title, 'FontName', fontName);
                
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca);

%% Fig 2:  single task performance over the course of multitask training 

cla;
plotSEM = 0;

fig = figure(1);
set(fig, 'Position', [100 100 900 250]);
x = 1:iterations_MultiTrain;

subplot(1,2,1);

legendText = {};

for selectedMultiCond = 1:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;

    y_MultiCon = squeeze(nanmean(Performance_MultitaskConTraining_Single_timeCourse(:,multiCond,:),1));
    y_MultiCon_sem = squeeze(nanstd(Performance_MultitaskConTraining_Single_timeCourse(:,multiCond,:),[],1)) / sqrt(numRepetitions);

    if(plotSEM)
        errorbar(x, y_MultiCon(x), y_MultiCon_sem(x), '-', 'LineWidth', lineWidth-2, 'Color', colorGradient(selectedMultiCond,:));
    else
        plot(x, y_MultiCon(x), '-', 'LineWidth', lineWidth+1, 'Color', colorGradientCon(selectedMultiCond,:));
    end
    hold on;
    
    legendText{selectedMultiCond} = [num2str(sglProp) '% Single Task / ' num2str(multiProp) '% Multitask Training'];

end
hold off;

ylabel('MSE For Single Tasks', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Training Iterations','FontSize', fontSize_xlabel, 'FontName', fontName);
                
l = legend(legendText,'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca);


subplot(1,2,2);

legendText = {};

for selectedMultiCond = 1:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;

    y_MultiCon = squeeze(nanmean(Performance_MultitaskIncTraining_Single_timeCourse(:,multiCond,:),1));
    y_MultiCon_sem = squeeze(nanstd(Performance_MultitaskIncTraining_Single_timeCourse(:,multiCond,:),[],1)) / sqrt(numRepetitions);

    if(plotSEM)
        errorbar(x, y_MultiCon(x), y_MultiCon_sem(x), '-', 'LineWidth', lineWidth-2, 'Color', colorGradientInc(selectedMultiCond,:));
    else
        plot(x, y_MultiCon(x), '-', 'LineWidth', lineWidth+1, 'Color', colorGradientInc(selectedMultiCond,:));
    end
    hold on;
    
    legendText{selectedMultiCond} = [num2str(sglProp) '% Single Task / ' num2str(multiProp) '% Multitask Training'];

end
hold off;

ylabel('MSE For Single Tasks', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Training Iterations','FontSize', fontSize_xlabel, 'FontName', fontName);
                
l = legend(legendText,'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca);

%% Fig 3:  multitasking performance over the course of multitask training 

cla;
plotSEM = 0;

fig = figure(1);
set(fig, 'Position', [100 100 900 250]);
x = 1:iterations_MultiTrain;

subplot(1,2,1);

legendText = {};

for selectedMultiCond = 1:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;

    y_MultiCon = squeeze(nanmean(Performance_MultitaskConTraining_Multi_timeCourse(:,multiCond,:),1));
    y_MultiCon_sem = squeeze(nanstd(Performance_MultitaskConTraining_Multi_timeCourse(:,multiCond,:),[],1)) / sqrt(numRepetitions);

    if(plotSEM)
        errorbar(x, y_MultiCon(x), y_MultiCon_sem(x), '-', 'LineWidth', lineWidth-2, 'Color', colorGradientInc(selectedMultiCond,:));
    else
        plot(x, y_MultiCon(x), '-', 'LineWidth', lineWidth+1, 'Color', colorGradientInc(selectedMultiCond,:));
    end
    hold on;
    
    legendText{selectedMultiCond} = [num2str(sglProp) '% Single Task / ' num2str(multiProp) '% Multitask Training'];

end
hold off;

ylabel('MSE For Multitasking', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Training Iterations','FontSize', fontSize_xlabel, 'FontName', fontName);
                
l = legend(legendText,'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca);


%% Fig 4: MDS plot
width = 1000; % 200
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

% plot single task training;
subplot(2,5,1);

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

% plot congruent task training
for selectedMultiCond = 2:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    subplot(2,5,selectedMultiCond);

    x = MultitaskConTraining_MDS_hidden_Phase1(repetitionIdx,multiCond,:,1);
    y = MultitaskConTraining_MDS_hidden_Phase1(repetitionIdx,multiCond,:,2);

    scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth);
    ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
    xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
    title({[num2str(multiProp) '% Congruent']; 'Multitask Training'},'FontSize', fontSize_title, 'FontName', fontName, 'Color', colors(cSingle,:));

    xlim([xlimit(1) xlimit(2)]) 
    ylim([ylimit(1) ylimit(2)]);

end

% plot congruent task training
for selectedMultiCond = 2:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    subplot(2,5,5+selectedMultiCond);

    x = MultitaskIncTraining_MDS_hidden_Phase1(repetitionIdx,multiCond,:,1);
    y = MultitaskIncTraining_MDS_hidden_Phase1(repetitionIdx,multiCond,:,2);

    scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth);
    ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
    xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
    title({[num2str(multiProp) '% Incongruent']; 'Multitask Training'},'FontSize', fontSize_title, 'FontName', fontName, 'Color', colors(cSingle,:));

    xlim([xlimit(1) xlimit(2)]) 
    ylim([ylimit(1) ylimit(2)]);

end

%% FIG 5: Learning iterations
width = 280;
height = 280;
fontSize_title = 11;
fontSize_xlabel = 16;
fontSize_ylabel = fontSize_xlabel;
fontSize_legend = 12;

fig = figure(1);
set(fig, 'Position', [100 100 width height]);

% plot learning iterations

x_multi = selectedMultiConditions(2:end)*100;
x_sgl = 0;
y_sgl = nanmean(singleTaskLearning_MultiCon(:,1),1);
y_sgl_sem = nanstd(singleTaskLearning_MultiCon(:,1),1)/sqrt(numRepetitions);
y_con = nanmean(singleTaskLearning_MultiCon(:,2:end),1);
y_con_sem = nanstd(singleTaskLearning_MultiCon(:,2:end),1)/sqrt(numRepetitions);
y_inc = nanmean(singleTaskLearning_MultiInc(:,2:end),1);
y_inc_sem = nanstd(singleTaskLearning_MultiInc(:,2:end),1)/sqrt(numRepetitions);

errorbar(x_sgl, y_sgl, y_sgl_sem, 'o', 'LineWidth', lineWidth-1, 'Color', [0.65 0.65 0.65]); hold on;
errorbar(x_multi, y_con, y_con_sem, '-', 'LineWidth', lineWidth-1, 'Color', [0.3 0.3 0.3]); hold on;
errorbar(x_multi, y_inc, y_inc_sem, '-', 'LineWidth', lineWidth-1, 'Color', [0 0 0]); hold off;

xlim([-5 100]);
ylim([0 500]);
leg = legend({'100 % Single Task Training', 'Congr. Multitask Training', 'Incongr. Multitask Training'},'Location', 'northwest');
set(leg, 'FontSize', fontSize_legend, 'FontName', fontName);
set(gca, 'FontSize', fontSize_ylabel, 'FontName', fontName);

xlabel('% Multitask Training','FontSize', fontSize_ylabel, 'FontName', fontName, 'Color', 'k');
ylabel('Iterations Required to Train','FontSize', fontSize_xlabel, 'FontName', fontName, 'Color', 'k');

%% FIG 5: Multitasking Accuracy
width = 280;
height = 280;
fontSize_title = 11;
fontSize_xlabel = 16;
fontSize_ylabel = fontSize_xlabel;
fontSize_legend = 12;

fig = figure(1);
set(fig, 'Position', [100 100 width height]);

% plot LCA multitasking accuracy

y_sgl = nanmean(MultitaskConTraining_LCA_Multi(:,1),1)*100;
y_sgl_sem = nanstd(MultitaskConTraining_LCA_Multi(:,1),1)*100 /sqrt(numRepetitions);
y_con = nanmean(MultitaskConTraining_LCA_Multi(:,2:end),1)*100 ;
y_con_sem = nanstd(MultitaskConTraining_LCA_Multi(:,2:end))*100 /sqrt(numRepetitions);
y_inc = nanmean(MultitaskIncTraining_LCA_Multi(:,2:end),1)*100;
y_inc_sem = nanstd(MultitaskIncTraining_LCA_Multi(:,2:end),1)*100 /sqrt(numRepetitions);

errorbar(x_sgl, y_sgl, y_sgl_sem, 'o', 'LineWidth', lineWidth-1, 'Color', [0.65 0.65 0.65]); hold on;
errorbar(x_multi, y_con, y_con_sem, '-', 'LineWidth', lineWidth-1, 'Color', [0.3 0.3 0.3]); hold on;
errorbar(x_multi, y_inc, y_inc_sem, '-', 'LineWidth', lineWidth-1, 'Color', [0 0 0]); hold off;

xlim([-5 100]);
ylim([0 100]);
leg = legend({'100 % Single Task Training', 'Congr. Multitask Training', 'Incongr. Multitask Training'},'Location', 'southeast');
set(leg, 'FontSize', fontSize_legend, 'FontName', fontName);
set(gca, 'FontSize', fontSize_ylabel, 'FontName', fontName);

xlabel('% Multitask Training','FontSize', fontSize_ylabel, 'FontName', fontName, 'Color', 'k');
ylabel('Multitasking Accuracy (%)','FontSize', fontSize_xlabel, 'FontName', fontName, 'Color', 'k');

%% plot mean correlations and sds for representational similarity
clc;

y_sgl = nanmean(MultitaskConTraining_similarTasksCorr_hidden(:,1,end),1);
y_sgl_std = nanstd(MultitaskConTraining_similarTasksCorr_hidden(:,1, end),1);
y_con = nanmean(MultitaskConTraining_similarTasksCorr_hidden(:,2:end, end),1);
y_con_std = nanstd(MultitaskConTraining_similarTasksCorr_hidden(:,2:end, end));
y_inc = nanmean(MultitaskIncTraining_similarTasksCorr_hidden(:,2:end, end),1);
y_inc_std = nanstd(MultitaskIncTraining_similarTasksCorr_hidden(:,2:end, end),1);

disp('LEARNING ITERATIONS');
disp(['0% con multi: M = ' num2str(round(y_sgl*100)/100) ', SD = ' num2str(round(y_sgl_std*100)/100)]);
disp('--');
for selectedMultiCond = 2:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    disp([num2str(multiProp) '% con multi: M = ' num2str(round(y_con(multiCond-1)*100)/100) ', SD = ' num2str(round(y_con_std(multiCond-1)*100)/100)]);
end
disp('--');
for selectedMultiCond = 2:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    disp([num2str(multiProp) '% inc multi: M = ' num2str(round(y_inc(multiCond-1)*100)/100) ', SD = ' num2str(round(y_inc_std(multiCond-1)*100)/100)]);
end

%% plot mean correlations and sds for final LCA multitasking accuracy
clc;

disp('------');
disp('MULTITASKING ACCURACY');

y_sgl = nanmean(MultitaskConTraining_LCA_Multi(:,1),1);
y_sgl_std = nanstd(MultitaskConTraining_LCA_Multi(:,1),1);
y_con = nanmean(MultitaskConTraining_LCA_Multi(:,2:end),1);
y_con_std = nanstd(MultitaskConTraining_LCA_Multi(:,2:end));
y_inc = nanmean(MultitaskIncTraining_LCA_Multi(:,2:end),1);
y_inc_std = nanstd(MultitaskIncTraining_LCA_Multi(:,2:end),1);

disp(['0% con multi: M = ' num2str(round(y_sgl*1000)/1000) ', SD = ' num2str(round(y_sgl_std*1000)/1000)]);
disp('--');
for selectedMultiCond = 2:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    disp([num2str(multiProp) '% con multi: M = ' num2str(round(y_con(multiCond-1)*1000)/1000) ', SD = ' num2str(round(y_con_std(multiCond-1)*1000)/1000)]);
end
disp('--');
for selectedMultiCond = 2:length(selectedMultiConditions)
    multiCond = find(multiTrainingProportion == selectedMultiConditions(selectedMultiCond));
    multiProp = multiTrainingProportion(multiCond)*100;
    sglProp = 100-multiProp;
    
    disp([num2str(multiProp) '% inc multi: M = ' num2str(round(y_inc(multiCond-1)*1000)/1000) ', SD = ' num2str(round(y_inc_std(multiCond-1)*1000)/1000)]);
end
