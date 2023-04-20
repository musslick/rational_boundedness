%% data settings

%  Manuscript FIG_DEPENDENCY_RESULTS

clear all;
clc;
close all;


datafilePrefix = 'PsychReview_Part1_Sim1_sharedRepresentation_Suppl_5P3F_10tasks_23_h100_r';

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

output_net_task = nan(numRepetitions, NPathways);

for rep = 1:numRepetitions
    
    avg_output_net_task = batch_log(numHiddenIdx, rep).avg_output_net_task;
   
    for i=1:length(avg_output_net_task)
        output_net_task(rep, i) = avg_output_net_task(i);
    end
end

plotSettings;
fontSize_gca = fontSize_gca-2;

fig1 = figure;
set(fig1, 'Position', [470 400 400 250]);
x = 1:size(output_net_task,2);
y_mean = mean(output_net_task);
y_sem = nanstd(output_net_task, [], 1) / sqrt(numRepetitions);

errorbar(x, y_mean, y_sem, '-k', 'LineWidth',3); hold on;
ylabel({'Average Net Input', 'from Task Layer'}, 'fontSize', fontSize_ylabel-2);
xlabel({'Number of Tasks Performed'}, 'fontSize', fontSize_xlabel-2);
set(gca, 'FontSize', fontSize_gca);