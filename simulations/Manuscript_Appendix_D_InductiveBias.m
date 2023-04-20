function Manuscript_Simulation_4_InductiveBias(inputArg)

inputArgDouble = str2double(inputArg);
rep = inputArgDouble;

tic

% meta simulation parameters
log_version = 3;

% set up network parameters
nHidden = 100;              % number of hidden units
hiddenPathSize = 1;         % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;         % group size of output units that receive the same weights from the task layer)
learningRate = 0.1;         % learning rate
thresh = 0.01;            % mean-squared error stopping criterion
decay = 0.0000;             % weight penalization parameter
bias = -4;                  % weight from bias units to hidden & output units
init_scale = 0.1;           % scales for initialized random weights
iterations_maxTrain = 5000;     % number of training iterations for multitask training
replications = 20;
MDSMetric = 'euclidean'; % 'euclidean', 'correlation', 'cosine'

% task environment parameters 
NPathways = 6;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
NFeatures = 3;                     % the number of features per feature dimension

sameStimuliAcrossTasks = 1;
sdScale = 0;
samplesPerTask = [];

% create 3 different initial training environments
tasksToPerform = [1];
[inputSgl_1shared1, tasksSgl_1shared1, trainSgl_1shared1] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);
tasksToPerform = [1 8];
[inputSgl_1shared2, tasksSgl_1shared2, trainSgl_1shared2] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);
tasksToPerform = [1 8 15];
[inputSgl_1shared3, tasksSgl_1shared3, trainSgl_1shared3] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);
% tasksToPerform = [1 2 8 9 13 15];
% [inputSgl_2shared, tasksSgl_2shared, trainSgl_2shared] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);
% tasksToPerform = [1 2 3 7 8 9 13 14 15];
% [inputSgl_3shared, tasksSgl_3shared, trainSgl_3shared] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);

% create test training environment 
tasksToPerform = [4 11 18];
[inputSgl_tested, tasksSgl_tested, trainSgl_tested] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);
nTasks_0shared = length(tasksToPerform);
tasksToPerform_0shared = tasksToPerform;

tasksToPerform_1 = [4];
[inputSgl_tested_1, tasksSgl_tested_1, trainSgl_tested_1] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform_1, 0);
tasksToPerform_1_shared_1 = tasksToPerform_1;

tasksToPerform_2 = [11];
[inputSgl_tested_2, tasksSgl_tested_2, trainSgl_tested_2] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform_2, 0);
tasksToPerform_1_shared_2 = tasksToPerform_2;

tasksToPerform_3 = [18];
[inputSgl_tested_3, tasksSgl_tested_3, trainSgl_tested_3] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform_3, 0);
tasksToPerform_1_shared_3 = tasksToPerform_3;

% create full environments
tasksToPerform = 1:(NPathways*NPathways);
[inputSgl, tasksSgl, trainSgl] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, 0);

tasksToPerform_1shared1 = [1 4 11 18];
nTasks_1shared1 = length(tasksToPerform_1shared1);

tasksToPerform_1shared2 = [1 8 4 11 18];
nTasks_1shared2 = length(tasksToPerform_1shared2);

tasksToPerform_1shared3 = [1 8 15 4 11 18];
nTasks_1shared3 = length(tasksToPerform_1shared3);

%%  TRAINING

for rep = 1:replications

    % initialize network
    taskNet = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
    taskNet.setData(inputSgl_1shared3, tasksSgl_1shared3, trainSgl_1shared3);
    taskNet.configure(); 
    taskNet.NPathways = NPathways;

    % generate 3 network conditions
    taskNet_0shared = NNmodel(taskNet);
    taskNet_1shared1 = NNmodel(taskNet);
    taskNet_1shared2 = NNmodel(taskNet);
    taskNet_1shared3 = NNmodel(taskNet);
%     taskNet_2shared = NNmodel(taskNet);
%     taskNet_3shared = NNmodel(taskNet);

    % pretrain all three networks until cirterion
    taskNet_1shared1.setData(inputSgl_1shared1, tasksSgl_1shared1, trainSgl_1shared1);
    taskNet_1shared2.setData(inputSgl_1shared2, tasksSgl_1shared2, trainSgl_1shared2);
    taskNet_1shared3.setData(inputSgl_1shared3, tasksSgl_1shared3, trainSgl_1shared3);
%     taskNet_2shared.setData(inputSgl_2shared, tasksSgl_2shared, trainSgl_2shared);
%     taskNet_3shared.setData(inputSgl_3shared, tasksSgl_3shared, trainSgl_3shared);

    disp('single task pretraining... 1 shared 1');
    taskNet_1shared1.trainOnline(iterations_maxTrain);
    disp('single task pretraining... 1 shared 2');
    taskNet_1shared2.trainOnline(iterations_maxTrain);
    disp('single task pretraining... 1 shared 3');
    taskNet_1shared3.trainOnline(iterations_maxTrain);
%     disp('single task pretraining... 2 shared');
%     taskNet_2shared.trainOnline(iterations_maxTrain);
%     disp('single task pretraining... 3 shared');
%     taskNet_3shared.trainOnline(iterations_maxTrain);

    % train networks on final training pattern
    taskNet_0shared.setData(inputSgl_tested, tasksSgl_tested, trainSgl_tested);
    taskNet_1shared1.setData(inputSgl_tested, tasksSgl_tested, trainSgl_tested);
    taskNet_1shared2.setData(inputSgl_tested, tasksSgl_tested, trainSgl_tested);
    taskNet_1shared3.setData(inputSgl_tested, tasksSgl_tested, trainSgl_tested);
%     taskNet_2shared.setData(inputSgl_tested, tasksSgl_tested, trainSgl_tested);
%     taskNet_3shared.setData(inputSgl_tested, tasksSgl_tested, trainSgl_tested);

    disp('single task training... 0 shared');
    %taskNet_0shared.trainOnline(iterations_maxTrain);
    MSE_log_taskNet_0_shared_1 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_0_shared_2 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_0_shared_3 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_0_shared = zeros(1, iterations_maxTrain);
    
    for i = 1:iterations_maxTrain
        taskNet_0shared.trainOnline(1);
        
        [~, ~, MSE] = taskNet_0shared.runSet(inputSgl_tested_1, tasksSgl_tested_1, trainSgl_tested_1);
        MSE_log_taskNet_0_shared_1(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_0shared.runSet(inputSgl_tested_2, tasksSgl_tested_2, trainSgl_tested_2);
        MSE_log_taskNet_0_shared_2(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_0shared.runSet(inputSgl_tested_3, tasksSgl_tested_3, trainSgl_tested_3);
        MSE_log_taskNet_0_shared_3(i) = mean(MSE);

        MSE_log_taskNet_0_shared(i) = taskNet_0shared.MSE_log;
        
        if taskNet_0shared.MSE_log < thresh
            break
        end
    end
    
    disp('single task training... 1 shared 1');
    %taskNet_1shared1.trainOnline(iterations_maxTrain);
    MSE_log_taskNet_1_shared_1 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_1_shared_2 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_1_shared_3 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_1_shared = zeros(1, iterations_maxTrain);
    
    for i = 1:iterations_maxTrain
        taskNet_1shared1.trainOnline(1);
        
        [~, ~, MSE] = taskNet_1shared1.runSet(inputSgl_tested_1, tasksSgl_tested_1, trainSgl_tested_1);
        MSE_log_taskNet_1_shared_1(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_1shared1.runSet(inputSgl_tested_2, tasksSgl_tested_2, trainSgl_tested_2);
        MSE_log_taskNet_1_shared_2(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_1shared1.runSet(inputSgl_tested_3, tasksSgl_tested_3, trainSgl_tested_3);
        MSE_log_taskNet_1_shared_3(i) = mean(MSE);

        MSE_log_taskNet_1_shared(i) = taskNet_1shared1.MSE_log;
        
        if taskNet_1shared1.MSE_log < thresh
            break
        end
    end
    
    disp('single task training... 1 shared 2');
    %taskNet_1shared2.trainOnline(iterations_maxTrain);
    MSE_log_taskNet_2_shared_1 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_2_shared_2 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_2_shared_3 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_2_shared = zeros(1, iterations_maxTrain);
    
    for i = 1:iterations_maxTrain
        taskNet_1shared2.trainOnline(1);
        
        [~, ~, MSE] = taskNet_1shared2.runSet(inputSgl_tested_1, tasksSgl_tested_1, trainSgl_tested_1);
        MSE_log_taskNet_2_shared_1(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_1shared2.runSet(inputSgl_tested_2, tasksSgl_tested_2, trainSgl_tested_2);
        MSE_log_taskNet_2_shared_2(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_1shared2.runSet(inputSgl_tested_3, tasksSgl_tested_3, trainSgl_tested_3);
        MSE_log_taskNet_2_shared_3(i) = mean(MSE);

        MSE_log_taskNet_2_shared(i) = taskNet_1shared2.MSE_log;
        
        if taskNet_1shared2.MSE_log < thresh
            break
        end
    end
    
    disp('single task training... 1 shared 3');
    %taskNet_1shared3.trainOnline(iterations_maxTrain);
    MSE_log_taskNet_3_shared_1 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_3_shared_2 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_3_shared_3 = zeros(1, iterations_maxTrain);
    MSE_log_taskNet_3_shared = zeros(1, iterations_maxTrain);
    
    for i = 1:iterations_maxTrain
        taskNet_1shared3.trainOnline(1);
        
        [~, ~, MSE] = taskNet_1shared3.runSet(inputSgl_tested_1, tasksSgl_tested_1, trainSgl_tested_1);
        MSE_log_taskNet_3_shared_1(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_1shared3.runSet(inputSgl_tested_2, tasksSgl_tested_2, trainSgl_tested_2);
        MSE_log_taskNet_3_shared_2(i) = mean(MSE);
        
        [~, ~, MSE] = taskNet_1shared3.runSet(inputSgl_tested_3, tasksSgl_tested_3, trainSgl_tested_3);
        MSE_log_taskNet_3_shared_3(i) = mean(MSE);

        MSE_log_taskNet_3_shared(i) = taskNet_1shared3.MSE_log;
        
        if taskNet_1shared3.MSE_log < thresh
            break
        end
    end
       
%     disp('single task training... 2 shared');
%     taskNet_2shared.trainOnline(iterations_maxTrain);
%     disp('single task training... 3 shared');
%     taskNet_3shared.trainOnline(iterations_maxTrain);

    batch_log.MSE_log_0shared(rep,:) = MSE_log_taskNet_0_shared;
    batch_log.MSE_log_1shared1(rep,:) = MSE_log_taskNet_1_shared;
    batch_log.MSE_log_1shared2(rep,:) = MSE_log_taskNet_2_shared;
    batch_log.MSE_log_1shared3(rep,:) = MSE_log_taskNet_3_shared;
    
    batch_log.MSE_log_taskNet_0_shared_1(rep, :) = MSE_log_taskNet_0_shared_1;
    batch_log.MSE_log_taskNet_1_shared_1(rep, :) = MSE_log_taskNet_1_shared_1;
    batch_log.MSE_log_taskNet_2_shared_1(rep, :) = MSE_log_taskNet_2_shared_1;
    batch_log.MSE_log_taskNet_3_shared_1(rep, :) = MSE_log_taskNet_3_shared_1;
    
    batch_log.MSE_log_taskNet_0_shared_2(rep, :) = MSE_log_taskNet_0_shared_2;
    batch_log.MSE_log_taskNet_1_shared_2(rep, :) = MSE_log_taskNet_1_shared_2;
    batch_log.MSE_log_taskNet_2_shared_2(rep, :) = MSE_log_taskNet_2_shared_2;
    batch_log.MSE_log_taskNet_3_shared_2(rep, :) = MSE_log_taskNet_3_shared_2;
    
    batch_log.MSE_log_taskNet_0_shared_3(rep, :) = MSE_log_taskNet_0_shared_3;
    batch_log.MSE_log_taskNet_1_shared_3(rep, :) = MSE_log_taskNet_1_shared_3;
    batch_log.MSE_log_taskNet_2_shared_3(rep, :) = MSE_log_taskNet_2_shared_3;
    batch_log.MSE_log_taskNet_3_shared_3(rep, :) = MSE_log_taskNet_3_shared_3;
%     batch_log.MSE_log_2shared(rep,:) = taskNet_2shared.MSE_log;
%     batch_log.MSE_log_3shared(rep,:) = taskNet_3shared.MSE_log;

    %% MDS ANALYSIS

    % perform MDS on hidden data 0 shared
    nTasks = nTasks_0shared;
    [~, ~, hiddenData] = computeTaskSimilarity(taskNet_0shared, inputSgl, tasksSgl, 'tasks', tasksToPerform_0shared);
    X = transpose(hiddenData);
    dist = MDSMetric;
    distances = pdist(X, dist);
    try
    Y = mdscale(distances,2);
    catch
       disp(input);
       Y = zeros(size(input,1), 2);
       warning('No MDS possible for average single task reps'); 
    end
    batch_log.MSE_hidden_0shared{rep} = Y(1:nTasks, :);
    
    % perform MDS on hidden data 1 shared 1
    nTasks = nTasks_1shared1;
    [~, ~, hiddenData] = computeTaskSimilarity(taskNet_1shared1, inputSgl, tasksSgl, 'tasks', tasksToPerform_1shared1);
    X = transpose(hiddenData);
    dist = MDSMetric;
    distances = pdist(X, dist);
    try
    Y = mdscale(distances,2);
    catch
       disp(input);
       Y = zeros(size(input,1), 2);
       warning('No MDS possible for average single task reps'); 
    end
    batch_log.MSE_hidden_1shared1{rep} = Y(1:nTasks, :);
    
    % perform MDS on hidden data 1 shared 2
    nTasks = nTasks_1shared2;
    [~, ~, hiddenData] = computeTaskSimilarity(taskNet_1shared2, inputSgl, tasksSgl, 'tasks', tasksToPerform_1shared2);
    X = transpose(hiddenData);
    dist = MDSMetric;
    distances = pdist(X, dist);
    try
    Y = mdscale(distances,2);
    catch
       disp(input);
       Y = zeros(size(input,1), 2);
       warning('No MDS possible for average single task reps'); 
    end
    batch_log.MSE_hidden_1shared2{rep} = Y(1:nTasks, :);
    
    % perform MDS on hidden data 1 shared 3
    nTasks = nTasks_1shared3;
    [~, ~, hiddenData] = computeTaskSimilarity(taskNet_1shared3, inputSgl, tasksSgl, 'tasks', tasksToPerform_1shared3);
    X = transpose(hiddenData);
    dist = MDSMetric;
    distances = pdist(X, dist);
    try
    Y = mdscale(distances,2);
    catch
       disp(input);
       Y = zeros(size(input,1), 2);
       warning('No MDS possible for average single task reps'); 
    end
    batch_log.MSE_hidden_1shared3{rep} = Y(1:nTasks, :);
    
    
end

%% save log file
toc

save(['logfiles/Part2/PychReview_Part2_Sim3_InductiveBias_V2_' num2str(NPathways) 'P' num2str(NFeatures) 'F_v' num2str(log_version) '_h' num2str(nHidden(1)) '_r' num2str(rep)], '-v7.3');

%% determine statistics

criterion = 0.01;

for rep = 1:size(batch_log.MSE_log_0shared,1)
    numIter_0shared(rep) = find(batch_log.MSE_log_0shared(rep,:) < criterion,1);
    numIter_1shared(rep) = find(batch_log.MSE_log_1shared1(rep,:) < criterion,1);
    numIter_2shared(rep) = find(batch_log.MSE_log_1shared2(rep,:) < criterion,1);
    numIter_3shared(rep) = find(batch_log.MSE_log_1shared3(rep,:) < criterion,1);
end

disp('----');
disp(['0 shared: M = ' num2str(mean(numIter_0shared)), ', SD = ' num2str(std(numIter_0shared))]);
disp(['1 shared: M = ' num2str(mean(numIter_1shared)), ', SD = ' num2str(std(numIter_1shared))]);
disp(['2 shared: M = ' num2str(mean(numIter_2shared)), ', SD = ' num2str(std(numIter_2shared))]);
disp(['3 shared: M = ' num2str(mean(numIter_3shared)), ', SD = ' num2str(std(numIter_3shared))]);

Y = [numIter_0shared'; numIter_1shared'; numIter_2shared'; numIter_3shared'];
X = [repmat(0,length(numIter_0shared),1); repmat(1,length(numIter_1shared),1); repmat(2,length(numIter_2shared),1); repmat(3,length(numIter_3shared),1)];
lm = fitlm(X,Y)

%% plot learning curves
load('logfiles/Part2/PychReview_Part2_Sim3_InductiveBias_V2_6P3F_v3_h100_r20');

plotSEM = 1;

fig1 = figure(1);
set(fig1, 'Position', [100 100 700 200]);

maxIterationsPlot = 30;

mean_0shared = mean(batch_log.MSE_log_0shared(:,1:maxIterationsPlot));
mean_1shared1 = mean(batch_log.MSE_log_1shared1(:,1:maxIterationsPlot));
mean_1shared2 = mean(batch_log.MSE_log_1shared2(:,1:maxIterationsPlot));
mean_1shared3 = mean(batch_log.MSE_log_1shared3(:,1:maxIterationsPlot));
% mean_2shared = mean(batch_log.MSE_log_2shared(:,1:maxIterationsPlot));
% mean_3shared = mean(batch_log.MSE_log_3shared(:,1:maxIterationsPlot));

sem_0shared = std(batch_log.MSE_log_0shared(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared1 = std(batch_log.MSE_log_1shared1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared2 = std(batch_log.MSE_log_1shared2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared3 = std(batch_log.MSE_log_1shared3(:,1:maxIterationsPlot)) / sqrt(replications);
% sem_2shared = std(batch_log.MSE_log_2shared(:,1:maxIterationsPlot)) / sqrt(replications);
% sem_3shared = std(batch_log.MSE_log_3shared(:,1:maxIterationsPlot)) / sqrt(replications);

plotSettings;

colors(5,:) = repmat(0.75,1,3);
colors(1,:) = repmat(0.5,1,3);
colors(7,:) = repmat(0.25,1,3);
colors(8,:) = repmat(0,1,3);

if(~plotSEM)
    plot(1:maxIterationsPlot, mean_0shared, 'LineWidth', 3, 'color', colors(5,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared2, 'LineWidth', 3, 'color', colors(7,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared3, 'LineWidth', 3, 'color', colors(8,:)); hold on;
%     plot(1:maxIterationsPlot, mean_2shared, 'LineWidth', 3, 'color', colors(2,:)); hold on;
%     plot(1:maxIterationsPlot, mean_3shared, 'LineWidth', 3, 'color', colors(3,:)); hold off;
else
    errorbar(1:maxIterationsPlot, mean_0shared, sem_0shared, 'LineWidth', 3, 'color', colors(5,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_1shared1, sem_1shared1,'LineWidth', 3, 'color', colors(1,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_1shared2, sem_1shared2, 'LineWidth', 3, 'color', colors(7,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_1shared3, sem_1shared3, 'LineWidth', 3, 'color', colors(8,:)); hold on;
%     errorbar(1:maxIterationsPlot, mean_2shared, sem_2shared, 'LineWidth', 3, 'color', colors(2,:)); hold on;
%     errorbar(1:maxIterationsPlot, mean_3shared, sem_3shared, 'LineWidth', 3, 'color', colors(3,:)); hold off;
end

xlabel('Training Iterations','FontSize', fontSize_xlabel-1);
ylabel({'MSE on', 'Target Single Tasks'}, 'FontSize', fontSize_ylabel-1, 'FontName', fontName);
l = legend('Pre-Trained on No Auxiliary Tasks', 'Pre-Trained on 1 Auxiliary Task', 'Pre-Trained on 2 Auxiliary Tasks', 'Pre-Trained on 3 Auxiliary Tasks','Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend-3);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca-1);

%% Individual Learning Curve plots
clf;
close all;
plotSEM = 1;

plotSettings;
width = 200;
height = 200;
fontSize_title = 11;
fontSize_xlabel = 12;
fontSize_ylabel = fontSize_xlabel;

maxIterationsPlot = 30;


mean_0shared1 = mean(batch_log.MSE_log_taskNet_0_shared_1(:,1:maxIterationsPlot));
mean_0shared2 = mean(batch_log.MSE_log_taskNet_0_shared_2(:,1:maxIterationsPlot));
mean_0shared3 = mean(batch_log.MSE_log_taskNet_0_shared_3(:,1:maxIterationsPlot));

mean_1shared1 = mean(batch_log.MSE_log_taskNet_1_shared_1(:,1:maxIterationsPlot));
mean_1shared2 = mean(batch_log.MSE_log_taskNet_1_shared_2(:,1:maxIterationsPlot));
mean_1shared3 = mean(batch_log.MSE_log_taskNet_1_shared_3(:,1:maxIterationsPlot));

mean_2shared1 = mean(batch_log.MSE_log_taskNet_2_shared_1(:,1:maxIterationsPlot));
mean_2shared2 = mean(batch_log.MSE_log_taskNet_2_shared_2(:,1:maxIterationsPlot));
mean_2shared3 = mean(batch_log.MSE_log_taskNet_2_shared_3(:,1:maxIterationsPlot));

mean_3shared1 = mean(batch_log.MSE_log_taskNet_3_shared_1(:,1:maxIterationsPlot));
mean_3shared2 = mean(batch_log.MSE_log_taskNet_3_shared_2(:,1:maxIterationsPlot));
mean_3shared3 = mean(batch_log.MSE_log_taskNet_3_shared_3(:,1:maxIterationsPlot));

sem_0shared1 = std(batch_log.MSE_log_taskNet_0_shared_1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_0shared2 = std(batch_log.MSE_log_taskNet_0_shared_2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_0shared3 = std(batch_log.MSE_log_taskNet_0_shared_3(:,1:maxIterationsPlot)) / sqrt(replications);

sem_1shared1 = std(batch_log.MSE_log_taskNet_1_shared_1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared2 = std(batch_log.MSE_log_taskNet_1_shared_2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared3 = std(batch_log.MSE_log_taskNet_1_shared_3(:,1:maxIterationsPlot)) / sqrt(replications);

sem_2shared1 = std(batch_log.MSE_log_taskNet_2_shared_1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_2shared2 = std(batch_log.MSE_log_taskNet_2_shared_2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_2shared3 = std(batch_log.MSE_log_taskNet_2_shared_3(:,1:maxIterationsPlot)) / sqrt(replications);

sem_3shared1 = std(batch_log.MSE_log_taskNet_3_shared_1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_3shared2 = std(batch_log.MSE_log_taskNet_3_shared_2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_3shared3 = std(batch_log.MSE_log_taskNet_3_shared_3(:,1:maxIterationsPlot)) / sqrt(replications);

mean_0shared = mean(batch_log.MSE_log_0shared(:,1:maxIterationsPlot));
mean_1shared = mean(batch_log.MSE_log_1shared1(:,1:maxIterationsPlot));
mean_2shared = mean(batch_log.MSE_log_1shared2(:,1:maxIterationsPlot));
mean_3shared = mean(batch_log.MSE_log_1shared3(:,1:maxIterationsPlot));

sem_0shared = std(batch_log.MSE_log_0shared(:,1:maxIterationsPlot)) / sqrt(replications);
sem_1shared = std(batch_log.MSE_log_1shared1(:,1:maxIterationsPlot)) / sqrt(replications);
sem_2shared = std(batch_log.MSE_log_1shared2(:,1:maxIterationsPlot)) / sqrt(replications);
sem_3shared = std(batch_log.MSE_log_1shared3(:,1:maxIterationsPlot)) / sqrt(replications);

color = nan(3,3);
color(1,:) = colors(cContrast1,:);
color(2,:) = colors(cContrast2,:);
color(3,:) = colors(cContrast3,:);

fig = figure(1);
set(fig, 'Position', [100 100 width height]);

if(~plotSEM)
    plot(1:maxIterationsPlot, mean_0shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    plot(1:maxIterationsPlot, mean_0shared2, 'LineWidth', 3, 'color', colors(2,:)); hold on;
    plot(1:maxIterationsPlot, mean_0shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
else
    errorbar(1:maxIterationsPlot, mean_0shared1, sem_0shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_0shared2, sem_0shared2,'LineWidth', 3, 'color', colors(2,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_0shared3, sem_0shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_0shared, sem_0shared, 'LineWidth', 3, 'color', colors(5,:)); 
end

xlabel('Training Iterations','FontSize', fontSize_xlabel-1);
ylabel({'MSE on', 'Single Targets Tasks'}, 'FontSize', fontSize_ylabel-1, 'FontName', fontName);
l = legend('Task A', 'Task B', 'Task C', 'Mean', 'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend-3);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca-1);

% 1 TASK TRAINED

fig = figure(2);
set(fig, 'Position', [300 100 width height]);

if(~plotSEM)
    plot(1:maxIterationsPlot, mean_1shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared2, 'LineWidth', 3, 'color', colors(2,:)); hold on;
    plot(1:maxIterationsPlot, mean_1shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
else
    errorbar(1:maxIterationsPlot, mean_1shared1, sem_1shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_1shared2, sem_1shared2,'LineWidth', 3, 'color', colors(2,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_1shared3, sem_1shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_1shared, sem_1shared, 'LineWidth', 3, 'color', colors(5,:)); 
end

xlabel('Training Iterations','FontSize', fontSize_xlabel-1);
ylabel({'MSE on', 'Single Targets Tasks'}, 'FontSize', fontSize_ylabel-1, 'FontName', fontName);
l = legend('Task A', 'Task B', 'Task C', 'Mean', 'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend-3);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca-1);

% 2 TASKS TRAINED

fig = figure(3);
set(fig, 'Position', [500 100 width height]);

if(~plotSEM)
    plot(1:maxIterationsPlot, mean_2shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    plot(1:maxIterationsPlot, mean_2shared2, 'LineWidth', 3, 'color', colors(2,:)); hold on;
    plot(1:maxIterationsPlot, mean_2shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
else
    errorbar(1:maxIterationsPlot, mean_2shared1, sem_2shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_2shared2, sem_2shared2,'LineWidth', 3, 'color', colors(2,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_2shared3, sem_2shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_2shared, sem_2shared, 'LineWidth', 3, 'color', colors(5,:)); 
end

xlabel('Training Iterations','FontSize', fontSize_xlabel-1);
ylabel({'MSE on', 'Single Targets Tasks'}, 'FontSize', fontSize_ylabel-1, 'FontName', fontName);
l = legend('Task A', 'Task B', 'Task C', 'Mean', 'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend-3);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca-1);

% 3 TASKS TRAINED

fig = figure(4);
set(fig, 'Position', [700 100 width height]);

if(~plotSEM)
    plot(1:maxIterationsPlot, mean_3shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    plot(1:maxIterationsPlot, mean_3shared2, 'LineWidth', 3, 'color', colors(2,:)); hold on;
    plot(1:maxIterationsPlot, mean_3shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
else
    errorbar(1:maxIterationsPlot, mean_3shared1, sem_3shared1, 'LineWidth', 3, 'color', colors(1,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_3shared2, sem_3shared2, 'LineWidth', 3, 'color', colors(2,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_3shared3, sem_3shared3, 'LineWidth', 3, 'color', colors(3,:)); hold on;
    errorbar(1:maxIterationsPlot, mean_3shared, sem_3shared, 'LineWidth', 3, 'color', colors(5,:)); 
end

xlabel('Training Iterations','FontSize', fontSize_xlabel-1);
ylabel({'MSE on', 'Single Targets Tasks'}, 'FontSize', fontSize_ylabel-1, 'FontName', fontName);
l = legend('Task A', 'Task B', 'Task C', 'Mean', 'Location', 'northeast');
set(l, 'FontName', fontName, 'fontSize', fontSize_legend-3);
set(gcf, 'Color', [1 1 1])
set(gca, 'FontSize', fontSize_gca-1);

%% MDS PLOTS
plotSettings;

width = 150;
height = 130;
fontSize_title = 11;
fontSize_xlabel = 12;
fontSize_ylabel = fontSize_xlabel;

repetitionIdx  = 2;
scalar_0 = 5;
scalar_1 = 5;
scalar_2 = scalar_1;
scalar_3 = scalar_1;

scalar_rest = 5;

markerSize = 50;
markerLineWidth = 2;
markerLineWidth_pretrained = 1;

fig = figure(1);
set(fig, 'Position', [100 100 width height]);

% plot single task training;
xlimit = [-1 1] * scalar_0;
ylimit = [-1 1] * scalar_0;

x = batch_log.MSE_hidden_0shared{repetitionIdx}(:,1);
y = batch_log.MSE_hidden_0shared{repetitionIdx}(:,2);

color = nan(numel(x,1),3);
color(1,:) = colors(cContrast1,:);
color(2,:) = colors(cContrast2,:);
color(3,:) = colors(cContrast3,:);

scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth);
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

% plot 1 shared
xlimit = [-1 1] * scalar_1;
ylimit = [-1 1] * scalar_1;

fig = figure(2);
set(fig, 'Position', [100 300 width height]);

x = batch_log.MSE_hidden_1shared1{repetitionIdx}(2:end,1);
y = batch_log.MSE_hidden_1shared1{repetitionIdx}(2:end,2);
x_pretrained = batch_log.MSE_hidden_1shared1{repetitionIdx}(1:1,1);
y_pretrained = batch_log.MSE_hidden_1shared1{repetitionIdx}(1:1,2);

color = zeros(numel(x,1),3);
color(1,:) = colors(cContrast1,:);
color(2,:) = colors(cContrast2,:);
color(3,:) = colors(cContrast3,:);
scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth); hold on;

color = zeros(numel(x_pretrained,1),3);
color(1,:) = colors(cContrast1,:);
scatter(x_pretrained, y_pretrained, markerSize, color, 'LineWidth', markerLineWidth_pretrained); hold off;
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

% plot 2 shared
xlimit = [-1 1] * scalar_2;
ylimit = [-1 1] * scalar_2;

fig = figure(3);
set(fig, 'Position', [100 500 width height]);

x = batch_log.MSE_hidden_1shared2{repetitionIdx}(3:end,1);
y = batch_log.MSE_hidden_1shared2{repetitionIdx}(3:end,2);
x_pretrained = batch_log.MSE_hidden_1shared2{repetitionIdx}(1:2,1);
y_pretrained = batch_log.MSE_hidden_1shared2{repetitionIdx}(1:2,2);

color = zeros(numel(x,1),3);
color(1,:) = colors(cContrast1,:);
color(2,:) = colors(cContrast2,:);
color(3,:) = colors(cContrast3,:);
scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth); hold on;

color = zeros(numel(x_pretrained,1),3);
color(1,:) = colors(cContrast1,:);
color(2,:) = colors(cContrast2,:);
scatter(x_pretrained, y_pretrained, markerSize, color, 'LineWidth', markerLineWidth_pretrained); hold off;
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

% plot 3 shared
xlimit = [-1 1] * scalar_3;
ylimit = [-1 1] * scalar_3;

fig = figure(4);
set(fig, 'Position', [100 700 width height]);

x = batch_log.MSE_hidden_1shared3{repetitionIdx}(4:end,1);
y = batch_log.MSE_hidden_1shared3{repetitionIdx}(4:end,2);
x_pretrained = batch_log.MSE_hidden_1shared3{repetitionIdx}(1:3,1);
y_pretrained = batch_log.MSE_hidden_1shared3{repetitionIdx}(1:3,2);

color = zeros(numel(x,1),3);
color(1,:) = colors(cContrast1,:);
color(2,:) = colors(cContrast2,:);
color(3,:) = colors(cContrast3,:);

scatter(x, y, markerSize, color, 'LineWidth', markerLineWidth); hold on;
scatter(x_pretrained, y_pretrained, markerSize, color, 'LineWidth', markerLineWidth_pretrained); hold off;
ylabel('Dimension 2', 'FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Dimension 1','FontSize', fontSize_xlabel, 'FontName', fontName);
           
xlim([xlimit(1) xlimit(2)]) 
ylim([ylimit(1) ylimit(2)]);

end