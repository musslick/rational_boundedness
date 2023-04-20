function fitVal = Manuscript_Simulation_2_TaskStrength_fun(inputArg)

inputArg = '1';
fitVal = 0;
inputArgDouble = str2double(inputArg);

%% initialization
% clc;
tic

log_version = 10;

% update meta simulation parameters
correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 
interferenceTaskStrength = [0:0.1:1.5];

% set up network parameters
numRepetitions = 10;
init_scale = 0.1;           % scales for initialized random weights
learningRate = 0.3;             % learning rate
decay = 0.0000;                 % weight penalization parameter
bias = -2;                          % weight from bias units to hidden & output units
iterations_Train = 50;         % maximum number of training iterations
thresh = 0.00;                  % mean-squared error stopping criterion

hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

% training environment parameters
NPathways = 2;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
NFeatures = 3;                      % number of feature units per stimulus input dimension
PathwayOverlap = 2;             % number of tasks per stimulus input dimension
numSglTasks = PathwayOverlap * NPathways;               % total number of tasks sampled
useCostomizedTasks = 1;     % set to 1 if using customized set of tasks (used for plotting purposes)
samplesPerTask_train = 100;   % stimuli per task for training (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
sdScale = 0;                        % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;     % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)

rng('default')
rng('shuffle');
 %% simulation loop
                      
rep = inputArgDouble;
                   
disp(strcat('repetition:',num2str(rep)));
numHiddenIdx = 1;

%% select subset of possible single tasks to train
interferenceTasks = [2 3];
primaryTasks = [1 4];

tasksToPerform = sort([interferenceTasks primaryTasks]);

samplesPerTask = [];
[inputSgl, tasksSgl, trainSgl, tasksIdxSgl, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);

for rep = 1:numRepetitions
    
    % generate task environment
    samplesPerTask = samplesPerTask_train;
    [input, tasks, train, ~, ~, ~, ~, ~, ~] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
  
    % set up model
    nHidden =  numHiddenUnits(numHiddenIdx);
    taskNet_org = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
    taskNet_org.NPathways = NPathways;

    % train on all single tasks
    taskNet_org.setData(input, tasks, train);
    taskNet_org.configure(); 

    for taskStrengthIdx = 1:length(interferenceTaskStrength)

        taskStrength = interferenceTaskStrength(taskStrengthIdx);

        % log task environment
        batch_log(taskStrengthIdx, rep).tasksToPerform(:) = tasksToPerform;
        
        % train neural network
        taskNet = NNmodel(taskNet_org);

        % if training set is complete then train on full set interleaved

        MSE_log = zeros(1, iterations_Train);
        for iter =1:iterations_Train

            % sample primary training patterns
            samplesPerTask = samplesPerTask_train;
            [input_primary, tasks_primary, train_primary] = createTaskPatterns(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, primaryTasks);

            % sub/over sample interference tasks
            samplesPerTask = round(samplesPerTask_train * taskStrength);
            [input_interference, tasks_interference, train_interference] = createTaskPatterns(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, interferenceTasks);

            % put it all together
            input = [input_primary; input_interference];
            tasks = [tasks_primary; tasks_interference];
            train = [train_primary; train_interference];

            % train network on samples
            taskNet.trainOnline(1, input, tasks, train);

            % is performance criterion reached?
            
            [~, ~, MSE] = taskNet.runSet(input_primary, tasks_primary, train_primary);
            MSE_log(iter) = mean(MSE);
            if (mean(MSE) <= taskNet.thresh)
                break
            end
            disp([iter mean(MSE)]);
        end

        % performance log
        batch_log(taskStrengthIdx, rep).MSE_log(:) = taskNet.MSE_log;
        batch_log(taskStrengthIdx, rep).CE_log(:) = taskNet.CE_log;
        batch_log(taskStrengthIdx, rep).CF_log(:) = taskNet.CF_log;
        batch_log(taskStrengthIdx, rep).DimCF_log(:) = taskNet.DimCF_log;

        disp('Completed single task training.');
        %% perform analysis for average task correlations

        % compute correlation
        [R_hidden, R_output] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform);

        % correlation matrix
        batch_log(taskStrengthIdx, rep).R_hidden(:,:) = R_hidden;       % stimulus wise correlations between tasks
        batch_log(taskStrengthIdx, rep).R_output(:,:) = R_output;         % stimulus wise correlations between tasks

        % validate extracted MIS

        for corrIdx = 1:length(correlation_thresholds)

            corr_threshold = correlation_thresholds(corrIdx);

            % COMPUTE DEPENDENCY GRAPH AND MIS
            NPathways = taskNet.NPathways;

            % find MIS
            [pathwayCapacities, maxCarryingCapacity,  BK_MIS, A_bipartite, A_tasksIdx, A_dual] = getMaxCarryingCapacity(R_hidden, R_output, corr_threshold);

            % formatting outputs
            BK_MIS = extend_BK_MIS(BK_MIS);
            % remove duplicates

            A_tasksIdx(A_tasksIdx > 0) = tasksToPerform(A_tasksIdx(A_tasksIdx > 0));
            A_tasksOrder = A_tasksIdx(:);
            A_tasksOrder(A_tasksOrder == 0) = [];
            pathwayCapacities = [A_tasksOrder pathwayCapacities];

            disp('Completed dependence graph computation.');

            disp(['Tested correlation threshold ' num2str(corrIdx) '/' num2str(length(correlation_thresholds))]);
        end

        %% COMPUTE MULTITASKING ACCURACY

        % FUNCTIONAL DEPENDENCE

        taskPair = [1 4];
        disp('TASK STRENGTH IDX');
        disp(taskStrengthIdx);
        thresholdIdx = 3;
        [CDFs, LCA_AB_Accuracy] = performTownsendAnalysis(taskNet, taskPair, inputSgl, tasksSgl, trainSgl, tasksIdxSgl, multiCap, thresholdIdx);

        batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy = LCA_AB_Accuracy;
        batch_log(taskStrengthIdx, rep).funcDependence.CDFs = CDFs;

    end
    
end

save(['logfiles/Part1/PsychReview_Part1_Sim2_TaskStrength_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_r' num2str(rep)]);


%% PREPARE DATA

numTaskStrengths = length(interferenceTaskStrength);

maxSize = 0;
for rep = 1:numRepetitions
    for taskStrengthIdx = 1:length(interferenceTaskStrength)
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.A_B)]);
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.A_B_1)]);
        maxSize = max([maxSize length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.min_A_B)]);
    end
end

numTimeSteps = maxSize;

LCA_funcDependence = nan(numRepetitions, numTaskStrengths);
CDFs_funcDependence.A_B = nan(numRepetitions, numTaskStrengths, numTimeSteps);
CDFs_funcDependence.A_B_1 = nan(numRepetitions, numTaskStrengths, numTimeSteps);
CDFs_funcDependence.min_A_B = nan(numRepetitions, numTaskStrengths, numTimeSteps);
CDFs_funcDependence.x = nan(numRepetitions, numTaskStrengths, numTimeSteps);
CDFs_funcDependence.C = nan(numRepetitions, numTaskStrengths, numTimeSteps);

for rep = 1:numRepetitions
    
   for taskStrengthIdx = 1:length(interferenceTaskStrength)
       LCA_funcDependence(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy;
       
       currentSize = length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.A_B);
       CDFs_funcDependence.A_B(rep, taskStrengthIdx, 1:currentSize) = batch_log(taskStrengthIdx, rep).funcDependence.CDFs.A_B;
       
       currentSize = length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.A_B_1);
       CDFs_funcDependence.A_B_1(rep, taskStrengthIdx, 1:currentSize) = batch_log(taskStrengthIdx, rep).funcDependence.CDFs.A_B_1;
       
       currentSize = length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.min_A_B);
       CDFs_funcDependence.min_A_B(rep, taskStrengthIdx, 1:currentSize) = batch_log(taskStrengthIdx, rep).funcDependence.CDFs.min_A_B;
       
       currentSize = length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.x);
       CDFs_funcDependence.x(rep, taskStrengthIdx, 1:currentSize) = batch_log(taskStrengthIdx, rep).funcDependence.CDFs.x;
       
       currentSize = length(batch_log(taskStrengthIdx, rep).funcDependence.CDFs.C);
       CDFs_funcDependence.C(rep, taskStrengthIdx, 1:currentSize) = batch_log(taskStrengthIdx, rep).funcDependence.CDFs.C;
   end
end

%% PLOT MULTITASKING ACCURACY

plotSettings;

plotBars = 0;

% set up figure
fig = figure(1);

if(plotBars)
    set(fig, 'Position', [100 100 250 250]);
else
    set(fig, 'Position', [100 100 500 200]);
end

% plot error bar data

grey = [0.8 0.8 0.8];

scale = 100;

% gather data
xLabels = {};
bardata_mean = nan(1, 0 + numTaskStrengths);
bardata_sem = nan(1, 0 + numTaskStrengths);
for i = 1:numTaskStrengths
    bardata_mean(0+i) = mean(LCA_funcDependence(:, i));
    bardata_sem(0+i) = std(LCA_funcDependence(:, i)) / sqrt(numRepetitions);
    xLabels{i} = num2str(round(interferenceTaskStrength(i)*100));
end
bardata_mean = bardata_mean * 100;
bardata_sem = bardata_sem * 100;

if(plotBars)
    % plot bars                      
    bar(bardata_mean, 'FaceColor', grey, 'EdgeColor', colors(cSingle,:)); hold on;
    % plot error bars
    errorbar(bardata_mean, bardata_sem, '.', 'LineWidth', lineWidth, 'Color', colors(cSingle,  :));
    hold off;
else
    errorbar(bardata_mean, bardata_sem, '-', 'LineWidth', lineWidth, 'Color', colors(cSingle,  :));
end

% set axes
ylim([0 1 * scale]);
if(plotBars)
    xlim([0.5 numTaskStrengths+0.5]);
end
if(plotBars)
    set(gca,'XTickLabel',{' '},'FontSize', fontSize_xlabel, 'FontName', fontName);
else
    set(gca, 'XTick', 1:length(interferenceTaskStrength));
    set(gca,'XTickLabel',xLabels,'FontSize', fontSize_xlabel-1, 'FontName', fontName);
end
ylabel('Multitasking Accuracy (%)','FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel({'% Training On Tasks C And D Compared To Tasks A And B'},'FontSize', fontSize_ylabel, 'FontName', fontName);

%% PLOT REACTION TIME INEQUALITY SERIIES (TOWNSEND & WENGER, 2004)

plotSEM = 1;
xlimit = [0.3 0.45];
ylimit = [-1 1];
plottedTaskStrengths = [0 0.5 1.5];

% plot

fig1 = figure(1);
set(fig1, 'Position', [100 100 900 200]);

% INDEPENDENCE
for i = 1:length(plottedTaskStrengths)

    relevantTaskStrengthIdx = find(interferenceTaskStrength == plottedTaskStrengths(i));

    x_data = squeeze(nanmean(CDFs_funcDependence.x(:,relevantTaskStrengthIdx,:)));
    y_funcDependence_A_B_mean = squeeze(nanmean(CDFs_funcDependence.A_B(:,relevantTaskStrengthIdx,:)));
    y_funcDependence_A_B_1_mean = squeeze(nanmean(CDFs_funcDependence.A_B_1(:,relevantTaskStrengthIdx,:)));
    y_funcDependence_min_A_B_mean = squeeze(nanmean(CDFs_funcDependence.min_A_B(:,relevantTaskStrengthIdx,:)));
    y_funcDependence_A_B_sem = squeeze(nanstd(CDFs_funcDependence.A_B(:,relevantTaskStrengthIdx,:)) / sqrt(numRepetitions));
    y_funcDependence_A_B_1_sem = squeeze(nanstd(CDFs_funcDependence.A_B_1(:,relevantTaskStrengthIdx,:)) / sqrt(numRepetitions));
    y_funcDependence_min_A_B_sem = squeeze(nanstd(CDFs_funcDependence.min_A_B(:,relevantTaskStrengthIdx,:)) / sqrt(numRepetitions));

    subplot(1, length(plottedTaskStrengths), i);

    if(plotSEM)
        errorbar(x_data, y_funcDependence_A_B_mean, y_funcDependence_A_B_sem, '-k', 'LineWidth',3); hold on;
        errorbar(x_data, y_funcDependence_A_B_1_mean, y_funcDependence_A_B_1_sem, '-k', 'LineWidth',1);
        errorbar(x_data, y_funcDependence_min_A_B_mean, y_funcDependence_min_A_B_sem, '--k', 'LineWidth',1); 
        hold off;
    else
        plot(x_data, y_funcDependence_A_B_mean, '-k', 'LineWidth',3); hold on;
        plot(x_data, y_funcDependence_A_B_1_mean, '-k', 'LineWidth',1);
        plot(x_data, y_funcDependence_min_A_B_mean, '--k', 'LineWidth',1); 
        hold off;    
    end

    ylabel({'Probability',  'Of Response Before t'}, 'fontSize', fontSize_ylabel-2);
    xlabel({'t (s)'}, 'fontSize', fontSize_xlabel-2);
%     title('Independence', 'FontSize', fontSize_title);
    % leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
    leg =  legend({'A AND B', 'A + B - 1', 'min(A, B)'}, 'Location', 'southeast');
    ylim(ylimit);
    xlim(xlimit);
    set(gca, 'FontSize', fontSize_gca);

end

%% FUNCTIONAL DEPENDENCE 1

relevantTaskStrengthIdx = find(interferenceTaskStrength == 1);

x_data = squeeze(nanmean(CDFs_funcDependence.x(:,relevantTaskStrengthIdx,:)));
y_funcDependence_A_B_mean = squeeze(nanmean(CDFs_funcDependence.A_B(:,relevantTaskStrengthIdx,:)));
y_funcDependence_A_B_1_mean = squeeze(nanmean(CDFs_funcDependence.A_B_1(:,relevantTaskStrengthIdx,:)));
y_funcDependence_min_A_B_mean = squeeze(nanmean(CDFs_funcDependence.min_A_B(:,relevantTaskStrengthIdx,:)));
y_funcDependence_A_B_sem = squeeze(nanstd(CDFs_funcDependence.A_B(:,relevantTaskStrengthIdx,:)) / sqrt(numRepetitions));
y_funcDependence_A_B_1_sem = squeeze(nanstd(CDFs_funcDependence.A_B_1(:,relevantTaskStrengthIdx,:)) / sqrt(numRepetitions));
y_funcDependence_min_A_B_sem = squeeze(nanstd(CDFs_funcDependence.min_A_B(:,relevantTaskStrengthIdx,:)) / sqrt(numRepetitions));

subplot(1, 3, 2);

if(plotSEM)
    errorbar(x_data, y_funcDependence_A_B_mean, y_funcDependence_A_B_sem, '-k', 'LineWidth',3); hold on;
    errorbar(x_data, y_funcDependence_A_B_1_mean, y_funcDependence_A_B_1_sem, '-k', 'LineWidth',1);
    errorbar(x_data, y_funcDependence_min_A_B_mean, y_funcDependence_min_A_B_sem, '--k', 'LineWidth',1); 
    hold off;
else
    plot(x_data, y_funcDependence_A_B_mean, '-k', 'LineWidth',3); hold on;
    plot(x_data, y_funcDependence_A_B_1_mean, '-k', 'LineWidth',1);
    plot(x_data, y_funcDependence_min_A_B_mean, '--k', 'LineWidth',1); 
    hold off;    
end

ylabel({'Probability',  'Of Response Before t'}, 'fontSize', fontSize_ylabel-2);
xlabel({'t (s)'}, 'fontSize', fontSize_xlabel-2);
title('Functional Dependence', 'FontSize', fontSize_title);
% leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
leg =  legend({'A AND B', 'A + B - 1', 'min(A, B)'}, 'Location', 'southeast');
ylim(ylimit);
xlim(xlimit);
set(gca, 'FontSize', fontSize_gca);

%% PLOT CAPACITY MEASURE BY TOWNSEND & WENGER (2004)


plotSEM = 0;
xlimit = [0.31 0.35];
ylimit = [0 1.1];

% plot

fig1 = figure(1);
set(fig1, 'Position', [100 100 400 200]);

colors = nan(numTaskStrengths+1, 3);
colors(:,1) = linspace(0.8, 0, numTaskStrengths+1);
colors(:,2) = linspace(0.8, 0, numTaskStrengths+1);
colors(:,3) = linspace(0.8, 0, numTaskStrengths+1);

xLabels = {};
for taskStrengthIdx = 1:numTaskStrengths
    strength = interferenceTaskStrength(taskStrengthIdx) * 100;
    xLabels{taskStrengthIdx} = ['Func. Dependence, ' num2str(strength) '% Task Strength'];
    
    x_data = squeeze(nanmean(CDFs_funcDependence.x(:,taskStrengthIdx,:)));
    y_funcDependence_C_mean = squeeze(nanmean(CDFs_funcDependence.C(:,taskStrengthIdx,:)));
    y_funcDependence_C_sem = squeeze(nanstd(CDFs_funcDependence.C(:,taskStrengthIdx,:)) / sqrt(numRepetitions)); 
    
    if(plotSEM)
        errorbar(x_data, y_funcDependence_C_mean, y_funcDependence_C_sem, '-', 'Color', colors(taskStrengthIdx*1,:), 'LineWidth',3); hold on;
    else
        plot(x_data, y_funcDependence_C_mean, '-', 'Color', colors(taskStrengthIdx+1,:), 'LineWidth',3); hold on     
    end
end
hold off;

ylabel({'Capacity C(t)'}, 'fontSize', fontSize_ylabel-2);
xlabel({'t (s)'}, 'fontSize', fontSize_xlabel-2);
title('Capacity based on Townsend & Wenger (2004)', 'FontSize', fontSize_title);
% leg =  legend({'$P_{AB}(T_A \leq t, T_B \leq t)$', '$P_{A}(T_A \leq t) + P_{B}(T_B \leq t) - 1$', '$min[P_{A}(T_A \leq t), P_{B}(T_B \leq t)]$'},'Interpreter','latex', 'Location', 'southoutside');
leg =  legend(xLabels, 'Location', 'northeast');
ylim(ylimit);
xlim(xlimit);
set(gca, 'FontSize', fontSize_gca);



