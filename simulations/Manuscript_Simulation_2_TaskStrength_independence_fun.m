function fitVal = Manuscript_Simulation_2_TaskStrength_independence_fun(inputArg)

inputArg = '1';
fitVal = 0;
inputArgDouble = str2double(inputArg);

init_taskCorr = 0;
fixOutputWeights = 1;
fixHiddenWeights = 1;
init_task_scale = 5;
outputWeightStrength = 2;
iterations_Train = 500;         % maximum number of training iterations
thresh = 0.001;                  % mean-squared error stopping criterion

NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)

if(NPathways == 3)
    interferenceTasks = [2 3 4 6 7 8];
    primaryTasks = [1 5 9];
    taskPair = [1 5]; % primaryTasks
elseif(NPathways == 2)
    interferenceTasks = [2 3];
    primaryTasks = [1 4];
    taskPair = [1 4]; % primaryTasks
end


%% initialization
% clc;
tic

log_version = 10;

% update meta simulation parameters
correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 
interferenceTaskStrength = [0 1]; %[0:0.1:1.5];

% set up network parameters
numRepetitions = 2;
init_scale = 0.1;           % scales for initialized random weights
learningRate = 0.3;             % learning rate
decay = 0.0000;                 % weight penalization parameter
bias = -2;                          % weight from bias units to hidden & output units

hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

% training environment parameters
NFeatures = 3;                      % number of feature units per stimulus input dimension
PathwayOverlap = 2;             % number of tasks per stimulus input dimension
numSglTasks = PathwayOverlap * NPathways;               % total number of tasks sampled
useCostomizedTasks = 1;     % set to 1 if using customized set of tasks (used for plotting purposes)
samplesPerTask_train = 100;   % stimuli per task for training (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
sdScale = 0;                        % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;     % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)
orthogonalizeReferenceVectors = 1;      % flag: set to 1 if reference vectors should be entirely orthogonal


rng('default')
rng('shuffle');
 %% simulation loop
                      
rep = inputArgDouble;
                   
disp(strcat('repetition:',num2str(rep)));
numHiddenIdx = 1;

%% select subset of possible single tasks to train

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

    %% set up output weights
    outputWeights = nan(size(taskNet_org.weights.W_TO));
    if(fixOutputWeights)
        if(NPathways == 2 & NFeatures == 3)
            vec = transpose([1 1 1 0 0 0]) * outputWeightStrength;
            outputWeights(:,1) = vec;
            outputWeights(:,3) = vec;
            outputWeights(:,2) = -vec;
            outputWeights(:,4) = -vec;
            taskNet_org.weights.W_TO = outputWeights;
        elseif(NPathways == 3 & NFeatures == 3)
%             vec1 = [1 1 1 -1 -1 -1 -1 -1 -1] * outPutWeightStrength;
%             vec2 = [-1 -1 -1 1 1 1 -1 -1 -1] * outPutWeightStrength;
%             vec3 = [-1 -1 -1 -1 -1 -1 1 1 1] * outPutWeightStrength;
            
            vec1 = [1 1 1 0 0 0 0 0 0] * outputWeightStrength;
            vec2 = [0 0 0 1 1 1 0 0 0] * outputWeightStrength;
            vec3 = [0 0 0 0 0 0 1 1 1] * outputWeightStrength;

            outputWeights(:,1) = vec1;
            outputWeights(:,2) = vec2;
            outputWeights(:,3) = vec3;
            outputWeights(:,4) = vec1;
            outputWeights(:,5) = vec2;
            outputWeights(:,6) = vec3;
            outputWeights(:,7) = vec1;
            outputWeights(:,8) = vec2;
            outputWeights(:,9) = vec3;
        else
            error('Could not fix output weights.');
        end
    end
    
    %% set up task weights according to max. cosine similarity

    num_points = NPathways;
    dim = nHidden(end);
    weightCorr = nan(1,NPathways);

    hiddenGroups = floor(dim/num_points);
    taskWeights = nan(dim, numSglTasks);

    for i = 1:num_points

        W = nan(dim, num_points);

        % create random initial vector
        x_ref_mask = zeros(dim, 1);
        x_ref_mask(((i-1)*hiddenGroups+1):(i*hiddenGroups)) = 1;
        if(~orthogonalizeReferenceVectors)
            x_ref = randn(dim, 1)*0.1 + x_ref_mask;
        else
            x_ref = x_ref_mask;
        end
        x_ref = x_ref./norm(x_ref);

        for p = 1:num_points

            r = rand*(1-init_taskCorr) + init_taskCorr;
            alpha = acos(r);

            % 1) sample uniformly in S^(n-2) in R^(n-1)
            l = sin(alpha);
            y = rand(dim-1, 1)-0.5;
            y = y./norm(y) * l;

            % 2) assume a ~= e_i for i = 1:dim; use Gram-Schmidt to get orthonormal
            % basis
            I = eye(dim);
            [Q,R] = qr([x_ref I(:,2:end)]);

            S = Q(:,2:end);

            % define T: R^(n-1) -> R^n

            Ty = S*y + r*x_ref;

            W(:,p) = Ty;

        end

        weightIdx = find(tril(ones(size(corr(W))),-1));
        W_corr = corr(W);
        weightCorr(i) = mean(W_corr(weightIdx));

        taskWeights(:,((i-1)*num_points+1):(i*num_points)) = W;

    end

    taskWeights = taskWeights .*init_task_scale(end);

    % set up network weights
    taskNet_org.weights.W_TH = taskWeights;

    % log initial weight correlation
    weightCorr_log(rep) = mean(weightCorr);
    init_taskCorr_log(rep) = init_taskCorr;

    %%
    for taskStrengthIdx = 1:length(interferenceTaskStrength)

        taskStrength = interferenceTaskStrength(taskStrengthIdx);

        % log task environment
        batch_log(taskStrengthIdx, rep).tasksToPerform(:) = tasksToPerform;
        
        % train neural network
        taskNet = NNmodel(taskNet_org);
        taskNet.setFixedWeights([NNmodel.W_TASK_HIDDEN NNmodel.W_TASK_OUTPUT]);

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
                disp('Learned primary tasks');
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

        disp('TASK STRENGTH IDX');
        disp(taskStrengthIdx);
        thresholdIdx = [];

        % identify task patterns
        taskProfile = zeros(1, size(tasksSgl, 2));
        taskProfile(taskPair) = 1;
        cap = length(taskPair);
        patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));

        input_LCA = multiCap{cap}.input(patternIdx, :);
        tasks_LCA = multiCap{cap}.tasks(patternIdx, :);
        train_LCA = multiCap{cap}.train(patternIdx, :);

        % LCA settings
        loadLCASettings;
        LCA_settings.responseThreshold = 0:0.05:2;
        LCA_settings.ITI = 1;

        % Multi LCA call
        [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
        = taskNet.runLCA(LCA_settings, input_LCA, tasks_LCA, train_LCA);
        LCA_Accuracy = nanmean(optAccuracy);

        batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy = LCA_Accuracy;

    end
    
end

save(['logfiles/Part1/PsychReview_Part1_Sim2_TaskStrength_independence_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_r' num2str(rep)]);


%% PREPARE DATA

numTaskStrengths = length(interferenceTaskStrength);

LCA_funcDependence = nan(numRepetitions, numTaskStrengths);

for rep = 1:numRepetitions
    
   for taskStrengthIdx = 1:length(interferenceTaskStrength)
       LCA_funcDependence(rep, taskStrengthIdx) = batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy;
       
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
