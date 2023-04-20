function fitVal = Manuscript_Simulation_ModifiedStroop(inputArg)

clc;
clear all;
inputArg = '1';

fitVal = 0;
inputArgDouble = str2double(inputArg);

initWeights = 0;
taskCorrs = [0 0.25 0.5 0.75 1];
init_taskCorr = taskCorrs(inputArgDouble);
fullOrthogonalization = 0;
trainMultitasking = 0;
fixOutputWeights = 0;
fixHiddenWeights = 0;
init_task_scale = 5;
iterations_Train = 500;         % maximum number of training iterations
thresh = 0.01;                  % mean-squared error stopping criterion

NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)

interferenceTasks = [4];
primaryTasks = [1 5 8];
taskPair_dependent = [1 5]; % primaryTasks
taskPair_independent = [1 8]; % primaryTasks

%% description
% 
% Simulation # 1
% 
% The simulation involves training a neural networks on a set of tasks and
% then predicting the multitasking performance of the network based on an
% interference graph that is extracted from the networks learned
% single-task representations
%
% author: Sebastian Musslick

%% initialization
% clc;
tic

log_version = 16;

% update meta simulation parameters
correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 
interferenceTaskStrength = [2]; %[0:0.1:1.5];

% set up network parameters
numRepetitions = 27;
init_scale = 0.1;           % scales for initialized random weights
learningRate = 0.3;             % learning rate
decay = 0.0000;                 % weight penalization parameter
bias = -2;                          % weight from bias units to hidden & output units

hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

% training environment parameters
NFeatures = 4;                      % number of feature units per stimulus input dimension
PathwayOverlap = 2;             % number of tasks per stimulus input dimension
numSglTasks = PathwayOverlap * NPathways;               % total number of tasks sampled
useCostomizedTasks = 1;     % set to 1 if using customized set of tasks (used for plotting purposes)
samplesPerTask_train = 100;   % stimuli per task for training (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
sdScale = 0;                        % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;     % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)
orthogonalizeReferenceVectors = 1;      % flag: set to 1 if reference vectors should be entirely orthogonal
MDSMetric = 'euclidean'; % 'euclidean', 'correlation', 'cosine'

rng('default')
rng('shuffle');
 
% simulation loop
                      
rep = inputArgDouble;
                   
disp(strcat('repetition:',num2str(rep)));
numHiddenIdx = 1;

% select subset of possible single tasks to train

tasksToPerform = sort([interferenceTasks primaryTasks]);

samplesPerTask = [];
[inputSgl, tasksSgl, trainSgl, tasksIdxSgl, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);

for rep = 1:numRepetitions
    
    % generate task environment
    samplesPerTask = samplesPerTask_train;
    [input, tasks, train, ~, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
  
    % set up model
    nHidden =  numHiddenUnits(numHiddenIdx);
    taskNet_org = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
    taskNet_org.NPathways = NPathways;

    % train on all single tasks
    taskNet_org.setData(input, tasks, train);
    taskNet_org.configure(); 

    %% set up output weights
    if(initWeights || fullOrthogonalization)
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
    orthogonalGroups = floor(dim/(NPathways*NPathways));
    taskWeights = nan(dim, NPathways*NPathways);

    if(fullOrthogonalization)
        taskWeights = zeros(dim, NPathways*NPathways);
        for i = 1:(NPathways*NPathways)
            taskWeights((i-1)*orthogonalGroups+[1:orthogonalGroups],i) = randn(orthogonalGroups,1);
        end
    else
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
    end
    taskWeights = taskWeights .*init_task_scale(end);

    % set up network weights
    taskNet_org.weights.W_TH = taskWeights;

    % log initial weight correlation
    weightCorr_log(rep) = mean(weightCorr);
    init_taskCorr_log(rep) = init_taskCorr;
    else
        taskNet_org.weights.W_TH  = taskNet_org.weights.W_TH; %* 2 * init_task_scale;
    end
    %%

    for taskStrengthIdx = 1:length(interferenceTaskStrength)

        taskStrength = interferenceTaskStrength(taskStrengthIdx);

        % log task environment
        batch_log(taskStrengthIdx, rep).tasksToPerform(:) = tasksToPerform;
        
        % train neural network
        taskNet = NNmodel(taskNet_org);
        if(fixOutputWeights)
            taskNet.setFixedWeights([NNmodel.W_TASK_OUTPUT]);
        end
        if(fixHiddenWeights)
            taskNet.setFixedWeights([NNmodel.W_TASK_HIDDEN]);
        end

        % if training set is complete then train on full set interleaved

        MSE_log = zeros(1, iterations_Train);
        for iter =1:iterations_Train

            % sample primary training patterns
            samplesPerTask = samplesPerTask_train;
            [input_primary, tasks_primary, train_primary] = createTaskPatterns(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, primaryTasks);

            % train on multitasking set
            if(trainMultitasking)
%                 % independent tasks
%                 taskProfile = zeros(1, size(tasksSgl, 2));
%                 taskProfile(taskPair_independent) = 1;
%                 cap = length(taskPair_independent);
%                 patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));
%                 patternIdx = patternIdx(randperm(length(patternIdx)));
%                 samplesPerTask = 1:samplesPerTask_train;
%                 sampledPatternIdx = patternIdx(samplesPerTask);
%                 
%                 input_primary = [input_primary; multiCap{2}.input(sampledPatternIdx,:)];
%                 tasks_primary = [tasks_primary; multiCap{2}.tasks(sampledPatternIdx,:)];
%                 train_primary = [train_primary; multiCap{2}.train(sampledPatternIdx,:)];
%                 
                % functionally dependent tasks
                taskProfile = zeros(1, size(tasksSgl, 2));
                taskProfile(taskPair_dependent) = 1;
                cap = length(taskPair_dependent);
                patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));
                patternIdx = patternIdx(randperm(length(patternIdx)));
                samplesPerTask = 1:samplesPerTask_train;
                sampledPatternIdx = patternIdx(samplesPerTask);
                
                input_primary = [input_primary; multiCap{2}.input(sampledPatternIdx,:)];
                tasks_primary = [tasks_primary; multiCap{2}.tasks(sampledPatternIdx,:)];
                train_primary = [train_primary; multiCap{2}.train(sampledPatternIdx,:)];
            end
            
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
            
            if(~trainMultitasking)
                [~, ~, MSE] = taskNet.runSet(input_primary, tasks_primary, train_primary);
            else
                [~, ~, MSE] = taskNet.runSet(input, tasks, train);
            end
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
        [R_hidden, R_output, hiddenData] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform);

        % correlation matrix
        batch_log(taskStrengthIdx, rep).R_hidden(:,:) = R_hidden;       % stimulus wise correlations between tasks
        batch_log(taskStrengthIdx, rep).R_output(:,:) = R_output;         % stimulus wise correlations between tasks

        % perform MDS analysis
        nTasks = length(tasksToPerform);
        X = transpose(hiddenData);
        dist = MDSMetric;
        distances = pdist(X, dist);
        try
        Y = mdscale(distances,2);
        catch
           disp(inputSgl);
           Y = zeros(size(inputSgl,1), 2);
           warning('No MDS possible for average single task reps'); 
        end
        batch_log(taskStrengthIdx, rep).MDS_data = Y(1:nTasks, :);
        
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

        %%
        relevantTasks = [primaryTasks interferenceTasks];
        hiddenCorr = corr(taskNet.weights.W_TH(:,relevantTasks));
        outputCorr = corr(taskNet.weights.W_TO(:,relevantTasks));
        batch_log(taskStrengthIdx, rep).hiddenCorr = hiddenCorr;
        batch_log(taskStrengthIdx, rep).outputCorr = outputCorr;

        %% COMPUTE MULTITASKING ACCURACY

        % FUNCTIONAL DEPENDENCE

        disp('TASK STRENGTH IDX');
        disp(taskStrengthIdx);

        % identify task patterns
        taskProfile = zeros(1, size(tasksSgl, 2));
        taskProfile(taskPair_dependent) = 1;
        cap = length(taskPair_dependent);
        patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));

        input_LCA = multiCap{cap}.input(patternIdx, :);
        tasks_LCA = multiCap{cap}.tasks(patternIdx, :);
        train_LCA = multiCap{cap}.train(patternIdx, :);

        % identify congruent & incongruent stimuli for the two DEPENDENT tasks
        congruentIdx = [];
        incongruentIdx = [];
        inputDimTaskA  = ceil(taskPair_dependent(1) / taskNet.NPathways);
        inputDimTaskB = ceil(taskPair_dependent(2) / taskNet.NPathways);
        for i = 1:size(input_LCA,1)
            A = reshape(input_LCA(i,:), NFeatures, NPathways);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end

        % LCA settings
        loadLCASettings;

        % Multi LCA call
        [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
        = taskNet.runLCA(LCA_settings, input_LCA, tasks_LCA, train_LCA);
        LCA_Accuracy = nanmean(optAccuracy);
        LCA_Accuracy_con = nanmean(optAccuracy(congruentIdx));
        LCA_Accuracy_inc = nanmean(optAccuracy(incongruentIdx));

        batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy = LCA_Accuracy;
        batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy_inc = LCA_Accuracy_inc;
        batch_log(taskStrengthIdx, rep).funcDependence.LCA_Accuracy_con = LCA_Accuracy_con;

        % INDEPENDENCE

        % identify task patterns
        taskProfile = zeros(1, size(tasksSgl, 2));
        taskProfile(taskPair_independent) = 1;
        cap = length(taskPair_independent);
        patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));

        input_LCA = multiCap{cap}.input(patternIdx, :);
        tasks_LCA = multiCap{cap}.tasks(patternIdx, :);
        train_LCA = multiCap{cap}.train(patternIdx, :);

        % identify congruent & incongruent stimuli for the two DEPENDENT tasks
        congruentIdx = [];
        incongruentIdx = [];
        inputDimTaskA  = ceil(taskPair_dependent(1) / taskNet.NPathways);
        inputDimTaskB = ceil(taskPair_dependent(2) / taskNet.NPathways);
        for i = 1:size(input_LCA,1)
            A = reshape(input_LCA(i,:), NFeatures, NPathways);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end

        % Multi LCA call
        [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
        = taskNet.runLCA(LCA_settings, input_LCA, tasks_LCA, train_LCA);
        LCA_Accuracy = nanmean(optAccuracy);
        LCA_Accuracy_con = nanmean(optAccuracy(congruentIdx));
        LCA_Accuracy_inc = nanmean(optAccuracy(incongruentIdx));

        batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy = LCA_Accuracy;
        batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy_inc = LCA_Accuracy_inc;
        batch_log(taskStrengthIdx, rep).Independence.LCA_Accuracy_con = LCA_Accuracy_con;

        % TASK A
        % identify task patterns
        taskProfile = zeros(1, size(tasksSgl, 2));
        taskProfile(taskPair_dependent(1)) = 1;
        cap = 1;
        patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));

        input_LCA = multiCap{cap}.input(patternIdx, :);
        tasks_LCA = multiCap{cap}.tasks(patternIdx, :);
        train_LCA = multiCap{cap}.train(patternIdx, :);
        
        % identify congruent & incongruent stimuli for the two DEPENDENT tasks
        congruentIdx = [];
        incongruentIdx = [];
        inputDimTaskA  = ceil(taskPair_dependent(1) / taskNet.NPathways);
        inputDimTaskB = ceil(taskPair_dependent(2) / taskNet.NPathways);
        for i = 1:size(input_LCA,1)
            A = reshape(input_LCA(i,:), NFeatures, NPathways);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end

        % Multi LCA call
        [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
        = taskNet.runLCA(LCA_settings, input_LCA, tasks_LCA, train_LCA);
        LCA_Accuracy = nanmean(optAccuracy);
        LCA_Accuracy_con = nanmean(optAccuracy(congruentIdx));
        LCA_Accuracy_inc = nanmean(optAccuracy(incongruentIdx));

        batch_log(taskStrengthIdx, rep).taskA.LCA_Accuracy = LCA_Accuracy;
        batch_log(taskStrengthIdx, rep).taskA.LCA_Accuracy_inc = LCA_Accuracy_inc;
        batch_log(taskStrengthIdx, rep).taskA.LCA_Accuracy_con = LCA_Accuracy_con;

         % TASK B
        % identify task patterns
        taskProfile = zeros(1, size(tasksSgl, 2));
        taskProfile(taskPair_dependent(1)) = 1;
        cap = 1;
        patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));

        input_LCA = multiCap{cap}.input(patternIdx, :);
        tasks_LCA = multiCap{cap}.tasks(patternIdx, :);
        train_LCA = multiCap{cap}.train(patternIdx, :);
        
        % identify congruent & incongruent stimuli for the two DEPENDENT tasks
        congruentIdx = [];
        incongruentIdx = [];
        inputDimTaskA  = ceil(taskPair_dependent(1) / taskNet.NPathways);
        inputDimTaskB = ceil(taskPair_dependent(2) / taskNet.NPathways);
        for i = 1:size(input_LCA,1)
            A = reshape(input_LCA(i,:), NFeatures, NPathways);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end

        % Multi LCA call
        [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
        = taskNet.runLCA(LCA_settings, input_LCA, tasks_LCA, train_LCA);
        LCA_Accuracy = nanmean(optAccuracy);
        LCA_Accuracy_con = nanmean(optAccuracy(congruentIdx));
        LCA_Accuracy_inc = nanmean(optAccuracy(incongruentIdx));

        batch_log(taskStrengthIdx, rep).taskB.LCA_Accuracy = LCA_Accuracy;
        batch_log(taskStrengthIdx, rep).taskB.LCA_Accuracy_inc = LCA_Accuracy_inc;
        batch_log(taskStrengthIdx, rep).taskB.LCA_Accuracy_con = LCA_Accuracy_con;


        % TASK C
        % identify task patterns
        taskProfile = zeros(1, size(tasksSgl, 2));
        taskProfile(taskPair_independent(2)) = 1;
        cap = 1;
        patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));

        input_LCA = multiCap{cap}.input(patternIdx, :);
        tasks_LCA = multiCap{cap}.tasks(patternIdx, :);
        train_LCA = multiCap{cap}.train(patternIdx, :);
        
        % identify congruent & incongruent stimuli for the two DEPENDENT tasks
        congruentIdx = [];
        incongruentIdx = [];
        inputDimTaskA  = ceil(taskPair_dependent(1) / taskNet.NPathways);
        inputDimTaskB = ceil(taskPair_dependent(2) / taskNet.NPathways);
        for i = 1:size(input_LCA,1)
            A = reshape(input_LCA(i,:), NFeatures, NPathways);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end

        % Multi LCA call
        [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
        = taskNet.runLCA(LCA_settings, input_LCA, tasks_LCA, train_LCA);
        LCA_Accuracy = nanmean(optAccuracy);
        LCA_Accuracy_con = nanmean(optAccuracy(congruentIdx));
        LCA_Accuracy_inc = nanmean(optAccuracy(incongruentIdx));

        batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy = LCA_Accuracy;
        batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy_inc = LCA_Accuracy_inc;
        batch_log(taskStrengthIdx, rep).taskC.LCA_Accuracy_con = LCA_Accuracy_con;

     %% TOWNSEND ANLAYIS

    % INDEPENDENCE

    % identify task patterns
    taskProfile = zeros(1, size(tasksSgl, 2));
    taskProfile(taskPair_independent) = 1;
    cap = length(taskPair_independent);
    patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));
    type0Idx = multiCap{cap}.taskIdx(patternIdx);
    assert(length(unique(type0Idx)) == 1);
    type0Idx = unique(type0Idx);

    [CDFs] = getTownsendCDF(taskNet, NPathways, type0Idx, multiCap);
    batch_log(taskStrengthIdx, rep).Independence.A_B = CDFs{1}.A_B; 
    batch_log(taskStrengthIdx, rep).Independence.A_B_1 = CDFs{1}.A_B_1; 
    batch_log(taskStrengthIdx, rep).Independence.min_A_B = CDFs{1}.min_A_B; 
    batch_log(taskStrengthIdx, rep).Independence.C = CDFs{1}.C; 

    % INDEPENDENCE

    % identify task patterns
    taskProfile = zeros(1, size(tasksSgl, 2));
    taskProfile(taskPair_dependent) = 1;
    cap = length(taskPair_dependent);
    patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));
    type2Idx = multiCap{cap}.taskIdx(patternIdx);
    assert(length(unique(type2Idx)) == 1);
    type2Idx = unique(type2Idx);

    [CDFs] = getTownsendCDF(taskNet, NPathways, type2Idx, multiCap);
    batch_log(taskStrengthIdx, rep).funcDependence.A_B = CDFs{1}.A_B; 
    batch_log(taskStrengthIdx, rep).funcDependence.A_B_1 = CDFs{1}.A_B_1; 
    batch_log(taskStrengthIdx, rep).funcDependence.min_A_B = CDFs{1}.min_A_B; 
    batch_log(taskStrengthIdx, rep).funcDependence.C = CDFs{1}.C; 

    end
    
end

if(~trainMultitasking)
    save(['logfiles/Part2/PsychReview_Part2_Exp1_modifiedStroop_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_ortho' num2str(fullOrthogonalization)]);
else
    save(['logfiles/Part2/PsychReview_Part2_Exp1_modifiedStroop_multi_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_ortho' num2str(fullOrthogonalization)]);
end

disp('DONE');
%% PREPARE DATA

numTaskStrengths = length(interferenceTaskStrength);

LCA_funcDependence = nan(numRepetitions, numTaskStrengths);

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
       
   end
end



%% plot corrleation matrix
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
set(fig2, 'Position', [100 100 230 230]);
imagesc(M);

map = colormap(copper);
map(1:end, 1) = linspace(0, 1, size(map,1));
map(1:end, 2) = linspace(0, 1, size(map,1));
map(1:end, 3) = linspace(0, 1, size(map,1));
map = map(fliplr(1:size(map,1)),:);
colormap(map);
% colorbar;
caxis([0 1]);
set(gca, 'FontSize', 12);
set(gca, 'XTick', 1:size(M,1));
set(gca, 'YTick', 1:size(M,2));
set(gca, 'XTickLabel', {'A', 'B', 'C', 'D', 'E'});
set(gca, 'YTickLabel', {'A', 'B', 'C', 'D', 'E'});
ylabel('Tasks','FontSize', fontSize_ylabel, 'FontName', fontName);
xlabel('Tasks','FontSize', fontSize_xlabel, 'FontName', fontName);
title({'Correlation of Task Representations' ' at Associative Layer'},'FontSize', fontSize_xlabel-2, 'FontName', fontName);

end
