%% train net and see how MSE looks like

% TODO: do only incongruent for multitaking
% identity matrix from input to hidden
% delete inhibitory weights for co-executed task

clear all;

% meta simulation parameters
log_version = 1;
replications = 100;                   % # replications per simulation

% set up network parameters
nHidden = 200;              % number of hidden units
hiddenPathSize = 1;         % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;         % group size of output units that receive the same weights from the task layer)
learningRate = 0.3;         % learning rate
thresh = 0.0001;            % mean-squared error stopping criterion
decay = 0.0000;             % weight penalization parameter
bias = -2;                  % weight from bias units to hidden & output units
init_scale = 0.1;           % scales for initialized random weights
iterations_Train = 500;     % number of training iterations
corr_thresholds = 0.1:0.05:0.9;      % threshold for detecting merged task representations in correlation matrix

% task environment parameters 
NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
Nactive = NPathways;               % how many tasks to perform simultaneously
NFeatures = 5;                     % the number of features per feature dimension
currFeaturesDims = 1:NPathways;    % lists all feature dimensions
nTasks = NPathways^2;
samplesPerTask = 400;
samplesPerTask_gen = 400;

% sameClassifierAcrossTasks = 1;
sameStimuliAcrossTasks = 1;
sdScale = 0.4;
samplesPerFeatureCombination = 1;

taskSimilarity = 1;
init_taskCorrelations = [1];
init_task_scale = 5;


switch_tau = 0.1;                   % rate constant for net input integration 
switch_iterations = 100;            % maximum switch time
switch_ddmp.nSimulations = 1;     % how many accumulator simulations
switch_ddmp.d = 0.01;                % acumulation drift rate
switch_ddmp.z = 0.1:0.1:10;          % thresholds
switch_ddmp.c = 0.0;                % noise
switch_ddmp.respOnset = 2;          % response onset (ignore all time points before until this time point)

init_scales = 0.1;%[0.1:0.2:1.0];        % magnitude of initial weights

% numCores = 8;               % specify number of available cores
% maxNumCompThreads(numCores);


%batch_log = repmat(struct('hiddenTaskRep',zeros(length(init_taskCorrelations), replications, nTasks, nTasks), ...
%                          'outputTaskRep',zeros(length(init_taskCorrelations), replications, nTasks, nTasks)),1, 1); % zeros(replications,nTasks,nTasks)

for rep = 1:replications
tic
    [inputSgl_cut tasksSgl_cut trainSgl_cut tasksIdxSgl_cut inputSgl_mask, tasksSgl_mask, trainSgl_mask, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerFeatureCombination, sdScale, sameStimuliAcrossTasks);
           
    taskNet1 = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
    taskNet1.setData(inputSgl_cut, tasksSgl_cut, trainSgl_cut);
    %taskNet1.init_task_scale = init_task_scale;
    taskNet1.fixedWeights.W_TH = [0];
    taskNet1.configure(); 

    %% set up task weights

    num_points = NPathways;
    dim = nHidden(end);

    hiddenGroups = floor(dim/num_points);
    taskWeights = nan(dim, nTasks);

    for i = 1:num_points

        W = nan(dim, num_points);

        % create random initial vector
        x_ref_mask = zeros(dim, 1);
        x_ref_mask(((i-1)*hiddenGroups+1):(i*hiddenGroups)) = 1;
        x_ref = randn(dim, 1)*0.1 + x_ref_mask;
        x_ref = x_ref./norm(x_ref);

        for p = 1:num_points

            init_taskCorr = -1 + rand*2;
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

        taskWeights(:,((i-1)*num_points+1):(i*num_points)) = W;

    end

    taskWeights = taskWeights .*init_task_scale(end);

    %%
    taskNet1.weights.W_TH = taskWeights;       

    taskNet1.trainOnline(iterations_Train,[],[],[],1);

    % get task representations
    taskIDs = unique(tasksIdxSgl_cut);

    [outData, hiddenData, MSE] = taskNet1.runSet(inputSgl_cut, tasksSgl_cut, trainSgl_cut);
    
    batch_log{rep}.activationData.outputLayer = outData;
    batch_log{rep}.activationData.hiddenLayer = hiddenData;
    batch_log{rep}.activationData.taskIndices = tasksIdxSgl_cut;
    
    hiddenTaskReps = zeros(nTasks, size(hiddenData,2));
    taskMSE = zeros(nTasks,1);
    for i = 1:length(taskIDs)
        task = taskIDs(i);
        hiddenTaskReps(i,:) = mean(hiddenData(tasksIdxSgl_cut == task,:),1);
        taskMSE(i) = mean(MSE(tasksIdxSgl_cut == task));
    end

    outputTaskReps = zeros(nTasks, size(outData,2));
    for i = 1:length(taskIDs)
        task = taskIDs(i);
        outputTaskReps(i,:) = mean(outData(tasksIdxSgl_cut == task,:),1);
    end

    R_hidden = corr(hiddenTaskReps');
    batch_log{rep}.hiddenTaskRep(:,:) = R_hidden;
    R_output = corr(outputTaskReps');
    batch_log{rep}.outputTaskRep(:,:) = R_output;
    
    batch_log{rep}.performance{1}.numberOfTasks = 1;
    batch_log{rep}.performance{1}.MSE = taskMSE;
    batch_log{rep}.performance{1}.taskIndicies = taskIDs;
    batch_log{rep}.performance{1}.nTaskCombinations = length(taskIDs);

    % test multitasking for each multitasking condition
    for cap = 2:length(multiCap) 
       
        batch_log{rep}.performance{cap}.numberOfTasks = cap;
        
        taskConditions = unique(multiCap{cap}.tasks,'rows');
        taskIndicies = nan(size(taskConditions,1), cap);
        multiMSE = nan(size(taskConditions,1), 1);
        
        [~, ~, multi_MSE_log] = taskNet1.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
        
        % determine congruency
        congruency = nan(1, size(multiCap{cap}.train,1));
        for row = 1:length(congruency)
            congrSum = sum(reshape(multiCap{cap}.train(row,:), NFeatures, NPathways),2);
            if(max(congrSum) <= 1)
                congruency(row) = 0;
            else
                congruency(row) = 1;
            end
        end
        
        for taskComb = 1:size(taskConditions,1)
            taskIndicies(taskComb,:) = find(taskConditions(taskComb,:) == 1);
            
            taskCond = taskConditions(taskComb,:);
            taskCondIdx = find(transpose(ismember(multiCap{cap}.tasks, taskCond, 'rows')) & congruency == 0);
            
            multiMSE(taskComb) = mean(multi_MSE_log(taskCondIdx));
        end
        
        batch_log{rep}.performance{cap}.nTaskCombinations = size(taskConditions,1);
        batch_log{rep}.performance{cap}.taskIndicies = taskIndicies;
        batch_log{rep}.performance{cap}.MSE = multiMSE;
        
    end

    progress = rep/replications;
    disp(['rep ' num2str(rep) '/' num2str(replications)]);
    disp(['progress: ' num2str(progress*100) '%']);
toc
    disp('---');


end


save(['logfiles/Nesreen_WeightedInterference_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(log_version)]);


