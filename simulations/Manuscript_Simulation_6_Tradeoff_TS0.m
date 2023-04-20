function Manuscript_Simulation_6_Tradeoff_TS0(inputArg)

inputArgDouble = str2double(inputArg);
rep = inputArgDouble;

% meta simulation parameters
log_version = 2;

% set up network parameters
nHidden = 100;              % number of hidden units
hiddenPathSize = 1;         % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;         % group size of output units that receive the same weights from the task layer)
learningRate = 0.3;         % learning rate
thresh = 0.01;            % mean-squared error stopping criterion
decay = 0.0000;             % weight penalization parameter
bias = -2;                  % weight from bias units to hidden & output units
init_scale = 0.1;           % scales for initialized random weights
iterations_Train = 1000;     % number of training iterations
corr_thresholds = 0.1:0.05:0.9;      % threshold for detecting merged task representations in correlation matrix
goodPerformanceThresholds = [0.9 0.95 0.99]; % performance threshold for determining optimal control policy

% task environment parameters 
NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
NFeatures = 5;                     % the number of features per feature dimension
nTasks = NPathways^2;

samplesPerFeatureCombination = 1;

sameClassifierAcrossTasks = 1;
sameStimuliAcrossTasks = 1;
taskSimilarity = 0;
% taskSimilarities = [0 0.2 0.4 0.6 0.8 1];
init_taskCorrelations = linspace(-1, 0.95, 40);  % 9*4    % magnitude of initial weights
init_task_scale = 5;
orthogonalizeReferenceVectors = 0;      % flag: set to 1 if reference vectors should be entirely orthogonal
runLCA = 0;

% samplesPerTask = NFeatures^(NPathways)*samplesPerFeatureCombination;
%samplesPerTask_gen = NFeatures^(NPathways)*samplesPerFeatureCombination;

samplesPerTask = 50; 
samplesPerTask_gen = 100;

batch_log = repmat(struct('MSE_log',nan(length(init_taskCorrelations), iterations_Train), ...
                          'CE_log',nan(length(init_taskCorrelations), iterations_Train), ...
                          'CF_log',nan(length(init_taskCorrelations), iterations_Train), ...
                          'DimCF_log',nan(length(init_taskCorrelations), iterations_Train), ...
                          'train_MSE',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ...  
                          'train_CE',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ...  
                          'train_CF',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ... 
                          'train_DimCF',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ... 
                          'train_Accuracy',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ... 
                          'train_Pcorrect',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ... 
                          'test_MSE',nan(length(init_taskCorrelations), samplesPerTask_gen*nTasks), ...  
                          'test_CE',nan(length(init_taskCorrelations), samplesPerTask_gen*nTasks), ... 
                          'test_CF',nan(length(init_taskCorrelations), samplesPerTask_gen*nTasks), ... 
                          'test_DimCF',nan(length(init_taskCorrelations), samplesPerTask_gen*nTasks), ... 
                          'test_Accuracy',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ... 
                          'test_Pcorrect',nan(length(init_taskCorrelations), samplesPerTask*nTasks), ... 
                          'hiddenTaskRep',nan(length(init_taskCorrelations), iterations_Train, nTasks, nTasks), ...
                          'outputTaskRep',nan(length(init_taskCorrelations), iterations_Train, nTasks, nTasks), ...
                          'init_taskCorr',nan(length(init_taskCorrelations), 1), ...
                          'weightCorr',nan(length(init_taskCorrelations), 1), ...
                          'maxCarryingCapacity', nan(length(init_taskCorrelations), iterations_Train, length(corr_thresholds)), ...
                          'taskCardinalityLog_meanAcc', nan(length(init_taskCorrelations), iterations_Train, length(goodPerformanceThresholds)), ...
                          'taskCardinalityLog_respProb', nan(length(init_taskCorrelations), iterations_Train, length(goodPerformanceThresholds))),1, 1); 

tic
% create trainign environment
[inputSgl_cut tasksSgl_cut trainSgl_cut tasksIdxSgl_cut multiCap classLog classificationFunction] = createTaskPatterns_GenN(NPathways, NFeatures, samplesPerTask, taskSimilarity, sameClassifierAcrossTasks, sameStimuliAcrossTasks, [], length(nHidden));
[multiCap_con, multiCap_inc] = splitTrainingPatternsByCongruency(multiCap, NFeatures, NPathways);

for init_taskCorr_idx = 1:length(init_taskCorrelations)

    % initialize network
    taskNet = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
    taskNet.setData(inputSgl_cut, tasksSgl_cut, trainSgl_cut);
    taskNet.configure(); 
    taskNet.NPathways = NPathways;

    %% set up task weights according to max. cosine similarity

    weightCorr = nan(1,NPathways);

    % get task weights as a function of initial weight vector correlation
    init_taskCorr = init_taskCorrelations(init_taskCorr_idx);
    taskWeights = generateWeightBasisSet(nHidden, nTasks, init_taskCorr);

    taskWeights = taskWeights .*init_task_scale(end);

    % set up network weights
    taskNet.weights.W_TH = taskWeights;

    % log initial weight correlation
    batch_log.init_hiddenTaskRepWeights = corr(taskNet.weights.W_TH);
    batch_log.weightCorr(init_taskCorr_idx) = mean(weightCorr);
    batch_log.init_taskCorr(init_taskCorr_idx) = init_taskCorr;

    %% training loop
    for iter = 1:iterations_Train

        % train network for one iteration
        taskNet.trainOnline(1);
        
        if(taskNet.MSE_log < thresh)
            break
        end

        % record error
        batch_log.MSE_log(init_taskCorr_idx,iter) = taskNet.MSE_log;

        % run network on traning and test set
        [outData, hiddenData, MSE_train, hidden_net, output_net, ceError_train, classError_train, classDimError_train]  = taskNet.runSet(inputSgl_cut, tasksSgl_cut, trainSgl_cut);

        % store training performance
        batch_log.train_MSE(init_taskCorr_idx, iter) = mean(MSE_train);


        %% compute MIS

        % get task representations
        [R_hidden, R_output] = computeTaskSimilarity(taskNet, inputSgl_cut, tasksSgl_cut);

        % store task representations
        batch_log.hiddenTaskRep(init_taskCorr_idx, iter, :,:) = R_hidden;
        batch_log.outputTaskRep(init_taskCorr_idx, iter, :,:) = R_output;
        batch_log.hiddenTaskRepWeights(init_taskCorr_idx, iter, :,:) = corr(taskNet.weights.W_TH);

        %% test multitasking performance

        for cap = 2:length(multiCap) 

            % incongruent
            [outData, ~, MSE_multi, ~, ~, classError_multi, ceError_multi] = taskNet.runSet(multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train);
            [~, PCorrect_mean_train] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  multiCap_inc{cap}.train, multiCap_inc{cap}.tasks);

            batch_log.MSE_multi_inc{init_taskCorr_idx, iter, cap} = mean(MSE_multi);
            batch_log.PCorrect_multi_inc{init_taskCorr_idx, iter, cap} = mean(PCorrect_mean_train);
            
            % congruent
            [outData, ~, MSE_multi, ~, ~, classError_multi, ceError_multi] = taskNet.runSet(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train);
            [~, PCorrect_mean_train] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  multiCap_con{cap}.train, multiCap_con{cap}.tasks);

            batch_log.MSE_multi_con{init_taskCorr_idx, iter, cap} = mean(MSE_multi);
            batch_log.PCorrect_multi_con{init_taskCorr_idx, iter, cap} = mean(PCorrect_mean_train);
            
        end

        % print current training iteration
        disp([iter taskNet.MSE_log]);
    end
    
    % compute final multitaskin performance
    
    % overall
    if runLCA == 1
        loadLCASettings;
        cap = length(multiCap);
        
        [outData, ~, MSE_multi, ~, ~, classError_multi, ceError_multi] = taskNet.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
        [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
        = taskNet.runLCA(LCA_settings, multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train, outData);
        batch_log.LCA_Multi(init_taskCorr_idx, cap) = mean(optAccuracy);
    end

    %% print progress
    progress = init_taskCorr_idx/length(init_taskCorrelations);
    disp([ 'init_scale ' num2str(init_taskCorr_idx) '/' num2str(length(init_taskCorrelations))]);
    disp(['progress: ' num2str(progress*100) '%']);
    disp('---');
end
toc




save(['logfiles/Part2/PychReview_Part2_Sim4_Tradeoff_TS' num2str(taskSimilarity*100) '_r' num2str(rep)], '-v7.4');

end