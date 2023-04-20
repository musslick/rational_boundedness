function Manuscript_Simulation_4_FeatureOverlap(inputArg)

inputArgDouble = str2double(inputArg);
rep = inputArgDouble;


% meta simulation parameters
log_version = 7;

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
tasksToPerform = 1:nTasks;

samplesPerFeatureCombination = 1;

sameClassifierAcrossTasks = 1;
sameStimuliAcrossTasks = 1;
taskSimilarities = [0 0.2 0.4 0.6 0.8 1];
init_task_scale = init_scale;

% samplesPerTask = NFeatures^(NPathways)*samplesPerFeatureCombination;
%samplesPerTask_gen = NFeatures^(NPathways)*samplesPerFeatureCombination;

samplesPerTask = 50; % 50
samplesPerTask_gen = 100;

batch_log = repmat(struct('MSE_log',nan(length(taskSimilarities),1), ...
                          'CE_log',nan(length(taskSimilarities),1), ...
                          'CF_log',nan(length(taskSimilarities),1), ...
                          'DimCF_log',nan(length(taskSimilarities),1), ...
                          'train_MSE',nan(length(taskSimilarities), samplesPerTask*nTasks), ...  
                          'train_CE',nan(length(taskSimilarities), samplesPerTask*nTasks), ...  
                          'train_CF',nan(length(taskSimilarities), samplesPerTask*nTasks), ... 
                          'train_DimCF',nan(length(taskSimilarities), samplesPerTask*nTasks), ... 
                          'train_Accuracy',nan(length(taskSimilarities), samplesPerTask*nTasks), ... 
                          'train_Pcorrect',nan(length(taskSimilarities), samplesPerTask*nTasks), ... 
                          'test_MSE',nan(length(taskSimilarities), samplesPerTask_gen*nTasks), ...  
                          'test_CE',nan(length(taskSimilarities), samplesPerTask_gen*nTasks), ... 
                          'test_CF',nan(length(taskSimilarities), samplesPerTask_gen*nTasks), ... 
                          'test_DimCF',nan(length(taskSimilarities), samplesPerTask_gen*nTasks), ... 
                          'test_Accuracy',nan(length(taskSimilarities), samplesPerTask*nTasks), ... 
                          'test_Pcorrect',nan(length(taskSimilarities), samplesPerTask*nTasks), ... 
                          'hiddenTaskRep',nan(length(taskSimilarities), nTasks, nTasks), ...
                          'outputTaskRep',nan(length(taskSimilarities), nTasks, nTasks), ...
                          'hiddenTaskRepWeights',nan(length(taskSimilarities), nHidden, nTasks), ...
                          'outputTaskRepWeights',nan(length(taskSimilarities), NFeatures*NPathways, nTasks), ...
                          'maxCarryingCapacity', nan(length(taskSimilarities), length(corr_thresholds)), ...
                          'taskCardinalityLog_meanAcc', nan(length(taskSimilarities), length(goodPerformanceThresholds)), ...
                          'taskCardinalityLog_respProb', nan(length(taskSimilarities), length(goodPerformanceThresholds))),1, 1); 

tic

for taskSimilarity_idx = 1:length(taskSimilarities)

    taskSimilarity = taskSimilarities(taskSimilarity_idx);
    batch_log(1).featureOverlap(taskSimilarity_idx) = taskSimilarity;
    
    % create trainign environment
    [inputSgl_cut, tasksSgl_cut, trainSgl_cut, ~, multiCap, ~, classificationFunction] = createTaskPatterns_GenN(NPathways, NFeatures, samplesPerTask, taskSimilarity, sameClassifierAcrossTasks, sameStimuliAcrossTasks, [], length(nHidden));
    [inputSgl_gen, tasksSgl_gen, trainSgl_gen] = createTaskPatterns_GenN(NPathways, NFeatures, samplesPerTask_gen, taskSimilarity, sameClassifierAcrossTasks, sameStimuliAcrossTasks, classificationFunction);
    [multiCap_con, multiCap_inc] = splitTrainingPatternsByCongruency(multiCap, NFeatures, NPathways);

    % initialize network
    taskNet = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
    taskNet.setData(inputSgl_cut, tasksSgl_cut, trainSgl_cut);
    taskNet.configure(); 
    taskNet.NPathways = NPathways;

    %% training loop
    for iter = 1:iterations_Train

        % train network for one iteration
        taskNet.trainOnline(1);
        
        if(taskNet.MSE_log < thresh)
            break
        end

        % record error
        batch_log.MSE_log(taskSimilarity_idx, iter) = taskNet.MSE_log(end);
        
        % print current training iteration
        disp([iter taskNet.MSE_log]);
    end
    

    % run network on traning and test set
    [outData, hiddenData, MSE_train, hidden_net, output_net, ceError_train, classError_train, classDimError_train]  = taskNet.runSet(inputSgl_cut, tasksSgl_cut, trainSgl_cut);
    [Accuracy_mean_train] = taskNet.calculateAbsoluteErrorTasks(outData,  trainSgl_cut, tasksSgl_cut);
    [~, PCorrect_mean_train] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  trainSgl_cut, tasksSgl_cut);

    [outData, hiddenData, MSE_test, hidden_net, output_net, ceError_test, classError_test, classDimError_test] = taskNet.runSet(inputSgl_gen, tasksSgl_gen, trainSgl_gen);
    [Accuracy_mean_test] = taskNet.calculateAbsoluteErrorTasks(outData, trainSgl_gen, tasksSgl_gen);
    [~, PCorrect_mean_test] = taskNet.calculateOutcomeProbabilitiesTasks(outData, trainSgl_gen, tasksSgl_gen);

    % store training performance
    batch_log.train_MSE(taskSimilarity_idx) = mean(MSE_train);
    batch_log.train_CE(taskSimilarity_idx) = mean(ceError_train);
    batch_log.train_CF(taskSimilarity_idx) = mean(classError_train);
    batch_log.train_DimCF(taskSimilarity_idx) = mean(classDimError_train);
    batch_log.train_Accuracy(taskSimilarity_idx) = mean(Accuracy_mean_train);
    batch_log.train_PCorrect(taskSimilarity_idx) = mean(PCorrect_mean_train);

    % store testing performance
    batch_log.test_MSE(taskSimilarity_idx) = mean(MSE_test);
    batch_log.test_CE(taskSimilarity_idx) = mean(ceError_test);
    batch_log.test_CF(taskSimilarity_idx) = mean(classError_test);
    batch_log.test_DimCF(taskSimilarity_idx) = mean(classDimError_test);
    batch_log.test_Accuracy(taskSimilarity_idx) = mean(Accuracy_mean_test);
    batch_log.test_PCorrect(taskSimilarity_idx) = mean(PCorrect_mean_test);

    %% compute MIS

    % get task representations
    [R_hidden, R_output] = computeTaskSimilarity(taskNet, inputSgl_cut, tasksSgl_cut);

    % store task representations
    batch_log.hiddenTaskRep(taskSimilarity_idx, :,:) = R_hidden;
    batch_log.outputTaskRep(taskSimilarity_idx, :,:) = R_output;
    batch_log.hiddenTaskRepWeights(taskSimilarity_idx, :,:) = taskNet.weights.W_TH;
    batch_log.outputTaskRepWeights(taskSimilarity_idx, :,:) = taskNet.weights.W_TO;

    %% test multitasking performance

    for cap = 2:length(multiCap) 

        % incongruent

        [outData, ~, MSE_multi, ~, ~, classError_multi, ceError_multi] = taskNet.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
        [Accuracy_mean_train] = taskNet.calculateAbsoluteErrorTasks(outData,  multiCap{cap}.train, multiCap{cap}.tasks);
        [~, PCorrect_mean_train] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  multiCap{cap}.train, multiCap{cap}.tasks);

        % LCA settings
        loadLCASettings;

        [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
        = taskNet.runLCA(LCA_settings, multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train, outData);

        batch_log.MSE_multi{taskSimilarity_idx, cap} = mean(MSE_multi);
        batch_log.LCA_multi{taskSimilarity_idx, cap} = mean(optAccuracy);

        % congruent

%         [outData, ~, MSE_multi, ~, ~, classError_multi, ceError_multi] = taskNet.runSet(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train);
%         [Accuracy_mean_train] = taskNet.calculateAbsoluteErrorTasks(outData,  multiCap_con{cap}.train, multiCap_con{cap}.tasks);
%         [~, PCorrect_mean_train] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  multiCap_con{cap}.train, multiCap_con{cap}.tasks);
% 
%         [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
%         = taskNet.runLCA(LCA_settings, multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train, outData);
% 
%         batch_log.MSE_multi_con{taskSimilarity_idx, cap} = mean(MSE_multi);
%         batch_log.LCA_multi_con{taskSimilarity_idx, cap} = mean(optAccuracy);
%         batch_log.classError_multi_con{taskSimilarity_idx, cap} = mean(classError_multi);
%         batch_log.ceError_multi_con{taskSimilarity_idx, cap} = mean(ceError_multi);
%         batch_log.Accuracy_multi_con{taskSimilarity_idx, cap} = mean(Accuracy_mean_train);
%         batch_log.PCorrect_multi_con{taskSimilarity_idx, cap} = mean(PCorrect_mean_train);

    end

    %% print progress
    progress = taskSimilarity_idx/length(taskSimilarities);
    disp([ 'feature overlap ' num2str(taskSimilarity_idx) '/' num2str(length(taskSimilarities))]);
    disp(['progress: ' num2str(progress*100) '%']);
    disp('---');
end
toc




save(['logfiles/Part2/PychReview_Part2_Sim1_FeatureOverlap_' num2str(NPathways) 'P' num2str(NFeatures) 'F_v' num2str(log_version) '_h' num2str(nHidden(1))  '_r' num2str(rep)], '-v7.3');

end