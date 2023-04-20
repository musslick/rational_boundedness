function Manuscript_Simulation_5_BasisVsTensorTraining(inputArg)

inputArgDouble = str2double(inputArg);
rep = inputArgDouble;

tic

% meta simulation parameters
log_version = 1;

% set up network parameters
nHidden = 100;              % number of hidden units
hiddenPathSize = 1;         % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;         % group size of output units that receive the same weights from the task layer)
learningRate = 0.3;         % learning rate
thresh = 0.001;            % mean-squared error stopping criterion
decay = 0.0000;             % weight penalization parameter
bias = -2;                  % weight from bias units to hidden & output units
init_scale = 0.1;           % scales for initialized random weights
iterations_MultiTrain = 1000;     % number of training iterations for multitask training
MDSMetric = 'euclidean'; % 'euclidean', 'correlation', 'cosine'
multiTrainingProportion = [0 0.2 0.4 0.6 0.8];
constantNumPatterns = 1;

% task environment parameters 
NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
NFeatures = NPathways;                     % the number of features per feature dimension

sameStimuliAcrossTasks = 1;
sdScale = 0;
samplesPerTask = [];
tasksToPerform = 1:(NPathways^2);
nOutput = NFeatures*NPathways;
if(constantNumPatterns)
    numSamplesPerTrainingIteration = 500;
else
    numSamplesPerTrainingIteration = NFeatures^NFeatures;
end

cap = NPathways;                               % number of tasks to be multitasked
nTasks = length(tasksToPerform);

% create full training environment
[inputSgl, tasksSgl, trainSgl, ~, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
[multiCap_con, multiCap_inc] = splitTrainingPatternsByCongruency(multiCap, NFeatures, NPathways);
nStimuliSingleTaskTraining = size(inputSgl, 1)/nTasks;

% initialize network
taskNet = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
taskNet.setData(inputSgl, tasksSgl, trainSgl);
taskNet.configure(); 
taskNet.NPathways = NPathways;

% training protocol
%
% train network on particular multitasking condition
% - generate 2 experimental conditions: 
% - network 2 continues to be trained on single + multitask (congruent)
% - network 3 continues to be trained on single + multitask (incongruent) 
% - track task representations over course of learning
% - track task similarities over course of learning
      
%%  TRAINING

for curr_multiTrainingProportion_idx = 1:length(multiTrainingProportion)

    curr_multiTrainingProportion = multiTrainingProportion(curr_multiTrainingProportion_idx);
    batch_log(1).multiTrainingProportion(curr_multiTrainingProportion_idx) = curr_multiTrainingProportion;

    % generate 3 network conditions
    taskNet_SingleTaskTraining = NNmodel(taskNet);
    taskNet_MultitaskConTraining = NNmodel(taskNet);
    taskNet_MultitaskIncTraining = NNmodel(taskNet);

    taskNet_SingleTaskTraining.thresh = 0;
    taskNet_MultitaskConTraining.thresh = 0;
    taskNet_MultitaskIncTraining.thresh = 0;

    disp('multitasking training...');
    % train networks on multitasking (& single tasks)
    for iter = 1:iterations_MultiTrain
        
        if(constantNumPatterns == 0)
            % vary number of patterns PER TASK
            
            numSglPatterns = round(numSamplesPerTrainingIteration * (1-curr_multiTrainingProportion));
            numMultiPatterns = numSamplesPerTrainingIteration - numSglPatterns;

            % subsample stimuli for single task training
            [inputSgl_Phase2_base, tasksSgl_Phase2_base, trainSgl_Phase2_base] = ...
                                                    subSampleStimuli(inputSgl, tasksSgl, trainSgl, numSglPatterns);

            % subsample stimuli for congruent multitasking training
            [inputMultiCon_Phase2, tasksMultiCon_Phase2, trainMultiCon_Phase2] = ...
                                                    subSampleStimuli(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train, numMultiPatterns);

            % subsample stimuli for incongruent multitasking training
            [inputMultiInc_Phase2, tasksMultiInc_Phase2, trainMultiInc_Phase2] = ...
                                                    subSampleStimuli(multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train, numMultiPatterns);

        else
           % vary number of pattern PER SINGLE VS. MULTI CONDITION
           numSglPatterns = round(numSamplesPerTrainingIteration * (1-curr_multiTrainingProportion));
           numMultiPatterns = numSamplesPerTrainingIteration - numSglPatterns;
            
           totalSglPatterns = size(inputSgl,1);
           totalMultiConPatterns = size(multiCap_con{cap}.input,1);
           totalMultiIncPatterns = size(multiCap_con{cap}.input,1);
           
           sglSamples = randsample(totalSglPatterns, numSglPatterns, 1);
           multiConSamples = randsample(totalMultiConPatterns, numMultiPatterns, 1);
           multiIncSamples = randsample(totalMultiIncPatterns, numMultiPatterns, 1);
           
           inputSgl_Phase2_base = inputSgl(sglSamples, :);
           tasksSgl_Phase2_base = tasksSgl(sglSamples, :);
           trainSgl_Phase2_base = trainSgl(sglSamples, :);
           
           inputMultiCon_Phase2 = multiCap_con{cap}.input(multiConSamples, :);
           tasksMultiCon_Phase2 = multiCap_con{cap}.tasks(multiConSamples, :);
           trainMultiCon_Phase2 = multiCap_con{cap}.train(multiConSamples, :);
           
           inputMultiInc_Phase2 = multiCap_inc{cap}.input(multiIncSamples, :);
           tasksMultiInc_Phase2 = multiCap_inc{cap}.tasks(multiIncSamples, :);
           trainMultiInc_Phase2 = multiCap_inc{cap}.train(multiIncSamples, :);
           
        end
        
        % set data for congruent multitasking training network
        input_MultitaskConTraining = [inputSgl_Phase2_base; inputMultiCon_Phase2];
        tasks_MultitaskConTraining = [tasksSgl_Phase2_base; tasksMultiCon_Phase2];
        train_MultitaskConTraining = [trainSgl_Phase2_base; trainMultiCon_Phase2];

        taskNet_MultitaskConTraining.setData(input_MultitaskConTraining, tasks_MultitaskConTraining, train_MultitaskConTraining);

        % set data for incongruent multitasking training network
        input_MultitaskIncTraining = [inputSgl_Phase2_base; inputMultiInc_Phase2];
        tasks_MultitaskIncTraining = [tasksSgl_Phase2_base; tasksMultiInc_Phase2];
        train_MultitaskIncTraining = [trainSgl_Phase2_base; trainMultiInc_Phase2];

        taskNet_MultitaskIncTraining.setData(input_MultitaskIncTraining, tasks_MultitaskIncTraining, train_MultitaskIncTraining);
        
        % compute similarity of task representations for congruent multitasking training network
        [R_hidden, R_output, hiddenTaskRep, outputTaskRep] = computeTaskSimilarity(taskNet_MultitaskConTraining, inputSgl, tasksSgl, 'tasks', tasksToPerform);

        batch_log.MultitaskConTraining_R_hidden_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = R_hidden;
        batch_log.MultitaskConTraining_R_output_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = R_output;
        batch_log.MultitaskConTraining_hiddenTaskRep_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = hiddenTaskRep;
        batch_log.MultitaskConTraining_outputTaskRep_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = outputTaskRep;

        % compute similarity of task representations for incongruent multitasking training network
        [R_hidden, R_output, hiddenTaskRep, outputTaskRep] = computeTaskSimilarity(taskNet_MultitaskIncTraining, inputSgl, tasksSgl, 'tasks', tasksToPerform);

        batch_log.MultitaskIncTraining_R_hidden_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = R_hidden;
        batch_log.MultitaskIncTraining_R_output_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = R_output;
        batch_log.MultitaskIncTraining_hiddenTaskRep_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = hiddenTaskRep;
        batch_log.MultitaskIncTraining_outputTaskRep_Phase2(curr_multiTrainingProportion_idx, :, :, iter) = outputTaskRep;

        % compute network performance for congruent multitask training
        [~,~,meanMSE] = taskNet_MultitaskConTraining.runSet(inputSgl, tasksSgl, trainSgl);
        batch_log.MultitaskConTraining_MSE_Single_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);
        
        [~,~,meanMSE] = taskNet_MultitaskConTraining.runSet(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train);
        batch_log.MultitaskConTraining_MSE_conMulti_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);

        [~,~,meanMSE] = taskNet_MultitaskConTraining.runSet(multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train);
        batch_log.MultitaskConTraining_MSE_incMulti_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);
        
        [~,~,meanMSE] = taskNet_MultitaskConTraining.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
        batch_log.MultitaskConTraining_MSE_Multi_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);

        % compute network performance for incongruent multitask training
        [~,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(inputSgl, tasksSgl, trainSgl);
        batch_log.MultitaskIncTraining_MSE_Single_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);
        
        [~,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train);
        batch_log.MultitaskIncTraining_MSE_conMulti_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);

        [~,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train);
        batch_log.MultitaskIncTraining_MSE_incMulti_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);

        [~,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
        batch_log.MultitaskIncTraining_MSE_Multi_timeCourse(curr_multiTrainingProportion_idx, iter) = mean(meanMSE);

        % train networks for one iteration
        taskNet_MultitaskConTraining.trainOnline(1);
        taskNet_MultitaskIncTraining.trainOnline(1);

        % print current training iteration
        disp(num2str(iter));
    end

    % compute network performance for congruent multitask training
    [~,~,meanMSE] = taskNet_MultitaskConTraining.runSet(inputSgl, tasksSgl, trainSgl);
    batch_log.MultitaskConTraining_MSE_Single(curr_multiTrainingProportion_idx) = mean(meanMSE);
    
    [outData_con,~,meanMSE] = taskNet_MultitaskConTraining.runSet(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train);
    batch_log.MultitaskConTraining_MSE_conMulti(curr_multiTrainingProportion_idx) = mean(meanMSE);

    [outData_inc,~,meanMSE] = taskNet_MultitaskConTraining.runSet(multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train);
    batch_log.MultitaskConTraining_MSE_incMulti(curr_multiTrainingProportion_idx) = mean(meanMSE);

    [outData,~,meanMSE] = taskNet_MultitaskConTraining.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
    batch_log.MultitaskConTraining_MSE_Multi(curr_multiTrainingProportion_idx) = mean(meanMSE);
    
    loadLCASettings;
    
%     [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
%     = taskNet_MultitaskConTraining.runLCA(LCA_settings, multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train, outData_inc);
%     batch_log.MultitaskConTraining_LCA_incMulti(curr_multiTrainingProportion_idx) = mean(optAccuracy);
%     
%     [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
%     = taskNet_MultitaskConTraining.runLCA(LCA_settings, multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train, outData_con);
%     batch_log.MultitaskConTraining_LCA_conMulti(curr_multiTrainingProportion_idx) = mean(optAccuracy);

    [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
    = taskNet_MultitaskConTraining.runLCA(LCA_settings, multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train, outData);
    batch_log.MultitaskConTraining_LCA_Multi(curr_multiTrainingProportion_idx) = mean(optAccuracy);

    % compute network performance for incongruent multitask training
    [~,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(inputSgl, tasksSgl, trainSgl);
    batch_log.MultitaskIncTraining_MSE_Single(curr_multiTrainingProportion_idx) = mean(meanMSE);
    
    [outData_con,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train);
    batch_log.MultitaskIncTraining_MSE_conMulti(curr_multiTrainingProportion_idx) = mean(meanMSE);

    [outData_inc,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train);
    batch_log.MultitaskIncTraining_MSE_incMulti(curr_multiTrainingProportion_idx) = mean(meanMSE);

    [outData,~,meanMSE] = taskNet_MultitaskIncTraining.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
    batch_log.MultitaskIncTraining_MSE_Multi(curr_multiTrainingProportion_idx) = mean(meanMSE);
    
    loadLCASettings;
    
%     [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
%     = taskNet_MultitaskIncTraining.runLCA(LCA_settings, multiCap_inc{cap}.input, multiCap_inc{cap}.tasks, multiCap_inc{cap}.train, outData_inc);
%     batch_log.MultitaskIncTraining_LCA_incMulti(curr_multiTrainingProportion_idx) = mean(optAccuracy);
%     
%     [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
%     = taskNet_MultitaskIncTraining.runLCA(LCA_settings, multiCap_con{cap}.input, multiCap_con{cap}.tasks, multiCap_con{cap}.train, outData_con);
%     batch_log.MultitaskIncTraining_LCA_conMulti(curr_multiTrainingProportion_idx) = mean(optAccuracy);

    [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
    = taskNet_MultitaskIncTraining.runLCA(LCA_settings, multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train, outData);
    batch_log.MultitaskIncTraining_LCA_Multi(curr_multiTrainingProportion_idx) = mean(optAccuracy);

    %% postprocess data

    % generate template to access similar t
    basis_template = eye(nTasks,nTasks);
    for row = 1:size(basis_template,1)
        basis_template(row, (ceil(row/NPathways)-1)*NPathways+(1:NPathways)) = 1;
        basis_template(row, row:end) = 0;
    end

    % compute mean task correlation among similar tasks for phase 1

    for  iter = 1:iterations_MultiTrain

        % MultitaskConTraining
        R_hidden = squeeze(batch_log.MultitaskConTraining_R_hidden_Phase2(curr_multiTrainingProportion_idx,:,:,iter));
        R_output = squeeze(batch_log.MultitaskConTraining_R_output_Phase2(curr_multiTrainingProportion_idx,:,:,iter));

        batch_log.MultitaskConTraining_similarTasksCorr_hidden(curr_multiTrainingProportion_idx,iter) = mean(R_hidden(basis_template==1));
        batch_log.MultitaskConTraining_similarTasksCorr_output(curr_multiTrainingProportion_idx,iter) = mean(R_output(basis_template==1));

        % MultitaskIncTraining
        R_hidden = squeeze(batch_log.MultitaskIncTraining_R_hidden_Phase2(curr_multiTrainingProportion_idx,:,:,iter));
        R_output = squeeze(batch_log.MultitaskIncTraining_R_output_Phase2(curr_multiTrainingProportion_idx,:,:,iter));

        batch_log.MultitaskIncTraining_similarTasksCorr_hidden(curr_multiTrainingProportion_idx,iter) = mean(R_hidden(basis_template==1));
        batch_log.MultitaskIncTraining_similarTasksCorr_output(curr_multiTrainingProportion_idx,iter) = mean(R_output(basis_template==1));

    end

    % perform MDS on task representations across all networks & training iterations
    % MDS for all inputs

    % phase 2 - congruent multitasking training
    MultitaskConTraining_hiddenData_Phase1 = nan(nTasks, nHidden);
    MultitaskConTraining_outputData_Phase1 = nan(nTasks, nOutput);

    idx = 1;
    MultitaskConTraining_taskIdx_Phase2 = nan(1, nTasks);
    MultitaskConTraining_iterationIdx_Phase2 = nan(1, nTasks);

    for task = 1:nTasks
        MultitaskConTraining_hiddenData_Phase1(idx, :) = batch_log.MultitaskConTraining_hiddenTaskRep_Phase2(curr_multiTrainingProportion_idx,:, task, end);
        MultitaskConTraining_outputData_Phase1(idx, :) = batch_log.MultitaskConTraining_outputTaskRep_Phase2(curr_multiTrainingProportion_idx,:, task, end);
        MultitaskConTraining_taskIdx_Phase2(idx) = tasksToPerform(task);
        MultitaskConTraining_iterationIdx_Phase2(idx) = iter;
        idx = idx + 1;
    end

    % phase 2 - incongruent multitasking training
    MultitaskIncTraining_hiddenData_Phase1 = nan(nTasks, nHidden);
    MultitaskIncTraining_outputData_Phase1 = nan(nTasks, nOutput);

    idx = 1;
    MultitaskIncTraining_taskIdx_Phase2 = nan(1, nTasks);
    MultitaskIncTraining_iterationIdx_Phase2 = nan(1, nTasks);

    for task = 1:nTasks
        MultitaskIncTraining_hiddenData_Phase1(idx, :) = batch_log.MultitaskIncTraining_hiddenTaskRep_Phase2(curr_multiTrainingProportion_idx, :, task, end);
        MultitaskIncTraining_outputData_Phase1(idx, :) = batch_log.MultitaskIncTraining_outputTaskRep_Phase2(curr_multiTrainingProportion_idx, :, task, end);
        MultitaskIncTraining_taskIdx_Phase2(idx) = tasksToPerform(task);
        MultitaskIncTraining_iterationIdx_Phase2(idx) = iter;
        idx = idx + 1;
    end

    % merge task representation data
    hiddenData = [
                          MultitaskConTraining_hiddenData_Phase1; ...
                          MultitaskIncTraining_hiddenData_Phase1];

    outputData = [
                          MultitaskConTraining_outputData_Phase1; ...
                          MultitaskIncTraining_outputData_Phase1];

    % perform MDS on hidden data
    X = hiddenData;
    dist = MDSMetric;
    distances = pdist(X, dist);
    try
    Y = mdscale(distances,2);
    catch
       disp(input);
       Y = zeros(size(input,1), 2);
       warning('No MDS possible for average single task reps'); 
    end
    batch_log.MultitaskConTraining_MDS_hidden_Phase1{curr_multiTrainingProportion_idx} = Y(1:nTasks, :);
    batch_log.MultitaskIncTraining_MDS_hidden_Phase1{curr_multiTrainingProportion_idx} = Y((1:nTasks) + nTasks, :);

    % perform MDS on output data
    X = outputData;
    dist = MDSMetric;
    distances = pdist(X, dist);
    try
    Y = mdscale(distances,2);
    catch
       Y = zeros(size(inputSgl,1), 2);
       warning('No MDS possible for average single task reps'); 
    end
    batch_log.MultitaskConTraining_MDS_output_Phase1{curr_multiTrainingProportion_idx} = Y(1:nTasks, :);
    batch_log.MultitaskIncTraining_MDS_output_Phase1{curr_multiTrainingProportion_idx} = Y(1:nTasks + nTasks, :);

end

%% save log file
toc

save(['logfiles/Part2/PychReview_Part2_Sim2_BasisVsTensorTraining_' num2str(NPathways) 'P' num2str(NFeatures) 'F_v' num2str(log_version) '_h' num2str(nHidden(1)) '_r' num2str(rep)], 'batch_log', '-v7.3');

end