function [inputSgl, tasksSgl, trainSgl, tasksIdxSgl, multiCap, classLog, classificationFunction] = createTaskPatterns_GenN(NPathways, NFeatures, samplesPerTask, taskSimilarity, sameClassifierAcrossTasks, sameStimuliAcrossTasks, varargin)

% createTrainingPatterns generate single tasking and multitasking training
% patterns
% 
% PARAMETERS:
%
% NPathways     ...Number of pathways (i.e. number of feature dimensions % output dimensions)
% Nactive       ...How many tasks to perform simultaneously
% NFeatures     ...lists all feature dimensions
%
% RETURN VALUES:
%
% inputSgl      ...Single-task training set for input stimuli 
% tasksSgl      ...Single-task training set for input tasks
% trainSgl      ...Single-task training set for correct network outputs
%
% inputMulti    ...Multi-task training set for input stimuli 
% tasksMulti    ...Multi-task training set for input tasks
% trainMulti    ...Multi-task training set for correct network outputs
%

numInputUnits = NPathways^2*NFeatures;
numTaskUnits = NPathways^2;
numOutputUnits = NFeatures*NPathways;

numPatterns = numTaskUnits*samplesPerTask;
numSharedFeatures = round(taskSimilarity*NFeatures);
numNonSharedFeatures = NFeatures - numSharedFeatures;

inputSgl = nan(numPatterns, numInputUnits);
tasksSgl = nan(numPatterns, numTaskUnits);
trainSgl = nan(numPatterns, numOutputUnits);
tasksIdxSgl = nan(numPatterns, 1);
classLog = nan(samplesPerTask,numTaskUnits);

if(~isempty(varargin))
    useExistingClassification = 0;
    if(~isempty(varargin{1}))
        useExistingClassification = 1;
        classificationFunction = varargin{1};
    end
else
    useExistingClassification = 0;
    classificationFunction = {};
end

if(length(varargin) > 1)
   nHiddenLayer = varargin{2}; 
else
   nHiddenLayer = 1;
end

%% create single tasking patterns

% for each task
for taskIdx = 1:numTaskUnits
    
    if(useExistingClassification)
        
        weights = classificationFunction(taskIdx).weights;
        bias = classificationFunction(taskIdx).bias;
        percentiles = classificationFunction(taskIdx).percentiles;
        
    else
        
        if(taskIdx == 1 || ~sameClassifierAcrossTasks)
            [weights, bias, percentiles] = createClassificationFunction(NFeatures, nHiddenLayer);
        end

        classificationFunction(taskIdx).weights = weights;
        classificationFunction(taskIdx).bias = bias;
        classificationFunction(taskIdx).percentiles = percentiles;

    end
    
    % for each sample per task
    for patternIdx = 1:samplesPerTask

        row = (taskIdx-1)*samplesPerTask + patternIdx;

        if(taskIdx == 1 || ~sameStimuliAcrossTasks)
            % specify input
            inputSgl(row,:) = rand(1, numInputUnits);
        
        else
            inputSgl(row,:) = inputSgl(patternIdx,:);
        end

        % specify task
        tasksSgl(row,:) = zeros(1, numTaskUnits);
        tasksSgl(row,taskIdx) = 1;

        tasksIdxSgl(row) = taskIdx;

        % specify output

        outputDimension = mod(taskIdx-1, NPathways)+1;
        taskRelevantDimension = ceil(taskIdx/NPathways);
        dimensionPointer = NPathways*NFeatures*(taskRelevantDimension-1);

        taskRelevantInput = nan(1, NFeatures);

        % select overlapping features
        targetUnits = dimensionPointer+(1:numSharedFeatures);
        taskRelevantInput(1:numSharedFeatures) = inputSgl(row, targetUnits);

        % select non-overlapping features
        targetUnits = (dimensionPointer+numSharedFeatures+(outputDimension-1)*numNonSharedFeatures+1):(dimensionPointer+numSharedFeatures+outputDimension*numNonSharedFeatures);
        taskRelevantInput((numSharedFeatures+1):(numSharedFeatures+numNonSharedFeatures)) =  inputSgl(row, targetUnits);

        class = classifyInput(taskRelevantInput, weights, bias, percentiles);
        taskRelevantOutput = zeros(1, NFeatures);
        taskRelevantOutput(class) = 1;

        classLog(patternIdx, taskIdx) = class;
        trainSgl(row,:) = zeros(1, numOutputUnits);
        trainSgl(row, (NFeatures*(outputDimension-1)+1):(NFeatures*outputDimension)) = taskRelevantOutput;

    end
    
end

%% create multitasking patterns

multiCap{1}.input = inputSgl;
multiCap{1}.tasks = tasksSgl;
multiCap{1}.train = trainSgl;
multiCap{1}.taskCombs = 1:(NPathways^2);

% Create matrix to keep track of which rows of tasks correspond to
% which rows of taskCombs.
relevantTasks = transpose(1:(NPathways^2));
taskIdx = repmat(1:length(relevantTasks), [samplesPerTask,1]); taskIdx = taskIdx(:);
multiCap{1}.taskIdx = taskIdx;


for Nactive = 2:NPathways
    
    % specify number of possible task combinations (this number blows up as we
    % increase NPathways) for multitasking condition (2 tasks at the same time)
    % also specify the number of possible input stimuli (only one feature per
    % feature dimension can be active)
    taskPaths = [zeros(1,NPathways-Nactive) 1:NPathways];
    taskCombs = taskPaths;

    for k = 2:1:NPathways
        taskCombs = combvec(taskCombs,taskPaths); % 1:NPathways
        for j = 1:(k-1)
            idicies = find(taskCombs(j,:) == taskCombs(k,:) & taskCombs(j,:) ~= 0);
            taskCombs(:,idicies) = [];
        end
    end

    % delete multitasking conditions that exceed specified number of active
    % tasks
    numNonZeroElements = taskCombs ~= 0;
    numNonZeroElements = sum(numNonZeroElements);
    taskCombs(:,numNonZeroElements ~= Nactive) = [];
    taskCombs = transpose(unique(transpose(taskCombs),'rows'));

    inputMulti = nan(samplesPerTask * size(taskCombs,2),numInputUnits);
    tasksMulti = nan(samplesPerTask * size(taskCombs,2),numTaskUnits);
    trainMulti = nan(samplesPerTask * size(taskCombs,2),numOutputUnits);

    taskCombs_log = [];
    NUM_MULTI_TASKS = size(taskCombs,2);
    
    for currTaskComb = 1:size(taskCombs,2)

        % build task input
        currTasksM = zeros(NPathways, NPathways);

        for k = 1:NPathways;
            if(taskCombs(k,currTaskComb) ~= 0)
                currTasksM(k,taskCombs(k,currTaskComb)) = 1;
            end
        end

        if(max(sum(currTasksM,1)) > 1) % hard constraint
            error('Overlapping tasks: Can''t use one output modality for two different feature dimensions');
        end
        if(max(sum(currTasksM,2)) > 1) % soft constraint (maybe remove soft constraint)
            error('Overlapping feature dimensions: Can''t perform two tasks on the same input features');
        end

        currTasks = currTasksM';
        currTasks = repmat(currTasks(:)', samplesPerTask, 1); % backtransform: reshape(currTasks(1,:),10,10)'

        idx = find(currTasks(1,:) == 1);
        taskCombs_log = [taskCombs_log; idx];
        
        tasksMulti(((currTaskComb-1)*samplesPerTask+1):(currTaskComb*samplesPerTask),:) = currTasks;

        % build stimulus input
        if(sameStimuliAcrossTasks)
            inputMulti(((currTaskComb-1)*samplesPerTask+1):(currTaskComb*samplesPerTask),:) = inputSgl(1:samplesPerTask,:);
        else
            inputMulti(((currTaskComb-1)*samplesPerTask+1):(currTaskComb*samplesPerTask),:) = rand(samplesPerTask, numInputUnits);
        end

        % build training output

        currTrain = zeros(samplesPerTask,NPathways*NFeatures);

        activeTasks = find(currTasks(1,:) == 1);
        for k = 1:length(activeTasks)
            
            task = activeTasks(k);
            
            if(sameClassifierAcrossTasks)
                taskTemplate = zeros(1, size(currTasks,2));
                taskTemplate(task) = 1;
                
                taskData = trainSgl(ismember(tasksSgl, taskTemplate, 'rows'),:);
                currTrain = currTrain + taskData;
            else
                inputData = inputMulti(((currTaskComb-1)*samplesPerTask+1):(currTaskComb*samplesPerTask),:);
                taskIdx = task;
                                    
                for patternIdx = 1:size(inputData,1)
                    
                    row = (taskIdx-1)*samplesPerTask + patternIdx;
                    
                    outputDimension = mod(taskIdx-1, NPathways)+1;
                    taskRelevantDimension = ceil(taskIdx/NPathways);
                    dimensionPointer = NPathways*NFeatures*(taskRelevantDimension-1);

                    taskRelevantInput = nan(1, NFeatures);

                    % select overlapping features
                    targetUnits = dimensionPointer+(1:numSharedFeatures);
                    taskRelevantInput(1:numSharedFeatures) = inputSgl(row, targetUnits);

                    % select non-overlapping features
                    targetUnits = (dimensionPointer+numSharedFeatures+(outputDimension-1)*numNonSharedFeatures+1):(dimensionPointer+numSharedFeatures+outputDimension*numNonSharedFeatures);
                    taskRelevantInput((numSharedFeatures+1):(numSharedFeatures+numNonSharedFeatures)) =  inputSgl(row, targetUnits);

                    class = classifyInput(taskRelevantInput, weights, bias, percentiles);
                    taskRelevantOutput = zeros(1, NFeatures);
                    taskRelevantOutput(class) = 1;

                    classLog(patternIdx, taskIdx) = class;
                    currTrain(row,:) = zeros(1, numOutputUnits);
                    currTrain(row, (NFeatures*(outputDimension-1)+1):(NFeatures*outputDimension)) = taskRelevantOutput;

                end
            end
        end

        trainMulti(((currTaskComb-1)*samplesPerTask+1):(currTaskComb*samplesPerTask),:) = currTrain;
    end
    
    % Create matrix to keep track of which rows of tasks correspond to
    % which rows of taskCombs.
    taskIdx = repmat(1:NUM_MULTI_TASKS, [samplesPerTask,1]); taskIdx = taskIdx(:);
        
    
    multiCap{Nactive}.input = inputMulti;
    multiCap{Nactive}.tasks = tasksMulti;
    multiCap{Nactive}.train = trainMulti;
    multiCap{Nactive}.taskCombs = taskCombs_log;
    multiCap{Nactive}.taskIdx = taskIdx;

end    
    

end


function [weights, bias, percentiles] = createClassificationFunction(NFeatures, nHiddenLayer)

    nInput = NFeatures;
    nOutput = 1;
    nHidden = 3;

    weights.W_IH = (rand(nInput, nHidden)-0.5)*30;
    
    weights.W_HH = [];
    for i  = 1:(nHiddenLayer-1)
       weights.W_HH{i} = rand(nHidden, nHidden);
    end
    
    weights.W_HO = rand(nHidden, nOutput);

    bias.hidden = rand(1, nHidden) * -2.5;
    bias.output = rand(1, nOutput) * -2.5;

    % generate data

    nSamples = 1000;

    input = rand(nSamples, nInput);

    for i = 1:nHiddenLayer
        if(i == 1)
            net_hidden = input*weights.W_IH + repmat(bias.hidden, size(input,1), 1);
            act_hidden = 1./(1+exp(-net_hidden));
        else
            net_hidden = act_hidden*weights.W_HH{i-1} + repmat(bias.hidden, size(input,1), 1);
            act_hidden = 1./(1+exp(-net_hidden));
        end
    end

    net_output = act_hidden*weights.W_HO + repmat(bias.output, size(input,1), 1);
    act_output = 1./(1+exp(-net_output));

    P = (0:1/NFeatures:1)*100;
    percentiles = prctile(act_output,P);
    percentiles(1) = 0;
    percentiles(end) = 100;
    
    hist(act_output);

end

function [class] = classifyInput(taskRelevantInput, weights, bias, percentiles)

    % run model

    net_hidden = taskRelevantInput*weights.W_IH + repmat(bias.hidden, size(taskRelevantInput,1), 1);
    act_hidden = 1./(1+exp(-net_hidden));
    
    for i = 1:length(weights.W_HH)
        net_hidden = act_hidden*weights.W_HH{i} + repmat(bias.hidden, size(taskRelevantInput,1), 1);
        act_hidden = 1./(1+exp(-net_hidden));
    end

    net_output = act_hidden*weights.W_HO + repmat(bias.output, size(taskRelevantInput,1), 1);
    act_output = 1./(1+exp(-net_output));
    
    idx = find(act_output > percentiles);
    class = max(idx);

end
