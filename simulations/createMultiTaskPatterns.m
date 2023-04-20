function [inputMulti tasksMulti trainMulti tasksIdxMulti] = createMultiTaskPatterns(NPathways, Nactive, NFeatures)
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

% specify number of possible task combinations (this number blows up as we
% increase NPathways) for multitasking condition (2 tasks at the same time)
% also specify the number of possible input stimuli (only one feature per
% feature dimension can be active)

% author: Sebastian Musslick

taskPaths = [zeros(1,NPathways-Nactive) 1:NPathways];
taskCombs = taskPaths;

for i = 2:1:NPathways
    taskCombs = combvec(taskCombs,taskPaths); % 1:NPathways
    for j = 1:(i-1)
        idicies = find(taskCombs(j,:) == taskCombs(i,:) & taskCombs(j,:) ~= 0);
        taskCombs(:,idicies) = [];
    end
end

% delete multitasking conditions that exceed specified number of active
% tasks
numNonZeroElements = taskCombs ~= 0;
numNonZeroElements = sum(numNonZeroElements);
taskCombs(:,numNonZeroElements ~= Nactive) = [];
taskCombs = transpose(unique(transpose(taskCombs),'rows'));

stimCombs = [1:NFeatures];
for i = 2:NPathways
    stimCombs = combvec(stimCombs,[1:NFeatures]); 
end
currFeaturesDims = 1:NPathways;%NPathways;

% create the task environment for performing 2 task at the same time (multitasking)
    
% loop through each new training set
input = [];
tasks = [];
train = [];
tasksIdx_Multi = [];

Nsets = size(taskCombs,2) * size(stimCombs,2);

for currTaskComb = 1:size(taskCombs,2)
    
    % build task input
    currTasksM = zeros(NPathways, NPathways);
    
    for i = 1:NPathways;
        if(taskCombs(i,currTaskComb) ~= 0)
            currTasksM(i,taskCombs(i,currTaskComb)) = 1;
        end
    end
    
    if(max(sum(currTasksM,1)) > 1) % hard constraint
        error('Overlapping tasks: Can''t use one output modality for two different feature dimensions');
    end
    if(max(sum(currTasksM,2)) > 1) % soft constraint (maybe remove soft constraint)
       error('Overlapping feature dimensions: Can''t perform two tasks on the same input features');
    end
    
    currTasks = currTasksM';
    currTasks = repmat(currTasks(:)', size(stimCombs,2), 1); % backtransform: reshape(currTasks(1,:),10,10)'

    % build feature input
    currInput = zeros(size(stimCombs,2),NPathways*NFeatures);
    
    for i = 1:size(currInput,1)
        currInput(i,(currFeaturesDims-1).*(NFeatures)+stimCombs(:,i)') = 1;
    end
    
    % build output
    currTrain = zeros(size(stimCombs,2),NPathways*NFeatures);
    
    for i = 1:size(currInput,1)
        noTaskMask = (taskCombs(:,currTaskComb) == 0)';
        onFeatures = (taskCombs(:,currTaskComb)'-1).*(NFeatures)+stimCombs(:,i)';
        onFeatures(noTaskMask) = [];
        
        currTrain(i,onFeatures) = 1;
    end

    % build full training set
    tasks = [tasks; currTasks];
    input = [input; currInput];
    train = [train; currTrain];
    tasksIdx_Multi = [tasksIdx_Multi; currTaskComb];
end

inputMulti = input;
tasksMulti = tasks;
trainMulti = train;
tasksIdxMulti = tasksIdx_Multi;

end