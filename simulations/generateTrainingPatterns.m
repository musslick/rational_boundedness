clear all;
clc;

%% training environment parameters

NPathways = 3;                                      % number of pathways (i.e. number of stimulus input dimensions & output dimensions)
NFeatures = 4;                                      % number of feature units per stimulus input dimension
samplesPerFeatureCombination = 1;       % how many istances of each stimulus should be in each task set
sdScale = 0;                                          % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;                 % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)
nTasks = NPathways^2;                       % specify number of tasks

% generate task environment
[inputSgl, tasksSgl, trainSgl, tasksIdxSgl, stimIdxSgl, inputSgl_mask, tasksSgl_mask, trainSgl_mask, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerFeatureCombination, sdScale, sameStimuliAcrossTasks);

%% single task patterns

% inputSgl and tasksSgl indicate all possible stimulus-task combinations.

% these are all the input stimuli patterns. Each row constitutes a stimulus
% which NPathways * Feature units. In this example the stimulus input is 
% grouped into 3 dimensions. Each dimension contains 4 feature units, 
% with one feature activated per input dimension.
disp('Example of all input stimuli patterns');
disp(inputSgl);

% these are all task patterns. Each row pattern denotes a one-hot vector
% indicating the current task to be executed. A task indicates a mapping
% from a particular input to a particular output dimension of the network.
% As there are 3 input and 3 output dimensions in total, there are 3 * 3 =
% 9 possible tasks. In the single-task training patterns only one of those
% tasks is activated (set to 1). 
disp('Example of all task input patterns');
disp(tasksSgl);

% these are all correct response patterns for each stimulus-task pair. 
% Each output pattern is encoded among R = NPathway * NFeatures 
% units. For a given single task, only one  response dimension is relevant 
% while all other response dimensions are supressed (set to 0). In for
% each relevant response dimension only one output  unit is allowed to 
% be active (set to 1), namely the one that corresponds to the active feature of the
% relevant stimulus input dimension (according to the current task).
disp('Example of all correct response patterns');
disp(trainSgl);

%% multitaksing patterns
% This script also generates multitasking patterns. This means that e.g. two
% tasks can be executed at the same time, by mapping two different input
% dimensions to two different output dimensions.
% note that multiCap only contains multitask conditions in which all simultaneously 
% executed tasks rely on different input and output patterns.
displayMultitaskingPatterns = 1;

% example multitasking patterns for performing 2 tasks at the same time
if(displayMultitaskingPatterns)
    disp('MULTITASKING 2 TASKS AT THE SAME TIME');
    disp('Input stimuli patterns:');
    disp(multiCap{2}.input);
    disp('Input task patterns:');
    disp(multiCap{2}.tasks);
    disp('Correct response patterns:');
    disp(multiCap{2}.train);
end

% example multitasking patterns for performing 3 tasks at the same time
if(displayMultitaskingPatterns)
    disp('MULTITASKING 3 TASKS AT THE SAME TIME');
    disp('Input stimuli patterns:');
    disp(multiCap{3}.input);
    disp('Input task patterns:');
    disp(multiCap{3}.tasks);
    disp('Correct response patterns:');
    disp(multiCap{3}.train);
end