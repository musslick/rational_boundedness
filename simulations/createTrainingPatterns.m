function [inputSgl tasksSgl trainSgl inputMulti tasksMulti trainMulti tasksIdxSgl tasksIdxMulti] = createTrainingPatterns(NPathways, Nactive, NFeatures)
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
% author: Sebastian Musslick

[inputSgl tasksSgl trainSgl tasksIdxSgl] = createSglTaskPatterns(NPathways, Nactive, NFeatures);
[inputMulti tasksMulti trainMulti tasksIdxMulti] = createMultiTaskPatterns(NPathways, Nactive, NFeatures);

end