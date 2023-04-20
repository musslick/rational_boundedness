function [CDFs, LCA_AB_Accuracy, LCA_A_Accuracy, LCA_B_Accuracy] = performTownsendAnalysis(taskNet, taskPair, inputSgl, tasksSgl, trainSgl, tasksIdxSgl, multiCap, varargin)

    thresholdIdx = [];
    if(~isempty(varargin))
       thresholdIdx = varargin{1};
    end

    taskA = taskPair(1);
    taskB = taskPair(2);
    taskProfile = zeros(1, size(tasksSgl, 2));
    taskProfile(taskPair) = 1;
    
    % find corresponding test patterns
    cap = length(taskPair);
    patternIdx = find(ismember(multiCap{cap}.tasks, taskProfile, 'rows'));
    input_AB = multiCap{cap}.input(patternIdx, :);
    tasks_AB = multiCap{cap}.tasks(patternIdx, :);
    train_AB = multiCap{cap}.train(patternIdx, :);
    
    input_A = inputSgl(tasksIdxSgl == taskA,:);
    tasks_A = tasksSgl(tasksIdxSgl == taskA,:);
    train_A = trainSgl(tasksIdxSgl == taskA,:);
    
    input_B = inputSgl(tasksIdxSgl == taskB,:);
    tasks_B = tasksSgl(tasksIdxSgl == taskB,:);
    train_B = trainSgl(tasksIdxSgl == taskB,:);

    % LCA settings
    loadLCASettings;

    % Multi LCA call
    [~, ~, optAccuracy, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_AB] ...
    = taskNet.runLCA(LCA_settings, input_AB, tasks_AB, train_AB);
    LCA_AB_Accuracy = nanmean(optAccuracy);
    
    % Single task calls
    [~, ~, optAccuracy, ~, optThreshIdx_A, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_A] ...
    = taskNet.runLCA(LCA_settings, input_A, tasks_A, train_A);
    LCA_A_Accuracy = nanmean(optAccuracy);

    [~, ~, optAccuracy, ~, optThreshIdx_B, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, taskRT_correct_B] ...
    = taskNet.runLCA(LCA_settings, input_B, tasks_B, train_B);
    LCA_B_Accuracy = nanmean(optAccuracy);
    
    % set threshold to single task threshold
    thresholdIdx = round(nanmean(nanmean([optThreshIdx_A(:) optThreshIdx_B(:)])));
    
    % perform Townsend analysis
    [A_B, A_B_1, min_A_B, x, C]  = computeTownsendWengerCDF(taskRT_correct_AB, taskRT_correct_A, taskRT_correct_B, taskA, taskB, optThreshIdx_A, optThreshIdx_B, optThreshIdx_AB, thresholdIdx);
    CDFs.A_B = A_B;
    CDFs.A_B_1 = A_B_1;
    CDFs.min_A_B = min_A_B;
    CDFs.x = x;
    CDFs.C = C;

end