function [transitionAccuracy, transitionRT] = Sim1c_transitionAnalysis(taskNet, inputA, inputB, taskA, taskB, trainA, trainB, nStimuli, tauNet, tauTask, iterations)

    % sample stimuli
    if(~isempty(nStimuli))
        stimuliA = randi(size(inputA,1),1,nStimuli);
        stimuliB = randi(size(inputB,1),1,nStimuli);
        inputA = inputA(stimuliA, :);
        inputB = inputB(stimuliB, :);
        trainA = trainA(stimuliA,:);
        trainB = trainB(stimuliB,:);
    end 
   
    % configure task
    if(max(size(taskA)) == 1)
        task = zeros(size(inputA,1), taskNet.NPathways^2);
        task(:, taskA) = 1;
        tasksA = task;
    end
    
    if(max(size(taskB)) == 1)
        task = zeros(size(inputB,1), taskNet.NPathways^2);
        task(:, taskB) = 1;
        tasksB = task;
    end
    
    % compute initial condition
   [outputB] = taskNet.switchTasks(tauNet, tauTask, iterations, inputA, inputB, tasksA, tasksB, trainA, trainB);
   outputB = permute(outputB, [1 3 2]);
   % figure(); imagesc(squeeze(outputB(1,:,:))); colorbar; caxis([0 1]);
    % LCA settings
    loadLCASettings;

    % LCA call
    LCA_settings.responseThreshold = [0:0.1:1];         % tested threshold
    LCA_settings.numSimulations = 100;
    LCA_settings.dt_tau = 0.1;
%     LCA_settings.dt_tau = 0.05;
    
    [~, ~, optAccuracy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, optRT_correct] ...
    = taskNet.runLCA(LCA_settings, inputB, tasksB, trainB, outputB);

    transitionAccuracy = optAccuracy;
    transitionRT = optRT_correct;
    
end