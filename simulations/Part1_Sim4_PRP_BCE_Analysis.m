function [task1RT, task1Accuracy, task2RT, task2Accuracy, task1Iterations, task2Iterations, task1RT_all, task2RT_all, task1RT_distribution, task2RT_distribution] = Part1_Sim4_PRP_BCE_Analysis(SOA, taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, varargin)

    if(~isempty(varargin))
        verbose = varargin{1};
    else
        verbose = 0;
    end
    
    makePlot = 0;
    minimumAccuracyTaskA = 0.8; % minimum accuracy required for Task A. This effects how much Task B will get delayed
    maxTimeSteps = 50;
    selectNonEmptyStimuliOnly = 0;

    %% subsample congruent / incongruent stimuli
    if(~isempty(congruent))
        
        if(verbose)
            if(congruent == 1)
                disp('Congruent stimuli');
            elseif(congruent == 0)
                disp('Incongruent stimuli');
            end
        end
        
        congruentIdx = [];
        incongruentIdx = [];
        
        % identify congruent & incongruent stimuli for the two tasks
        inputDimTaskA  = ceil(taskA / taskNet.NPathways);
        inputDimTaskB = ceil(taskB / taskNet.NPathways);
        NFeatures = size(input, 2)/taskNet.NPathways;
        for i = 1:size(input,1)
            A = reshape(input(i,:), NFeatures, taskNet.NPathways);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end
        
    % prune
        if(congruent == 1)

            input(incongruentIdx,:) = [];
            train(incongruentIdx,:) = [];
            tasksIdxSgl(incongruentIdx,:) = [];
            
        elseif(congruent == 0)
            
            input(congruentIdx,:) = [];
            train(congruentIdx,:) = [];
            tasksIdxSgl(congruentIdx,:) = [];
            
        end
    end
    
    %% sample stimuli and zero out stimulus dimension of task 2 until SOA
    
    NPathways = taskNet.NPathways;
    NFeatures = size(input,2)/NPathways;
    
    if(selectNonEmptyStimuliOnly)
        % only select non-empty stimuli
        nonEmptyStimIdx = nan(1, size(input,1));
        for i = 1:size(input,1)
           stim = input(i,:);
           if(sum(stim) == NPathways)
               nonEmptyStimIdx(i) = 1;
           else
               nonEmptyStimIdx(i) = 0;
           end
        end

        input = input(nonEmptyStimIdx == 1, :);
        train = train(nonEmptyStimIdx == 1, :);
        tasksIdxSgl = tasksIdxSgl(nonEmptyStimIdx == 1, :);
    end
    
    % now create copy of input and train where dimension of Task B is empty
    input_taskA = input(tasksIdxSgl == taskA,:);
    train_taskA = train(tasksIdxSgl == taskA,:);
    input_taskB = input(tasksIdxSgl == taskB,:);
    train_taskB = train(tasksIdxSgl == taskB,:);
    
    % check if the stimuli for task A and task B are the same
    assert(isequal(input_taskA, input_taskB));

    inputDimTaskA = floor((taskA-1) / NPathways)+1;
    inputDimTaskAFeatures = (inputDimTaskA-1) * NFeatures + [1:NFeatures];
    inputDimNotTaskAFeatures = ones(1, NFeatures*NPathways);
    inputDimNotTaskAFeatures(inputDimTaskAFeatures) = 0;

    outputDimTaskB = mod(taskB-1, NPathways)+1;
    inputDimTaskB = floor((taskB-1) / NPathways)+1;
    outputDimTaskBFeatures = (outputDimTaskB-1) * NFeatures + [1:NFeatures];
    inputDimTaskBFeatures = (inputDimTaskB-1) * NFeatures + [1:NFeatures];

    inputDimNotTaskABFeatures = ones(1, NFeatures*NPathways);
    inputDimNotTaskABFeatures(inputDimTaskAFeatures) = 0;
    inputDimNotTaskABFeatures(inputDimTaskBFeatures) = 0;

    for i = 1:size(input_taskA,1)

        input_taskA(i,inputDimNotTaskAFeatures == 1) = 0;
        train_taskA(i,outputDimTaskBFeatures) = 0;
    end

    for i = 1:size(input_taskB,1)

        input_taskB(i,inputDimNotTaskABFeatures == 1) = 0;

    end

    
    % finally, sample stimuli
    
    % sample stimuli
    if(~isempty(nStimuli))
        stimuli = randi(size(train_taskA,1),1,nStimuli);
        inputA = input_taskA(stimuli, :);
        inputB = input_taskB(stimuli, :);
        trainA = train_taskA(stimuli,:);
        trainB = train_taskB(stimuli,:);
    else
        inputA = input_taskA;
        inputB = input_taskB;
        trainA = train_taskA;
        trainB = train_taskB;
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
    
    if(max(size(taskB)) == 1)
        tasksAB = zeros(size(inputB,1), taskNet.NPathways^2);
        task(:, taskA) = 1;
        task(:, taskB) = 1;
        tasksAB = task;
    end
    
probedOnsets = SOA:1:(SOA+min(max(50, SOA*5), maxTimeSteps-SOA));
    
totalRewardRate = nan(1, length(probedOnsets));
outcomes = {};

for onsetTaskB = probedOnsets
    
    disp(['Testing time step ' num2str(onsetTaskB) '/' num2str(length(probedOnsets))]);

    loadLCASettings;
    LCA_settings.maxTimeSteps = maxTimeSteps;
    LCA_settings.T0 = 0; % 0.2 % set T0 to zero to allow for better interpretability of RTs
%     LCA_settings.ITI = 3.250;
    
    %% generate full pattern for simulating RT of task A

    input_phase1 = nan([size(inputA) LCA_settings.maxTimeSteps]);
    tasks_phase1 = nan([size(tasksA) LCA_settings.maxTimeSteps]);

    for t = 1:size(input_phase1, 3)

        % before SOA, set input corresponds to input of task A
        if(t < SOA)
            input_phase1(:,:,t) = inputA;
            tasks_phase1(:,:,t) = tasksA;
        end

        if(t >= SOA)
            input_phase1(:,:,t) = inputB;
        end

        if(t >= onsetTaskB)
            tasks_phase1(:,:,t) = tasksAB;
        else
            tasks_phase1(:,:,t) = tasksA;
        end

    end

    % start from initial condition at first time step
    input_phase1(:,:,1) = 0;
    tasks_phase1(:,:,1) = 0;

    %% compute RT of first task

    % get output patterns
    input_series = input_phase1;
    tasks_series = tasks_phase1;
    [output_act_phase1] = taskNet.integrate(tauNet, tauTask, input_series, tasks_series);
    output_act_phase1 = permute(output_act_phase1, [1 3 2]);

    if(makePlot)
       patternIdx = 2;
       subplot(1,3,1);
       imagesc(squeeze(input_series(patternIdx,:,:)));
       subplot(1,3,2);
       imagesc(squeeze(tasks_series(patternIdx,:,:)));
       subplot(1,3,3);
       imagesc(squeeze(output_act_phase1(patternIdx,:,:)));
    end

    % run LCA analysis to determine number of time steps for first task

    % [optTaskAccuracy, optTaskRT] ...
    [optTaskAccuracy, optTaskRT, optAccuracy, optRT, optThreshIdx, RR, taskAccuracy, taskRT, taskRT_all, meanAccuracy, meanRT, outputProbabilities, RT_Simulations, noThresholdReachedPatternTaskPairs, minThresholdReachedPatternTaskPairs, maxThresholdReachedPatternTaskPairs, optRT_correct, optRT_incorrect, taskRT_correct, taskRT_incorrect] ...
    = taskNet.runLCA(LCA_settings, [], tasksA, trainA, output_act_phase1);
    taskA_RT = nanmean(optTaskRT(:,taskA));
    taskA_RT_all = optTaskRT(:,taskA);
    
    taskA_RT_distribution = nan(size(taskRT_all,2), size(taskRT_all,4));
    for pattern = 1:size(optThreshIdx,1)
        optimalThresh = optThreshIdx(pattern, taskA);
        taskA_RT_distribution(pattern,:) = transpose(squeeze(taskRT_all(optimalThresh, pattern, taskA,:)));
    end
    
    taskA_Accuracy = nanmean(optTaskAccuracy(:,taskA));

    % determine optimal number of iterations for that task A
    tau_default = 0.1;
    numIterationsTask1 = round((taskA_RT - LCA_settings.T0)/(LCA_settings.dt_tau * tau_default));
    task1Iterations = numIterationsTask1;


    %% now generate input based on response of first task
    
    input_phase2 = nan([size(inputA) LCA_settings.maxTimeSteps]);
    tasks_phase2 = nan([size(tasksA) LCA_settings.maxTimeSteps]);
    
    for t = 1:size(input_phase2, 3)
        
        % before SOA, set input corresponds to input of task A
        if(t < SOA)
            input_phase2(:,:,t) = inputA;
            tasks_phase2(:,:,t) = tasksA;
        end
        
        % after SOA, set input corresponds to input of task B
        if(t >= SOA)
            input_phase2(:,:,t) = inputB;
        end
        
        % turn off task unit A after response on Task A
        if(t > numIterationsTask1)
                tasks_phase2(:,:,t) = tasksB;   
        else
            % turn on task B at onset
            if(t >= onsetTaskB)
                tasks_phase2(:,:,t) = tasksAB;
            else % otherwise keep just task A
                tasks_phase2(:,:,t) = tasksA;
            end
        end

        
    end
    
    % start from initial condition at first time step
    input_phase2(:,:,1) = 0;
    tasks_phase2(:,:,1) = 0;
    
    %% compute RT for second task

    % get output patterns
    input_series = input_phase2;
    tasks_series = tasks_phase2;
    [output_act_phase2] = taskNet.integrate(tauNet, tauTask, input_series, tasks_series);
    output_act_phase2 = permute(output_act_phase2, [1 3 2]);
    output_act_phase2_org = output_act_phase2;
    
    % make sure to start with onset of task 2 stimulus (at SOA)
    output_act_phase2(:,:,1:SOA) = [];
    LCA_settings.maxTimeSteps = size(output_act_phase2,3);
    
    if(makePlot)
       patternIdx = 1;
       subplot(1,4,1);
       imagesc(squeeze(input_series(patternIdx,:,:)));
       subplot(1,4,2);
       imagesc(squeeze(tasks_series(patternIdx,:,:)));
       subplot(1,4,3);
       imagesc(squeeze(output_act_phase2_org(patternIdx,:,:)));
       subplot(1,4,4);
       imagesc(squeeze(output_act_phase2(patternIdx,:,:)));
    end
    
    % run LCA analysis to determine RT for second task
    [optTaskAccuracy, optTaskRT, optAccuracy, optRT, optThreshIdx, RR, taskAccuracy, taskRT, taskRT_all, meanAccuracy, meanRT, outputProbabilities, RT_Simulations, noThresholdReachedPatternTaskPairs, minThresholdReachedPatternTaskPairs, maxThresholdReachedPatternTaskPairs, optRT_correct, optRT_incorrect, taskRT_correct, taskRT_incorrect] ...
    = taskNet.runLCA(LCA_settings, [], tasksB, trainB, output_act_phase2);
    taskB_RT = nanmean(optTaskRT(:,taskB));
    taskB_RT_all = optTaskRT(:,taskB);
    
    taskB_RT_distribution = nan(size(taskRT_all,2), size(taskRT_all,4));
    for pattern = 1:size(optThreshIdx,1)
        optimalThresh = optThreshIdx(pattern, taskB);
        taskB_RT_distribution(pattern,:) = transpose(squeeze(taskRT_all(optimalThresh, pattern, taskB,:)));
    end
    
    taskB_Accuracy = nanmean(optTaskAccuracy(:,taskB));
    
    % determine optimal number of iterations for that task A
    tau_default = 0.1;
    numIterationsTask2 = round((taskB_RT - LCA_settings.T0)/(LCA_settings.dt_tau * tau_default));
    task2Iterations = numIterationsTask2;
    
     if(verbose)
         disp(['Accuracy Task A: ' num2str(mean(taskA_Accuracy))]);
         disp(['Accuracy Task B: ' num2str(mean(taskB_Accuracy))]);
         disp(['RT Task A: ' num2str(mean(taskA_RT))]);
         disp(['RT Task B: ' num2str(mean(taskB_RT))]);
     end

     task1RT = taskA_RT;
     task1RT_all = taskA_RT_all;
     task1RT_distribution = taskA_RT_distribution;
     task1Accuracy = taskA_Accuracy;
     task2RT = taskB_RT;
     task2RT_all = taskB_RT_all;
     task2RT_distribution = taskB_RT_distribution;
     task2Accuracy = taskB_Accuracy;
     
     % compute reward rate
     currentIdx = onsetTaskB-SOA+1;
     totalRewardRate(currentIdx) = task1Accuracy * task2Accuracy / (task2RT + LCA_settings.ITI + SOA * (LCA_settings.dt_tau * tau_default));
     outcome{currentIdx}.task1RT = task1RT;
     outcome{currentIdx}.task1RT_all = task1RT_all;
     outcome{currentIdx}.task1RT_distribution = task1RT_distribution;
     outcome{currentIdx}.task1Accuracy = task1Accuracy;
     outcome{currentIdx}.task2RT = task2RT;
     outcome{currentIdx}.task2RT_all = task2RT_all;
     outcome{currentIdx}.task2RT_distribution = task2RT_distribution;
     outcome{currentIdx}.task2Accuracy = task2Accuracy;
     
end

bestRRIdx = find(totalRewardRate == max(totalRewardRate));
if(bestRRIdx == length(totalRewardRate))
    totalRewardRate
    task1RT
    task1Accuracy
    task2RT
    task2Accuracy
    SOA
    onsetTaskB
    tauNet
    warning(['tested onset range not large enough: ' num2str(bestRRIdx)]);
end
disp(['best onset index: ' num2str(bestRRIdx)]);

task1RT = outcome{bestRRIdx}.task1RT;
task1RT_all = outcome{bestRRIdx}.task1RT_all;
task1RT_distribution = outcome{bestRRIdx}.task1RT_distribution;
task1Accuracy = outcome{bestRRIdx}.task1Accuracy;
task2RT = outcome{bestRRIdx}.task2RT;
task2RT_all = outcome{bestRRIdx}.task2RT_all;
task2RT_distribution = outcome{bestRRIdx}.task2RT_distribution;
task2Accuracy = outcome{bestRRIdx}.task2Accuracy;


