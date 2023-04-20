function [task1RT, task1Accuracy, task2RT, task2Accuracy, task1Iterations, task2Iterations] = Part1_Sim4_PRP_Analysis(SOA, taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, varargin)

    if(~isempty(varargin))
        verbose = varargin{1};
    else
        verbose = 0;
    end
    
    plot = 0;
    useRegularStimuliForTask1RT = 0; % set to 0 in order to only present stimulus features that are task-relevant (all others will be turned off)
    findOptimalTaskBOnset = 1;  % set to 1 if the onset of Task B should be optimized
    minimumAccuracyTaskA = 0.8; % minimum accuracy required for Task A. This effects how much Task B will get delayed
    
    loadLCASettings;
    LCA_settings.maxTimeSteps = 100;
    LCA_settings.T0 = 0; % 0.2 % set T0 to zero to allow for better interpretability of RTs
    LCA_settings.responseThreshold = [0:0.1:1.5];

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
        for i = 1:size(input,1)
            A = reshape(input(i,:), 3, 3);
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
    
    % only select non-empty stimuli
    NPathways = taskNet.NPathways;
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
    
    % now create copy of input and train where dimension of Task B is empty
    input_taskA = input(tasksIdxSgl == taskA,:);
    train_taskA = train(tasksIdxSgl == taskA,:);
    input_taskB = input(tasksIdxSgl == taskB,:);
    train_taskB = train(tasksIdxSgl == taskB,:);
    
    % check if the stimuli for task A and task B are the same
    assert(isequal(input_taskA, input_taskB));
    
    if(~useRegularStimuliForTask1RT)
        
        NFeatures = size(input,2)/NPathways;
        
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
    
if(~findOptimalTaskBOnset) % don't optimize onset of Task B
    %% generate full pattern for simulating RT of task A
    
    input_phase1 = nan([size(inputA) LCA_settings.maxTimeSteps]);
    tasks_phase1 = nan([size(tasksA) LCA_settings.maxTimeSteps]);
    
    for t = 1:size(input_phase1, 3)
        
        if(~useRegularStimuliForTask1RT)
            % before SOA, set input corresponds to input of task A
            if(t < SOA)
                input_phase1(:,:,t) = inputA;
                tasks_phase1(:,:,t) = tasksA;
            end

            if(t >= SOA)
                input_phase1(:,:,t) = inputB;
                tasks_phase1(:,:,t) = tasksAB;
            end
       else
            input_phase1(:,:,t) = inputA;
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
    
    if(plot)
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
    taskA_Accuracy = nanmean(optTaskAccuracy(:,taskA));
    
    % determine optimal number of iterations for that task A
    tau_default = 0.1;
    numIterationsTask1 = round((taskA_RT - LCA_settings.T0)/(LCA_settings.dt_tau * tau_default));
    task1Iterations = numIterationsTask1;
    
else % OPTIMIZE ONSET OF TASK B
    
    if(useRegularStimuliForTask1RT)
        error('useRegularStimuliForTask1RT needs to be 0 when optimizing for onset of Task B.');
    end
    bestRR = nan;
    bestOnsetTaskB = SOA;
    
    for onsetTaskB = SOA:1:LCA_settings.maxTimeSteps
        
        %% generate full pattern for simulating RT of task A

        input_phase1 = nan([size(inputA) LCA_settings.maxTimeSteps]);
        tasks_phase1 = nan([size(tasksA) LCA_settings.maxTimeSteps]);

        for t = 1:size(input_phase1, 3)

            if(~useRegularStimuliForTask1RT)
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
           else
                input_phase1(:,:,t) = inputA;
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

        if(plot)
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
        taskA_Accuracy = nanmean(optTaskAccuracy(:,taskA));

        if(taskA_Accuracy > minimumAccuracyTaskA)
            bestOnsetTaskB = onsetTaskB;
            % determine optimal number of iterations for that task A
            tau_default = 0.1;
            numIterationsTask1 = round((taskA_RT - LCA_settings.T0)/(LCA_settings.dt_tau * tau_default));
            task1Iterations = numIterationsTask1;
            break;
        end
        

    end
    
end
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
        
        if(~findOptimalTaskBOnset)
            % turn off task unit A after response on Task A
            if(t > numIterationsTask1)
                    tasks_phase2(:,:,t) = tasksB;   
            else
                    tasks_phase2(:,:,t) = tasksAB;
            end
        else
            % turn off task unit A after response on Task A
            if(t > numIterationsTask1)
                    tasks_phase2(:,:,t) = tasksB;   
            else
                % turn on task B at onset
                if(t >= bestOnsetTaskB)
                    tasks_phase2(:,:,t) = tasksAB;
                else % otherwise keep just task A
                    tasks_phase2(:,:,t) = tasksA;
                end
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
    
    if(plot)
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
    [optTaskAccuracy, optTaskRT] ...
    = taskNet.runLCA(LCA_settings, [], tasksB, trainB, output_act_phase2);
    taskB_RT = nanmean(optTaskRT(:,taskB));
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
     task1Accuracy = taskA_Accuracy;
     task2RT = taskB_RT;
     task2Accuracy = taskB_Accuracy;
         
end


