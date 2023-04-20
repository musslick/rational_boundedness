function [repeatAccuracy, switchAccuracy, repeatRT, switchRT, inputB_switch, inputB_repeat] = Part1_Sim3_TaskSwitching_Analysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, varargin)

    univalentStimuli = 0;
    
    if(~isempty(varargin))
        verbose = varargin{1};
        
        if(length(varargin) >= 2)
            univalentStimuli = varargin{2};
        end
        
    else
        verbose = 0;
    end
    
    NFeatures = size(input,2) / taskNet.NPathways;
    
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

    % adjust number of time steps
    loadLCASettings;
    iterations = LCA_settings.maxTimeSteps;
    
     % SWITCH
    
     % generate network outputs as a function of time
     inputA = input( tasksIdxSgl == taskA,:);
%      inputA = zeros(size(inputA));
     trainA = train( tasksIdxSgl == taskA,:);

     inputB = input( tasksIdxSgl == taskB,:);
     trainB = train( tasksIdxSgl == taskB,:);
     [switchAccuracy, switchRT, ~, ~, inputB_switch] = Part1_Sim3_Transition_Analysis(taskNet, inputA, inputB, taskA, taskB, trainA, trainB, nStimuli, tauNet, tauTask, iterations);
     % figure(); imagesc(squeeze(outputB(1,:,:))); colorbar; caxis([0 1]);
     
     if(univalentStimuli)
        
        inputDimTaskA = floor((taskA-1) / taskNet.NPathways)+1;
        inputDimTaskAFeatures = (inputDimTaskA-1) * NFeatures + [1:NFeatures];
        inputDimNotTaskAFeatures = ones(1, NFeatures*taskNet.NPathways);
        inputDimNotTaskAFeatures(inputDimTaskAFeatures) = 0;
        
        inputDimTaskB = floor((taskB-1) / taskNet.NPathways)+1;
        inputDimTaskBFeatures = (inputDimTaskB-1) * NFeatures + [1:NFeatures];
        inputDimNotTaskBFeatures = ones(1, NFeatures*taskNet.NPathways);
        inputDimNotTaskBFeatures(inputDimTaskBFeatures) = 0;
        
        inputA(:, inputDimNotTaskAFeatures == 1) = 0;
        inputB(:, inputDimNotTaskBFeatures == 1) = 0;
    end
     
     if(verbose)
         disp('=========');
         disp('switch:');
         mean(switchAccuracy)
         mean(switchRT)
     end

     % REPEAT
     inputA = input( tasksIdxSgl == taskB,:);
     inputA = zeros(size(inputA));
     trainA = train( tasksIdxSgl == taskB,:);

     inputB = input( tasksIdxSgl == taskB,:);
     trainB = train( tasksIdxSgl == taskB,:);
     
     if(univalentStimuli)
        
        inputDimTaskA = floor((taskA-1) / taskNet.NPathways)+1;
        inputDimTaskAFeatures = (inputDimTaskA-1) * NFeatures + [1:NFeatures];
        inputDimNotTaskAFeatures = ones(1, NFeatures*taskNet.NPathways);
        inputDimNotTaskAFeatures(inputDimTaskAFeatures) = 0;
        
        inputDimTaskB = floor((taskB-1) / taskNet.NPathways)+1;
        inputDimTaskBFeatures = (inputDimTaskB-1) * NFeatures + [1:NFeatures];
        inputDimNotTaskBFeatures = ones(1, NFeatures*taskNet.NPathways);
        inputDimNotTaskBFeatures(inputDimTaskBFeatures) = 0;
        
        inputA(:, inputDimNotTaskAFeatures == 1) = 0;
        inputB(:, inputDimNotTaskBFeatures == 1) = 0;
    end
     
     taskA = taskB;
     [repeatAccuracy, repeatRT, ~, ~, inputB_repeat] = Part1_Sim3_Transition_Analysis(taskNet, inputA, inputB, taskA, taskB, trainA, trainB, nStimuli, tauNet, tauTask, iterations);
     % figure(); imagesc(squeeze(outputB(1,:,:))); colorbar; caxis([0 1]);
     
     if(verbose)
         disp('repeat:');
         mean(repeatAccuracy)
         mean(repeatRT)
         
         disp(['switch cost (accuracy): ' num2str(mean(repeatAccuracy-switchAccuracy))]);
         disp(['switch cost (RT): ' num2str(mean(switchRT-repeatRT))]);
     end

         
end


