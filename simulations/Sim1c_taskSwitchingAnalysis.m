function [repeatAccuracy, switchAccuracy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, varargin)

    if(~isempty(varargin))
        verbose = varargin{1};
    else
        verbose = 0;
    end
    
    
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
    
     % generate network outputs as a function of time
     inputA = input( tasksIdxSgl == taskA,:);
     inputA = zeros(size(inputA));
     trainA = train( tasksIdxSgl == taskA,:);

     % Switch
     inputB = input( tasksIdxSgl == taskB,:);
     trainB = train( tasksIdxSgl == taskB,:);
     [switchAccuracy, switchRT] = Sim1c_transitionAnalysis(taskNet, inputA, inputB, taskA, taskB, trainA, trainB, nStimuli, tauNet, tauTask, iterations);

     if(verbose)
         disp('=========');
         disp('switch:');
         mean(switchAccuracy)
         mean(switchRT)
     end

     % Repeat
     taskB = taskA;
     inputB = input( tasksIdxSgl == taskB,:);
     trainB = train(tasksIdxSgl == taskB,:);
     [repeatAccuracy, repeatRT] = Sim1c_transitionAnalysis(taskNet, inputA, inputB, taskA, taskB, trainA, trainB, nStimuli, tauNet, tauTask, iterations);

     if(verbose)
         disp('repeat:');
         mean(repeatAccuracy)
         mean(repeatRT)
         
         disp(['switch cost (accuracy): ' num2str(mean(repeatAccuracy-switchAccuracy))]);
         disp(['switch cost (RT): ' num2str(mean(switchRT-repeatRT))]);
     end

         
end


