function [dualTaskCombs, taskCombData] = validateDualTaskDependencies(taskNet, NPathways, NFeatures, taskAccuracies, A_tasksIdx, inputSgl, tasksSgl, multiCap, CDFs)

            % test MSE of individual tasks depending on dependency type 3 (indirect, asummetric)
            
            % dependency types
            %
            %  O    O
            %  |      |
            %  O    O   type 0 (independent)
            %
            %
            %   O       type 1 (fan-in)
            %  / \   
            %
            %  \ /
            %   O       type 2 (fan-out)
            %
            %  O    O
            %  |  \  |
            %  O    O   type 3 (indirect, asymmetric)
            %
            %  O    O
            %  |  /\ |
            %  O    O   type 4 (indirect, symmetric)
            %
            
            % combinations of all dual task conditions
            dualTaskCombs = transpose(multiCap{2}.taskCombs);
            dualTaskCombs = [dualTaskCombs; zeros(1, size(dualTaskCombs,2))];
            
            % for each condition, check type of dependency
            
            for dualCondIdx = 1:size(dualTaskCombs,2)
                
                task1 = dualTaskCombs(1, dualCondIdx);
                task2 = dualTaskCombs(2, dualCondIdx);
                
                % determine dependency type
                [rowTask1, colTask1] = find(A_tasksIdx == task1);
                [rowTask2, colTask2] = find(A_tasksIdx == task2);
                
                if(colTask1 == colTask2 && rowTask1 ~= rowTask2) % fan-in dependency
                    dualTaskCombs(3, dualCondIdx) = 1;
                end
                
                if(rowTask1 == rowTask2 && colTask1 ~= colTask2) % fan-out dependency
                    dualTaskCombs(3, dualCondIdx) = 2;
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && ((A_tasksIdx(rowTask2, colTask1) ~= 0 ...
                              && A_tasksIdx(rowTask1, colTask2) == 0) )) % asymmetric indirect dependency
                    dualTaskCombs(3, dualCondIdx) = 3;
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && A_tasksIdx(rowTask2, colTask1) == 0 ...
                        && A_tasksIdx(rowTask1, colTask2) ~= 0) % asymmetric indirect dependency (mirrored)
                    dualTaskCombs(3, dualCondIdx) = 3;
                    dualTaskCombs(1:2, dualCondIdx) = flipud(dualTaskCombs(1:2, dualCondIdx));
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && A_tasksIdx(rowTask2, colTask1) ~= 0 ...
                        && A_tasksIdx(rowTask1, colTask2) ~= 0) % symmetric indirect dependency
                    dualTaskCombs(3, dualCondIdx) = 4;
                end
                
                if(rowTask1 ~= rowTask2 && colTask1 ~= colTask2 ...
                        && A_tasksIdx(rowTask2, colTask1) == 0 ...
                        && A_tasksIdx(rowTask1, colTask2) == 0) % independent tasks
                    dualTaskCombs(3, dualCondIdx) = 0;
                    dualTaskCombs(1:2, dualCondIdx) = dualTaskCombs(randperm(2), dualCondIdx);  % randomly assign task 1 & task 2
                end
                
            end

            % test multitasking performance for critical conditions (type 3 & 4)
            
            relevantTaskCombs = dualTaskCombs(:, dualTaskCombs(3,:) == 3 | dualTaskCombs(3,:) == 4  | dualTaskCombs(3,:) == 0);
            
            taskCombData.responseActivation_con = nan(size(relevantTaskCombs, 2), 3 * 2);   % mean activation of response units for both tasks on congruent trials
            taskCombData.responseActivation_inc = nan(size(relevantTaskCombs, 2), 3 * 2);   % mean activation of response units for both tasks on incongruent trials
            taskCombData.similarity = nan(size(relevantTaskCombs, 2), 2);                   %  similarity for between dual-tasks & potentially interfering tasks
            
            for combIdx = 1:size(relevantTaskCombs,2);
                
                tasks = relevantTaskCombs(1:2, combIdx);
                
                % identify relevant output dimension
                currentTaskComb = zeros(1, NPathways^2);
                currentTaskComb(tasks(1)) = tasks(1);
                currentTaskComb(tasks(2)) = tasks(2);
                taskM = reshape(currentTaskComb, NPathways, NPathways);
                [relevantOutputDims(1) relevantInputDim(1)] = find(taskM == tasks(1));
                [relevantOutputDims(2) relevantInputDim(2)] = find(taskM == tasks(2));
                
                 % generate mask to extract relevant output dimension
                relevantOutputMask = zeros(NFeatures, NPathways);
                relevantOutputMask(:, relevantOutputDims) = repmat(tasks', NFeatures, 1);
                relevantOutputMask = relevantOutputMask(:)';
                
                % generate mask to extract relevant input dimension
                relevantInputMask = zeros(NFeatures, NPathways);
                relevantInputMask(:, relevantInputDim) = repmat(tasks', NFeatures, 1);
                relevantInputMask = relevantInputMask(:)';

                % find corresponding test patterns
                currentTaskComb(currentTaskComb > 0) = 1;
                patternIdx = find(ismember(multiCap{2}.tasks, currentTaskComb, 'rows'));
                inputPatterns = multiCap{2}.input(patternIdx, :);
                tasksPatterns = multiCap{2}.tasks(patternIdx, :);
                trainPatterns = multiCap{2}.train(patternIdx, :);

                % compute output for all task patterns of a task
                [outputPatterns] = taskNet.runSet(inputPatterns, tasksPatterns, trainPatterns);
                
                % compute MSE for executed tasks
                
                
                % find relevant features
                task1InputMask = relevantInputMask == tasks(1);
                task2InputMask = relevantInputMask == tasks(2);
                task1OutputMask = relevantOutputMask == tasks(1);
                task2OutputMask = relevantOutputMask == tasks(2);
                
                % compute raw feature activation
                congruentCounter = 0;
                incongruentCounter = 0;
                responseActivations_con = zeros(1, 3 * 2);
                responseActivations_inc = zeros(1, 3 * 2);
                congruency = nan(1, size(inputPatterns, 2));
                
                for stimIdx = 1:size(inputPatterns, 1)

                    % store features of relevant input dimension
                    task1CorrectFeature = find(inputPatterns(stimIdx, task1InputMask) ==1);
                    task2CorrectFeature = find(inputPatterns(stimIdx, task2InputMask) ==1);
                    irrelevantFeatures = find(inputPatterns(stimIdx, task1InputMask) ~=1 & inputPatterns(stimIdx, task2InputMask) ~=1);
                    
                    % store activation values for output dimensions of each task
                    task1ResponseActivation = outputPatterns(stimIdx, task1OutputMask);
                    task2ResponseActivation = outputPatterns(stimIdx, task2OutputMask);
                    
                    % check congruency
                    if(task1CorrectFeature == task2CorrectFeature)
                        % congruent
                        congruency(stimIdx) = 1;
                        responseActivations_con_curr = nan(1, 3 * 2);
                        
                        % task 1 features
                        responseActivations_con_curr(1) = task1ResponseActivation(task1CorrectFeature); % 1st feature is relevant feature == feature of irrelevant task
                        responseActivations_con_curr(2) = mean(task1ResponseActivation(irrelevantFeatures)); % 2nd feature is irrelevant
                        responseActivations_con_curr(3) = mean(task1ResponseActivation(irrelevantFeatures)); % 3rd feature is irrelevant
                        
                        % task 2 features
                        responseActivations_con_curr(4) = task2ResponseActivation(task1CorrectFeature); % 1st feature is relevant feature
                        responseActivations_con_curr(5) = mean(task2ResponseActivation(irrelevantFeatures)); % 2nd feature is irrelevant
                        responseActivations_con_curr(6) = mean(task2ResponseActivation(irrelevantFeatures)); % 3rd feature is irrelevant
                        
                        % add to existing data
                        responseActivations_con = responseActivations_con + responseActivations_con_curr;
                        congruentCounter = congruentCounter + 1;
                    else
                        % incongruent
                        congruency(stimIdx) = 0;
                        responseActivations_inc_curr = nan(1, 3 * 2);
                        
                        % task 1 features
                        responseActivations_inc_curr(1) = task1ResponseActivation(task1CorrectFeature); % 1st feature is relevant feature
                        responseActivations_inc_curr(2) = task1ResponseActivation(task2CorrectFeature); % 2nd feature is incongruent feature (from competing task)
                        responseActivations_inc_curr(3) = mean(task1ResponseActivation(irrelevantFeatures)); % 3rd feature is irrelevant
                        
                        % task 2 features
                        responseActivations_inc_curr(4) = task2ResponseActivation(task1CorrectFeature); % 1st feature is relevant feature
                        responseActivations_inc_curr(5) = task2ResponseActivation(task2CorrectFeature); % 2nd feature is incongruent feature (from competing task)
                        responseActivations_inc_curr(6) = mean(task2ResponseActivation(irrelevantFeatures)); % 3rd feature is irrelevant
                        
                        % add to existing data
                        responseActivations_inc = responseActivations_inc + responseActivations_inc_curr;
                        incongruentCounter = incongruentCounter + 1;
                    end
                    
                end
                
%                 %% compute task performance for congruent trials
%                 inputPatterns_con = inputPatterns(congruency == 1, :);
%                 tasksPatterns_con = tasksPatterns(congruency == 1, :);
%                 trainPatterns_con = trainPatterns(congruency == 1, :);
% 
%                  % compute output for all task patterns of a task
%                 [outputPatterns_con] = taskNet.runSet(inputPatterns_con, tasksPatterns_con, trainPatterns_con);
% 
%                 %% compute task performance for incongruent trials
%                 inputPatterns_inc = inputPatterns(congruency == 0, :);
%                 tasksPatterns_inc = tasksPatterns(congruency == 0, :);
%                 trainPatterns_inc = trainPatterns(congruency == 0, :);
% 
%                 % compute output for all task patterns of a task
%                 [outputPatterns_inc] = taskNet.runSet(inputPatterns_inc, tasksPatterns_inc, trainPatterns_inc);

                %% test similarity between dual-tasks & potentially interfering tasks

                if(relevantTaskCombs(3, combIdx) == 0) % tasks are independent (no shared representations)
                    
                    task1Similarity = nan;
                    task2Similarity = nan;
                    
                else    
                    
                    % find interfering task
                    [rowTask1, colTask1] = find(A_tasksIdx == tasks(1));
                    [rowTask2, colTask2] = find(A_tasksIdx == tasks(2));
                    task3 = A_tasksIdx(rowTask2, colTask1);

                    % get relevant task patterns
                    currentTasks = zeros(2, NPathways^2);
                    currentTasks(1, tasks(2)) = 1;
                    currentTasks(2, task3) = 1;

                    % find corresponding test patterns
                    patternIdx = find(ismember(tasksSgl, currentTasks, 'rows'));
                    inputPatterns = inputSgl(patternIdx, :);
                    tasksPatterns = tasksSgl(patternIdx, :);

                    % compute (default) task pattern similarity
                    [hiddenSimilarity] = computeTaskSimilarity(taskNet, inputPatterns, tasksPatterns);
                    task1Similarity = hiddenSimilarity(2,1);

                    % add similarity for symmetric case
                    if(relevantTaskCombs(3, combIdx) == 3)

                        task2Similarity = nan;

                    elseif(relevantTaskCombs(3, combIdx) == 4) % add similarity for symmetric case

                        task4 = A_tasksIdx(rowTask2, colTask1);

                        % get relevant task patterns
                        currentTasks = zeros(2, NPathways^2);
                        currentTasks(1, tasks(1)) = 1;
                        currentTasks(2, task4) = 1;

                        % find corresponding test patterns
                        patternIdx = find(ismember(tasksSgl, currentTasks, 'rows'));
                        inputPatterns = inputSgl(patternIdx, :);
                        tasksPatterns = tasksSgl(patternIdx, :);

                        [hiddenSimilarity] = computeTaskSimilarity(taskNet, inputPatterns, tasksPatterns);
                        task2Similarity = hiddenSimilarity(2,1);

                    end
                
                end
                
                % average response activations across stimuli
                responseActivations_con = responseActivations_con / congruentCounter; % mean activation of response units for both tasks on congruent trials
                responseActivations_inc = responseActivations_inc / incongruentCounter; % mean activation of response units for both tasks on incongruent trials
                
                % assign to data set
                task1Idx = find(tasks(1) == taskAccuracies{2}.allTaskCombs(combIdx,:));
                task2Idx = find(tasks(1) == taskAccuracies{2}.allTaskCombs(combIdx,:));
                taskCombData.task1Accuracy_LCA{combIdx} = taskAccuracies{2}.LCA(combIdx,task1Idx);
                taskCombData.task2Accuracy_LCA{combIdx} = taskAccuracies{2}.LCA(combIdx,task2Idx);
                taskCombData.taskAccuracies_LCA{combIdx} = taskAccuracies{2}.LCA(combIdx,:);
                taskCombData.A_B{combIdx} = CDFs{combIdx}.A_B;                                      % CDF for tasks A & B
                taskCombData.A_B_1{combIdx} = CDFs{combIdx}.A_B_1;                                      % lower bound for CDF
                taskCombData.min_A_B{combIdx} = CDFs{combIdx}.min_A_B;                        % upper bound for CDF

                taskCombData.responseActivation_con(combIdx,:) = responseActivations_con;             % mean activation of response units for both tasks on congruent trials
                taskCombData.responseActivation_inc(combIdx,:) = responseActivations_inc;               % mean activation of response units for both tasks on incongruent trials
                taskCombData.similarity(combIdx,:) = [task1Similarity task2Similarity];       % similarity for between dual-tasks & potentially interfering tasks

            end
            
            % assign to data set
            taskCombData.taskDependencies = relevantTaskCombs;                          % info about task combinations and their dependency types
            
end