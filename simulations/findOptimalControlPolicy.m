function [maximumPerformanceCurve, taskAccuracies, CDFs] = findOptimalControlPolicy(taskNet, NPathways, NFeatures, multiCap)

  CDFs = {};
  
    for cap = 1:length(multiCap)
        
        % find all available tasks
        allTaskCombs = multiCap{cap}.taskCombs;
        
        if(~isempty(allTaskCombs))
            
            taskAccuracies{cap}.LCA = zeros(size(allTaskCombs, 1), cap);
            taskAccuracies{cap}.allTaskCombs = allTaskCombs;

            for combIdx = 1:size(allTaskCombs, 1)

                % find corresponding test patterns
                patternIdx = multiCap{cap}.taskIdx == combIdx;
                input = multiCap{cap}.input(patternIdx, :);
                tasks = multiCap{cap}.tasks(patternIdx, :);
                train = multiCap{cap}.train(patternIdx, :);

                % Get the encoding for the task combination
                currentTaskComb = tasks(1,:);

                % check task cardinality
                if(cap ~= sum(currentTaskComb))
                    warning(['Number of tasks does not match assigned cardinality. Cardinality is ' num2str(cap) ' and number of tasks is ' num2str(sum(currenTasks))]);
                end

                % compute output for all task patterns of a task
                [outputPatterns] = taskNet.runSet(input, tasks, train);

                %% COMPUTE OPTIMIZATION CRITERION

                % identify relevant output dimension
                taskM = reshape(currentTaskComb, NPathways, NPathways);
                [relevantOutputDims relevantInputDim] = find(taskM == 1);

                % test if two tasks afford a response at the same output dimension
                if(length(unique(relevantOutputDims)) ~= length(relevantOutputDims))
                    warning('Tested multitasking pair affords response at the same output dimension.');
                end

                %% optimal LCA accuracy as performance criterion

                % COMPUTE LCA ACCURACY

                % LCA settings
                loadLCASettings;
                    
                % LCA call
                tic
                [optTaskAccuracy, ~, ~, ~, optThreshIdx, ~, ~, ~, taskRT_all] ...
                = taskNet.runLCA(LCA_settings, input, tasks, train);
                toc

                disp(['max thresh: ' num2str(max(max(optThreshIdx)))]);

                taskIDs = find(currentTaskComb == 1);

                taskAccuracies{cap}.LCA(combIdx, :) = nanmean(optTaskAccuracy(:, taskIDs),1);
                
                % perform Townsend analysis
                
                if(cap == 2)
                    taskA = allTaskCombs(combIdx,1);
                    taskB = allTaskCombs(combIdx,2);
                    
                    tasksIdxSgl = multiCap{1}.tasksIdxSgl;
                    input_A = multiCap{1}.input(tasksIdxSgl == taskA,:);
                    tasks_A = multiCap{1}.tasks(tasksIdxSgl == taskA,:);
                    train_A = multiCap{1}.train(tasksIdxSgl == taskA,:);

                    input_B = multiCap{1}.input(tasksIdxSgl == taskB,:);
                    tasks_B = multiCap{1}.tasks(tasksIdxSgl == taskB,:);
                    train_B = multiCap{1}.train(tasksIdxSgl == taskB,:);

                    % Multi LCA call
                    [optTaskAccuracy, ~, ~, ~, optThreshIdx_AB, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,~,taskRT_correct_AB_all] ...
                    = taskNet.runLCA(LCA_settings, input, tasks, train);

                    disp(['max thresh AB: ' num2str(max(max(optThreshIdx_AB)))]);

                    % Single task calls
                    [~, ~, ~, ~, optThreshIdx_A, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,~,taskRT_correct_A_all] ...
                    = taskNet.runLCA(LCA_settings, input_A, tasks_A, train_A);

                    disp(['max thresh A: ' num2str(max(max(optThreshIdx_A)))]);

                    [~, ~, ~, ~, optThreshIdx_B, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,~,taskRT_correct_B_all] ...
                    = taskNet.runLCA(LCA_settings, input_B, tasks_B, train_B);

                    disp(['max thresh B: ' num2str(max(max(optThreshIdx_B)))]);
                    
                    thresholdIdx = round(nanmean([optThreshIdx_A(:) optThreshIdx_B(:)]));
                    
                    [A_B, A_B_1, min_A_B, x, C]  = computeTownsendWengerCDF(taskRT_correct_AB_all, taskRT_correct_A_all, taskRT_correct_B_all, taskA, taskB, optThreshIdx_A, optThreshIdx_B, optThreshIdx_AB, thresholdIdx);
                    CDFs{combIdx}.A_B = A_B;
                    CDFs{combIdx}.A_B_1 = A_B_1;
                    CDFs{combIdx}.min_A_B = min_A_B;
                    CDFs{combIdx}.x = x;
                    CDFs{combIdx}.C = C;
                end

                disp(['cap: ' num2str(cap) ', comb ' num2str(combIdx) '/' num2str(size(allTaskCombs, 1))]);
            end

        end
    end
    
    %% determine maximum 
    
    maximumPerformanceCurve.LCA = nan(1, length(multiCap));
    
    % for all capacities, look for maximum performance of the worst task in a set 
    for cap = 1:length(taskAccuracies)

            maximumPerformanceCurve.LCA(cap) = max(min(taskAccuracies{cap}.LCA, [], 2));

    end

    
end