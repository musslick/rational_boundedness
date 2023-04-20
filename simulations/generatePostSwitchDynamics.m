function [MSE_2, MSE_1, accuracy_2, accuracy_1, pcorrect_2, pcorrect_1, correctResponseAct_2, correctResponseAct_1, interferenceResponseAct_2, interferenceResponseAct_1, correctResponseActSum_2, correctResponseActSum_1, interferenceResponseActSum_2, interferenceResponseActSum_1] = generatePostSwitchDynamics(task1, task2, taskNet, NPathways, NFeatures, multiCap, tau, iterations_taskSwitching)
        
        task2_pattern = zeros(1, NPathways^2);
        task1_pattern = zeros(1, NPathways^2);
        task2_pattern(task2) = 1;
        task1_pattern(task1) = 1;
        multitask_pattern = task2_pattern + task1_pattern;
        
        % generate full task & input streams
        relevantInputPatterns = multiCap{2}.input(ismember(multiCap{2}.tasks, multitask_pattern, 'rows'),:);
        relevantTrainPatterns = multiCap{2}.train(ismember(multiCap{2}.tasks, multitask_pattern, 'rows'),:);
        
        inputStream = repmat(relevantInputPatterns, 1, 1, iterations_taskSwitching);
        
        taskStream = repmat(task2_pattern, size(relevantInputPatterns,1), 1, iterations_taskSwitching);
        taskStream(:,:,1) = repmat(task1_pattern, size(relevantInputPatterns,1), 1);
        
        trainStream = repmat(relevantTrainPatterns, 1, 1, iterations_taskSwitching);
        
        % remove irrelevant training dimension depending on task
        outputDim_2 = mod(task2-1,taskNet.NPathways)+1;
        inputDim_2 = ceil(task2/taskNet.NPathways); 
        relevantOutputFeatures_2 = (NFeatures*(outputDim_2-1)+1) : (NFeatures*outputDim_2);
        relevantInputFeatures_2 = (NFeatures*(inputDim_2-1)+1) : (NFeatures*inputDim_2);

        outputDim_1 = mod(task1-1,taskNet.NPathways)+1;
        inputDim_1 = ceil(task1/taskNet.NPathways); 
        relevantOutputFeatures_1 = (NFeatures*(outputDim_1-1)+1) : (NFeatures*outputDim_1);
        relevantInputFeatures_1 = (NFeatures*(inputDim_1-1)+1) : (NFeatures*inputDim_1);
        
        % mod 03/23/17:
        %trainStream(:,relevantOutputFeatures_2,1) = 0;
        %trainStream(:,relevantOutputFeatures_1,1:end) = 0; % 2:end
        % mod 06/09/17
        % trainStream(:,relevantOutputFeatures_1,1:end) = 0; % 2:end

        % run task switch simulation
        [output_act_log, ~, ~, MSE_tasks, ~, AE_tasks, ~, ~, PCorrect_tasks] = taskNet.integrateNetInput(tau, inputStream, taskStream, trainStream);
        
        % mod 06/12/17
        % cut off initial condition
        output_act_log = output_act_log(:,2:end,:);
        MSE_tasks = MSE_tasks(:,:,2:end);
        AE_tasks = AE_tasks(:,:, 2:end);
        PCorrect_tasks = PCorrect_tasks(:,:,2:end);
        
        MSE_2 = squeeze(mean(MSE_tasks(:,task2,:),1));
        MSE_1 = squeeze(mean(MSE_tasks(:,task1,:),1));
        
        accuracy_2 = squeeze(mean(AE_tasks(:,task2,:),1));
        accuracy_1 = squeeze(mean(AE_tasks(:,task1,:),1));
        
        pcorrect_2 = squeeze(mean(PCorrect_tasks(:,task2,:),1));
        pcorrect_1 = squeeze(mean(PCorrect_tasks(:,task1,:),1));
        
        % compute activity of correct response unit for tasks A & B
        correctResponseAct_2 = nan(size(inputStream,1), size(inputStream,3));
        correctResponseAct_1 = nan(size(inputStream,1), size(inputStream,3));
        interferenceResponseAct_2 = nan(size(inputStream,1), size(inputStream,3));
        interferenceResponseAct_1 = nan(size(inputStream,1), size(inputStream,3));
        
        for patternIdx = 1:size(inputStream,1)
            
            % track activation values of corrept response units for each task
             % mod 06/09/17 (changed to -1)
            for iter = 1:(size(inputStream,3)-1)
                % determine relevant input features
                correctResponse_2 = relevantOutputFeatures_2(find(inputStream(patternIdx, relevantInputFeatures_2,iter) == 1,1));
                correctResponse_1 = relevantOutputFeatures_1(find(inputStream(patternIdx, relevantInputFeatures_1,iter) == 1,1));
                interferenceResponse_2 = relevantOutputFeatures_2(find(inputStream(patternIdx, relevantInputFeatures_1,iter) == 1,1));
                interferenceResponse_1 = relevantOutputFeatures_1(find(inputStream(patternIdx, relevantInputFeatures_2,iter) == 1,1));
                
                % store relevant response unit activations
                correctResponseAct_2(patternIdx, iter) = output_act_log(patternIdx, iter, correctResponse_2);
                correctResponseAct_1(patternIdx, iter) = output_act_log(patternIdx, iter, correctResponse_1);
                interferenceResponseAct_2(patternIdx, iter) = output_act_log(patternIdx, iter, interferenceResponse_2);
                interferenceResponseAct_1(patternIdx, iter) = output_act_log(patternIdx, iter, interferenceResponse_1);
            end
            
        end
        
        correctResponseAct_2 = mean(correctResponseAct_2);
        correctResponseAct_1 = mean(correctResponseAct_1);
        correctResponseActSum_2 = cumsum(correctResponseAct_2);
        correctResponseActSum_1 = cumsum(correctResponseAct_1);
        interferenceResponseAct_2 = mean(interferenceResponseAct_2);
        interferenceResponseAct_1 = mean(interferenceResponseAct_1);
        interferenceResponseActSum_2 = cumsum(interferenceResponseAct_2);
        interferenceResponseActSum_1 = cumsum(interferenceResponseAct_1);

end

