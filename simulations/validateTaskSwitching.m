function [taskSwitchData] = validateTaskSwitching(taskNet, NPathways, NFeatures, taskDependencies, multiCap, tau, iterations_taskSwitching, onsets_taskSwitching)

    % consider two types of switching for the following type of
    % interference:
    %
    %  O    O
    %  |  \  |
    %  A  \ B
    %  |    \|
    %  O    O   type 3 (indirect, asymmetric)
    %
    % 1) switch to A after execution of B
    % 2) execute A and then then try to activate B at any time
    
    % TODO: validate that for same task A (in comparison, you get different
    % switch dynamics depending on what task it is paired with!
    
    taskSwitchData.taskDependencies = taskDependencies;
    taskSwitchData.onsets = onsets_taskSwitching;
            
    %% 1) switch to A after execution of B
    
    %  for each task pair
    for taskCombIdx = 1:size(taskDependencies,2)
        
        % set up task patterns
        taskA = taskDependencies(1,taskCombIdx);
        taskB = taskDependencies(2,taskCombIdx);
        
        % simulate A after B
        task1 = taskB;
        task2 = taskA;
        [MSE_1, MSE_2, accuracy_1, accuracy_2, pcorrect_1, pcorrect_2, ...
            correctResponseAct_1, correctResponseAct_2, interferenceResponseAct_1, interferenceResponseAct_2, ...
            correctResponseActSum_1, correctResponseActSum_2, interferenceResponseActSum_1, interferenceResponseActSum_2] = generatePostSwitchDynamics(task1, task2, taskNet, NPathways, NFeatures, multiCap, tau, iterations_taskSwitching);
        
        taskSwitchData.A_after_B.MSE_1(taskCombIdx, :) = MSE_1;
        taskSwitchData.A_after_B.MSE_2(taskCombIdx, :) = MSE_2;
        taskSwitchData.A_after_B.accuracy_1(taskCombIdx, :) = accuracy_1;
        taskSwitchData.A_after_B.accuracy_2(taskCombIdx, :) = accuracy_2;
        taskSwitchData.A_after_B.pcorrect_1(taskCombIdx, :) = pcorrect_1;
        taskSwitchData.A_after_B.pcorrect_2(taskCombIdx, :) = pcorrect_2;
        taskSwitchData.A_after_B.correctResponseAct_1(taskCombIdx, :) = correctResponseAct_1;
        taskSwitchData.A_after_B.correctResponseAct_2(taskCombIdx, :) = correctResponseAct_2;
        taskSwitchData.A_after_B.interferenceResponseAct_1(taskCombIdx, :) = interferenceResponseAct_1;
        taskSwitchData.A_after_B.interferenceResponseAct_2(taskCombIdx, :) = interferenceResponseAct_2;
        taskSwitchData.A_after_B.correctResponseActSum_1(taskCombIdx, :) = correctResponseActSum_1;
        taskSwitchData.A_after_B.correctResponseActSum_2(taskCombIdx, :) = correctResponseActSum_2;
        taskSwitchData.A_after_B.interferenceResponseActSum_1(taskCombIdx, :) = interferenceResponseActSum_1;
        taskSwitchData.A_after_B.interferenceResponseActSum_2(taskCombIdx, :) = interferenceResponseActSum_2;
        
        % simulate B after A
        task1 = taskA;
        task2 = taskB;
        [MSE_1, MSE_2, accuracy_1, accuracy_2, pcorrect_1, pcorrect_2, ...
            correctResponseAct_1, correctResponseAct_2, interferenceResponseAct_1, interferenceResponseAct_2, ...
            correctResponseActSum_1, correctResponseActSum_2, interferenceResponseActSum_1, interferenceResponseActSum_2] = generatePostSwitchDynamics(task1, task2, taskNet, NPathways, NFeatures, multiCap, tau, iterations_taskSwitching);
        
        taskSwitchData.B_after_A.MSE_1(taskCombIdx, :) = MSE_1;
        taskSwitchData.B_after_A.MSE_2(taskCombIdx, :) = MSE_2;
        taskSwitchData.B_after_A.accuracy_1(taskCombIdx, :) = accuracy_1;
        taskSwitchData.B_after_A.accuracy_2(taskCombIdx, :) = accuracy_2;
        taskSwitchData.B_after_A.pcorrect_1(taskCombIdx, :) = pcorrect_1;
        taskSwitchData.B_after_A.pcorrect_2(taskCombIdx, :) = pcorrect_2;
        taskSwitchData.B_after_A.correctResponseAct_1(taskCombIdx, :) = correctResponseAct_1;
        taskSwitchData.B_after_A.correctResponseAct_2(taskCombIdx, :) = correctResponseAct_2;
        taskSwitchData.B_after_A.interferenceResponseAct_1(taskCombIdx, :) = interferenceResponseAct_1;
        taskSwitchData.B_after_A.interferenceResponseAct_2(taskCombIdx, :) = interferenceResponseAct_2;
        taskSwitchData.B_after_A.correctResponseActSum_1(taskCombIdx, :) = correctResponseActSum_1;
        taskSwitchData.B_after_A.correctResponseActSum_2(taskCombIdx, :) = correctResponseActSum_2;
        taskSwitchData.B_after_A.interferenceResponseActSum_1(taskCombIdx, :) = interferenceResponseActSum_1;
        taskSwitchData.B_after_A.interferenceResponseActSum_2(taskCombIdx, :) = interferenceResponseActSum_2;
        
        % simulate onset of B after initiation of A
        for onsetIdx = 1:length(onsets_taskSwitching)
            onset = onsets_taskSwitching(onsetIdx);
            task1 = taskA;
            task2 = taskB;
            [MSE_1, MSE_2, accuracy_1, accuracy_2, pcorrect_1, pcorrect_2, ...
                correctResponseAct_1, correctResponseAct_2, interferenceResponseAct_1, interferenceResponseAct_2, ...
                correctResponseActSum_1, correctResponseActSum_2, interferenceResponseActSum_1, interferenceResponseActSum_2] = generateOnsetSwitchDynamics(task1, task2, taskNet, NPathways, NFeatures, multiCap, tau, iterations_taskSwitching, onset);

            taskSwitchData.B_onsetAfter_A(onsetIdx).onset = onset;
            taskSwitchData.B_onsetAfter_A(onsetIdx).MSE_1(taskCombIdx, :) = MSE_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).MSE_2(taskCombIdx, :) = MSE_2;
            taskSwitchData.B_onsetAfter_A(onsetIdx).accuracy_1(taskCombIdx, :) = accuracy_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).accuracy_2(taskCombIdx, :) = accuracy_2;
            taskSwitchData.B_onsetAfter_A(onsetIdx).pcorrect_1(taskCombIdx, :) = pcorrect_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).pcorrect_2(taskCombIdx, :) = pcorrect_2;
            taskSwitchData.B_onsetAfter_A(onsetIdx).correctResponseAct_1(taskCombIdx, :) = correctResponseAct_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).correctResponseAct_2(taskCombIdx, :) = correctResponseAct_2;
            taskSwitchData.B_onsetAfter_A(onsetIdx).interferenceResponseAct_1(taskCombIdx, :) = interferenceResponseAct_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).interferenceResponseAct_2(taskCombIdx, :) = interferenceResponseAct_2;
            taskSwitchData.B_onsetAfter_A(onsetIdx).correctResponseActSum_1(taskCombIdx, :) = correctResponseActSum_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).correctResponseActSum_2(taskCombIdx, :) = correctResponseActSum_2;
            taskSwitchData.B_onsetAfter_A(onsetIdx).interferenceResponseActSum_1(taskCombIdx, :) = interferenceResponseActSum_1;
            taskSwitchData.B_onsetAfter_A(onsetIdx).interferenceResponseActSum_2(taskCombIdx, :) = interferenceResponseActSum_2;
        end
        
        disp([num2str(taskCombIdx) ' / ' num2str(size(taskDependencies,2)) ' switch condition']);
    end
    
    figure(1);
    plotPostSwitchDynamics(taskSwitchData, 'A_after_B');
    figure(2);
    plotPostSwitchDynamics(taskSwitchData, 'B_after_A');
    figure(3);
    plotOnsetSwitchDynamics(taskSwitchData, 'B_onsetAfter_A', 2);
    
end
