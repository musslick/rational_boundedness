function plotOnsetSwitchDynamics(taskSwitchData, field, onsetIdx)
    
    if(isfield(taskSwitchData, field))
            
        taskSwitchField = eval(strcat('taskSwitchData.', field));
        onset = taskSwitchField(onsetIdx).onset;
                
        x = 1:50;

        type3Idx = find(taskSwitchData.taskDependencies(3, :) == 3);
        type0Idx = find(taskSwitchData.taskDependencies(3, :) == 0);
        type4Idx = find(taskSwitchData.taskDependencies(3, :) == 4);

        MSE_2_t3 = mean(taskSwitchField(onsetIdx).MSE_2(type3Idx, :));
        MSE_2_t0 = mean(taskSwitchField(onsetIdx).MSE_2(type0Idx, :));
        MSE_2_t4 = mean(taskSwitchField(onsetIdx).MSE_2(type4Idx, :));
        MSE_2_t3_sem = std(taskSwitchField(onsetIdx).MSE_2(type3Idx, :)) / sqrt(length(type3Idx));
        MSE_2_t0_sem = std(taskSwitchField(onsetIdx).MSE_2(type0Idx, :)) / sqrt(length(type0Idx));
        MSE_2_t4_sem = std(taskSwitchField(onsetIdx).MSE_2(type4Idx, :)) / sqrt(length(type4Idx));

        accuracy_2_t3 = mean(taskSwitchField(onsetIdx).accuracy_2(type3Idx, :));
        accuracy_2_t0 = mean(taskSwitchField(onsetIdx).accuracy_2(type0Idx, :));
        accuracy_2_t4 = mean(taskSwitchField(onsetIdx).accuracy_2(type4Idx, :));
        accuracy_2_t3_sem = std(taskSwitchField(onsetIdx).accuracy_2(type3Idx, :)) / sqrt(length(type3Idx));
        accuracy_2_t0_sem = std(taskSwitchField(onsetIdx).accuracy_2(type0Idx, :)) / sqrt(length(type0Idx));
        accuracy_2_t4_sem = std(taskSwitchField(onsetIdx).accuracy_2(type4Idx, :)) / sqrt(length(type4Idx));

        pcorrect_2_t3 = mean(taskSwitchField(onsetIdx).pcorrect_2(type3Idx, :));
        pcorrect_2_t0 = mean(taskSwitchField(onsetIdx).pcorrect_2(type0Idx, :));
        pcorrect_2_t4 = mean(taskSwitchField(onsetIdx).pcorrect_2(type4Idx, :));
        pcorrect_2_t3_sem = std(taskSwitchField(onsetIdx).pcorrect_2(type3Idx, :)) / sqrt(length(type3Idx));
        pcorrect_2_t0_sem = std(taskSwitchField(onsetIdx).pcorrect_2(type0Idx, :)) / sqrt(length(type0Idx));
        pcorrect_2_t4_sem = std(taskSwitchField(onsetIdx).pcorrect_2(type4Idx, :)) / sqrt(length(type4Idx));
        
        correctResponseAct_2_t3 = mean(taskSwitchField(onsetIdx).correctResponseAct_2(type3Idx, :));
        correctResponseAct_2_t0 = mean(taskSwitchField(onsetIdx).correctResponseAct_2(type0Idx, :));
        correctResponseAct_2_t4 = mean(taskSwitchField(onsetIdx).correctResponseAct_2(type4Idx, :));
        correctResponseAct_2_t3_sem = std(taskSwitchField(onsetIdx).correctResponseAct_2(type3Idx, :)) / sqrt(length(type3Idx));
        correctResponseAct_2_t0_sem = std(taskSwitchField(onsetIdx).correctResponseAct_2(type0Idx, :)) / sqrt(length(type0Idx));
        correctResponseAct_2_t4_sem = std(taskSwitchField(onsetIdx).correctResponseAct_2(type4Idx, :)) / sqrt(length(type4Idx));
        
        interferenceResponseAct_2_t3 = mean(taskSwitchField(onsetIdx).interferenceResponseAct_2(type3Idx, :));
        interferenceResponseAct_2_t0 = mean(taskSwitchField(onsetIdx).interferenceResponseAct_2(type0Idx, :));
        interferenceResponseAct_2_t4 = mean(taskSwitchField(onsetIdx).interferenceResponseAct_2(type4Idx, :));
        interferenceResponseAct_2_t3_sem = std(taskSwitchField(onsetIdx).interferenceResponseAct_2(type3Idx, :)) / sqrt(length(type3Idx));
        interferenceResponseAct_2_t0_sem = std(taskSwitchField(onsetIdx).interferenceResponseAct_2(type0Idx, :)) / sqrt(length(type0Idx));
        interferenceResponseAct_2_t4_sem = std(taskSwitchField(onsetIdx).interferenceResponseAct_2(type4Idx, :)) / sqrt(length(type4Idx));
        
        correctResponseActSum_2_t3 = mean(taskSwitchField(onsetIdx).correctResponseActSum_2(type3Idx, :));
        correctResponseActSum_2_t0 = mean(taskSwitchField(onsetIdx).correctResponseActSum_2(type0Idx, :));
        correctResponseActSum_2_t4 = mean(taskSwitchField(onsetIdx).correctResponseActSum_2(type4Idx, :));
        correctResponseActSum_2_t3_sem = std(taskSwitchField(onsetIdx).correctResponseActSum_2(type3Idx, :)) / sqrt(length(type3Idx));
        correctResponseActSum_2_t0_sem = std(taskSwitchField(onsetIdx).correctResponseActSum_2(type0Idx, :)) / sqrt(length(type0Idx));
        correctResponseActSum_2_t4_sem = std(taskSwitchField(onsetIdx).correctResponseActSum_2(type4Idx, :)) / sqrt(length(type4Idx));
        
        interferenceResponseActSum_2_t3 = mean(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type3Idx, :));
        interferenceResponseActSum_2_t0 = mean(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type0Idx, :));
        interferenceResponseActSum_2_t4 = mean(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type4Idx, :));
        interferenceResponseActSum_2_t3_sem = std(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type3Idx, :)) / sqrt(length(type3Idx));
        interferenceResponseActSum_2_t0_sem = std(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type0Idx, :)) / sqrt(length(type0Idx));
        interferenceResponseActSum_2_t4_sem = std(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type4Idx, :)) / sqrt(length(type4Idx));
        
        set(gcf, 'Position', [100 100 900 600]);
        subplot(3,3,1);
        errorbar(x, MSE_2_t3(x), MSE_2_t3_sem(x), 'r'); hold on;
        errorbar(x, MSE_2_t0(x), MSE_2_t0_sem(x), 'b'); hold on;
        errorbar(x, MSE_2_t4(x), MSE_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([MSE_2_t3 MSE_2_t0 MSE_2_t4]) max([MSE_2_t3 MSE_2_t0 MSE_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('MSE (task 2)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,2);
        errorbar(x, accuracy_2_t3(x), accuracy_2_t3_sem(x), 'r'); hold on;
        errorbar(x, accuracy_2_t0(x), accuracy_2_t0_sem(x), 'b'); hold on;
        errorbar(x, accuracy_2_t4(x), accuracy_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([accuracy_2_t3 accuracy_2_t0 accuracy_2_t4]) max([accuracy_2_t3 accuracy_2_t0 accuracy_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Accuracy (task 2)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,3);
        errorbar(x, pcorrect_2_t3(x), pcorrect_2_t3_sem(x), 'r'); hold on;
        errorbar(x, pcorrect_2_t0(x), pcorrect_2_t0_sem(x), 'b'); hold on;
        errorbar(x, pcorrect_2_t4(x), pcorrect_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([pcorrect_2_t3 pcorrect_2_t0 pcorrect_2_t4]) max([pcorrect_2_t3 pcorrect_2_t0 pcorrect_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('P(correct|task 2)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,4);
        errorbar(x, correctResponseAct_2_t3(x), correctResponseAct_2_t3_sem(x), 'r'); hold on;
        errorbar(x, correctResponseAct_2_t0(x), correctResponseAct_2_t0_sem(x), 'b'); hold on;
        errorbar(x, correctResponseAct_2_t4(x), correctResponseAct_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([correctResponseAct_2_t3 correctResponseAct_2_t0 correctResponseAct_2_t4]) max([correctResponseAct_2_t3 correctResponseAct_2_t0 correctResponseAct_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Activation of Correct Response of Task 2');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,5);
        errorbar(x, interferenceResponseAct_2_t3(x), interferenceResponseAct_2_t3_sem(x), 'r'); hold on;
        errorbar(x, interferenceResponseAct_2_t0(x), interferenceResponseAct_2_t0_sem(x), 'b'); hold on;
        errorbar(x, interferenceResponseAct_2_t4(x), interferenceResponseAct_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([interferenceResponseAct_2_t3 interferenceResponseAct_2_t0 interferenceResponseAct_2_t4]) max([interferenceResponseAct_2_t3 interferenceResponseAct_2_t0 interferenceResponseAct_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Activation of Incorrect Response of Task 2');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,7);
        errorbar(x, correctResponseActSum_2_t3(x), correctResponseActSum_2_t3_sem(x), 'r'); hold on;
        errorbar(x, correctResponseActSum_2_t0(x), correctResponseActSum_2_t0_sem(x), 'b'); hold on;
        errorbar(x, correctResponseActSum_2_t4(x), correctResponseActSum_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([correctResponseActSum_2_t3 correctResponseActSum_2_t0 correctResponseActSum_2_t4]) max([correctResponseActSum_2_t3 correctResponseActSum_2_t0 correctResponseActSum_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Cumulative Activation of Correct Resp of Task 2');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,8);
        errorbar(x, interferenceResponseActSum_2_t3(x), interferenceResponseActSum_2_t3_sem(x), 'r'); hold on;
        errorbar(x, interferenceResponseActSum_2_t0(x), interferenceResponseActSum_2_t0_sem(x), 'b'); hold on;
        errorbar(x, interferenceResponseActSum_2_t4(x), interferenceResponseActSum_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([interferenceResponseActSum_2_t3 interferenceResponseActSum_2_t0 interferenceResponseActSum_2_t4]) max([interferenceResponseActSum_2_t3 interferenceResponseActSum_2_t0 interferenceResponseActSum_2_t4])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Cumulative Activation of Incorrect Resp of Task 2');
        legend('type 3', 'type 0', 'type 4');
    end

end