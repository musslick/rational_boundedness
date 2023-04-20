function plotPostSwitchDynamics(taskSwitchData, field)
    
    if(isfield(taskSwitchData, field))
            
        taskSwitchField = eval(strcat('taskSwitchData.', field));
                
        x = 1:50;

        type3Idx = find(taskSwitchData.taskDependencies(3, :) == 3);
        type0Idx = find(taskSwitchData.taskDependencies(3, :) == 0);
        type4Idx = find(taskSwitchData.taskDependencies(3, :) == 4);

        MSE_1_t3 = mean(taskSwitchField.MSE_1(type3Idx, :));
        MSE_1_t0 = mean(taskSwitchField.MSE_1(type0Idx, :));
        MSE_1_t4 = mean(taskSwitchField.MSE_1(type4Idx, :));
        MSE_1_t3_sem = std(taskSwitchField.MSE_1(type3Idx, :)) / sqrt(length(type3Idx));
        MSE_1_t0_sem = std(taskSwitchField.MSE_1(type0Idx, :)) / sqrt(length(type0Idx));
        MSE_1_t4_sem = std(taskSwitchField.MSE_1(type4Idx, :)) / sqrt(length(type4Idx));

        accuracy_1_t3 = mean(taskSwitchField.accuracy_1(type3Idx, :));
        accuracy_1_t0 = mean(taskSwitchField.accuracy_1(type0Idx, :));
        accuracy_1_t4 = mean(taskSwitchField.accuracy_1(type4Idx, :));
        accuracy_1_t3_sem = std(taskSwitchField.accuracy_1(type3Idx, :)) / sqrt(length(type3Idx));
        accuracy_1_t0_sem = std(taskSwitchField.accuracy_1(type0Idx, :)) / sqrt(length(type0Idx));
        accuracy_1_t4_sem = std(taskSwitchField.accuracy_1(type4Idx, :)) / sqrt(length(type4Idx));

        pcorrect_1_t3 = mean(taskSwitchField.pcorrect_1(type3Idx, :));
        pcorrect_1_t0 = mean(taskSwitchField.pcorrect_1(type0Idx, :));
        pcorrect_1_t4 = mean(taskSwitchField.pcorrect_1(type4Idx, :));
        pcorrect_1_t3_sem = std(taskSwitchField.pcorrect_1(type3Idx, :)) / sqrt(length(type3Idx));
        pcorrect_1_t0_sem = std(taskSwitchField.pcorrect_1(type0Idx, :)) / sqrt(length(type0Idx));
        pcorrect_1_t4_sem = std(taskSwitchField.pcorrect_1(type4Idx, :)) / sqrt(length(type4Idx));
        
        correctResponseAct_1_t3 = mean(taskSwitchField.correctResponseAct_1(type3Idx, :));
        correctResponseAct_1_t0 = mean(taskSwitchField.correctResponseAct_1(type0Idx, :));
        correctResponseAct_1_t4 = mean(taskSwitchField.correctResponseAct_1(type4Idx, :));
        correctResponseAct_1_t3_sem = std(taskSwitchField.correctResponseAct_1(type3Idx, :)) / sqrt(length(type3Idx));
        correctResponseAct_1_t0_sem = std(taskSwitchField.correctResponseAct_1(type0Idx, :)) / sqrt(length(type0Idx));
        correctResponseAct_1_t4_sem = std(taskSwitchField.correctResponseAct_1(type4Idx, :)) / sqrt(length(type4Idx));
        
        interferenceResponseAct_1_t3 = mean(taskSwitchField.interferenceResponseAct_1(type3Idx, :));
        interferenceResponseAct_1_t0 = mean(taskSwitchField.interferenceResponseAct_1(type0Idx, :));
        interferenceResponseAct_1_t4 = mean(taskSwitchField.interferenceResponseAct_1(type4Idx, :));
        interferenceResponseAct_1_t3_sem = std(taskSwitchField.interferenceResponseAct_1(type3Idx, :)) / sqrt(length(type3Idx));
        interferenceResponseAct_1_t0_sem = std(taskSwitchField.interferenceResponseAct_1(type0Idx, :)) / sqrt(length(type0Idx));
        interferenceResponseAct_1_t4_sem = std(taskSwitchField.interferenceResponseAct_1(type4Idx, :)) / sqrt(length(type4Idx));
        
        correctResponseActSum_1_t3 = mean(taskSwitchField.correctResponseActSum_1(type3Idx, :));
        correctResponseActSum_1_t0 = mean(taskSwitchField.correctResponseActSum_1(type0Idx, :));
        correctResponseActSum_1_t4 = mean(taskSwitchField.correctResponseActSum_1(type4Idx, :));
        correctResponseActSum_1_t3_sem = std(taskSwitchField.correctResponseActSum_1(type3Idx, :)) / sqrt(length(type3Idx));
        correctResponseActSum_1_t0_sem = std(taskSwitchField.correctResponseActSum_1(type0Idx, :)) / sqrt(length(type0Idx));
        correctResponseActSum_1_t4_sem = std(taskSwitchField.correctResponseActSum_1(type4Idx, :)) / sqrt(length(type4Idx));
        
        interferenceResponseActSum_1_t3 = mean(taskSwitchField.interferenceResponseActSum_1(type3Idx, :));
        interferenceResponseActSum_1_t0 = mean(taskSwitchField.interferenceResponseActSum_1(type0Idx, :));
        interferenceResponseActSum_1_t4 = mean(taskSwitchField.interferenceResponseActSum_1(type4Idx, :));
        interferenceResponseActSum_1_t3_sem = std(taskSwitchField.interferenceResponseActSum_1(type3Idx, :)) / sqrt(length(type3Idx));
        interferenceResponseActSum_1_t0_sem = std(taskSwitchField.interferenceResponseActSum_1(type0Idx, :)) / sqrt(length(type0Idx));
        interferenceResponseActSum_1_t4_sem = std(taskSwitchField.interferenceResponseActSum_1(type4Idx, :)) / sqrt(length(type4Idx));
        
        set(gcf, 'Position', [100 100 900 600]);
        subplot(3,3,1);
        errorbar(x, MSE_1_t3(x), MSE_1_t3_sem(x), 'r'); hold on;
        errorbar(x, MSE_1_t0(x), MSE_1_t0_sem(x), 'b'); hold on;
        errorbar(x, MSE_1_t4(x), MSE_1_t4_sem(x), 'g'); hold off;
        title('MSE (task 1)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,2);
        errorbar(x, accuracy_1_t3(x), accuracy_1_t3_sem(x), 'r'); hold on;
        errorbar(x, accuracy_1_t0(x), accuracy_1_t0_sem(x), 'b'); hold on;
        errorbar(x, accuracy_1_t4(x), accuracy_1_t4_sem(x), 'g'); hold off;
        title('Accuracy (task 1)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,3);
        errorbar(x, pcorrect_1_t3(x), pcorrect_1_t3_sem(x), 'r'); hold on;
        errorbar(x, pcorrect_1_t0(x), pcorrect_1_t0_sem(x), 'b'); hold on;
        errorbar(x, pcorrect_1_t4(x), pcorrect_1_t4_sem(x), 'g'); hold off;
        title('P(correct|task 1)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,4);
        errorbar(x, correctResponseAct_1_t3(x), correctResponseAct_1_t3_sem(x), 'r'); hold on;
        errorbar(x, correctResponseAct_1_t0(x), correctResponseAct_1_t0_sem(x), 'b'); hold on;
        errorbar(x, correctResponseAct_1_t4(x), correctResponseAct_1_t4_sem(x), 'g'); hold off;
        title('Activation of Correct Response of Task 1');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,5);
        errorbar(x, interferenceResponseAct_1_t3(x), interferenceResponseAct_1_t3_sem(x), 'r'); hold on;
        errorbar(x, interferenceResponseAct_1_t0(x), interferenceResponseAct_1_t0_sem(x), 'b'); hold on;
        errorbar(x, interferenceResponseAct_1_t4(x), interferenceResponseAct_1_t4_sem(x), 'g'); hold off;
        title('Activation of Incorrect Response of Task 1');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,7);
        errorbar(x, correctResponseActSum_1_t3(x), correctResponseActSum_1_t3_sem(x), 'r'); hold on;
        errorbar(x, correctResponseActSum_1_t0(x), correctResponseActSum_1_t0_sem(x), 'b'); hold on;
        errorbar(x, correctResponseActSum_1_t4(x), correctResponseActSum_1_t4_sem(x), 'g'); hold off;
        title('Cumulative Activation of Correct Resp of Task 1');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,8);
        errorbar(x, interferenceResponseActSum_1_t3(x), interferenceResponseActSum_1_t3_sem(x), 'r'); hold on;
        errorbar(x, interferenceResponseActSum_1_t0(x), interferenceResponseActSum_1_t0_sem(x), 'b'); hold on;
        errorbar(x, interferenceResponseActSum_1_t4(x), interferenceResponseActSum_1_t4_sem(x), 'g'); hold off;
        title('Cumulative Activation of Incorrect Resp of Task 1');
        legend('type 3', 'type 0', 'type 4');
    end

end