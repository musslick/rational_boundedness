function plotOnsetSwitchDynamics_batch(x, iterations_taskSwitching, onsets_taskSwitching, batch_log, onsetIdx, correlation_thresholds, varargin)
    
field = 'B_onsetAfter_A';

numRepetitions =  size(batch_log,2);
numCorrThresholds =  size(batch_log(1,1).taskSwitchData,2);
onset = onsets_taskSwitching(onsetIdx);

plotSEM = 1;

% plot across all correlation thresholds?
if(length(varargin) >= 3)
    performanceMeasure = varargin{1};
    corrThreshIdx = varargin{2};
    plotSettings = varargin{3};
else
    performanceMeasure = 'MSE';
    corrThreshIdx = -1;
    plotSettings = [];
    
    fontSize_title = 14;
    fontSize_gca = 14;
    fontSize_xlabel = 14;
    fontSize_ylabel = 14;
    fontSize_legend = 14;

    fontName = 'Helvetica';

    lineWidth = 3;

colors = [253 120 21; ... % orange
              31 104 172; ... % blue
              44 155 37; ... % green
              0     0   0  ; ... % black
            142 142 142; ... % grey 
            255 255 255; ... % white 
            128 190 207] / 255; % blue 2 
        
cContrast1 = 1;
cContrast2 = 2;
cContrast3 = 3;
cContrast4 = 7;
cSingle = 4;
cWeak = 5;
cWhite = 6;
cBlack = 4;

% for export
plotSettings.fontSize_title = fontSize_title;
plotSettings.fontSize_gca = fontSize_gca;
plotSettings.fontSize_xlabel = fontSize_xlabel;
plotSettings.fontSize_ylabel = fontSize_ylabel;
plotSettings.fontSize_legend = fontSize_legend;
plotSettings.fontName = fontName;
plotSettings.lineWidth = lineWidth;
plotSettings.colors = colors;
plotSettings.cContrast1 = cContrast1;
plotSettings.cContrast2 = cContrast2;
plotSettings.cContrast3 = cContrast3;
plotSettings.cContrast4 = cContrast4;
plotSettings.cSingle = cSingle;
plotSettings.cWeak = cWeak;
plotSettings.cWhite = cWhite;

end

colors = plotSettings.colors;
cContrast1 = plotSettings.cContrast1;
cContrast2 = plotSettings.cContrast2;
cContrast3 = plotSettings.cContrast3;
cContrast4 = plotSettings.cContrast4;
cSingle = plotSettings.cSingle;
cWhite = plotSettings.cWhite;

for corrThreshIdx_tmp = 1:numCorrThresholds
    
    if(isempty(varargin) | corrThreshIdx == corrThreshIdx_tmp)
        
        % initialize data
        
        MSE_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        MSE_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        MSE_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        accuracy_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        accuracy_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        accuracy_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        pcorrect_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        pcorrect_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        pcorrect_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        correctResponseAct_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        correctResponseAct_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        correctResponseAct_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        interferenceResponseAct_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        interferenceResponseAct_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        interferenceResponseAct_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        correctResponseActSum_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        correctResponseActSum_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        correctResponseActSum_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        interferenceResponseActSum_2_t3 = nan(numRepetitions, iterations_taskSwitching);
        interferenceResponseActSum_2_t0 = nan(numRepetitions, iterations_taskSwitching);
        interferenceResponseActSum_2_t4 = nan(numRepetitions, iterations_taskSwitching);

        % collect data from each repetition

        for rep = 1:numRepetitions
            
            taskSwitchData = batch_log(1,rep).taskSwitchData{corrThreshIdx_tmp};

            if(isfield(taskSwitchData, field))

                taskSwitchField = eval(strcat('taskSwitchData.', field));
                onset = taskSwitchField(onsetIdx).onset;
                    
                type3Idx = find(taskSwitchData.taskDependencies(3, :) == 3);
                type0Idx = find(taskSwitchData.taskDependencies(3, :) == 0);
                type4Idx = find(taskSwitchData.taskDependencies(3, :) == 4);

                MSE_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).MSE_2(type3Idx, :));
                MSE_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).MSE_2(type0Idx, :));
                MSE_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).MSE_2(type4Idx, :));

                accuracy_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).accuracy_2(type3Idx, :));
                accuracy_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).accuracy_2(type0Idx, :));
                accuracy_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).accuracy_2(type4Idx, :));

                pcorrect_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).pcorrect_2(type3Idx, :));
                pcorrect_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).pcorrect_2(type0Idx, :));
                pcorrect_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).pcorrect_2(type4Idx, :));

                correctResponseAct_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).correctResponseAct_2(type3Idx, :));
                correctResponseAct_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).correctResponseAct_2(type0Idx, :));
                correctResponseAct_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).correctResponseAct_2(type4Idx, :));

                interferenceResponseAct_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).interferenceResponseAct_2(type3Idx, :));
                interferenceResponseAct_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).interferenceResponseAct_2(type0Idx, :));
                interferenceResponseAct_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).interferenceResponseAct_2(type4Idx, :));
  
                correctResponseActSum_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).correctResponseActSum_2(type3Idx, :));
                correctResponseActSum_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).correctResponseActSum_2(type0Idx, :));
                correctResponseActSum_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).correctResponseActSum_2(type4Idx, :));
 
                interferenceResponseActSum_2_t3(rep,:) = mean(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type3Idx, :));
                interferenceResponseActSum_2_t0(rep,:) = mean(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type0Idx, :));
                interferenceResponseActSum_2_t4(rep,:) = mean(taskSwitchField(onsetIdx).interferenceResponseActSum_2(type4Idx, :));
                
            end
            
        end
        
        % compute mean and sem across replications
        
        MSE_2_t3_mean = nanmean(MSE_2_t3, 1);
        MSE_2_t0_mean = nanmean(MSE_2_t0, 1);
        MSE_2_t4_mean = nanmean(MSE_2_t4, 1);
        MSE_2_t3_sem = nanstd(MSE_2_t3, [], 1) / sqrt(sum(isnan(MSE_2_t3(:,end))));
        MSE_2_t0_sem = nanstd(MSE_2_t0, [], 1) / sqrt(sum(isnan(MSE_2_t0(:,end))));
        MSE_2_t4_sem = nanstd(MSE_2_t4, [], 1) / sqrt(sum(isnan(MSE_2_t4(:,end))));

        accuracy_2_t3_mean = nanmean(accuracy_2_t3, 1);
        accuracy_2_t0_mean = nanmean(accuracy_2_t0, 1);
        accuracy_2_t4_mean = nanmean(accuracy_2_t4, 1);
        accuracy_2_t3_sem = nanstd(accuracy_2_t3, [], 1) / sqrt(sum(isnan(accuracy_2_t3(:,end))));
        accuracy_2_t0_sem = nanstd(accuracy_2_t0, [], 1) / sqrt(sum(isnan(accuracy_2_t0(:,end))));
        accuracy_2_t4_sem = nanstd(accuracy_2_t4, [], 1) / sqrt(sum(isnan(accuracy_2_t4(:,end))));

        pcorrect_2_t3_mean = nanmean(pcorrect_2_t3, 1);
        pcorrect_2_t0_mean = nanmean(pcorrect_2_t0, 1);
        pcorrect_2_t4_mean = nanmean(pcorrect_2_t4, 1);
        pcorrect_2_t3_sem = nanstd(pcorrect_2_t3, [], 1) / sqrt(sum(isnan(pcorrect_2_t3(:,end))));
        pcorrect_2_t0_sem = nanstd(pcorrect_2_t0, [], 1) / sqrt(sum(isnan(pcorrect_2_t0(:,end))));
        pcorrect_2_t4_sem = nanstd(pcorrect_2_t4, [], 1) / sqrt(sum(isnan(pcorrect_2_t4(:,end))));

        correctResponseAct_2_t3_mean = nanmean(correctResponseAct_2_t3, 1);
        correctResponseAct_2_t0_mean = nanmean(correctResponseAct_2_t0, 1);
        correctResponseAct_2_t4_mean = nanmean(correctResponseAct_2_t4, 1);
        correctResponseAct_2_t3_sem = nanstd(correctResponseAct_2_t3, [], 1) / sqrt(sum(isnan(correctResponseAct_2_t3(:,end))));
        correctResponseAct_2_t0_sem = nanstd(correctResponseAct_2_t0, [], 1) / sqrt(sum(isnan(correctResponseAct_2_t0(:,end))));
        correctResponseAct_2_t4_sem = nanstd(correctResponseAct_2_t4, [], 1) / sqrt(sum(isnan(correctResponseAct_2_t4(:,end))));

        interferenceResponseAct_2_t3_mean = nanmean(interferenceResponseAct_2_t3, 1);
        interferenceResponseAct_2_t0_mean = nanmean(interferenceResponseAct_2_t0, 1);
        interferenceResponseAct_2_t4_mean = nanmean(interferenceResponseAct_2_t4, 1);
        interferenceResponseAct_2_t3_sem = nanstd(interferenceResponseAct_2_t3, [], 1) / sqrt(sum(isnan(interferenceResponseAct_2_t3(:,end))));
        interferenceResponseAct_2_t0_sem = nanstd(interferenceResponseAct_2_t0, [], 1) / sqrt(sum(isnan(interferenceResponseAct_2_t0(:,end))));
        interferenceResponseAct_2_t4_sem = nanstd(interferenceResponseAct_2_t4, [], 1) / sqrt(sum(isnan(interferenceResponseAct_2_t4(:,end))));

        correctResponseActSum_2_t3_mean = nanmean(correctResponseActSum_2_t3, 1);
        correctResponseActSum_2_t0_mean = nanmean(correctResponseActSum_2_t0, 1);
        correctResponseActSum_2_t4_mean = nanmean(correctResponseActSum_2_t4, 1);
        correctResponseActSum_2_t3_sem = nanstd(correctResponseActSum_2_t3, [], 1) / sqrt(sum(isnan(correctResponseActSum_2_t3(:,end))));
        correctResponseActSum_2_t0_sem = nanstd(correctResponseActSum_2_t0, [], 1) / sqrt(sum(isnan(correctResponseActSum_2_t0(:,end))));
        correctResponseActSum_2_t4_sem = nanstd(correctResponseActSum_2_t4, [], 1) / sqrt(sum(isnan(correctResponseActSum_2_t4(:,end))));

        interferenceResponseActSum_2_t3_mean = nanmean(interferenceResponseActSum_2_t3, 1);
        interferenceResponseActSum_2_t0_mean = nanmean(interferenceResponseActSum_2_t0, 1);
        interferenceResponseActSum_2_t4_mean = nanmean(interferenceResponseActSum_2_t4, 1);
        interferenceResponseActSum_2_t3_sem = nanstd(interferenceResponseActSum_2_t3, [], 1) / sqrt(sum(isnan(interferenceResponseActSum_2_t3(:,end))));
        interferenceResponseActSum_2_t0_sem = nanstd(interferenceResponseActSum_2_t0, [], 1) / sqrt(sum(isnan(interferenceResponseActSum_2_t0(:,end))));
        interferenceResponseActSum_2_t4_sem = nanstd(interferenceResponseActSum_2_t4, [], 1) / sqrt(sum(isnan(interferenceResponseActSum_2_t4(:,end))));

        figure
        set(gcf, 'Position', [100 100 900 600]);
        subplot(3,3,1);
        errorbar(x, MSE_2_t3_mean(x), MSE_2_t3_sem(x), 'r'); hold on;
        errorbar(x, MSE_2_t0_mean(x), MSE_2_t0_sem(x), 'b'); hold on;
        errorbar(x, MSE_2_t4_mean(x), MSE_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([MSE_2_t3_mean MSE_2_t0_mean MSE_2_t4_mean]) max([MSE_2_t3_mean MSE_2_t0_mean MSE_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('MSE (task 2)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,2);
        errorbar(x, accuracy_2_t3_mean(x), accuracy_2_t3_sem(x), 'r'); hold on;
        errorbar(x, accuracy_2_t0_mean(x), accuracy_2_t0_sem(x), 'b'); hold on;
        errorbar(x, accuracy_2_t4_mean(x), accuracy_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([accuracy_2_t3_mean accuracy_2_t0_mean accuracy_2_t4_mean]) max([accuracy_2_t3_mean accuracy_2_t0_mean accuracy_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Accuracy (task 2)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,3);
        errorbar(x, pcorrect_2_t3_mean(x), pcorrect_2_t3_sem(x), 'r'); hold on;
        errorbar(x, pcorrect_2_t0_mean(x), pcorrect_2_t0_sem(x), 'b'); hold on;
        errorbar(x, pcorrect_2_t4_mean(x), pcorrect_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([pcorrect_2_t3_mean pcorrect_2_t0_mean pcorrect_2_t4_mean]) max([pcorrect_2_t3_mean pcorrect_2_t0_mean pcorrect_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('P(correct|task 2)');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,4);
        errorbar(x, correctResponseAct_2_t3_mean(x), correctResponseAct_2_t3_sem(x), 'r'); hold on;
        errorbar(x, correctResponseAct_2_t0_mean(x), correctResponseAct_2_t0_sem(x), 'b'); hold on;
        errorbar(x, correctResponseAct_2_t4_mean(x), correctResponseAct_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([correctResponseAct_2_t3_mean correctResponseAct_2_t0_mean correctResponseAct_2_t4_mean]) max([correctResponseAct_2_t3_mean correctResponseAct_2_t0_mean correctResponseAct_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Activation of Correct Response of Task 2');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,5);
        errorbar(x, interferenceResponseAct_2_t3_mean(x), interferenceResponseAct_2_t3_sem(x), 'r'); hold on;
        errorbar(x, interferenceResponseAct_2_t0_mean(x), interferenceResponseAct_2_t0_sem(x), 'b'); hold on;
        errorbar(x, interferenceResponseAct_2_t4_mean(x), interferenceResponseAct_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([interferenceResponseAct_2_t3_mean interferenceResponseAct_2_t0_mean interferenceResponseAct_2_t4_mean]) max([interferenceResponseAct_2_t3_mean interferenceResponseAct_2_t0_mean interferenceResponseAct_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Activation of Incorrect Response of Task 2');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,7);
        errorbar(x, correctResponseActSum_2_t3_mean(x), correctResponseActSum_2_t3_sem(x), 'r'); hold on;
        errorbar(x, correctResponseActSum_2_t0_mean(x), correctResponseActSum_2_t0_sem(x), 'b'); hold on;
        errorbar(x, correctResponseActSum_2_t4_mean(x), correctResponseActSum_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([correctResponseActSum_2_t3_mean correctResponseActSum_2_t0_mean correctResponseActSum_2_t4_mean]) max([correctResponseActSum_2_t3_mean correctResponseActSum_2_t0_mean correctResponseActSum_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Cumulative Activation of Correct Resp of Task 2');
        legend('type 3', 'type 0', 'type 4');
        subplot(3,3,8);
        errorbar(x, interferenceResponseActSum_2_t3_mean(x), interferenceResponseActSum_2_t3_sem(x), 'r'); hold on;
        errorbar(x, interferenceResponseActSum_2_t0_mean(x), interferenceResponseActSum_2_t0_sem(x), 'b'); hold on;
        errorbar(x, interferenceResponseActSum_2_t4_mean(x), interferenceResponseActSum_2_t4_sem(x), 'g'); hold on;
        plotLimits = [min([interferenceResponseActSum_2_t3_mean interferenceResponseActSum_2_t0_mean interferenceResponseActSum_2_t4_mean]) max([interferenceResponseActSum_2_t3_mean interferenceResponseActSum_2_t0_mean interferenceResponseActSum_2_t4_mean])];
        plot([onset onset], [plotLimits(1) plotLimits(end)], '--k'); hold off;
        title('Cumulative Activation of Incorrect Resp of Task 2');
        legend('type 3', 'type 0', 'type 4');
        text(.030,-12.0, ['corrThresh = ' num2str(correlation_thresholds(corrThreshIdx_tmp))], 'fontSize', 16);
        
         %% plot for particular correlation threshold
        
        if(corrThreshIdx == corrThreshIdx_tmp)
            
            % check performance measure
            disp(performanceMeasure)
            switch performanceMeasure
                case 'MSE'
                    y_t3_mean = MSE_2_t3_mean;
                    y_t0_mean = MSE_2_t0_mean;
                    y_t4_mean = MSE_2_t4_mean;
                    y_t3_sem = MSE_2_t3_sem;
                    y_t0_sem = MSE_2_t0_sem;
                    y_t4_sem = MSE_2_t4_sem;
                  case 'accuracy'
                    y_t3_mean = accuracy_2_t3_mean;
                    y_t0_mean = accuracy_2_t0_mean;
                    y_t4_mean = accuracy_2_t4_mean;
                    y_t3_sem = accuracy_2_t3_sem;
                    y_t0_sem = accuracy_2_t0_sem;
                    y_t4_sem = accuracy_2_t4_sem;
                  case 'pcorrect'
                    y_t3_mean = pcorrect_2_t3_mean;
                    y_t0_mean = pcorrect_2_t0_mean;
                    y_t4_mean = pcorrect_2_t4_mean;
                    y_t3_sem = pcorrect_2_t3_sem;
                    y_t0_sem = pcorrect_2_t0_sem;
                    y_t4_sem = pcorrect_2_t4_sem;
                  case 'correctResponseAct'
                    y_t3_mean = correctResponseAct_2_t3_mean;
                    y_t0_mean = correctResponseAct_2_t0_mean;
                    y_t4_mean = correctResponseAct_2_t4_mean;
                    y_t3_sem = correctResponseAct_2_t3_sem;
                    y_t0_sem = correctResponseAct_2_t0_sem;
                    y_t4_sem = correctResponseAct_2_t4_sem;
                  case 'interferenceResponseAct'
                    y_t3_mean = interferenceResponseAct_2_t3_mean;
                    y_t0_mean = interferenceResponseAct_2_t0_mean;
                    y_t4_mean = interferenceResponseAct_2_t4_mean;
                    y_t3_sem = interferenceResponseAct_2_t3_sem;
                    y_t0_sem = interferenceResponseAct_2_t0_sem;
                    y_t4_sem = interferenceResponseAct_2_t4_sem;
                  case 'correctResponseActSum'
                    y_t3_mean = correctResponseActSum_2_t3_mean;
                    y_t0_mean = correctResponseActSum_2_t0_mean;
                    y_t4_mean = correctResponseActSum_2_t4_mean;
                    y_t3_sem = correctResponseActSum_2_t3_sem;
                    y_t0_sem = correctResponseActSum_2_t0_sem;
                    y_t4_sem = correctResponseActSum_2_t4_sem;
                  case 'interferenceResponseActSum'
                    y_t3_mean = interferenceResponseActSum_2_t3_mean;
                    y_t0_mean = interferenceResponseActSum_2_t0_mean;
                    y_t4_mean = interferenceResponseActSum_2_t4_mean;
                    y_t3_sem = interferenceResponseActSum_2_t3_sem;
                    y_t0_sem = interferenceResponseActSum_2_t0_sem;
                    y_t4_sem = interferenceResponseActSum_2_t4_sem;
            end
            
            % actual plot
            fig1 = figure;
            set(fig1, 'Position', [470 200 250 220]);
            if(~any(isnan(y_t3_mean(2:end))) && ~any(isnan(y_t0_mean(2:end))) && ~any(isnan(y_t4_mean(2:end)))) 
                if(plotSEM)
                    errorbar(x, y_t0_mean(x), y_t0_sem(x), 'LineWidth', plotSettings.lineWidth*0.5, 'Color', colors(cContrast4,:)); hold on;
                    errorbar(x, y_t3_mean(x), y_t3_sem(x), 'LineWidth', plotSettings.lineWidth*0.5, 'Color', colors(cContrast1,:)); hold on;
%                     errorbar(x, y_t4_mean(x), y_t4_sem(x), 'LineWidth', plotSettings.lineWidth*0.5, 'Color', colors(cContrast3,:)); hold on;
                else
                    plot(x, y_t0_mean(x), 'LineWidth', plotSettings.lineWidth, 'Color', colors(cContrast4,:)); hold on;
                    plot(x, y_t3_mean(x), 'LineWidth', plotSettings.lineWidth, 'Color', colors(cContrast1,:)); hold on;
%                     plot(x, y_t4_mean(x), 'LineWidth', plotSettings.lineWidth, 'Color', colors(cContrast3,:)); hold on;
                end
                switch performanceMeasure
                    case 'MSE'
                        plotLimits = [0 max([y_t0_mean(x) y_t3_mean(x) y_t4_mean(x)])];
                    case 'accuracy'
                        plotLimits = [min([y_t0_mean(x) y_t3_mean(x) y_t4_mean(x)])  1];
                end
                plot([onset onset], [plotLimits(1) plotLimits(end)], '--', 'Color', colors(cWhite,:), 'LineWidth', plotSettings.lineWidth); hold off;
                performanceLabel = performanceMeasure;
                performanceLabel(1) = upper(performanceLabel(1));
                if(strcmp(performanceLabel, 'MSE'))
                    performanceLabel = 'Error'
                end
                ylabel([performanceLabel ' Of a'], 'FontSize', plotSettings.fontSize_ylabel, 'FontName', plotSettings.fontName, 'Color', 'w');
                xlabel('Time Steps','FontSize', plotSettings.fontSize_xlabel, 'Color', 'w');
                field_label = 'Effect Of B Onset on A';
                
%                 title({[field_label ', corrThresh=' num2str(correlation_thresholds(corrThreshIdx)) ]},'FontSize', plotSettings.fontSize_title, 'FontName', plotSettings.fontName, , 'Color', 'w');
                xlim([0 x(end)]);
                ylim([plotLimits(1) plotLimits(end)]);
%                 l = legend('type 0','type 3','type 4','Onset Of Task B','Location', 'northeast');
%                 set(l, 'FontName', plotSettings.fontName, 'fontSize', plotSettings.fontSize_legend-1, 'Color', 'w');
                set(gca, 'FontSize', plotSettings.fontSize_gca);
                
                set(gcf, 'Color', 'k');
                set(gca, 'Color', 'k');
                set(gca, 'xColor', 'w');
                set(gca, 'yColor', 'w');
                set(gca, 'zColor', 'w');
            end
        end
    
    end
    
end

end