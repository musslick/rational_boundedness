function [A_B, A_B_1, min_A_B, x, C]  = computeTownsendWengerCDF(taskRT_all_AB, taskRT_all_A, taskRT_all_B, taskA, taskB, optThreshIdx_A, optThreshIdx_B, optThreshIdx_AB, thresholdIdx)

loadLCASettings;

%% compute Townsend capacity C_a(t)
if(isempty(thresholdIdx))
    pickOptimalThreshold = 1;
else
    pickOptimalThreshold = 0;
end

assert(isequal(size(taskRT_all_AB), size(taskRT_all_A)));
assert(isequal(size(taskRT_all_AB), size(taskRT_all_B)));

numSimulations = size(taskRT_all_AB, 4);
maxRT = max([ taskRT_all_AB(:); taskRT_all_A(:); taskRT_all_B(:)]);
maxRT = maxRT(1);
dt = 0.01;

numPatterns = size(taskRT_all_AB, 2);

A_B_1 = nan(numPatterns, round(maxRT/dt));
A_B = nan(numPatterns, round(maxRT/dt));
min_A_B = nan(numPatterns, round(maxRT/dt));

for patternIdx = 1:numPatterns
    
    for t = dt:dt:maxRT

        t_index = round(t/dt);

        if(pickOptimalThreshold)
            thresholdIdx = optThreshIdx_A(patternIdx,taskA);
        end
        nonNaN_taskA = squeeze(~isnan(taskRT_all_A(thresholdIdx(1), patternIdx, taskA, :)));
        F_A = sum(taskRT_all_A(thresholdIdx(1), patternIdx, taskA, nonNaN_taskA) < t) / sum(nonNaN_taskA);
        if(pickOptimalThreshold)
            thresholdIdx = [nan optThreshIdx_B(patternIdx,taskB)];
        end
        nonNaN_taskB = squeeze(~isnan(taskRT_all_B(thresholdIdx(2), patternIdx, taskB, :)));
        F_B = sum(taskRT_all_B(thresholdIdx(2), patternIdx, taskB, nonNaN_taskB) < t) / sum(nonNaN_taskB);
        if(pickOptimalThreshold)
            thresholdIdx1 = optThreshIdx_AB(patternIdx,taskA);
            thresholdIdx2 = optThreshIdx_AB(patternIdx,taskB);
        else
            thresholdIdx1 = thresholdIdx(1);
            thresholdIdx2 = thresholdIdx(2);
        end
        
        AND_thresholdIdx = max([thresholdIdx1 thresholdIdx2]);
        
        nonNaN_taskA = squeeze(~isnan(taskRT_all_AB(AND_thresholdIdx, patternIdx, taskA, :)));
        nonNaN_taskB = squeeze(~isnan(taskRT_all_AB(AND_thresholdIdx, patternIdx, taskB, :)));
        nonNaN_both = nonNaN_taskA & nonNaN_taskB;
        F_AB = sum(taskRT_all_AB(AND_thresholdIdx, patternIdx, taskA, nonNaN_both') < t & taskRT_all_AB(AND_thresholdIdx, patternIdx, taskB, nonNaN_both') < t) / sum(nonNaN_both);

        A_B_1(patternIdx, t_index) = F_A + F_B -1;
        A_B(patternIdx, t_index) = F_AB;
        min_A_B(patternIdx, t_index) = min([F_A, F_B]);
        
        K_A = log(F_A);
        K_B = log(F_B);
        K_AB = log(F_AB);
        C(patternIdx, t_index) = (K_A + K_B) / K_AB;

    end

end

%%

% pat = 6;
% figure(1);
% plot(min_A_B(pat,:), '-k'); hold on;
% plot(A_B_1(pat,:), '-r');
% plot(A_B(pat,:), '-b'); hold off;

% figure(2); 
% plot(C(pat,:));

%%
% figure(1);
% plot(mean(min_A_B), '-k'); hold on;
% plot(mean(A_B_1), '-r');
% plot(mean(A_B), '-b'); hold off;
% 
% nanmean(optThreshIdx_AB)
% figure(2); 
% plot(mean(C));
%%

A_B = nanmean(A_B, 1);
A_B_1 = nanmean(A_B_1, 1);
min_A_B = nanmean(min_A_B, 1);
C = nanmean(C, 1);
x = (dt:dt:maxRT) + LCA_settings.T0;


end

