% author: Sebastian Musslick

% multiPerformance_mean        % first row indicates capacity
                                              % second row indicates good MSE performance
                                              % third row indicates bad MSE performance

function [pathwayCapacities, maxCarryingCapacity, BK_MIS, A_bipartite, A_tasksIdx, multiPerformance_mean, multiPerformance_sem, MSEdata, A_dual] = validateMIS(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, multiCap, varargin)

    NPathways = taskNet.NPathways;
    
    % subsampling
    if~isempty(varargin)
        subSample = varargin{1};
    else
        subSample = 0;
    end

    % find MIS
    [pathwayCapacities, maxCarryingCapacity,  BK_MIS, A_bipartite, A_tasksIdx, A_dual] = getMaxCarryingCapacity(R_hidden, R_output, corr_threshold);
    
    % formatting outputs
    BK_MIS = extend_BK_MIS(BK_MIS);
    % remove duplicates
    
    A_tasksIdx(A_tasksIdx > 0) = tasksToPerform(A_tasksIdx(A_tasksIdx > 0));
    A_tasksOrder = A_tasksIdx(:);
    A_tasksOrder(A_tasksOrder == 0) = [];
    pathwayCapacities = [A_tasksOrder pathwayCapacities];

    % extract all possible multitasking conditions
    good_multiTaskConditions = zeros(size(BK_MIS,2), NPathways*NPathways); % holds task input patterns for good multitasking conditions (that constitute an independent set)
    good_maximumCapacity = sum(BK_MIS);

    for multiCase = 1:size(BK_MIS,2)

        multiTaskIdx = pathwayCapacities(BK_MIS(:,multiCase) == 1,1); % extract which tasks can be performed at the same time

        good_multiTaskConditions(multiCase,multiTaskIdx) = 1; % activate task input pattern units for corresponding set of tasks

    end
    

    % extract all capacities (larger than just single task)
    allCapacities = unique(good_maximumCapacity);
    allCapacities(allCapacities == 1) = [];
    allCapacities(allCapacities == 0) = [];

    % for each multitasking capacity, calculate performance for
    % multitaskable tasks vs. non-multitaskable tasks
    multiPerformance_mean = zeros(3, length(allCapacities));        % first row indicates capacity
                                                                                               % second row indicates good MSE performance
                                                                                               % third row indicates bad MSE performance
    multiPerformance_sem = zeros(size(multiPerformance_mean));
    
    multiPerformance_mean(1,:) = allCapacities;
    multiPerformance_sem(1,:) = allCapacities;
    
    MSEdata = {};

    for capIdx = 1:length(allCapacities)

        cap = allCapacities(capIdx);
        
        input_MultiCap = multiCap{cap}.input;
        tasks_MultiCap = multiCap{cap}.tasks;
        train_MultiCap = multiCap{cap}.train;
        
        % extract task inputs for good multitasking conditions
        good_multiTaskConditions_cap = good_multiTaskConditions(good_maximumCapacity == cap,:);
        removeInvalidIdx = [];
        for i = 1:size(good_multiTaskConditions_cap,1)
            all_taskConditionsMatrix = reshape(good_multiTaskConditions_cap(i,:),NPathways, NPathways);
            % select only allowed tasks
            currentTaskIDs = find(good_multiTaskConditions_cap(i,:));
            if(mean(ismember(currentTaskIDs, tasksToPerform)) ~= 1)
                removeInvalidIdx = [removeInvalidIdx i];
            else
                % select only type 3 interference
                if (max(sum(all_taskConditionsMatrix,1)) > 1 || max(sum(all_taskConditionsMatrix,2)) > 1)
                    removeInvalidIdx = [removeInvalidIdx i];
                end
            end
        end
        good_multiTaskConditions_cap(removeInvalidIdx,:) = [];
        
        % extract all valid multitasking conditions
        all_multiTaskConditions_cap = unique(tasks_MultiCap, 'rows');
        removeInvalidIdx = [];
        for i = 1:size(all_multiTaskConditions_cap,1)
            all_taskConditionsMatrix = reshape(all_multiTaskConditions_cap(i,:),NPathways, NPathways);
            % select only allowed tasks
            currentTaskIDs = find(all_multiTaskConditions_cap(i,:));
            if(mean(ismember(currentTaskIDs, tasksToPerform)) ~= 1)
                removeInvalidIdx = [removeInvalidIdx i];
            else
                % select only type 3 interference
                if (max(sum(all_taskConditionsMatrix,1)) > 1 || max(sum(all_taskConditionsMatrix,2)) > 1)
                    removeInvalidIdx = [removeInvalidIdx i];
                end
            end
        end
        all_multiTaskConditions_cap(removeInvalidIdx,:) = [];
        
        % extract all bad multitasking conditons
        bad_multiTaskConditions_cap = all_multiTaskConditions_cap(~ismember(all_multiTaskConditions_cap, good_multiTaskConditions_cap, 'rows'),:);
        
        % determine number of available samples in both groups
        nGoodConditions = size(unique(good_multiTaskConditions_cap, 'rows'),1);
        nBadConditions = size(unique(bad_multiTaskConditions_cap, 'rows'),1);

        % subsample to ensure same number of samples per group
        nSamples = min(nGoodConditions, nBadConditions);
        
        % if subsampling is enabled
        if(subSample)
            good_multiTaskConditions_cap_subsampled = good_multiTaskConditions_cap(randsample(nGoodConditions, nSamples), :);
            bad_multiTaskConditions_cap_subsampled = bad_multiTaskConditions_cap(randsample(nBadConditions, nSamples), :);
        else
            good_multiTaskConditions_cap_subsampled = good_multiTaskConditions_cap;
            bad_multiTaskConditions_cap_subsampled = bad_multiTaskConditions_cap;
        end
        
        % exception check
        if(sum(ismember(good_multiTaskConditions_cap_subsampled, bad_multiTaskConditions_cap_subsampled, 'rows')) > 0)
            error('Good and bad multitasking conditions overlap. This should not be the case.');
        end
        
        if(nSamples > 0) 
            % all available indices
            goodMultiIdx = find(ismember(tasks_MultiCap, good_multiTaskConditions_cap_subsampled, 'rows'));
            badMultiIdx = find(ismember(tasks_MultiCap, bad_multiTaskConditions_cap_subsampled, 'rows'));

            % test good multitask performance
            [outData_good, ~, MSE_multi_good] = taskNet.runSet(input_MultiCap(goodMultiIdx,:), tasks_MultiCap(goodMultiIdx,:), train_MultiCap(goodMultiIdx,:));
            [~,~,multi_good_tasksIdx] = unique(tasks_MultiCap(goodMultiIdx,:),'rows');
            
            % store mean performance on all good multitasking conditions
            [GroupId,~,index_j]=unique(multi_good_tasksIdx);
            MSE_multi_good_GroupMean=arrayfun(@(k) mean(MSE_multi_good(index_j==k)),1:length(GroupId));
            MSEdata(capIdx).cap = cap;
            MSEdata(capIdx).goodMSEData = MSE_multi_good_GroupMean;

            % test bad multitask performance

            [outData_bad, ~, MSE_multi_bad] = taskNet.runSet(input_MultiCap(badMultiIdx,:), tasks_MultiCap(badMultiIdx,:), train_MultiCap(badMultiIdx,:));
            [~,~,multi_bad_tasksIdx] = unique(tasks_MultiCap(badMultiIdx,:),'rows');
            
            % store mean performance on all good multitasking conditions
            [GroupId,~,index_j]=unique(multi_bad_tasksIdx);
            MSE_multi_bad_GroupMean=arrayfun(@(k) mean(MSE_multi_bad(index_j==k)),1:length(GroupId));
            MSEdata(capIdx).badMSEData = MSE_multi_bad_GroupMean;

            multiPerformance_mean(2, capIdx) = mean(MSE_multi_good);
            multiPerformance_sem(2,capIdx) = std(MSE_multi_good)/sqrt(length(MSE_multi_good));

            multiPerformance_mean(3, capIdx) = mean(MSE_multi_bad);
            multiPerformance_sem(3,capIdx) = std(MSE_multi_bad)/sqrt(length(MSE_multi_bad));
            
        else
            % no samples available
            MSEdata(capIdx).cap = cap;

        end
                                                    
    end
    

end