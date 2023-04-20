
function [pathwayCapacities, maxCarryingCapacity, BK_MIS, A_bipartite, A_tasksIdx, A_dual, depencencyLogLCA_mean, depencencyAvgLogLCA_mean] = validateTaskSets(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, taskAccuracies)

    NPathways = taskNet.NPathways;

    % find MIS
    [pathwayCapacities, maxCarryingCapacity,  BK_MIS, A_bipartite, A_tasksIdx, A_dual] = getMaxCarryingCapacity(R_hidden, R_output, corr_threshold);
    
    % formatting outputs
    BK_MIS = extend_BK_MIS(BK_MIS);
    % remove duplicates
    
    A_tasksIdx(A_tasksIdx > 0) = tasksToPerform(A_tasksIdx(A_tasksIdx > 0));
    A_tasksOrder = A_tasksIdx(:);
    A_tasksOrder(A_tasksOrder == 0) = [];
    pathwayCapacities = [A_tasksOrder pathwayCapacities];
    
    maxSetSize = length(taskAccuracies);
    % initialize log
    for cap = 1:length(taskAccuracies)
        for dependency = 1:( (maxSetSize^2 - maxSetSize)/2 + 1) % + 1 to account for 0 dependencies
            numDependenciesLog{cap, dependency}.meanPerformanceLCA = [];
        end
    end
    
    for dependency = 1:( (maxSetSize^2 - maxSetSize)/2 + 1)  % + 1 to account for 0 dependencies
            numDependenciesAvgLog{dependency}.meanPerformanceLCA = [];
    end

    % identify for each multitasking combination the number of interfering tasks in the dependency graph
    for cap = 1:length(taskAccuracies)
        
        for taskComb = 1:size(taskAccuracies{cap}.allTaskCombs, 1)
            currentTasks = taskAccuracies{cap}.allTaskCombs(taskComb, :);
            % build sub dependency graph
            subGraph = A_dual;
            toRemove = find(~ismember(tasksToPerform, currentTasks));
            subGraph(toRemove,:) = [];
            subGraph(:,toRemove) = [];
            numDependencies = numel(find(subGraph > 0))/2;
            numDependenciesLog{cap, numDependencies+1}.meanPerformanceLCA = [numDependenciesLog{cap, numDependencies+1}.meanPerformanceLCA ...
                                                                                                                                  mean(taskAccuracies{cap}.LCA(taskComb, :))];
            numDependenciesAvgLog{numDependencies+1}.meanPerformanceLCA = [numDependenciesAvgLog{numDependencies+1}.meanPerformanceLCA ...
                                                                                                                                  mean(taskAccuracies{cap}.LCA(taskComb, :))];
        end
    end
    
    % condense num dependencies log
    depencencyLogLCA_mean = zeros(NPathways, NPathways+1);
    depencencyAvgLogLCA_mean = zeros(1, NPathways+1);
    
    for cap = 1:NPathways
        for numDependency = 1:(NPathways+1)
            depencencyLogLCA_mean(cap, numDependency) = mean(numDependenciesLog{cap, numDependency}.meanPerformanceLCA);
        end
    end
    
   for numDependency = 1:(NPathways+1)
            depencencyAvgLogLCA_mean(numDependency) = mean(numDependenciesAvgLog{numDependency}.meanPerformanceLCA);
    end

end