function [meanMSE, meanAccuracy, meanPCorrect] = getMeanNetworkPerformance(taskNet, input, tasks, train)

    [outData, ~, meanMSE]  = taskNet.runSet(input, tasks, train);
    [meanAccuracy] = taskNet.calculateAbsoluteErrorTasks(outData,  train, tasks);
    [~, meanPCorrect] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  train, tasks);

end



