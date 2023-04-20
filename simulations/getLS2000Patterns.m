function [train] = getLS2000Patterns(input, tasks, NFeatures, NPathways)

    train = zeros(size(input));
    
    for patternIdx = 1:size(input,1)
         
        taskA = find(tasks(patternIdx, :) == 1);
        inputDim  = ceil(taskA / NPathways);
        outputDim = mod(taskA-1, NPathways)+1;
        feature = find(input(patternIdx, (1:NFeatures) + (inputDim-1)*NFeatures)==1);
        response = mod(feature,2)+1 + (outputDim-1)*NFeatures;
        train(patternIdx, response) = 1;
        
    end

end