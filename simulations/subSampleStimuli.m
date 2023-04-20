function [input_subsampled, tasks_subsampled, train_subsampled] = subSampleStimuli(input, tasks, train, nSampledStimuli, varargin)
    % optional argument: indices of subsampled stimuli

    nTasks = size(unique(tasks, 'rows'),1);
    nStim = size(input,1)/nTasks;
    
    if(~isempty(varargin))
        stimuliSubsampled = varargin{1};
        nSampledStimuli = length(stimuliSubsampled);
    else
        if(nStim >= nSampledStimuli)
            stimuliSubsampled = randsample(nStim, nSampledStimuli);
        else
            replace = 1;
            stimuliSubsampled = randsample(nStim, nSampledStimuli, replace);
            warning('Not enough stimuli per task available. Sampling stimuli with replacement.');
        end
    end
    
%     if(nSampledStimuli > nStim)
%         error('Number of sampled stimuli per task exceeds number of available stimuli per task.')
%     end

    input_subsampled = nan(nSampledStimuli * nTasks, size(input,2));
    tasks_subsampled = nan(nSampledStimuli * nTasks, size(tasks,2));
    train_subsampled = nan(nSampledStimuli * nTasks, size(train,2));

    for taskComb = 1:nTasks
        startRow = (taskComb-1)*nStim;
        startRow_subSampled = (taskComb-1)*nSampledStimuli;
        input_subsampled(startRow_subSampled+(1:nSampledStimuli), :) = input(startRow + stimuliSubsampled,:);
        tasks_subsampled(startRow_subSampled+(1:nSampledStimuli), :) = tasks(startRow + stimuliSubsampled,:);
        train_subsampled(startRow_subSampled+(1:nSampledStimuli), :) = train(startRow + stimuliSubsampled,:);
    end

end