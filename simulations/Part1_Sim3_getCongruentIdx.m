function [congruentIdx, incongruentIdx, inputA] = Part1_Sim3_getCongruentIdx(input, NPathways, taskA, taskB)

        NFeatures = size(input,2)/NPathways;

        % identify congruent & incongruent stimuli for the two tasks
        congruentIdx = [];
        incongruentIdx = [];
        
        inputDimTaskA  = ceil(taskA / NPathways);
        inputDimTaskB = ceil(taskB / NPathways);
        for i = 1:size(input,1)
            A = reshape(input(i,:), 3, 3);
            if(isequal(A(:, inputDimTaskA), A(:, inputDimTaskB)))
                congruentIdx = [congruentIdx i];
            else
                incongruentIdx = [incongruentIdx i];
            end
        end
        

end