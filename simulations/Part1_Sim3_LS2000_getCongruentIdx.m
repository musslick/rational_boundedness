function [congruentIdx, incongruentIdx] = Part1_Sim3_LS2000_getCongruentIdx(input, tasksIdxSgl, NPathways, taskA, taskB)

        NFeatures = size(input,2)/NPathways;

        % identify congruent & incongruent stimuli for the two tasks
        congruentIdx = [];
        incongruentIdx = [];
        
        inputA = input( tasksIdxSgl == taskA,:);
        
        inputDimA  = ceil(taskA / NPathways);
        inputDimB  = ceil(taskB / NPathways);
        
        for patternIdx = 1:size(inputA,1)
            
            featureA = find(inputA(patternIdx, (1:NFeatures) + (inputDimA-1)*NFeatures)==1);
            responseA = mod(featureA,2)+1;
            
            featureB = find(inputA(patternIdx, (1:NFeatures) + (inputDimB-1)*NFeatures)==1);
            responseB = mod(featureB,2)+1;
            
            if(responseA == responseB)
                congruentIdx = [congruentIdx patternIdx];
            else
                incongruentIdx = [incongruentIdx patternIdx];
            end
        end
        

end