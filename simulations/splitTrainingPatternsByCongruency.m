function [multiCap_con, multiCap_inc] = splitTrainingPatternsByCongruency(multiCap, NFeatures, NPathways)

        multiCap_con{1} = multiCap{1};
        multiCap_inc{1} = multiCap{1};
        
        for cap = 2:length(multiCap)

            congruency = nan(1, size(multiCap{cap}.train,1));
            for row = 1:length(congruency)
                congrSum = sum(reshape(multiCap{cap}.train(row,:), NFeatures, NPathways),2);
                numCorrectOutputUnits = sum(multiCap{cap}.train(row,:));
                if(max(congrSum) == numCorrectOutputUnits)
                    congruency(row) = 1;
                else
                    congruency(row) = 0;
                end
            end

            multiCap_con{cap}.input = multiCap{cap}.input(congruency == 1,:);
            multiCap_con{cap}.tasks = multiCap{cap}.tasks(congruency == 1,:);
            multiCap_con{cap}.train = multiCap{cap}.train(congruency == 1,:);
            multiCap_inc{cap}.input = multiCap{cap}.input(congruency == 0,:);
            multiCap_inc{cap}.tasks = multiCap{cap}.tasks(congruency == 0,:);
            multiCap_inc{cap}.train = multiCap{cap}.train(congruency == 0,:);

            disp(['cap ' num2str(cap) '/' num2str(length(multiCap))]);
        end
        
end