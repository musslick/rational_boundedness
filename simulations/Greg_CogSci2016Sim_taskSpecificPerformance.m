%% description
% Simulation study as reported in
% Musslick S., Dey B., Oezcimder K., Patwary, M., Willke, T. L., Cohen, J. D. (submitted). Controlled vs. Automatic Processing: A Graph-Theoretic Approach to the Analysis of Serial vs. Parallel Processing in Neural Network Architectures. Proceedings of the 38th Annual Meeting of the Cognitive Science Society.

%% initialization
clc;
clear all;
tic

% set up network parameters
hiddenPathSize = 1;         % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;         % group size of output units that receive the same weights from the task layer)
learningRate = 0.3;         % learning rate
decay = 0.0000;             % weight penalization parameter
bias = -2;                  % weight from bias units to hidden & output units
iterations_Train = 500;     % number of training iterations

%% load existing training set

load('trainingSets/MIS_simulationDemo_6P2F.mat');

% update after loading full training set

% update meta simulation parameters
log_version = 1;
replications = 100;                 % # replications per simulation

% update network parameters
NPathways = 6;                     % number of pathways (i.e. number of stimulus dimensions and output dimensions)
NFeatures = 2; 			 % number of units per stimulus dimension and output dimension
thresh = 0.0001;                  % mean-squared error stopping criterion
PathwayOverlap = 2;
corr_threshold = 0.8;      % threshold for detecting merged task representations in correlation matrix
useCostomizedTasks = 3;
numSglTasks = 12;
init_scale = 0.02;           % scales for initialized random weights


nHidden = 200;              % number of hidden units


%% simulation loop

NTasks = PathwayOverlap * NPathways;

batch_log = repmat(struct('taskRepsHidden',zeros(replications, NTasks, nHidden), ...  
                          'taskRepsOutput',zeros(replications, NTasks, NPathways*NFeatures), ...  
                          'R_hidden',zeros(replications, NTasks, NTasks), ...
                          'R_output',zeros(replications, NTasks, NTasks), ...
                          'taskMatrix',zeros(replications, NPathways, NPathways), ...
                          'tasksToPerform',zeros(replications, NTasks), ...
                          'MSE_log',zeros(replications, iterations_Train), ...
                          'maxCarryingCapacity',zeros(replications, 1), ...
                          'inDegree_distrComplexity',zeros(replications, 1), ...
                          'true_inDegree_distrComplexity',zeros(replications, 1), ...
                          'outDegree_distrComplexity',zeros(replications, 1), ...
                          'NComponents_input',zeros(replications, 1), ...
                          'NComponents_output',zeros(replications, 1), ...
                          'taskNet',NNmodel.empty(replications,0)),1, 1); % zeros(replications,nTasks,nTasks)



for rep = 1:replications
    
%% select subset of possible single tasks to train

input_sglT = input;
tasks_sglT = tasks;
train_sglT = train;

if(useCostomizedTasks == 3)
    % pick random example with fixed out-degree distribution complexity
    A = eye(NPathways);
    for col = 1:size(A,2)
        
        add = randsample(NPathways-1,(PathwayOverlap-1));
        
        column = zeros((size(A,1)-1),1);
        column(add) = 1;
        A(col, [1:(col-1) (col+1):end]) = column;
    end
    
    % create task index mask
    taskMask = reshape(1:(NPathways*NPathways), NPathways, NPathways)';
    
    % get task indicies
    tasksToPerform = taskMask(A == 1)';

elseif(useCostomizedTasks == 2)
    
    % select example tasks with given distribution complexity
    allTasksSets = find(distrComplexity-mean(distrComplexity) == min(abs(distrComplexity-mean(distrComplexity))));
    %allTasksSets = find(distrComplexity-mean(distrComplexity) == mean(abs(distrComplexity-mean(distrComplexity))));
    
    % pick one example set of tasks
    taskSet = allTasksSets(randsample(length(allTasksSets),1));
    
    % create bipartite graph of task structure
    A = eye(N);
    
    for inputNode = 1:N
        A(inputNode,[1:(inputNode-1) (inputNode+1):end]) = U(V(taskSet, inputNode),:);
    end
    
    % create task index mask
    taskMask = reshape(1:(NPathways*NPathways), NPathways, NPathways)';
    
    % get task indicies
    tasksToPerform = taskMask(A == 1)';
    
    
elseif (useCostomizedTasks == 1)
        %tasksToPerform = [6*[0:5]+[1:6]];   % train only on the 6 tasks that each specify a separate pathway
        %tasksToPerform = 1:(NPathways*NPathways);
        tasksToPerform = [1,2,8,9,14,15,20,22,27,29,33,36];
else
    % randomly pick tasks if no tasks are specified
    tasksToPerform = randsample(nTasks, numSglTasks)';
end

% delete all tasks from training set that the network is not trained on
idxDelete = [];
for i = 1:size(tasks,1)
   taskIdx = find( tasks(i,:) == 1);
   if(~ismember(taskIdx, tasksToPerform))
      idxDelete = [idxDelete i]; 
   end
end

input_sglT(idxDelete,:) = [];
tasks_sglT(idxDelete,:) = [];
train_sglT(idxDelete,:) = [];

batch_log(1,1).taskMatrix(rep,:,:) = A;
batch_log(1,1).tasksToPerform(rep,:) = tasksToPerform;

% calculate similarity to imposed task structure
% in degree distribution complexity
in_Degree_dist = sum(A,1);
in_Degree_PROB = in_Degree_dist/sum(in_Degree_dist);
%
batch_log(1,1).true_inDegree_distrComplexity(rep) =  ...
    -1*sum(in_Degree_PROB.*log2(in_Degree_PROB));


%%

            disp(strcat('repetition:',num2str(rep)));
            
            % set up model
            taskNet = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
            taskNet.NPathways = NPathways;
            
            % train on all single tasks
            taskNet.setData(input_sglT, tasks_sglT, train_sglT);
            taskNet.configure(); 
            taskNet.trainOnline(iterations_Train);

            % MSE log
            batch_log(1,1).MSE_log(rep,:) = taskNet.MSE_log;
            batch_log(1,1).taskNet(rep) = taskNet;

            % track unit activations (for hidden & output layer) for all input
            % patterns in a single task environment
            [outData, hiddenData, ~] = taskNet.runSet(input, tasks, train);

            totalAct_sglT = [hiddenData outData];

            % average hidden unit activations over all inputs under a curresponding
            % task, do that for each task
            hiddenAct_pure = zeros(NPathways*NPathways, taskNet.Nhidden);
            setSize = size(input,1)/(NPathways*NPathways);
            for j = 1:size(hiddenAct_pure,1)
                hiddenAct_pure(j,:) = mean(totalAct_sglT((setSize*(j-1)+1):(setSize*j),1:taskNet.Nhidden));
            end

            % average output unit activations over all inputs under a curresponding
            % task, do that for each task
            outputAct_pure = zeros(NPathways*NPathways, taskNet.Noutput);
            for j = 1:size(outputAct_pure,1)
                outputAct_pure(j,:) = mean(totalAct_sglT((setSize*(j-1)+1):(setSize*j),(taskNet.Nhidden+1):end));
            end

            hiddenAct_pure_selected = hiddenAct_pure(tasksToPerform,:);
            outputAct_pure_selected = outputAct_pure(tasksToPerform,:);
            
            batch_log(1,1).taskRepsHidden(rep,:,:) = hiddenAct_pure_selected;
            batch_log(1,1).taskRepsOutput(rep,:,:) = outputAct_pure_selected;

            R_hidden = corr(hiddenAct_pure_selected');
            R_output = corr(outputAct_pure_selected');

            % correlation matrix
            batch_log(1,1).R_hidden(rep,:,:) = R_hidden;
            batch_log(1,1).R_output(rep,:,:) = R_output;

            % validate MIS
            [pathwayCapacities maxCarryingCapacity  BK_MIS A_bipartite A_tasksIdx multiPerformance_mean multiPerformance_sem MSEdata A_dual A_dual_scores] = CogSci_validateMIS(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, multiCap);
            
            batch_log(1,1).maxCarryingCapacity(rep) = maxCarryingCapacity;
            batch_log(1,1).pathwayCapacities{rep} = pathwayCapacities;
            batch_log(1,1).BK_MIS{rep} = BK_MIS;
            batch_log(1,1).A_bipartite{rep} = A_bipartite;
            batch_log(1,1).A_dual{rep} = A_dual;
            batch_log(1,1).A_tasksIdx{rep} = A_tasksIdx;
            batch_log(1,1).multiPerformance_mean{rep} = multiPerformance_mean;
            batch_log(1,1).multiPerformance_sem{rep} = multiPerformance_sem;
            batch_log(1,1).MSEdata{rep} = MSEdata;
            batch_log(1,1).A_dual_scores{rep} = A_dual_scores;
            
            % calculate similarity to imposed task structure
            % in degree distribution complexity
            in_Degree_dist = sum(A_dual,1);
            in_Degree_PROB = in_Degree_dist/sum(in_Degree_dist);
            %
            batch_log(1,1).inDegree_distrComplexity(rep) =  ...
                -1*sum(in_Degree_PROB.*log2(in_Degree_PROB));
            
            % in degree distribution complexity
            out_Degree_dist = sum(A_dual,2);
            out_Degree_PROB = out_Degree_dist/sum(out_Degree_dist);
            %
            batch_log(1,1).outDegree_distrComplexity(rep) =  ...
                -1*sum(out_Degree_PROB.*log2(out_Degree_PROB));
            
            % # output components
            batch_log(1,1).NComponents_input(rep) = size(A_dual,1);
            % # input components
            batch_log(1,1).NComponents_input(rep) = size(A_dual,2);
            
            
            % TEST TASK SPECIFIC PERFORMANCE
            
            % test multitasking for each multitasking condition
            for cap = 2

                batch_log(1,1).performance{rep, cap}.numberOfTasks = cap;

                taskConditions = unique(multiCap{cap}.tasks,'rows');
                taskIndicies = nan(size(taskConditions,1), cap);
                multiMSE = nan(size(taskConditions,1), 1);
                multiPCorrect = nan(size(taskConditions,1), 1);

                [outData, ~, multi_MSE_log] = taskNet.runSet(multiCap{cap}.input, multiCap{cap}.tasks, multiCap{cap}.train);
                 [~, PCorrect_tasks_mean, PCorrect_tasks] = taskNet.calculateOutcomeProbabilitiesTasks(outData,  multiCap{cap}.train, multiCap{cap}.tasks);

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

                for taskComb = 1:size(taskConditions,1)
                    taskIndicies(taskComb,:) = find(taskConditions(taskComb,:) == 1);

                    taskCond = taskConditions(taskComb,:);
                    taskCondIdx = find(transpose(ismember(multiCap{cap}.tasks, taskCond, 'rows')) & congruency == 0);

                    multiMSE(taskComb) = mean(multi_MSE_log(taskCondIdx));
                    multiPCorrect(taskComb) = mean(PCorrect_tasks_mean(taskCondIdx));
                end

                batch_log(1,1).performance{rep, cap}.nTaskCombinations = size(taskConditions,1);
                batch_log(1,1).performance{rep, cap}.taskIndicies = taskIndicies;
                batch_log(1,1).performance{rep, cap}.MSE = multiMSE;
                batch_log(1,1).performance{rep, cap}.PCorrect = multiPCorrect;

            end

end


save(['logfiles/CogSci_simulation_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version)]);

toc