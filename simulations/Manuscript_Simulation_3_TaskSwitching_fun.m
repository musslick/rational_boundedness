function fitVal = Manuscript_Simulation_3_TaskSwitching_fun(inputArg)

fitVal = 0;
inputArgDouble = str2double(inputArg);

%% initialization
% clc;
tic

log_version = 11;

% update meta simulation parameters
correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 

% set up network parameters
init_scale = 0.1;           % scales for initialized random weights
learningRate = 0.3;             % learning rate
decay = 0.0000;                 % weight penalization parameter
bias = -2;                          % weight from bias units to hidden & output units
iterations_Train = 5000;         % maximum number of training iterations
thresh = 0.01;                  % mean-squared error stopping criterion
tau = 0.15; 
tauTaskOnly = 0;

hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

% training environment parameters
NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
NFeatures = 3;                      % number of feature units per stimulus input dimension
PathwayOverlap = 2;             % number of tasks per stimulus input dimension
numSglTasks = PathwayOverlap * NPathways;               % total number of tasks sampled
useCostomizedTasks = 1;     % set to 1 if using customized set of tasks (used for plotting purposes)
samplesPerTask_train = [];   % stimuli per task for training (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
samplesPerTask_test = []; % stimuli per task for testing (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
sdScale = 0;                        % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;     % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)

rng('default')
rng('shuffle');
 %% simulation loop
                      
rep = inputArgDouble;
                   
    disp(strcat('repetition:',num2str(rep)));
                    
    for numHiddenIdx = 1:length(numHiddenUnits)

        %% select subset of possible single tasks to train
        
        % set seed for random generator, make sure to use different seed
        % for different replication
%         rng(rep);

        if(useCostomizedTasks == 0)
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

            % keep in-degree the same, rather than out-degree (just flip A)
            A =  transpose(A);
            
            % get task indicies
            tasksToPerform = taskMask(A == 1)';

        elseif (useCostomizedTasks == 1)
                % pick tasks to train (for plotting)
                tasksToPerform = [1,2,4,5,9];
        elseif (useCostomizedTasks == 2)
                % pick tasks to train (for plotting)
                 A = eye(NPathways);
                 nTasks = PathwayOverlap * NPathways;
                implementedTasks = NPathways;
                while (implementedTasks < nTasks)
                    % sample random task
                    sample = randi(NPathways*NPathways);
                    if(A(sample) ~= 1)
                        A(sample) = 1;
                        implementedTasks = implementedTasks + 1;
                    end
                end
                % create task index mask
                taskMask = reshape(1:(NPathways*NPathways), NPathways, NPathways)';

                % keep in-degree the same, rather than out-degree (just flip A)
                A =  transpose(A);

                % get task indicies
                tasksToPerform = taskMask(A == 1)';
        elseif (useCostomizedTasks == 3)
                % randomly pick tasks if no tasks are specified
                tasksToPerform = randsample(nTasks, numSglTasks)';
        end
disp('training on the following tasks:');        
tasksToPerform
        %% generate task environment
        samplesPerTask = samplesPerTask_train;
        [input, tasks, train, tasksIdxSgl, stimIdxSgl, inputSgl_mask, tasksSgl_mask, trainSgl_mask, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
        
        [multiCap_con, multiCap_inc] = splitTrainingPatternsByCongruency(multiCap, NFeatures, NPathways);
        
        % log task environment
        batch_log(numHiddenIdx, rep).tasksToPerform(:) = tasksToPerform;

        %% train neural network
        
        % set up model
        nHidden =  numHiddenUnits(numHiddenIdx);
        taskNet = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
        taskNet.NPathways = NPathways;

        % train on all single tasks
        taskNet.setData(input, tasks, train);
        taskNet.configure(); 
        
        % if training set is complete then train on full set interleaved
        if(isempty(samplesPerTask_train))
            
            taskNet.trainOnline(iterations_Train);
            
        else % else training set is sampled for each training interation
            
            MSE_log = zeros(1, iterations_Train);
            for iter =1:iterations_Train

                % sample training patterns
                samplesPerTask = samplesPerTask_train;
                [input, tasks, train] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
                
                % train network on samples
                taskNet.trainOnline(1, input, tasks, train);

                % is performance criterion reached?
                MSE_log(iter) = taskNet.MSE_log(end);
                if (taskNet.MSE_log(end) <= taskNet.thresh)
                    break
                end
                disp([iter taskNet.MSE_log(end)]);
            end

        end

        % performance log
        batch_log(numHiddenIdx, rep).MSE_log(:) = taskNet.MSE_log;
        batch_log(numHiddenIdx, rep).CE_log(:) = taskNet.CE_log;
        batch_log(numHiddenIdx, rep).CF_log(:) = taskNet.CF_log;
        batch_log(numHiddenIdx, rep).DimCF_log(:) = taskNet.DimCF_log;
        
        
        % sample testing data
        samplesPerTask = samplesPerTask_test;
        [inputSgl, tasksSgl, ~, ~, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
        
        %  mean hidden unit activity for a given tasks (averaged across all stimuli)
        [R_hidden_CorrTaskMean, R_output_CorrTaskMean] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform, ...
                                                                                         'similarityMetric', 'CorrTaskAvg', 'param', 'Mean', 'Pearson');

        %  mean correlation for the activations of different tasks under the same stimulus
        [R_hidden_MeanTaskCorr, R_output_MeanTaskCorr] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform, ...
                                                                                         'similarityMetric', 'AvgTaskCorr', 'param', 'Mean', 'Pearson');

        %  mean hidden unit activity for a given tasks (averaged across all stimuli)
        [R_hidden_CorrTaskSpearman, R_output_CorrTaskSpearman] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform, ...
                                                                                         'similarityMetric', 'CorrTaskAvg', 'param', 'Mean', 'Spearman');
        
        % correlation matrix
        batch_log(numHiddenIdx, rep).R_hidden(:,:) = R_hidden_MeanTaskCorr;       % stimulus wise correlations between tasks
        batch_log(numHiddenIdx, rep).R_output(:,:) = R_output_MeanTaskCorr;         % stimulus wise correlations between tasks
        batch_log(numHiddenIdx, rep).R_hidden_avg(:,:) = R_hidden_CorrTaskMean; % task correlations based on mean task activities
        batch_log(numHiddenIdx, rep).R_output_avg(:,:) = R_output_CorrTaskMean;   % task correlations based on mean task activities
        batch_log(numHiddenIdx, rep).R_hidden_Spearman(:,:) = R_hidden_CorrTaskSpearman; % spearman task rank correlations based on median task activities
        batch_log(numHiddenIdx, rep).R_output_Spearman(:,:) = R_output_CorrTaskSpearman;   % spearman task rank correlations based on median task activities

        disp('Completed single task training.');
        %% perform analysis for average task correlations
        
        % perform graph-theoretic analysis on default similarity measure
        [R_hidden, R_output] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform);

        % validate extracted MIS
        
        for corrIdx = 1:length(correlation_thresholds)
            
            corr_threshold = correlation_thresholds(corrIdx);
            
            % COMPUTE DEPENDENCY GRAPH AND MIS
            NPathways = taskNet.NPathways;

            % find MIS
            [pathwayCapacities, maxCarryingCapacity,  BK_MIS, A_bipartite, A_tasksIdx, A_dual] = getMaxCarryingCapacity(R_hidden, R_output, corr_threshold);

            batch_log(numHiddenIdx, rep).A_bipartite = A_bipartite;
            batch_log(numHiddenIdx, rep).A_dual = A_dual;
            batch_log(numHiddenIdx, rep).maxCarryingCapacity = maxCarryingCapacity;
            % formatting outputs
            BK_MIS = extend_BK_MIS(BK_MIS);
            % remove duplicates

            A_tasksIdx(A_tasksIdx > 0) = tasksToPerform(A_tasksIdx(A_tasksIdx > 0));
            A_tasksOrder = A_tasksIdx(:);
            A_tasksOrder(A_tasksOrder == 0) = [];
            pathwayCapacities = [A_tasksOrder pathwayCapacities];
    
            disp('Completed dependence graph computation.');
            
            disp(['Tested correlation threshold ' num2str(corrIdx) '/' num2str(length(correlation_thresholds))]);
        end
        
        %% TASK SWITCHING
        verbose = 0;
        tauTask = tau;
        if(tauTaskOnly)
            tauNet = 1;
        else
            tauNet = tauTask;
        end
        nStimuli = 100;
        
        % structural dependence
        taskA = 4;
        taskB = 1;
        
        congruent = 1;
        [repeatAccuarcy, switchAccuarcy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
        batch_log(numHiddenIdx, rep).structDependence.ER_repCon = 1-repeatAccuarcy;
        batch_log(numHiddenIdx, rep).structDependence.ER_switchCon = 1-switchAccuarcy;
        batch_log(numHiddenIdx, rep).structDependence.RT_repCon = repeatRT;
        batch_log(numHiddenIdx, rep).structDependence.RT_switchCon = switchRT;
         
        congruent = 0;
        [repeatAccuarcy, switchAccuarcy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
        batch_log(numHiddenIdx, rep).structDependence.ER_repInc = 1-repeatAccuarcy;
        batch_log(numHiddenIdx, rep).structDependence.ER_switchInc = 1-switchAccuarcy;
        batch_log(numHiddenIdx, rep).structDependence.RT_repInc = repeatRT;
        batch_log(numHiddenIdx, rep).structDependence.RT_switchInc = switchRT;
        
        % functional dependence
        taskA = 5;
        taskB = 1;
        
        congruent = 1;
        [repeatAccuarcy, switchAccuarcy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
        batch_log(numHiddenIdx, rep).funcDependence.ER_repCon = 1-repeatAccuarcy;
        batch_log(numHiddenIdx, rep).funcDependence.ER_switchCon = 1-switchAccuarcy;
        batch_log(numHiddenIdx, rep).funcDependence.RT_repCon = repeatRT;
        batch_log(numHiddenIdx, rep).funcDependence.RT_switchCon = switchRT;
         
        congruent = 0;
        [repeatAccuarcy, switchAccuarcy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
        batch_log(numHiddenIdx, rep).funcDependence.ER_repInc = 1-repeatAccuarcy;
        batch_log(numHiddenIdx, rep).funcDependence.ER_switchInc = 1-switchAccuarcy;
        batch_log(numHiddenIdx, rep).funcDependence.RT_repInc = repeatRT;
        batch_log(numHiddenIdx, rep).funcDependence.RT_switchInc = switchRT;
        
        % independence
        taskA = 9;
        taskB = 1;
        
        congruent = 1;
        [repeatAccuarcy, switchAccuarcy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
        batch_log(numHiddenIdx, rep).independence.ER_repCon = 1-repeatAccuarcy;
        batch_log(numHiddenIdx, rep).independence.ER_switchCon = 1-switchAccuarcy;
        batch_log(numHiddenIdx, rep).independence.RT_repCon = repeatRT;
        batch_log(numHiddenIdx, rep).independence.RT_switchCon = switchRT;
         
        congruent = 0;
        [repeatAccuarcy, switchAccuarcy, repeatRT, switchRT] = Sim1c_taskSwitchingAnalysis(taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
        batch_log(numHiddenIdx, rep).independence.ER_repInc = 1-repeatAccuarcy;
        batch_log(numHiddenIdx, rep).independence.ER_switchInc = 1-switchAccuarcy;
        batch_log(numHiddenIdx, rep).independence.RT_repInc = repeatRT;
        batch_log(numHiddenIdx, rep).independence.RT_switchInc = switchRT;
        
        
    end



save(['logfiles/Part1/PsychReview_Part1_Sim3_TaskSwitching_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_r' num2str(rep)]);

