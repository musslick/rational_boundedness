function fitVal = Manuscript_Simulation_1_FigSim_11_12(inputArg)

% Manuscript FIG_DEPENDENCY_RESULTS

inputArgDouble = str2double(inputArg);

%% initialization
% clc;
tic

log_version = 23;

% update meta simulation parameters
correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 

% set up network parameters
init_scale = 0.1;           % scales for initialized random weights
learningRate = 0.3;             % learning rate
decay = 0.0000;                 % weight penalization parameter
bias = -2;                          % weight from bias units to hidden & output units
iterations_Train = 5000;         % maximum number of training iterations
thresh = 0.001;                  % mean-squared error stopping criterion

hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

% training environment parameters
NPathways = 5;                     % number of pathways (i.e. number of feature dimensions % output dimensions)
NFeatures = 3;                      % number of feature units per stimulus input dimension
PathwayOverlap = 2;             % number of tasks per stimulus input dimension
numSglTasks = PathwayOverlap * NPathways;               % total number of tasks sampled
useCostomizedTasks = 2;     % set to 1 if using customized set of tasks (used for plotting purposes)
samplesPerTask_train = [];   % stimuli per task for training (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
samplesPerTask_test = []; % stimuli per task for testing (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
sdScale = 0;                        % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;     % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)

rng('default');
rng('shuffle');


 %% simulation loop
                      
rep = inputArgDouble;
                   
    disp(strcat('repetition:',num2str(rep)));
                    
    for numHiddenIdx = 1:length(numHiddenUnits)

        %% select subset of possible single tasks to train
        
        % set seed for random generator, make sure to use different seed
        % for different replication
        % rng(rep);
	time = fix(clock);
	rng(time(end));

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
                tasksToPerform = [1,2,8,9,14,15,20,22,27,29,33,36];
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
        
        % log task environment
        batch_log(numHiddenIdx, rep).taskMatrix(:,:) = A;
        batch_log(numHiddenIdx, rep).tasksToPerform(:) = tasksToPerform;

        % calculate similarity to imposed task structure
        % in degree distribution complexity
        in_Degree_dist = sum(A,1);
        in_Degree_PROB = in_Degree_dist/sum(in_Degree_dist);
        batch_log(numHiddenIdx, rep).inDegree_distrComplexity(:) =  ...
            -1*sum(in_Degree_PROB.*log2(in_Degree_PROB));

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
                [input, tasks, train] = createTaskPatterns(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);

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
        
        %% Evaluate inhibition from task layer as a function of set size
        
        % for each multitasking set size
        avg_output_net_task = nan(1, length(multiCap));
        for set_size = 1:length(multiCap)
           
            taskCombs = multiCap{set_size}.taskCombs;
            output_net_task = nan(size(taskCombs));
            
            % for each combination of tasks
            for comb = 1:size(taskCombs,1)
                
                rel_tasks = taskCombs(comb, :);
                % for each relevant output dimension
                for rel_output_dim = 1:size(taskCombs,2)
                    
                    task = taskCombs(comb, rel_output_dim);
                    output_dim = mod(task-1, NPathways)+1;
                    rel_features = ((NFeatures * (output_dim-1))+1):(NFeatures * output_dim);
                    
                    output_net_task(comb, rel_output_dim) = mean(sum(taskNet.weights.W_TO(rel_features,rel_tasks),2));
                end
                
            end
            avg_output_net_task(set_size) = mean(mean(output_net_task, 2), 1);
        end
        
        batch_log(numHiddenIdx, rep).avg_output_net_task = avg_output_net_task;
        
        %% perform analysis for average task correlations
        
        % perform graph-theoretic analysis on default similarity measure
        [R_hidden, R_output] = computeTaskSimilarity(taskNet, inputSgl, tasksSgl, 'tasks', tasksToPerform);

        % validate extracted MIS
        
        for corrIdx = 1:length(correlation_thresholds)
            
            corr_threshold = correlation_thresholds(corrIdx);
            
            % MIS ANALYSIS
            
            % compute optimal control policy
            [maximumPerformanceCurve, taskAccuracies, CDFs] = findOptimalControlPolicy(taskNet, NPathways, NFeatures, multiCap);
            
            disp('Completed maximal performance curve analysis.');
            
            % log data
            batch_log(numHiddenIdx, rep).maximumPerformanceCurveLog.LCA = maximumPerformanceCurve.LCA;
            batch_log(numHiddenIdx, rep).CDFs = CDFs;

            % ASSESS MULTITASKING ACCURACY AS A FUNCTION OF NUMBER OF DEPENDENT TASKS
            [pathwayCapacities, maxCarryingCapacity, BK_MIS, A_bipartite, A_tasksIdx, A_dual, depencencyLogLCA_mean, depencencyLogAvgLCA_mean] = validateTaskSets(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, taskAccuracies);
            
            disp('Completed set dependence analysis.');
            
            % only proceed if no tasks got lost after extracting bipartite
            % graph structure (i.e.. neglect situations in which two different tasks
            % have identical edges in bipartite graph)
            if(all(ismember(tasksToPerform, A_tasksIdx))) 
                
                batch_log(numHiddenIdx, rep).depencencyLogLCA_mean{corrIdx} = depencencyLogLCA_mean;
                batch_log(numHiddenIdx, rep).depencencyLogAvgLCA_mean{corrIdx} = depencencyLogAvgLCA_mean;
            
                batch_log(numHiddenIdx, rep).maxCarryingCapacity(corrIdx) = maxCarryingCapacity;
                batch_log(numHiddenIdx, rep).pathwayCapacities{corrIdx} = pathwayCapacities;
                batch_log(numHiddenIdx, rep).BK_MIS{corrIdx} = BK_MIS;
                batch_log(numHiddenIdx, rep).A_bipartite{corrIdx} = A_bipartite;
                batch_log(numHiddenIdx, rep).A_dual{corrIdx} = A_dual;
                batch_log(numHiddenIdx, rep).A_tasksIdx{corrIdx} = A_tasksIdx;

                % calculate similarity to imposed task structure
                % in degree distribution complexity
                in_Degree_dist = sum(A_dual,1);
                in_Degree_PROB = in_Degree_dist/sum(in_Degree_dist);
                %
                batch_log(numHiddenIdx, rep).inDegree_distrComplexity(corrIdx) =  ...
                    -1*sum(in_Degree_PROB.*log2(in_Degree_PROB));

                % in degree distribution complexity
                out_Degree_dist = sum(A_dual,2);
                out_Degree_PROB = out_Degree_dist/sum(out_Degree_dist);
                %
                batch_log(numHiddenIdx, rep).outDegree_distrComplexity(corrIdx) =  ...
                    -1*sum(out_Degree_PROB.*log2(out_Degree_PROB));

                % # output components
                batch_log(numHiddenIdx, rep).NComponents_input(corrIdx) = size(A_dual,1);
                % # input components
                batch_log(numHiddenIdx, rep).NComponents_input(corrIdx) = size(A_dual,2);

                % DEPENDENCY TYPE ANALYSIS: validate task dependencies
                [dualTaskCombs, taskCombData] = validateDualTaskDependencies(taskNet, NPathways, NFeatures, taskAccuracies, A_tasksIdx, inputSgl, tasksSgl, multiCap, CDFs);
                batch_log(numHiddenIdx, rep).dualTaskCombs{corrIdx} = dualTaskCombs;
                batch_log(numHiddenIdx, rep).taskCombData{corrIdx} = taskCombData;
                disp('Completed set dependency constellation analysis.');
                
            end

            disp(['Tested correlation threshold ' num2str(corrIdx) '/' num2str(length(correlation_thresholds))]);
        end
        

        
    end



save(['logfiles/Part1/PsychReview_Part1_Sim1_sharedRepresentation_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_r' num2str(rep)]);

toc


%% PLOT: Correct Output vs. Network Output

outData = taskNet.runSet(input, tasks, train);

figure(1);
subplot(1,2,1);
imagesc(outData);
title('output');
colorbar;
caxis([0 1]);

subplot(1,2,2)
imagesc(train)
title('correct');
colorbar;
caxis([0 1]);

%% STATISTICS
