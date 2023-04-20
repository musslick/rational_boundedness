function fitVal = Manuscript_Simulation_3_PRP_LS2000_fun(inputArg)
 
% clear all;
% inputArg = '1';

fitVal = 0;
inputArgDouble = str2double(inputArg);

initWeights = 0;
taskCorrs = [0 0.25 0.5 0.75 1];
init_taskCorr = taskCorrs(1);
fullOrthogonalization = 0;
fixOutputWeights = 0;
fixHiddenWeights = 0;
init_task_scale = 5;
outputWeightStrength = 2;
iterations_Train = 500;         % maximum number of training iterations
thresh = 0.001;                  % mean-squared error stopping criterion
incongruentOnly = 1;
skipPRPAnalysis = 0;

NPathways = 2;                     % number of pathways (i.e. number of feature dimensions % output dimensions)

interferenceTasks = [2 3];
primaryTasks = [1 4];
taskPair_dependent = primaryTasks; % primaryTasks

tasksToPerform = sort([primaryTasks, interferenceTasks]);

%% initialization
tic

log_version = 21;

% update meta simulation parameters
tauTaskOnly = 0;
SOA_range = 10:10:80;
tau_range =  [0.1];
correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 

% set up network parameters
init_scale = 0.1;           % scales for initialized random weights
learningRate = 0.3;             % learning rate
decay = 0.0000;                 % weight penalization parameter
bias = -2;                          % weight from bias units to hidden & output units

hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

% training environment parameters
NFeatures = 8;                      % number of feature units per stimulus input dimension
PathwayOverlap = 2;             % number of tasks per stimulus input dimension
numSglTasks = PathwayOverlap * NPathways;               % total number of tasks sampled
useCostomizedTasks = 1;     % set to 1 if using customized set of tasks (used for plotting purposes)
samplesPerTask_train = [];   % stimuli per task for training (if samplesPerTask is empty then training patterns include entire space of stimuli per task)
sdScale = 0;                        % max. variance for around each stimulus input unit activation
sameStimuliAcrossTasks = 1;     % use same stimuli across tasks? (this parameter is only relevant if sdScale > 0)
orthogonalizeReferenceVectors = 1;      % flag: set to 1 if reference vectors should be entirely orthogonal


rng('default')
rng('shuffle');


 %% simulation loop
                      
rep = inputArgDouble;
                   
    disp(strcat('repetition:',num2str(rep)));
                    
    for numHiddenIdx = 1:length(numHiddenUnits)

         % generate task environment
        samplesPerTask = samplesPerTask_train;
        [input, tasks, train, ~, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
        [train] = getLS2000Patterns(input, tasks, NFeatures, NPathways);
        
        % set up model
        nHidden =  numHiddenUnits(numHiddenIdx);
        taskNet_org = NNmodel(nHidden, learningRate, bias, init_scale, thresh, decay, hiddenPathSize, outputPathSize);
        taskNet_org.NPathways = NPathways;

        % train on all single tasks
        taskNet_org.setData(input, tasks, train);
        taskNet_org.configure(); 

    %% set up output weights
    if(initWeights || fullOrthogonalization)
        outputWeights = nan(size(taskNet_org.weights.W_TO));
        if(fixOutputWeights)
            if(NPathways == 2 & NFeatures == 3)
                vec = transpose([1 1 1 0 0 0]) * outputWeightStrength;
                outputWeights(:,1) = vec;
                outputWeights(:,3) = vec;
                outputWeights(:,2) = -vec;
                outputWeights(:,4) = -vec;
                taskNet_org.weights.W_TO = outputWeights;
            elseif(NPathways == 3 & NFeatures == 3)

                vec1 = [1 1 1 0 0 0 0 0 0] * outputWeightStrength;
                vec2 = [0 0 0 1 1 1 0 0 0] * outputWeightStrength;
                vec3 = [0 0 0 0 0 0 1 1 1] * outputWeightStrength;

                outputWeights(:,1) = vec1;
                outputWeights(:,2) = vec2;
                outputWeights(:,3) = vec3;
                outputWeights(:,4) = vec1;
                outputWeights(:,5) = vec2;
                outputWeights(:,6) = vec3;
                outputWeights(:,7) = vec1;
                outputWeights(:,8) = vec2;
                outputWeights(:,9) = vec3;
            else
                error('Could not fix output weights.');
            end
        end

        %% set up task weights according to max. cosine similarity

        num_points = NPathways;
        dim = nHidden(end);
        weightCorr = nan(1,NPathways);

        hiddenGroups = floor(dim/num_points);
        orthogonalGroups = floor(dim/(NPathways*NPathways));
        taskWeights = nan(dim, NPathways*NPathways);

        if(fullOrthogonalization)
            taskWeights = zeros(dim, NPathways*NPathways);
            for i = 1:(NPathways*NPathways)
                taskWeights((i-1)*orthogonalGroups+[1:orthogonalGroups],i) = randn(orthogonalGroups,1);
            end
        else
            for i = 1:num_points

                W = nan(dim, num_points);

                % create random initial vector
                x_ref_mask = zeros(dim, 1);
                x_ref_mask(((i-1)*hiddenGroups+1):(i*hiddenGroups)) = 1;
                if(~orthogonalizeReferenceVectors)
                    x_ref = randn(dim, 1)*0.1 + x_ref_mask;
                else
                    x_ref = x_ref_mask;
                end
                x_ref = x_ref./norm(x_ref);

                for p = 1:num_points

                    r = rand*(1-init_taskCorr) + init_taskCorr;
                    alpha = acos(r);

                    % 1) sample uniformly in S^(n-2) in R^(n-1)
                    l = sin(alpha);
                    y = rand(dim-1, 1)-0.5;
                    y = y./norm(y) * l;

                    % 2) assume a ~= e_i for i = 1:dim; use Gram-Schmidt to get orthonormal
                    % basis
                    I = eye(dim);
                    [Q,R] = qr([x_ref I(:,2:end)]);

                    S = Q(:,2:end);

                    % define T: R^(n-1) -> R^n

                    Ty = S*y + r*x_ref;

                    W(:,p) = Ty;

                end

                weightIdx = find(tril(ones(size(corr(W))),-1));
                W_corr = corr(W);
                weightCorr(i) = mean(W_corr(weightIdx));

                taskWeights(:,((i-1)*num_points+1):(i*num_points)) = W;

            end
        end
        taskWeights = taskWeights .*init_task_scale(end);

        % set up network weights
        taskNet_org.weights.W_TH = taskWeights;

        % log initial weight correlation
        weightCorr_log(rep) = mean(weightCorr);
        init_taskCorr_log(rep) = init_taskCorr;
    else
        taskNet_org.weights.W_TH  = taskNet_org.weights.W_TH; %* 2 * init_task_scale;
    end
    
    %% NETWORK TRAINING
        
        % log task environment
        
        % train neural network
        taskNet = NNmodel(taskNet_org);
        if(fixOutputWeights)
            taskNet.setFixedWeights([NNmodel.W_TASK_OUTPUT]);
        end
        if(fixHiddenWeights)
            taskNet.setFixedWeights([NNmodel.W_TASK_HIDDEN]);
        end

        % if training set is complete then train on full set interleaved

        MSE_log = zeros(1, iterations_Train);
        for iter =1:iterations_Train

            % sample primary training patterns
            samplesPerTask = samplesPerTask_train;
            [input_single, tasks_single, train_single] = createTaskPatterns(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform);
            [train] = getLS2000Patterns(input, tasks, NFeatures, NPathways);
            
            % put it all together
            input = [input_single];
            tasks = [tasks_single];

            % train network on samples
            taskNet.trainOnline(1, input, tasks, train);

            % is performance criterion reached?
            
            [~, ~, MSE] = taskNet.runSet(input, tasks, train);
            MSE_log(iter) = mean(MSE);
            if (mean(MSE) <= taskNet.thresh)
                disp('Learned primary tasks');
                break
            end
            disp([iter mean(MSE)]);
        end

        % performance log
        batch_log(numHiddenIdx, rep).MSE_log(:) = taskNet.MSE_log;
        batch_log(numHiddenIdx, rep).CE_log(:) = taskNet.CE_log;
        batch_log(numHiddenIdx, rep).CF_log(:) = taskNet.CF_log;
        batch_log(numHiddenIdx, rep).DimCF_log(:) = taskNet.DimCF_log;

        disp('Completed single task training.');
        
        %% perform analysis for average task correlations
        
        % perform graph-theoretic analysis on default similarity measure
        [inputSgl, tasksSgl] = createTaskPatterns(NPathways, NFeatures, [], sdScale, sameStimuliAcrossTasks, tasksToPerform);
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
        
        relevantTasks = [primaryTasks interferenceTasks];
        hiddenCorr = corr(taskNet.weights.W_TH(:,relevantTasks));
        outputCorr = corr(taskNet.weights.W_TO(:,relevantTasks));
        batch_log(numHiddenIdx, rep).hiddenCorr = hiddenCorr;
        batch_log(numHiddenIdx, rep).outputCorr = outputCorr;
        
        %% PRP ANALYSIS
        if(~skipPRPAnalysis)
            
            [input, tasks, train, tasksIdxSgl] = createTaskPatterns(NPathways, NFeatures, [], sdScale, sameStimuliAcrossTasks, tasksToPerform);
            [train] = getLS2000Patterns(input, tasks, NFeatures, NPathways);
            
            for tau_idx = 1:length(tau_range)

                verbose = 0;
                tauTask = tau_range(tau_idx);
                if(tauTaskOnly)
                    tauNet = 1;
                else
                    tauNet = tauTask;
                end

                nStimuli = [];

                % functional dependence
                taskA = taskPair_dependent(2);
                taskB = taskPair_dependent(1);
                [congruentIdx, incongruentIdx] = Part1_Sim3_LS2000_getCongruentIdx(input, tasksIdxSgl, NPathways, taskA, taskB);
                congruent = [];

                task1RT_log = size(SOA_range);
                task1RT_all_log = {};
                task1Accuracy_log = size(SOA_range);
                task1Iterations_log = size(SOA_range);
                task2RT_log = size(SOA_range);
                task2RT_all_log = {};
                task2Accuracy_log = size(SOA_range);
                task2Iterations_log = size(SOA_range);

                disp('PRP dependent...');
                for SOA_idx = 1:length(SOA_range)
                    disp(['SOA ' num2str(SOA_idx) '/' num2str(length(SOA_range))]);
                    SOA = SOA_range(SOA_idx);
                    [task1RT, task1Accuracy, task2RT, task2Accuracy, task1Iterations, task2Iterations, task1RT_all, task2RT_all] = Part1_Sim4_PRP_Analysis_V2(SOA, taskNet, input, train, tasksIdxSgl, taskA, taskB, congruent, tauNet, tauTask, nStimuli, verbose);
                    
                    task1RT_log(SOA_idx) = task1RT;
                    task1RT_all_log{SOA_idx} = task1RT_all;
                    task1Accuracy_log(SOA_idx) = task1Accuracy;
                    task1Iterations_log(SOA_idx) = task1Iterations;
                    task2RT_log(SOA_idx) = task2RT;
                    task2RT_all_log{SOA_idx} = task2RT_all;
                    task2Accuracy_log(SOA_idx) = task2Accuracy;
                    task2Iterations_log(SOA_idx) = task2Iterations;
                    
                    task1RT_inc_log(SOA_idx) = mean(task1RT_all(incongruentIdx));
                    task1RT_all_inc_log{SOA_idx} = task1RT_all(incongruentIdx);
                    task2RT_inc_log(SOA_idx) = mean(task2RT_all(incongruentIdx));
                    task2RT_all_inc_log{SOA_idx} = task2RT_all(incongruentIdx);
                    
                    task1RT_con_log(SOA_idx) = mean(task1RT_all(congruentIdx));
                    task1RT_all_con_log{SOA_idx} = task1RT_all(congruentIdx);
                    task2RT_con_log(SOA_idx) = mean(task2RT_all(congruentIdx));
                    task2RT_all_con_log{SOA_idx} = task2RT_all(congruentIdx);
                end
                task2RT_log
                
                batch_log(numHiddenIdx, rep).funcDependence(tau_idx).task1Accuracy_log = task1Accuracy_log;
                batch_log(numHiddenIdx, rep).funcDependence(tau_idx).task2Accuracy_log = task2Accuracy_log;
                
                batch_log(numHiddenIdx, rep).funcDependence_inc(tau_idx).task1RT_log = task1RT_inc_log;
                batch_log(numHiddenIdx, rep).funcDependence_inc(tau_idx).task1RT_all_log = task1RT_all_inc_log;
                batch_log(numHiddenIdx, rep).funcDependence_inc(tau_idx).task2RT_log = task2RT_inc_log;
                batch_log(numHiddenIdx, rep).funcDependence_inc(tau_idx).task2RT_all_log = task2RT_all_inc_log;
                
                batch_log(numHiddenIdx, rep).funcDependence_con(tau_idx).task1RT_log = task1RT_con_log;
                batch_log(numHiddenIdx, rep).funcDependence_con(tau_idx).task1RT_all_log = task1RT_all_con_log;
                batch_log(numHiddenIdx, rep).funcDependence_con(tau_idx).task2RT_log = task2RT_con_log;
                batch_log(numHiddenIdx, rep).funcDependence_con(tau_idx).task2RT_all_log = task2RT_all_con_log;
  
            end
        
        end
        
    end


save(['logfiles/Part1/PsychReview_Part1_Sim4_PRP_LS2000_' num2str(NPathways) 'P' num2str(NFeatures) 'F_' num2str(length(tasksToPerform)) 'tasks_' num2str(log_version) '_h' num2str(numHiddenUnits(1)) '_r' num2str(rep)]);

