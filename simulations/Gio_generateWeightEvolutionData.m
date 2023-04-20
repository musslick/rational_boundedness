function fitVal = Gio_generateWeightEvolutionData(inputArg)

%     clear all;
%     inputArg = '1';

    fitVal = 0;
    inputArgDouble = str2double(inputArg);

    initWeights = 0;
    taskCorrs = [0 0.25 0.5 0.75 1];
    init_taskCorr = taskCorrs(1);
    fullOrthogonalization = 0;
    trainMultitasking = 1;
    fixOutputWeights = 0;
    fixHiddenWeights = 0;
    init_task_scale = 5;
    outputWeightStrength = 2;
    iterations_Train = 500;         % maximum number of training iterations
    thresh = 0.001;                  % mean-squared error stopping criterion
    incongruentOnly = 0;

    NPathways = 3;                     % number of pathways (i.e. number of feature dimensions % output dimensions)

    interferenceTasks = [];
    primaryTasks = [1:9];
    taskPair_structDependent = [1 4]; % primaryTasks
    taskPair_funcDependent = [1 5]; % primaryTasks
    taskPair_independent = [1 9]; % primaryTasks

    tasksToPerform = sort([primaryTasks, interferenceTasks]);
    %% description
    % 
    % Simulation # 1
    % 
    % The simulation involves training a neural networks on a set of tasks and
    % then predicting the multitasking performance of the network based on an
    % interference graph that is extracted from the networks learned
    % single-task representations
    %
    % author: Sebastian Musslick

    %% initialization
    tic

    log_version = 15;

    % update meta simulation parameters
    tauTaskOnly = 0;
    SOA_range = 5:5:60;
    tau_range = 0.1; %[0.1 0.2 0.5 1];
    correlation_thresholds =[0.5]; % thresholds picked for extracting interference graph from task similarity matrix
    numHiddenUnits = 100;%[200:100:1000];    % number of hidden units 

    % set up network parameters
    numRepetitions = 20;
    init_scale = 0.1;           % scales for initialized random weights
    learningRate = 0.3;             % learning rate
    decay = 0.0000;                 % weight penalization parameter
    bias = -2;                          % weight from bias units to hidden & output units

    hiddenPathSize = 1;             % group size of hidden units that receive the same weights from the task layer)
    outputPathSize = 1;             % group size of output units that receive the same weights from the task layer)

    % training environment parameters
    NFeatures = 3;                      % number of feature units per stimulus input dimension
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
        generateMultitaskingPatterns = 1;
        univalentStimuli = 0;
        [input, tasks, train, tasksIdxSgl, ~, ~, ~, ~, multiCap] = createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, generateMultitaskingPatterns, univalentStimuli);

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
        % log weights
        W_input_hidden = zeros([size(taskNet_org.weights.W_IH), iterations_Train]);
        W_task_hidden = zeros([size(taskNet_org.weights.W_TH), iterations_Train]);
        W_hidden_output = zeros([size(taskNet_org.weights.W_HO), iterations_Train]);
        W_task_output = zeros([size(taskNet_org.weights.W_TO), iterations_Train]);
        
        
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
            generateMultitaskingPatterns = 1;
            univalentStimuli = 0;
            [input_single, tasks_single, train_single] =createTaskPatterns_GenG(NPathways, NFeatures, samplesPerTask, sdScale, sameStimuliAcrossTasks, tasksToPerform, generateMultitaskingPatterns, univalentStimuli);

            % train on multitasking set
            if(trainMultitasking)

                input = [input_single; multiCap{3}.input];
                tasks = [tasks_single; multiCap{3}.tasks];
                train = [train_single; multiCap{3}.train];
                
            else
            
                % put it all together
                input = [input_single];
                tasks = [tasks_single];
                train = [train_single];
            end

            % train network on samples
            taskNet.trainOnline(1, input, tasks, train);

            % log weights
            W_input_hidden(:,:,iter) = taskNet.weights.W_IH;
            W_task_hidden(:,:,iter) = taskNet.weights.W_TH;
            W_hidden_output(:,:,iter) = taskNet.weights.W_HO;
            W_task_output(:,:,iter) = taskNet.weights.W_TO;
            
            % is performance criterion reached?
            
            [~, ~, MSE] = taskNet.runSet(input, tasks, train);
            MSE_log(iter) = mean(MSE);
            if (mean(MSE) <= taskNet.thresh)
                disp('Learned primary tasks');
                W_input_hidden(:,:,(iter+1):end) = [];
                W_task_hidden(:,:,(iter+1):end) = [];
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
        
    end
    
    if(~trainMultitasking)
        save(['gioLog/WeightEvolution_' num2str(rep)], 'W_input_hidden', 'W_task_hidden', 'W_hidden_output', 'W_task_output');
    else
        save(['gioLog/WeightEvolution_multi_' num2str(rep)], 'W_input_hidden', 'W_task_hidden', 'W_hidden_output', 'W_task_output');
    end
    
end

