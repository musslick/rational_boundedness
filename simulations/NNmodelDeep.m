%% description
% this class implements a neural network model with two input layers, one
% hidden layer and an output layer. One input layer represents the current
% stimulus, the other input layer represents the current task. The input
% layer projects to the hidden layer and the hidden layer to the output
% layer. The task input projects to both the hidden layer and the output
% layer. Learning is implemented using backpropagation with optional weight
% penalization and optional weight constraint over the weights coming from
% the task input layer.
%
% author: Sebastian Musslick

%%

classdef NNmodelDeep < handle
    
    properties(SetAccess = public)
        trainSet;           % training data
        inputSet;           % input data
        taskSet;            % task data
        
        Ninput;             % number of input units
        Ntask;              % number of task units
        Nhidden;            % number of hidden units
        Noutput;            % number of output units
        
        hiddenPathSize      % size of one path (group of units) in hidden layer
        outputPathSize      % size of one path (group of units) in output layer
        NPathways           % optional parameter: number of processing pathways (can be used to infer input and output dimensions)
        
        bias_weight;        % default bias weight
        hidden_bias;        % hidden bias weight
        output_bias;        % output bias weight
        
        coeff;              % learning rate
        thresh;             % stopping criterion
        decay;              % weight penalization parameter
        weights;            % network weights
        fixedWeights;      % ENUM vector specifying which weights are fixed
        
        init_scale;         % scales for initialized weights
        init_task_scale     % scales for initialized weights from task to hidden layer
        
        hidden_act;         % current hidden activity
        output_act;         % current output activity
        
        MSE_log;            % full MSE curve
        MSE_patterns_log    % full MSE curve for all patterns
        CE_log;                 % full cross entropy error curve
        CE_patterns_log;    % full cross entropy error curve for all patterns
        CF_log;                 % full classification error curve
        CF_patterns_log;    % full classification error curve for all patterns
        DimCF_log;          % full dimension-wise classification error curve
        DimCF_patterns_log; % full dimension-wise classification error curve for all patterns
        
        hidden_log;         % full hidden activity for input set
        output_log;         % full output for input set
               
    end
    
    properties (Constant)
        
        W_INPUT_HIDDEN = 1;
        W_TASK_HIDDEN = 2;
        W_HIDDEN_OUTPUT = 3;
        W_TASK_OUTPUT = 4;
        
    end
    
    methods
        
        % constructor
        function this = NNmodelDeep(varargin)
            
            % make a copy of existing network object
            if(isa(varargin{1},'NNmodel'))
                copy = varargin{1};
                this.trainSet       = copy.trainSet;
                this.inputSet       = copy.inputSet;
                this.taskSet        = copy.taskSet;
                this.Ninput         = copy.Ninput;
                this.Ntask          = copy.Ntask;
                this.Nhidden        = copy.Nhidden;
                this.Noutput        = copy.Noutput;
                this.bias_weight    = copy.bias_weight;
                this.hidden_bias    = copy.hidden_bias;
                this.output_bias    = copy.output_bias;
                this.coeff          = copy.coeff;
                this.weights        = copy.weights;
                this.fixedWeights = copy.fixedWeights;
                this.init_scale     = copy.init_scale;
                this.init_task_scale        = copy.init_task_scale;
                this.hidden_act     = copy.hidden_act;
                this.output_act     = copy.output_act;
                this.MSE_log        = copy.MSE_log;
                this.CE_log        = copy.CE_log;
                this.CF_log        = copy.CF_log;
                this.MSE_patterns_log        = copy.MSE_patterns_log;
                this.hidden_log     = copy.hidden_log;
                this.output_log     = copy.output_log;
                this.hiddenPathSize = copy.hiddenPathSize;
                this.outputPathSize = copy.outputPathSize;
                this.thresh         = copy.thresh;
                this.decay          = copy.decay;
                this.NPathways          = copy.NPathways;
                return;
                
                % parse arguments
            else
               % number of hidden layer units
               this.Nhidden = varargin{1};  
            end
            
            % fixed weights only for bias units by default
               this.fixedWeights.W_IH = 0;
               this.fixedWeights.W_TH = zeros(1, length(this.Nhidden));
               this.fixedWeights.W_HH = zeros(1, length(this.Nhidden)-1);
               this.fixedWeights.W_HO = 0;
               this.fixedWeights.W_TO = 0;
               this.fixedWeights.W_BH = ones(1, length(this.Nhidden));
               this.fixedWeights.W_BO = 1;
            
            % learning rate
            if(length(varargin)>=2)
               this.coeff = varargin{2};  
            else
               this.coeff = 0.3; 
            end
            
            % weight from bias units to hidden and output units
            if(length(varargin)>=3)
               this.bias_weight = varargin{3};  
            else
               this.bias_weight = -1; 
            end
            
            % maximum absolute magnitude of initial weights
            if(length(varargin)>=4)
               this.init_scale = varargin{4};  
            else
               this.init_scale = 1; 
            end
            
            % mean-squared error stopping criterion for learning
            if(length(varargin)>=5)
               this.thresh = varargin{5};  
            else
               this.thresh = 0.01; 
            end
            
            % weight penalization parameter
            if(length(varargin)>=6)
               this.decay = varargin{6};  
            else
               this.decay = 0.02; 
            end
            
            % size of one path (group of units) in hidden layer
            if(length(varargin)>=7)
               this.hiddenPathSize = varargin{7};  
            else
               this.hiddenPathSize = 1; 
            end
            
            % size of one path (group of units) in output layer
            if(length(varargin)>=8)
               this.outputPathSize = varargin{8};  
            else
               this.outputPathSize = 1; 
            end
            
            % initialization noise from task to hidden layer
            this.init_task_scale = this.init_scale;
            
            % assign optional parameters
            this.NPathways = [];
            
        end
        
        % configure net: set up weights and network size depending on
        % trainig patterns
        function configure(this, varargin)
            
            % evaluate input arguments
            if(length(varargin)==1)
                
                % use configuration of existing net
                if(isa(varargin{1},'NNmodel'))
                   netObj = varargin{1};
                   this.weights = netObj.weights;
                   this.Ninput = netObj.Ninput;
                   this.Ntask = netObj.Ntask;
                   this.Noutput = netObj.Noutput;
                   this.Nhidden = netObj.Nhidden;
                   this.hidden_bias = netObj.hidden_bias;
                   this.output_bias = netObj.output_bias;
                   this.bias_weight = netObj.bias_weight;
                   this.hiddenPathSize = netObj.hiddenPathSize;
                   this.outputPathSize = netObj.outputPathSize;
                   this.thresh         = netObj.thresh;
                   this.decay          = netObj.decay;
                end
            else
                
                % set input patterns if provided by arguments
                if(length(varargin)>=3)
                   this.inputSet = varargin{1};
                   this.taskSet =  varargin{2};
                   this.trainSet = varargin{3};
                   
                
                % check if network has inputs, tasks and output patterns 
                else
                   if(isempty(this.inputSet) || isempty(this.taskSet) || isempty(this.trainSet))
                       error('Input set and training set need to be specified in order to configure network.');
                   end  
                end
                
                % set number of units for each layer
                this.Ninput = size(this.inputSet,2);
                this.Ntask = size(this.taskSet,2);
                this.Noutput = size(this.trainSet,2);
                if(isempty(this.Nhidden))
                    this.Nhidden = size(this.inputSet,2);
                end
                
                % set bias inputs for hidden & output layers
                if(isempty(this.bias_weight))
                    this.bias_weight = -1;           
                end
                this.hidden_bias = repmat(this.bias_weight,1,1);    % bias is the same for all hidden units
                this.output_bias = repmat(this.bias_weight,1,1);    % bias is the same for all output units 
                            
                % weight initialization (random using seed)
                rand('state',sum(100*clock));

                % set up weight matrices
                this.weights.W_IH = (-1 +2.*rand(this.Nhidden(1),this.Ninput))*this.init_scale;      % input-to-hidden weights

                for hiddenLayer = 1:length(this.Nhidden)
                    if(length(this.init_task_scale) < hiddenLayer)
                            % if no init task scale for hidden layer is specified then don't connect hidden layer to task layer
                            this.weights.W_TH{hiddenLayer} = zeros(this.Nhidden(hiddenLayer),this.Ntask);
                            this.fixedWeights.W_TH(hiddenLayer) = 1;
                    else
                            % standard initialization by weight range
                             this.weights.W_TH{hiddenLayer} = (-1 +2.*rand(this.Nhidden(hiddenLayer),this.Ntask))*this.init_task_scale(hiddenLayer);
                    end
                    if(hiddenLayer > 1)
                        this.weights.W_HH{hiddenLayer-1} = (-1 +2.*rand(this.Nhidden(hiddenLayer),this.Nhidden(hiddenLayer-1)))*this.init_scale;
                    end
                    
                    this.weights.W_BH{hiddenLayer} = ones(this.Nhidden(hiddenLayer),1);                                         % bias-to-hidden weights
                end
                
                this.weights.W_TO = (-1 +2.*rand(this.Noutput,this.Ntask))*this.init_scale;              
                this.weights.W_HO = (-1 +2.*rand(this.Noutput,this.Nhidden(end)))*this.init_scale;     % output-to-hidden weights
                this.weights.W_BO = ones(this.Noutput,1);                                         % bias-to-output weights
            end
        end
        
        % train the network on all patterns
        function [] = trainOnline(this, iterations, varargin)
            
            % parse arguments: input patterns, task pattenrs, output patterns
            if(length(varargin)>=2)
               inputData =  varargin{1};
               taskData = varargin{2};
               trainData =  varargin{3};
            else
               inputData = this.inputSet;
               taskData = this.taskSet;
               trainData = this.trainSet;
            end
            
            % check if input and task datasets have equal number of patterns (rows)
            if(size(inputData,1) ~= size(taskData,1))
                error('Task data has to have same number of rows as input data.');
            end
            
            % check if input and training datasets have equal number of patterns (rows)
            if(size(inputData,1) ~= size(trainData,1))
                error('Training data has to have same number of rows as input data.');
            end
            
            Ndatasets = size(inputData,1);              % total number of datasets
            this.MSE_log = zeros(1,iterations);         % log mean-squared error (MSE)
            this.CE_log = zeros(1,iterations);         % log cross-entropy error
            this.CF_log = zeros(1,iterations);         % log classification error
            this.DimCF_log = zeros(1,iterations);         % log dimension classification error
            this.MSE_patterns_log = zeros(Ndatasets, iterations);  % log MSE for all patterns
            this.CE_patterns_log = zeros(Ndatasets, iterations);  % log cross-entropy error for all patterns
            this.CF_patterns_log = zeros(Ndatasets, iterations);  % log classification error for all patterns
            this.DimCF_patterns_log = zeros(Ndatasets, iterations);  % log dimension classification error for all patterns
            
            % for each learning iteration
            for i = 1:iterations
               
               % randomize training set for each learning iteration
               order = randperm(size(inputData,1));
               inputData = inputData(order,:);
               taskData = taskData(order,:);
               trainData = trainData(order,:);
                
               MSE = zeros(1,Ndatasets);                            % current mean-squared error for all patterns (datasets)
               this.hidden_log = zeros(Ndatasets,this.Nhidden(1));     % current log activity of hidden units for each dataset
               this.output_log = zeros(Ndatasets,this.Noutput);     % current log activity of output units for each dataset
               
               % loop through all the patterns (online training)
               for dset = 1:Ndatasets
                  [MSE(dset)] = trainTrial(this, inputData(dset,:),  taskData(dset,:), trainData(dset,:)); % trains weights on current pattern
                  this.hidden_log(dset,:) = this.hidden_act{end}';       % log hidden unit activity for this pattern
                  this.output_log(dset,:) = this.output_act';       % log output unit activity for this pattern
               end
               
               [this.MSE_log(i),  this.MSE_patterns_log(:,i)] = this.calculateMeanSquaredError(this.output_log, trainData);
               [this.CE_log(i), this.CE_patterns_log(:,i)] = this.calculateCrossEntropyError(this.output_log, trainData);
               [this.CF_log(:,i), this.CF_patterns_log(:,i)] = this.calculateClassificationError(this.output_log, trainData);
               [this.DimCF_log(i), this.DimCF_patterns_log(:,i)] = this.calculateDimClassificationError(this.output_log, trainData);
               
               % stop learning if the mean-squared error reaches a certain
               % threshold
               if(this.MSE_log(i)) < this.thresh
                  break; 
               end
               
               disp(['iteration:' num2str(i)]);
               
            end
            
        end
        
        % train a trial
        function [MSE] = trainTrial(this, input, task, train)
            
            
               % simulate trial, retrieve activation values for hidden and
               % output layer
               this.runTrial(input, task);
               
               % weight update (backpropagation):
               % delta_w = -coeff * delta * x_i
               % delta_w      ...weight adjustment
               % coeff        ...learning rate
               % delta        ...represents backpropagated error
               % x_i          ...activation of sending unit

               % calculate delta's for output layer: delta_output = (output_act - train) * f_act'(netj)
               % delta_output ...backpropagated error for output units
               % output_act   ...output activity
               % train        ...correct output
               % f_act'(netj) ...first derivative of activation function of output units with respect to the net input
               error_term = (this.output_act - transpose(train));
               error_term(isnan(error_term)) = 0;                   % if target value is not specified (NaN), then error should be 0
               delta_output = error_term .* this.output_act .* (1 - this.output_act);

               % calculate delta's for hidden layer: delta_hidden = sum(delta_output * W_HO) * f_act'(netj)
               % delta_hidden ...backpropagated error for hidden units
               % delta_output ...backpropagated error for output units
               % W_HO         ...weights from hidden (columns) to output layer (rows)
               % f_act'(netj) ...first derivative of activation function of hidden units with respect to the net input
               for hiddenLayerIdx = 1:length(this.Nhidden)
                   
                   hiddenLayer = length(this.Nhidden)-hiddenLayerIdx+1;

                   % hidden layer that connects to output layer
                   if(hiddenLayer == length(this.Nhidden))
                        delta_hidden{hiddenLayer} = sum(repmat(delta_output,1,size(this.weights.W_HO,2)) .* this.weights.W_HO,1)' .* this.hidden_act{hiddenLayer} .* (1 - this.hidden_act{hiddenLayer});
                   else
                   % preceding hidden layers
                        delta_hidden{hiddenLayer} = sum(repmat(delta_hidden{hiddenLayer+1} ,1,size(this.weights.W_HH{hiddenLayer},2)) .* this.weights.W_HH{hiddenLayer},1)' .* this.hidden_act{hiddenLayer} .* (1 - this.hidden_act{hiddenLayer});
                   end
                   
               end
               
               % if a pathway size for the hidden unit layer is specified, 
               % then the deltas for groups of hidden layer units that 
               % receive input from the task input layer will be averaged. 
               % The path size specifies the number of hidden units in a
               % group. Each task unit projects to all groups of hidden
               % units with the constraint that the projecting weights will 
               % be the same for the units within a group
               delta_hiddenTask = delta_hidden;
               for hiddenLayer = 1:length(this.Nhidden)
                   if(this.hiddenPathSize(min(hiddenLayer, length(this.hiddenPathSize))) > 1) % no need to do averaging if hidden pathway size is 1
                       % average paths for hidden-to-task backprop
                       Npaths = this.Nhidden(hiddenLayer)/this.hiddenPathSize(hiddenLayer);
                       refVec = repmat(1:Npaths,this.hiddenPathSize(hiddenLayer),1);
                       refVec = refVec(:);

                       for i = 1:Npaths
                           delta_hiddenTask{hiddenLayer}(refVec==i) = mean(delta_hidden{hiddenLayer}(refVec==i));
                       end
                   end 
               end
               
               % if a pathway size for the output unit layer is specified, 
               % then the deltas for groups of output layer units that 
               % receive input from the task input layer will be averaged. 
               % The path size specifies the number of output units in a
               % group. Each task unit projects to all groups of output
               % units with the constraint that the projecting weights will 
               % be the same for the units within a group
               delta_outputTask = delta_output;
               if(this.outputPathSize > 1)
                   % average paths for output-to-task backprop
                   Npaths = this.Noutput/this.outputPathSize;
                   refVec = repmat(1:Npaths,this.outputPathSize,1);
                   refVec = refVec(:);
                   for i = 1:Npaths
                       delta_outputTask(refVec==i) = mean(delta_output(refVec==i));
                   end
               end

               % adjust weights from input to hidden units
               if(~this.fixedWeights.W_IH)
                    this.weights.W_IH = this.weights.W_IH - (this.coeff * delta_hidden{1} * input) - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_IH));
                    weightLog.dW_IH = - this.coeff * delta_hidden{1} * input - this.coeff * this.decay * sign(this.weights.W_IH) - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_IH));
               else
                   weightLog.dW_IH = 0;
               end
               
               % adjust weights to hidden layers
               for hiddenLayer = 1:length(this.Nhidden)
                   
                   % adjust weights from hidden to hidden layer
                   if(hiddenLayer < length(this.Nhidden))
                       if(~this.fixedWeights.W_HH(hiddenLayer))
                            this.weights.W_HH{hiddenLayer} = this.weights.W_HH{hiddenLayer} - (this.coeff * delta_hidden{hiddenLayer+1} * this.hidden_act{hiddenLayer}') - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_HH{hiddenLayer}));
                            weightLog.dW_HH{hiddenLayer} = - (this.coeff * delta_hidden{hiddenLayer+1} * this.hidden_act{hiddenLayer}') - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_HH{hiddenLayer}));
                       else
                            weightLog.dW_HH{hiddenLayer} = 0;
                       end
                   end
               
                   % adjust weights from task to hidden units
                   if(~this.fixedWeights.W_TH(hiddenLayer))
                        this.weights.W_TH{hiddenLayer} = this.weights.W_TH{hiddenLayer} - (this.coeff * delta_hiddenTask{hiddenLayer} * task .* this.dg(this.weights.W_TH{hiddenLayer})) - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_TH{hiddenLayer}));
                        weightLog.dW_TH{hiddenLayer} = - this.coeff * delta_hiddenTask{hiddenLayer} * task * this.dg(this.weights.W_TH{hiddenLayer}) - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_TH{hiddenLayer}));
                   else
                        weightLog.dW_TH{hiddenLayer} = 0;
                   end
                   
                   if(~this.fixedWeights.W_BH(hiddenLayer))
                        this.weights.W_BH{hiddenLayer} = this.weights.W_BH{hiddenLayer} - this.coeff * delta_hidden{hiddenLayer} * this.hidden_bias - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_BH{hiddenLayer})); 
                   end
               
               end
               
               % adjust weights from task to output units
               if(~this.fixedWeights.W_TO)
                   this.weights.W_TO = this.weights.W_TO - (this.coeff * delta_outputTask * task) .* this.dg(this.weights.W_TO) - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_TO));
                   weightLog.dW_TO = - (this.coeff * delta_outputTask * task) .* this.dg(this.weights.W_TO) - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_TO));
               else
                   weightLog.dW_TO = 0;
               end
               
               % adjust weights from hidden to output units
               if(~this.fixedWeights.W_HO)
                   this.weights.W_HO = this.weights.W_HO - (this.coeff * delta_output * this.hidden_act{end}') - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_HO));
                   weightLog.dW_HO = - (this.coeff * delta_output * this.hidden_act{end}') - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_HO));
               else
                   weightLog.dW_HO = 0;
               end
               
               % this.decay * sign(this.weights.W_IH) ...penalize weights
               % this.dg(this.weights.W_TO)           ...derivative of transformation function of the weights
               
               
               % learning from bias units to output units
               if(~this.fixedWeights.W_BO)
                    this.weights.W_BO = this.weights.W_BO - this.coeff * delta_output * this.output_bias - (this.coeff * this.decay / size(this.inputSet,1) * sign(this.weights.W_BO));
               end
               
               % calculate mean-squared error
               [~,  MSE] = this.calculateMeanSquaredError(this.output_act', train);

        end
        
         % train the network on all patterns
        function [] = trainBatch(this, iterations, varargin)
            
            error('trainBatch function not yet implemented for NNmodelDeep.');
            
            % parse arguments: input patterns, task pattenrs, output patterns
            if(length(varargin)>=2)
               inputData =  varargin{1};
               taskData = varargin{2};
               trainData =  varargin{3};
            else
               inputData = this.inputSet;
               taskData = this.taskSet;
               trainData = this.trainSet;
            end
            
            % check if input and task datasets have equal number of patterns (rows)
            if(size(inputData,1) ~= size(taskData,1))
                error('Task data has to have same number of rows as input data.');
            end
            
            % check if input and training datasets have equal number of patterns (rows)
            if(size(inputData,1) ~= size(trainData,1))
                error('Training data has to have same number of rows as input data.');
            end
            
            Ndatasets = size(inputData,1);              % total number of datasets
            this.MSE_log = zeros(1,iterations);         % log mean-squared error (MSE)
            this.MSE_patterns_log = zeros(Ndatasets, iterations); % log MSE for all patterns
            
            % calculate groups of hidden units, if a pathway size from task to hidden layer
            % is specified (for details, see below)
            if(this.hiddenPathSize > 1) % no need to do averaging if hidden pathway size is 1
               % average paths for hidden-to-task backprop
               Npaths_hidden = this.Nhidden/this.hiddenPathSize;
               refVec_hidden = repmat(1:Npaths_hidden,this.hiddenPathSize,1);
               refVec_hidden = refVec_hidden(:);
            end
            
            % calculate groups of output units, if a pathway size from task to output layer
            % is specified (for details, see below)
            if(this.outputPathSize > 1) % no need to do averaging if output pathway size is 1
               % average paths for output-to-task backprop
               Npaths_output = this.Noutput/this.outputPathSize;
               refVec_output = repmat(1:Npaths_output,this.outputPathSize,1);
               refVec_output = refVec_output(:);
            end
            
            
            % for each learning iteration
            for i = 1:iterations
                
               this.hidden_log = zeros(Ndatasets,this.Nhidden);     % current log activity of hidden units for each dataset
               this.output_log = zeros(Ndatasets,this.Noutput);     % current log activity of output units for each dataset
               
               % simulate trial, retrieve activation values for hidden and
               % output layer for each dataset
               [outData, hiddenData, MSE] = this.runSet(inputData, taskData, trainData);
               
               % weight update (backpropagation):
               % delta_w = -coeff * delta * x_i
               % delta_w      ...weight adjustment
               % coeff        ...learning rate
               % delta        ...represents backpropagated error
               % x_i          ...activation of sending unit

               % calculate delta's for output layer: delta_output = (output_act - train) * f_act'(netj)
               % delta_output ...backpropagated error for output units
               % output_act   ...output activity
               % train        ...correct output
               % f_act'(netj) ...first derivative of activation function of output units with respect to the net input
               error_term = transpose(outData - trainData);
               error_term(isnan(error_term)) = 0;                   % if target value is not specified (NaN), then error should be 0
               delta_output = error_term .* this.output_act .* (1 - this.output_act);
               
               % calculate delta's for hidden layer: delta_hidden = sum(delta_output * W_HO) * f_act'(netj)
               % delta_hidden ...backpropagated error for hidden units
               % delta_output ...backpropagated error for output units
               % W_HO         ...weights from hidden (columns) to output layer (rows)
               
               % f_act'(netj) ...first derivative of activation function of hidden units with respect to the net input
               % delta_hidden = sum(repmat(delta_output,1,size(this.weights.W_HO,2)) .* this.weights.W_HO,1)' .* this.hidden_act .* (1 - this.hidden_act);
               delta_hidden = delta_output' * this.weights.W_HO .* hiddenData .* (1 - hiddenData);
               
               % if a pathway size for the hidden unit layer is specified, 
               % then the deltas for groups of hidden layer units that 
               % receive input from the task input layer will be averaged. 
               % The path size specifies the number of hidden units in a
               % group. Each task unit projects to all groups of hidden
               % units with the constraint that the projecting weights will 
               % be the same for the units within a group
               delta_hiddenTask = delta_hidden;
               if(this.hiddenPathSize > 1) % no need to do averaging if hidden pathway size is 1
                   % average paths for hidden-to-task backprop
                   [GroupId,~,index_j]=unique(refVec_hidden);
                   GroupMean=arrayfun(@(k) mean(delta_hidden(:,index_j==k),2),1:length(GroupId), 'UniformOutput', 0);
                   delta_hiddenTask=[GroupMean{index_j}];
               end
               
               % if a pathway size for the output unit layer is specified, 
               % then the deltas for groups of output layer units that 
               % receive input from the task input layer will be averaged. 
               % The path size specifies the number of output units in a
               % group. Each task unit projects to all groups of output
               % units with the constraint that the projecting weights will 
               % be the same for the units within a group
               delta_output = delta_output';
               delta_outputTask = delta_output;
               if(this.outputPathSize > 1) % no need to do averaging if hidden pathway size is 1
                   % average paths for hidden-to-task backprop
                   [GroupId,~,index_j]=unique(refVec_output);
                   GroupMean=arrayfun(@(k) mean(delta_output(:,index_j==k),2),1:length(GroupId), 'UniformOutput', 0);
                   delta_outputTask=[GroupMean{index_j}];
               end
               
               % adjust weights from input to hidden units
               if(~ismember(this.W_INPUT_HIDDEN, this.fixedWeights))
                    this.weights.W_IH = this.weights.W_IH - this.coeff * delta_hidden' * inputData - Ndatasets * this.coeff * this.decay * sign(this.weights.W_IH);
               end
               % adjust weights from task to hidden units
               if(~ismember(this.W_TASK_HIDDEN, this.fixedWeights))
                    this.weights.W_TH = this.weights.W_TH - this.coeff * delta_hiddenTask' * taskData .* this.dg(this.weights.W_TH);
               end

               % adjust weights from task to output units
               if(~ismember(this.W_TASK_OUTPUT, this.fixedWeights))
                    this.weights.W_TO = this.weights.W_TO - this.coeff * delta_outputTask' * taskData .* this.dg(this.weights.W_TO);
               end
               % adjust weights from hidden to output units
               if(~ismember(this.W_HIDDEN_OUTPUT, this.fixedWeights))
                    this.weights.W_HO = this.weights.W_HO - this.coeff * delta_output' * hiddenData - this.coeff * this.decay * sign(this.weights.W_HO);
               end
               
               % this.decay * sign(this.weights.W_IH) ...penalize weights
               % this.dg(this.weights.W_TO)           ...derivative of transformation function of the weights
               
               % learning of bias weights is turned off for now
               % this.weights.W_BH = this.weights.W_BH - this.coeff * sum(delta_hidden,2) * this.hidden_bias; 
               % this.weights.W_BO = this.weights.W_BO - this.coeff * sum(delta_output,2) * this.output_bias;

               % calculate mean-squared error
               MSE = sum((outData - trainData).^2,2);
               
               [this.MSE_log(i),  this.MSE_patterns_log(:,i)] = this.calculateMeanSquaredError(outData, trainData);
               
               this.hidden_log = outData;                   % log hidden unit activity for this pattern
               this.output_log = hiddenData;                % log output unit activity for this pattern
               [this.MSE_log(i),  this.MSE_patterns_log(:,i)] = this.calculateMeanSquaredError(outData, trainData);
               [this.CE_log(i), this.CE_patterns_log(:,i)] = this.calculateCrossEntropyError(outData, trainData);
               [this.CF_log(i), this.CF_patterns_log(:,i)] = this.calculateClassificationError(outData, trainData);
               [this.DimCE_log(i), this.DimCE_patterns_log(:,i)] = this.calculateDimClassificationError(outData, trainData);
               
               this.MSE_patterns_log(:,i) = MSE;            % calculate mean-squared error for whole set of patterns
               
               % stop learning if the mean-squared error reaches a certain
               % threshold
               if(this.MSE_log(i)) < this.thresh
                  break; 
               end
               
            end
            
        end

        % run through a data set (no training)
        function [outData, hiddenData, MSE, hidden_net, output_net, ceError, classError, classDimError] = runSet(this, varargin)
            
            % parse arguments: input patterns, task pattenrs, output patterns
            if(length(varargin) >=2)
                inputData = varargin{1};
                taskData = varargin{2};
            else
                inputData = this.inputSet;
                taskData = this.taskSet;
            end
            
            % check if input and task datasets have equal number of patterns (rows)
            if(size(inputData,1) ~= size(taskData,1))
                error('Task data has to have same number of rows as input data.');
            end
            
            Ndatasets = size(inputData,1);                  % total number of datasets
            
            % for each hidden layer 
            for hiddenLayer = 1:length(this.Nhidden)
    
                % calculate net inputs for hidden layer     
                hidden_net_task = this.g(this.weights.W_TH{hiddenLayer}) * transpose(taskData);                              % input from task layer (task cue)
                hidden_net_bias  = this.weights.W_BH{hiddenLayer} * (this.hidden_bias * ones(1,Ndatasets));   % input from hidden bias units
                
                if(hiddenLayer == 1)
                    hidden_net_input = this.weights.W_IH * transpose(inputData);                                    % input from input layer (stimulus)
                else
                    hidden_net_input = this.weights.W_HH{hiddenLayer-1} * this.hidden_act{hiddenLayer-1}; 
                end
                hidden_net{hiddenLayer} = hidden_net_input + hidden_net_task + hidden_net_bias;                              % total net input to hidden layer

                % calculate activation for hidden units
                this.hidden_act{hiddenLayer} = 1./(1+exp(-hidden_net{hiddenLayer}));                                                      % use sigmoid activation function

            end

           % calculate net input for output units
            output_net_task = this.g(this.weights.W_TO) * transpose(taskData);                              % input from task layer (task cue)
            output_net_hidden = this.weights.W_HO * this.hidden_act{end};                                        % input from hidden layer
            output_net_bias   = this.weights.W_BO * (this.output_bias * ones(1,size(hidden_net_input,2)));  % input from output bias units
            output_net = output_net_hidden + output_net_task + output_net_bias;                             % total net input to output layer

            % calculate activation of output units
            this.output_act = 1./(1+exp(-output_net)); 

            % final network output
            for hiddenLayer = 1:length(this.hidden_act)
                hiddenData{hiddenLayer} = this.hidden_act{hiddenLayer}';  % log activity of hidden units for each dataset
            end
            % hiddenData = this.hidden_act{end}';                      % log activity of hidden units for each dataset
            outData = this.output_act';                         % log activity of output units for each dataset
            for hiddenLayer = 1:length(this.Nhidden)
                hidden_net{hiddenLayer} = hidden_net{hiddenLayer}';
            end
            output_net = output_net';
            
            % calculate MSE if correct output provided (train)
            MSE = -1*ones(1,Ndatasets);
            if(length(varargin)>=3)
                trainData = varargin{3};
                if(size(trainData,2) == size(outData,2))
                    % MSE
                    [~, MSE] = this.calculateMeanSquaredError(outData, trainData);
                    
                    % classification error
                    [~, classError] = this.calculateClassificationError(outData, trainData);
                    
                    % dimension classification error
                    [~, classDimError] = this.calculateDimClassificationError(outData, trainData);
                    
                    % cross-entropy error
                    [~, ceError] = this.calculateCrossEntropyError(outData, trainData);
                else
                    warning('Training data has to have same number of rows as input data. Can''t calculate MSE for each dataset.');
                end
            end
            
            % log activities for hidden and output units (for all patterns)
            this.hidden_log = hiddenData;
            this.output_log = outData;
            
        end
        
        
        % run a trial (feedforward step, no training)
        function [output_act, hidden_act, MSE, hidden_net, output_net] = runTrial(this, input, task, varargin)
            
            % initialize output activity
            output_act = zeros(1,this.Noutput);
            
            % for each hidden layer 
            for hiddenLayer = 1:length(this.Nhidden)
                
                hidden_net_task = this.g(this.weights.W_TH{hiddenLayer}) * transpose(task);      % input from task layer (task cue)
                hidden_net_bias  = this.weights.W_BH{hiddenLayer} * this.hidden_bias;            % input from hidden bias units
                
                % calculate net inputs for hidden layers
                if(hiddenLayer == 1)
                    hidden_net_input = this.weights.W_IH * transpose(input);            % input from input layer (stimulus)       
                else
                    hidden_net_input = this.weights.W_HH{hiddenLayer-1} * this.hidden_act{hiddenLayer-1}; % input from previous hidden layer  
                end
                hidden_net = hidden_net_input + hidden_net_task + hidden_net_bias;  % total net input to hidden layer
                
                % calculate activation for hidden units
                this.hidden_act{hiddenLayer} = 1./(1+exp(-hidden_net));                          % use sigmoid activation function
            end


            % calculate net input for output units
            output_net_task = this.g(this.weights.W_TO) * transpose(task);      % input from task layer (task cue)
            output_net_hidden = this.weights.W_HO * this.hidden_act{end};            % input from hidden layer
            output_net_bias   = this.weights.W_BO * this.output_bias;           % input from output bias units
            output_net = output_net_hidden + output_net_task + output_net_bias; % total net input to output layer

            % calculate activation of output units
            this.output_act = 1./(1+exp(-output_net)); 

            % final network output
            output_act(:) = this.output_act';
            hidden_act = this.hidden_act;
            
            % calculate MSE if correct output provided (train)
            MSE = -1;
            if(~isempty(varargin))
                train = varargin{1};
                if(length(train) == length(output_act))
                    [~,  MSE] = this.calculateMeanSquaredError(this.output_act', train);
                end
            end
            
        end
        
        function [MSE_total, MSE_patterns] = calculateMeanSquaredError(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            MSE_patterns = sum((outData - trainData).^2, 2)./this.Noutput';
            MSE_total = mean(MSE_patterns);
            
        end
        
        function [CE_error_total, CE_error_patterns] = calculateCrossEntropyError(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            outData_normalized = outData./repmat(sum(outData,2),1,size(outData,2));
            CE_error_patterns = -sum(log(outData_normalized).*trainData,2);
            
            CE_error_total = mean(CE_error_patterns);
        end
        
        function [CF_error_total, CF_error_patterns] = calculateDimClassificationError(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            CF_error_patterns = nan(size(trainData,1),1);
            
            if(isempty(this.NPathways))
                warning('NPathways needs to be specified in order to calculate classificaiton error.');
            else
                
                % calculate number of features per output dimension
                NFeatures  = size(trainData,2)/this.NPathways;
                
                % loop through each pattern
                for i = 1:size(trainData,1)
                
                    groupedTrainPattern = reshape(trainData(i,:), NFeatures, this.NPathways);
                    groupedOutputPattern = reshape(outData(i,:), NFeatures, this.NPathways);
                    
                    CF_error_patterns(i) = 0;
                    
                    for col = 1:this.NPathways
                        
                        correctIdx = find(groupedTrainPattern(:, col) == 1);
                        
                        if(isempty(correctIdx))
                            continue;
                        else
                            maxIdx = find(groupedOutputPattern(:, col) == max(groupedOutputPattern(:, col)));
                            if(maxIdx ~= correctIdx)
                                CF_error_patterns(i) = CF_error_patterns(i)  + 1;
                            end
                        end
                        
                    end
                    
                end
            end
            
            CF_error_total = mean(CF_error_patterns);
        end
        
        function [CF_error_total, CF_error_patterns] = calculateClassificationError(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
        
            maxOutputs = max(outData,[],2);
            classes = zeros(size(outData));
            for i = 1:length(maxOutputs)
               classes(i,  outData(i,:) == maxOutputs(i)) = 1;
            end
            classError = sum(abs(classes-trainData),2);
            CF_error_patterns = zeros(length(classError),1);
            CF_error_patterns(classError > 0) = 1;
            
            CF_error_total = mean(CF_error_patterns);
        end
        
        % computes MSE only for tasks that are about to be executed
        % (ignores output dimensions that are not relevant)
        function [MSE_tasks_mean, MSE_tasks] = calculateMeanSquaredErrorTasks(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            taskData = this.taskSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(length(varargin) > 2)
                 taskData = varargin{3};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            if(size(trainData,1) ~= size(taskData,1))
                warning('Output data and task data has to have same number of patterns (rows) order to calculate classification error.');
            end
            
            % compute task mask for MSE computation
            
            NFeatures = size(outData,2) / this.NPathways;
            
            MSE_tasks= nan(size(taskData));
            
            % compute mean MSE for invididual task being performed at a
            % given time
            for patternIdx = 1:size(taskData,1)
                
                tasks = find(taskData(patternIdx, :));
                for taskIdx = 1:length(tasks)
                    
                    task = tasks(taskIdx);
                    outputDim = mod(task-1,this.NPathways)+1;
                    relevantFeatures = (NFeatures*(outputDim-1)+1) : (NFeatures*outputDim);
                    MSE_tasks(patternIdx,task) = mean((outData(patternIdx, relevantFeatures) - trainData(patternIdx, relevantFeatures)).^2);

                end
                
            end
            
            % compute mean MSE for performed tasks
            MSE_tasks_mean = nanmean(MSE_tasks,2);
            
        end
      
        % computes absolute error only for tasks that are about to be executed
        % (ignores output dimensions that are not relevant)
        function [AE_tasks_mean, AE_tasks] = calculateAbsoluteErrorTasks(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            taskData = this.taskSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(length(varargin) > 2)
                 taskData = varargin{3};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            if(size(trainData,1) ~= size(taskData,1))
                warning('Output data and task data has to have same number of patterns (rows) order to calculate classification error.');
            end
            
            % compute task mask for MSE computation
            
            NFeatures = size(outData,2) / this.NPathways;
            
            AE_tasks= nan(size(taskData));
            
            % compute mean MSE for invididual task being performed at a
            % given time
            for patternIdx = 1:size(taskData,1)
                
                tasks = find(taskData(patternIdx, :));
                for taskIdx = 1:length(tasks)
                    
                    task = tasks(taskIdx);
                    outputDim = mod(task-1,this.NPathways)+1;
                    relevantFeatures = (NFeatures*(outputDim-1)+1) : (NFeatures*outputDim);
                    AE_tasks(patternIdx,task) = 1 - mean(abs(outData(patternIdx, relevantFeatures) - trainData(patternIdx, relevantFeatures)));

                end
                
            end
            
            % compute mean MSE for performed tasks
            AE_tasks_mean = nanmean(AE_tasks,2);
            
        end
      
        % computes likelihood for each response for tasks that are about to
        % be executed by appying the softmax function to the output patterns
        % (ignores output dimensions that are not relevant). Also computes
        % general likelihood of correct response
        function [outcomeProbs, PCorrect_tasks_mean, PCorrect_tasks] = calculateSoftmaxProbabilitiesTasks(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            taskData = this.taskSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(length(varargin) > 2)
                 taskData = varargin{3};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            if(size(trainData,1) ~= size(taskData,1))
                warning('Output data and task data has to have same number of patterns (rows) order to calculate classification error.');
            end
            
            % compute task mask for MSE computation
            
            NFeatures = size(outData,2) / this.NPathways;
            
            outcomeProbs =nan(size(outData));
            PCorrect_tasks= nan(size(taskData));
            PCorrect_tasks_mean= nan(size(PCorrect_tasks,1),1);
            
            % compute mean MSE for invididual task being performed at a
            % given time
            for patternIdx = 1:size(taskData,1)
                
                tasks = find(taskData(patternIdx, :));
                for taskIdx = 1:length(tasks)
                    
                    task = tasks(taskIdx);
                    
                    outputDim = mod(task-1,this.NPathways)+1;
                    relevantFeatures = (NFeatures*(outputDim-1)+1) : (NFeatures*outputDim);
                    correctResponse = relevantFeatures(find(trainData(patternIdx, relevantFeatures) == 1,1));
                    
                    % apply softmax
                    outcomeProbs(patternIdx,relevantFeatures) = exp(outData(patternIdx, relevantFeatures))./sum(exp(outData(patternIdx, relevantFeatures)));
                    
                    % compute correct response probability
                    PCorrect_tasks(patternIdx, task) = outcomeProbs(patternIdx, correctResponse);

                end
                
            end
            
            % compute mean MSE for performed tasks
            PCorrect_tasks_mean = nanmean(PCorrect_tasks,2);
            
        end
      
        % computes likelihood for each response for tasks that are about to
        % be executed.
        % (ignores output dimensions that are not relevant). Also computes
        % general likelihood of correct response
        function [outcomeProbs, PCorrect_tasks_mean, PCorrect_tasks] = calculateOutcomeProbabilitiesTasks(this, varargin)
            
            outData = this.output_log;
            trainData = this.trainSet;
            taskData = this.taskSet;
            
            if(~isempty(varargin))
                 outData = varargin{1};
            end
            
            if(length(varargin) > 1)
                 trainData = varargin{2};
            end
            
            if(length(varargin) > 2)
                 taskData = varargin{3};
            end
            
            if(~isequal(size(trainData), size(outData)))
                warning('Output data and training data has to have same dimensions in order to calculate classification error.');
            end
            
            if(size(trainData,1) ~= size(taskData,1))
                warning('Output data and task data has to have same number of patterns (rows) order to calculate classification error.');
            end
            
            % compute task mask for MSE computation
            
            NFeatures = size(outData,2) / this.NPathways;
            
            outcomeProbs =nan(size(outData));
            PCorrect_tasks= nan(size(taskData));
            PCorrect_tasks_mean= nan(size(PCorrect_tasks,1),1);
            
            % compute mean MSE for invididual task being performed at a
            % given time
            for patternIdx = 1:size(taskData,1)
                
                tasks = find(taskData(patternIdx, :));
                for taskIdx = 1:length(tasks)
                    
                    task = tasks(taskIdx);
                    
                    outputDim = mod(task-1,this.NPathways)+1;
                    relevantFeatures = (NFeatures*(outputDim-1)+1) : (NFeatures*outputDim);
                    correctResponse = relevantFeatures(find(trainData(patternIdx, relevantFeatures) == 1,1));
                    
                    % normalize
                    outcomeProbs(patternIdx,relevantFeatures) = outData(patternIdx, relevantFeatures)./sum(outData(patternIdx, relevantFeatures));
                    
                    % compute correct response probability
                    PCorrect_tasks(patternIdx, task) = outcomeProbs(patternIdx, correctResponse);

                end
                
            end
            
            % compute mean MSE for performed tasks
            PCorrect_tasks_mean = nanmean(PCorrect_tasks,2);
            
        end
      
        function [output_act_log, hidden_act_log, MSE, MSE_tasks, AE_tasks_mean, AE_tasks, outcomeProbs, PCorrect_tasks_mean, PCorrect_tasks] = integrateNetInput(this, tau, inputStream, taskStream, varargin)
            
            error('The method integrateNetInput has not been implemented yet for class NNmodelDeep');
            % check inputs
            
            if(length(varargin) >=1)
                trainStream = varargin{1};
                
                if(~isequal(size(inputStream,1), size(trainStream,1)) || ~isequal(size(inputStream,3), size(trainStream,3)))
                    error('Number of input and training patterns need to match.'); 
                end
            else
                trainStream = [];
            end
             
            if(~isequal(size(inputStream,1), size(taskStream,1)) || ~isequal(size(inputStream,3), size(taskStream,3)))
               error('Number of input and task patterns need to match.'); 
            end
            
            % prepare log data
            Nsets = size(inputStream,1);
            iterations = size(inputStream,3);
            output_net_log = nan(Nsets, iterations, this.Noutput);
            output_act_log = nan(Nsets, iterations, this.Noutput);
            hidden_net_log = nan(Nsets, iterations, this.Nhidden);
            hidden_act_log = nan(Nsets, iterations, this.Nhidden);
            MSE = nan(Nsets, iterations);
            MSE_tasks = nan(Nsets, size(taskStream,2), iterations);
            AE_tasks_mean= nan(Nsets, iterations);
            AE_tasks= nan(Nsets, size(taskStream,2), iterations);
            outcomeProbs = nan(Nsets, size(trainStream,2), iterations);
            PCorrect_tasks_mean= nan(Nsets, iterations);
            PCorrect_tasks= nan(Nsets, size(taskStream,2), iterations);
            
            % set up initial condition (corresponds to first patterns in the set)
            [output_act_curr, hidden_act_curr, ~, hidden_net_curr, output_net_curr] = runSet(this, inputStream(:,:,1), taskStream(:,:,1), trainStream(:,:,1));

            output_net_log(:,1,:) = output_net_curr;
            output_act_log(:,1,:) = output_act_curr;
            hidden_net_log(:,1,:) = hidden_net_curr;
            hidden_act_log(:,1,:) = hidden_act_curr;
            
            % compute error measures if training patterns are available
            if(~isempty(trainStream))
                
                this.taskSet = taskStream(:,:,1);
                this.trainSet = trainStream(:,:,1);
            
                % compute error measures
                [MSE_tasks_mean_curr, MSE_tasks_curr] = this.calculateMeanSquaredErrorTasks();
                [AE_tasks_mean_curr, AE_tasks_curr] = this.calculateAbsoluteErrorTasks();
                [outcomeProbs_curr, PCorrect_tasks_mean_curr, PCorrect_tasks_curr] = this.calculateOutcomeProbabilitiesTasks();
                
                % store error measures
                MSE(:,1) = MSE_tasks_mean_curr;
                MSE_tasks(:,:,1) = MSE_tasks_curr;
                AE_tasks_mean(:,1) = AE_tasks_mean_curr;
                AE_tasks(:,:,1) = AE_tasks_curr;
                outcomeProbs(:,:,1) = outcomeProbs_curr;
                PCorrect_tasks_mean(:,1) = PCorrect_tasks_mean_curr;
                PCorrect_tasks(:,:,1) = PCorrect_tasks_curr;
                
            end
            
            for i = 2:iterations
                
                inputData = inputStream(:,:,i);
                taskData = taskStream(:,:,i);
                
                % calculate net inputs for hidden layer
                hidden_net_input = this.weights.W_IH * transpose(inputData);                                    % input from input layer (stimulus)       
                hidden_net_task = this.g(this.weights.W_TH) * transpose(taskData);                              % input from task layer (task cue)
                hidden_net_bias  = this.weights.W_BH * (this.hidden_bias * ones(1,size(hidden_net_input,2)));   % input from hidden bias units
                hidden_net = hidden_net_input + hidden_net_task + hidden_net_bias;                              % total net input to hidden layer
                hidden_net_log(:,i,:) = tau .* transpose(hidden_net) + (1-tau) .* squeeze(hidden_net_log(:,i-1,:)); % integrate input from previous time step and current time step
                hidden_net = squeeze(hidden_net_log(:,i,:));

                % calculate activation for hidden units
                this.hidden_act = 1./(1+exp(-hidden_net));                          % use sigmoid activation function
                hidden_act = this.hidden_act;
                hidden_act_log(:,i,:) = this.hidden_act; 

                % calculate net input for output units
                output_net_task = this.g(this.weights.W_TO) * transpose(taskData);                              % input from task layer (task cue)
                output_net_hidden = this.weights.W_HO * transpose(this.hidden_act);                                        % input from hidden layer
                output_net_bias   = this.weights.W_BO * (this.output_bias * ones(1,size(hidden_net_input,2)));  % input from output bias units
                output_net = output_net_hidden + output_net_task + output_net_bias;                             % total net input to output layer
                output_net_log(:,i,:) = tau*transpose(output_net) + (1-tau) * squeeze(output_net_log(:,i-1,:)); % integrate input from previous time step and current time step
                output_net = squeeze(output_net_log(:,i,:));
                    
                % calculate activation of output units
                this.output_act = 1./(1+exp(-output_net)); 
                output_act = this.output_act;
                output_act_log(:,i,:) = this.output_act; 
                

                 % compute error measures if training patterns are available
                    if(~isempty(trainStream))
                        this.output_log = this.output_act; 
                        this.taskSet = taskStream(:,:,i);
                        this.trainSet = trainStream(:,:,i);

                        % compute error measures
                        [MSE_tasks_mean_curr, MSE_tasks_curr] = this.calculateMeanSquaredErrorTasks();
                        [AE_tasks_mean_curr, AE_tasks_curr] = this.calculateAbsoluteErrorTasks();
                        [outcomeProbs_curr, PCorrect_tasks_mean_curr, PCorrect_tasks_curr] = this.calculateOutcomeProbabilitiesTasks();

                        % store error measures
                        MSE(:,i) = MSE_tasks_mean_curr;
                        MSE_tasks(:,:,i) = MSE_tasks_curr;
                        AE_tasks_mean(:,i) = AE_tasks_mean_curr;
                        AE_tasks(:,:,i) = AE_tasks_curr;
                        outcomeProbs(:,:,i) = outcomeProbs_curr;
                        PCorrect_tasks_mean(:,i) = PCorrect_tasks_mean_curr;
                        PCorrect_tasks(:,:,i) = PCorrect_tasks_curr;

                    end

            end
            
            
        end
        
        % calculate switch dynamics between inputs based on delayed net
        % input
        function [output_act_log, hidden_act_log, MSE, DDM_RTs, DDM_ERs, E, A] = switchTasks(this, tau, iterations, inputA, inputB, tasksA, tasksB, varargin)
            
            if(length(varargin) >=1)
                trainA = varargin{1};
            else
                trainA = [];
            end
            
            if(length(varargin) >=2)
                trainB = varargin{2};
            else
                trainB = [];
            end
            
            if(length(varargin) >=3)
                ddmp = varargin{3};
                runAccumulator = 1;
            else
                ddmp = [];
                runAccumulator = 0;
            end
            
            if(iterations < 2)
                warning('Number of iterations should be at least 2 in order to simulate a task switch.');
            end
            
            if(~isequal(size(inputA), size(inputB)))
               error('Dimensions of input patterns need to match.'); 
            end
            
            if(~isequal(size(tasksA), size(tasksB)))
               error('Dimensions of task patterns need to match.'); 
            end
            
            if(~isequal(size(inputA,1), size(tasksA,1)))
               error('Number of input and task patterns need to match.'); 
            end
            
            if(~isequal(size(inputA,1), size(trainB,1)))
               error('Number of input patterns and correct output patterns need to match.'); 
            end
            
            % prepare log data
            Nsets = size(inputA,1);
            output_net_log = nan(Nsets, iterations, this.Noutput);
            output_act_log = nan(Nsets, iterations, this.Noutput);
            hidden_net_log = nan(Nsets, iterations, this.Nhidden);
            hidden_act_log = nan(Nsets, iterations, this.Nhidden);
            MSE = nan(Nsets, iterations);
            
            % get first data
            [output_act_A, hidden_act_A, MSE_A, hidden_net_A, output_net_A] = runSet(this, inputA, tasksA, trainA);

            output_net_log(:,1,:) = output_net_A;
            output_act_log(:,1,:) = output_act_A;
            hidden_net_log(:,1,:) = hidden_net_A;
            hidden_act_log(:,1,:) = hidden_act_A;

            MSE(:,1) = MSE_A;
                
            for set = 1:Nsets
                
                input = inputB(set,:);
                task = tasksB(set,:);
                train = trainB(set,:);
                
                % loop through time steps
                for i = 2:iterations

                    % calculate net inputs for hidden layer
                    hidden_net_input = this.weights.W_IH * transpose(input);            % input from input layer (stimulus)       
                    hidden_net_task = this.g(this.weights.W_TH) * transpose(task);      % input from task layer (task cue)
                    hidden_net_bias  = this.weights.W_BH * this.hidden_bias;            % input from hidden bias units
                    hidden_net = hidden_net_input + hidden_net_task + hidden_net_bias;  % total net input to hidden layer
                    hidden_net_log(set,i,:) = tau*hidden_net + (1-tau) * squeeze(hidden_net_log(set,i-1,:)); % integrate input from previous time step and current time step
                    hidden_net = squeeze(hidden_net_log(set,i,:));
                    
                    % calculate activation for hidden units
                    this.hidden_act = 1./(1+exp(-hidden_net));                          % use sigmoid activation function
                    hidden_act = this.hidden_act;
                    hidden_act_log(set,i,:) = this.hidden_act; 
                    
                    % calculate net input for output units
                    output_net_task = this.g(this.weights.W_TO) * transpose(task);      % input from task layer (task cue)
                    output_net_hidden = this.weights.W_HO * this.hidden_act;            % input from hidden layer
                    output_net_bias   = this.weights.W_BO * this.output_bias;           % input from output bias units
                    output_net = output_net_hidden + output_net_task + output_net_bias; % total net input to output layer
                    output_net_log(set,i,:) = tau*output_net + (1-tau) * squeeze(output_net_log(set,i-1,:)); % integrate input from previous time step and current time step
                    output_net = squeeze(output_net_log(set,i,:));
                    
                    % calculate activation of output units
                    this.output_act = 1./(1+exp(-output_net)); 
                    output_act = this.output_act;
                    output_act_log(set,i,:) = this.output_act'; 
                    
                    
                    % calculate MSE
                    if(~isempty(train))
                        MSE(set,i) = sum((output_act' - train).^2)./this.Noutput;
                    else
                        MSE(set,i) = 0; 
                    end
                
                end

            end
            
            % evidence accumulator 
            if(runAccumulator)
                  
                % drift, e.g. d = 0.1
                if(isfield(ddmp, 'd'))
                    d = ddmp.d;
                else
                    d = 0.1;
                end

                % threshold, e.g. z = 1
                if(isfield(ddmp, 'z'))
                    z = ddmp.z;
                else
                    z = 1;
                end

                % sigma, e.g. c = 0.1
                if(isfield(ddmp, 'c'))
                    c = ddmp.c;
                else
                    c = 0.1;
                end

                % wait time until recording
                if(isfield(ddmp, 'runs'))
                    respOnset = ddmp.respOnset;
                else
                    respOnset = 2;
                end

                % number of simulations
                if(isfield(ddmp, 'runs'))
                    nSimulations = ddmp.nSimulations;
                else
                    nSimulations = 10;
                end
                
                % create evidence matrix
                [output_sorted output_sorted_idx] = sort(output_act_log,3,'descend');
                
                % currently largest output value
                output_max = squeeze(output_sorted(:,:,1));
                output_max_idx = squeeze(output_sorted_idx(:,:,1));
                
                % currently second largest output value
                output_max2 = squeeze(output_sorted(:,:,2));
                output_max_idx2 = squeeze(output_sorted_idx(:,:,2));
                
                % build evidence matrix
                E = output_act_log;
%                 for outputIdx = 1:size(E, 3)
%                    maxMask = output_max;
%                    maxMask(output_max_idx == outputIdx) = output_max2(output_max_idx == outputIdx);
%                    E(:,:,outputIdx) = E(:,:,outputIdx)  - maxMask;
%                 end
                
                % create accumulation matrix
                A = repmat(E, [ 1, 1, 1, nSimulations]);
                A = normrnd(A.*d, c, size(A));      % scale evidence by drift rate and add noise
                A = cumsum(A, 2);   % accumulate evidence over time
            
                % performance matrices for reaction time and error rate
                DDM_RTs = nan(Nsets, nSimulations, length(z));
                DDM_ERs = nan(Nsets, nSimulations, length(z));
                
                % if training data is provided figure out correct response
                % unit for each pattern
                if(~isempty(trainA) && ~isempty(trainB))
                    C = repmat(trainB, [1, 1, nSimulations]);
                end
                
                start_z = 1; % threshold to search from
                % look at each time point
                for t = respOnset:iterations 
                    
                    % get snapshot
                    S = squeeze(A(:,t,:,:));
                    
                    for z_idx = start_z:length(z)
                        
                        % mark accumulators that are past the current threshold
                        threshMask = S > z(z_idx);
                        S_thresholded = squeeze(sum(threshMask,2)) >= 1;
                        
                        % calculate ER
                        accuracy = threshMask == C;
                        S_accuracy = squeeze(sum(accuracy, 2) == size(accuracy,2));
                        DDM_ERs_zIdx = DDM_ERs(:,:,z_idx);
                        pastThreshold = isnan(DDM_ERs(:,:,z_idx)) & S_thresholded;
                        DDM_ERs_zIdx(pastThreshold) = S_accuracy(pastThreshold);
                        DDM_ERs(:,:,z_idx) = DDM_ERs_zIdx;
                        
                        
                        % if accumulator hasn't passed the current threshold yet but just crossed threshold then assign time index as RT
                        DDM_RTs_zIdx = DDM_RTs(:,:,z_idx);
                        pastThreshold = isnan(DDM_RTs(:,:,z_idx)) & S_thresholded;
                        
                        
                        DDM_RTs_zIdx(pastThreshold & S_accuracy) = t;
                        DDM_RTs_zIdx(pastThreshold & (1-S_accuracy)) = -t;
                        DDM_RTs(:,:,z_idx) = DDM_RTs_zIdx;                    
                        

                        % if all accumulators passed the threshold, then remove threshold from list to improve performance
                        if(isequal(S_thresholded, ones(size(S_thresholded))))
                            start_z = z_idx + 1;
                        end
                        
                        
                    end
                    
                    
                end 
                
            end
            
        end
        
        % generate output given hidden layer activity and task input to
        % output layer
        function [output_act, hidden_act, MSE] = generateOutput(this, hidden_act, varargin)
            
            % parse arguments: task patterns
            if(length(varargin) >=1)
                task = varargin{1};
            else
                task = this.taskSet;
            end
            
            % calculate net input for output units
            output_net_task = this.g(this.weights.W_TO) * transpose(task);      % input from task layer (task cue)
            output_net_hidden = this.weights.W_HO * transpose(hidden_act);            % input from hidden layer
            output_net_bias   = this.weights.W_BO * this.output_bias;           % input from output bias units
            output_net = output_net_hidden + output_net_task + repmat(output_net_bias,1,size(output_net_hidden,2)); % total net input to output layer

            % calculate activation of output units
            output_act = 1./(1+exp(-output_net)); 

            % final network output
            output_act = output_act';
            
            % calculate MSE if correct output provided (train)
            MSE = -1;
            if(length(varargin) >=2)
                train = varargin{2};
                if(length(train) == length(output_act))
                    disp(output_act);
                    disp(train);
                    MSE = sum((output_act - train).^2)./this.Noutput;
                end
            end
            
        end
        
        % weight transformation function g(w) - used to calculate net input
        function weightsMod = g(this, weights)
            % weightsMod = abs(weights); % this weight transformation function ensures only positive weights
            weightsMod = weights;
        end
        
        % derivative of weight transformation function dg(w) - used to
        % calculate deltas for backpropagation
        function dWeights = dg(this, weights)
           %dWeights = sign(weights); % % this weight transformation function ensures only positive weights
           dWeights = 1;
        end
        
        function setData(this, inputData, taskData, trainData)
           this.inputSet = inputData;
           this.taskSet = taskData;
           this.trainSet = trainData;
        end

    end
    
end

