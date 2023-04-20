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

classdef NNModelLinear < handle
    
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
    
    methods
        
        % constructor
        function this = NNModelLinear(varargin)
            
            % make a copy of existing network object
            if(isa(varargin{1},'NNModelLinear'))
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
                if(isa(varargin{1},'NNModelLinear'))
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
                this.weights.W_IH = (-1 +2.*rand(this.Nhidden,this.Ninput))*0.1;%*this.init_scale;      % input-to-hidden weights
                this.weights.W_TH = (-1 +2.*rand(this.Nhidden,this.Ntask))*this.init_task_scale;
                this.weights.W_BH = ones(this.Nhidden,1);                                         % bias-to-hidden weights
                this.weights.W_TO = (-1 +2.*rand(this.Noutput,this.Ntask))*0.1;%*this.init_scale;              
                this.weights.W_HO = (-1 +2.*rand(this.Noutput,this.Nhidden))*0.1;%*this.init_scale;     % output-to-hidden weights
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
               this.hidden_log = zeros(Ndatasets,this.Nhidden);     % current log activity of hidden units for each dataset
               this.output_log = zeros(Ndatasets,this.Noutput);     % current log activity of output units for each dataset
               
               % loop through all the patterns (online training)
               for dset = 1:Ndatasets
                  [MSE(dset)] = trainTrial(this, inputData(dset,:),  taskData(dset,:), trainData(dset,:)); % trains weights on current pattern
                  this.hidden_log(dset,:) = this.hidden_act';       % log hidden unit activity for this pattern
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
               delta_output = error_term .* 1;

               % calculate delta's for hidden layer: delta_hidden = sum(delta_output * W_HO) * f_act'(netj)
               % delta_hidden ...backpropagated error for hidden units
               % delta_output ...backpropagated error for output units
               % W_HO         ...weights from hidden (columns) to output layer (rows)
               % f_act'(netj) ...first derivative of activation function of hidden units with respect to the net input
               delta_hidden = sum(repmat(delta_output,1,size(this.weights.W_HO,2)) .* this.weights.W_HO,1)' .* 1;
               
               delta_hiddenTask = delta_hidden;
               
               % if a pathway size for the hidden unit layer is specified, 
               % then the deltas for groups of hidden layer units that 
               % receive input from the task input layer will be averaged. 
               % The path size specifies the number of hidden units in a
               % group. Each task unit projects to all groups of hidden
               % units with the constraint that the projecting weights will 
               % be the same for the units within a group
               if(this.hiddenPathSize > 1) % no need to do averaging if hidden pathway size is 1
                   % average paths for hidden-to-task backprop
                   Npaths = this.Nhidden/this.hiddenPathSize;
                   refVec = repmat(1:Npaths,this.hiddenPathSize,1);
                   refVec = refVec(:);
                   
                   for i = 1:Npaths
                       delta_hiddenTask(refVec==i) = mean(delta_hidden(refVec==i));
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
               this.weights.W_IH = (this.weights.W_IH - this.coeff * delta_hidden * input) - this.coeff * this.decay * sign(this.weights.W_IH);
               % dW_IH = dW_IH +  this.coeff * delta_hidden * input - this.coeff * this.decay * sign(this.weights.W_IH);
               % adjust weights from task to hidden units
               this.weights.W_TH = (this.weights.W_TH - this.coeff * delta_hiddenTask * task) .* this.dg(this.weights.W_TH);
               % dW_TH = dW_TH + this.coeff * delta_hiddenTask * task .* this.dg(this.weights.W_TH);
               % adjust weights from task to output units
               this.weights.W_TO = this.weights.W_TO - (this.coeff * delta_outputTask * task) .* this.dg(this.weights.W_TO);
               % dW_TO = dW_TO + (this.coeff * delta_outputTask * task) .* this.dg(this.weights.W_TO);
               % adjust weights from hidden to output units
               this.weights.W_HO = this.weights.W_HO - (this.coeff * delta_output * this.hidden_act') - this.coeff * this.decay * sign(this.weights.W_HO);
               % dW_HO = dW_HO + (this.coeff * delta_output * this.hidden_act') - this.coeff * this.decay * sign(this.weights.W_HO);
               
               % this.decay * sign(this.weights.W_IH) ...penalize weights
               % this.dg(this.weights.W_TO)           ...derivative of transformation function of the weights
               
               
               % learning of bias weights is turned off for now
               % this.weights.W_BH = this.weights.W_BH - this.coeff * delta_hidden * this.hidden_bias; 
               % this.weights.W_BO = this.weights.W_BO - this.coeff * delta_output * this.output_bias;
               
               % calculate mean-squared error
               [~,  MSE] = this.calculateMeanSquaredError(this.output_act', train);

        end
        
         % train the network on all patterns
        function [] = trainBatch(this, iterations, varargin)
            
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
               delta_output = error_term .* 1;
               
               % calculate delta's for hidden layer: delta_hidden = sum(delta_output * W_HO) * f_act'(netj)
               % delta_hidden ...backpropagated error for hidden units
               % delta_output ...backpropagated error for output units
               % W_HO         ...weights from hidden (columns) to output layer (rows)
               
               % f_act'(netj) ...first derivative of activation function of hidden units with respect to the net input
               % delta_hidden = sum(repmat(delta_output,1,size(this.weights.W_HO,2)) .* this.weights.W_HO,1)' .* this.hidden_act .* (1 - this.hidden_act);
               delta_hidden = delta_output' * this.weights.W_HO .* 1;
               
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
               this.weights.W_IH = this.weights.W_IH - this.coeff * delta_hidden' * inputData - Ndatasets * this.coeff * this.decay * sign(this.weights.W_IH);
               % adjust weights from task to hidden units
               this.weights.W_TH = this.weights.W_TH - this.coeff * delta_hiddenTask' * taskData .* this.dg(this.weights.W_TH);

               % adjust weights from task to output units
               this.weights.W_TO = this.weights.W_TO - this.coeff * delta_outputTask' * taskData .* this.dg(this.weights.W_TO);
               % adjust weights from hidden to output units
               this.weights.W_HO = this.weights.W_HO - this.coeff * delta_output' * hiddenData - this.coeff * this.decay * sign(this.weights.W_HO);
               
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
            
            % calculate net inputs for hidden layer
            hidden_net_input = this.weights.W_IH * transpose(inputData);                                    % input from input layer (stimulus)       
            hidden_net_task = this.g(this.weights.W_TH) * transpose(taskData);                              % input from task layer (task cue)
            hidden_net_bias  = this.weights.W_BH * (this.hidden_bias * ones(1,size(hidden_net_input,2)));   % input from hidden bias units
            hidden_net = hidden_net_input + hidden_net_task + hidden_net_bias;                              % total net input to hidden layer

            % calculate activation for hidden units
            this.hidden_act = hidden_net;                                                      % use sigmoid activation function

            % calculate net input for output units
            output_net_task = this.g(this.weights.W_TO) * transpose(taskData);                              % input from task layer (task cue)
            output_net_hidden = this.weights.W_HO * this.hidden_act;                                        % input from hidden layer
            output_net_bias   = this.weights.W_BO * (this.output_bias * ones(1,size(hidden_net_input,2)));  % input from output bias units
            output_net = output_net_hidden + output_net_task + output_net_bias;                             % total net input to output layer

            % calculate activation of output units
            this.output_act = output_net; 

            % final network output
            hiddenData = this.hidden_act';                      % log activity of hidden units for each dataset
            outData = this.output_act';                         % log activity of output units for each dataset
            hidden_net = hidden_net';
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
            
            % calculate net inputs for hidden layer
            hidden_net_input = this.weights.W_IH * transpose(input);            % input from input layer (stimulus)       
            hidden_net_task = this.g(this.weights.W_TH) * transpose(task);      % input from task layer (task cue)
            hidden_net_bias  = this.weights.W_BH * this.hidden_bias;            % input from hidden bias units
            hidden_net = hidden_net_input + hidden_net_task + hidden_net_bias;  % total net input to hidden layer

            % calculate activation for hidden units
            this.hidden_act = hidden_net;                          % use sigmoid activation function

            % calculate net input for output units
            output_net_task = this.g(this.weights.W_TO) * transpose(task);      % input from task layer (task cue)
            output_net_hidden = this.weights.W_HO * this.hidden_act;            % input from hidden layer
            output_net_bias   = this.weights.W_BO * this.output_bias;           % input from output bias units
            output_net = output_net_hidden + output_net_task + output_net_bias; % total net input to output layer

            % calculate activation of output units
            this.output_act = output_net; 

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
        
        function [MSE_error_total, MSE_error_patterns] = calculateMeanSquaredError(this, varargin)
            
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
            
            MSE_error_patterns = sum((outData - trainData).^2, 2)./this.Noutput';
            MSE_error_total = mean(MSE_error_patterns);
            
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

