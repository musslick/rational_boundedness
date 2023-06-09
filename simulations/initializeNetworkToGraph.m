function network = initializeNetworkToGraph(A, network, varargin)
% initializeNetworkToGraph(A, NNModel)
% 
% A ...adjacency matrix of bipartite graph (A)
% rows of A refer to input nodes
% columns of A refer to output nodes
%
% network ...instance of NNmodel
%
% initializeNetworkToGraph(A, NNModel, fixInputHiddenWeights, fixTaskHiddenWeights)
%
% optional boolean inputs fixInputHiddenWeights and fixTaskHiddenWeights
% specify if weights from input layer to hidden layer or task layer to
% hidden layer shall be fixed. E.g. for initializeNetworkToGraph(A, NNModel, 0, 1)
% the network will initialize and fix only weights from the task to the
% hidden layer while weights from the input to the hidden layer will be
% randmly initialized and can be modified through learning
%
%% preliminary steps

% parse input arguments
if(~isempty(varargin))
    fixInputHiddenWeights = varargin{1};
else
    fixInputHiddenWeights = 1;
end

if(length(varargin) >= 2)
    fixTaskHiddenWeights = varargin{2};
else
    fixTaskHiddenWeights = 1;
end

% determine number of input & output components
[inputComponents, outputComponents] = size(A);

NPathways = max(size(A));

% set network bias
network.bias_weight = -5;
network.hidden_bias = network.bias_weight;
network.output_bias = network.bias_weight;

% maximum activation value of a hidden unit
maxActivationValueByTask = 0.1;
maxInputValue = 1;

% calculate required weight magnitude to activate a hidden unit from a
% single input unit
% w_IH = (- log(1/maxActivationValue - 1) - network.bias_weight) / maxInputValue;
w_IH = 2;

% calculate required weight from task units to activate hidden unit that
% receives input from stimulus
w_TH = (- log(1/maxActivationValueByTask - 1) - network.bias_weight) / maxInputValue;

% quick check if we have the right amount of input units
NFeatures = network.Ninput / inputComponents;
if(round(NFeatures) ~= NFeatures)
    error('Number of input units does not match required number of input components (row number of A).');
end

% check if we do have enoug task units
if(network.Ntask < NPathways^2)
    error(['Network does not have enough task units for specified number of input and output dimensions. At least ' num2str(NPathways^2) ' task units are required.']);
end

% check if we do have enoug hidden units
if(network.Nhidden < NFeatures * inputComponents)
    error(['Network does not have enough hidden units for specified number of input components and features per input dimension. At least ' num2str(NFeatures*inputComponents) ' hidden units are required.']);
end

%% determine input-hidden weights

if(fixInputHiddenWeights)

    % set input weights
    subWeights = eye(NFeatures * inputComponents, NFeatures * inputComponents) * w_IH;

    % overwrite network weights
    network.weights.W_IH = zeros(network.Nhidden,network.Ninput);
    network.weights.W_IH(1:size(subWeights,1), 1:size(subWeights,2)) = subWeights;

    % fix weights
    network.setFixedWeights([network.getFixedWeights()  network.W_INPUT_HIDDEN]);

end

%% determine task-hidden weights

if(fixTaskHiddenWeights)

    network.weights.W_TH = zeros(network.Nhidden, network.Ntask);

    for hiddenCIdx = 1:inputComponents

        % determine which hidden units a task projects to based on common
        % hidden (input) component
        hiddenTaskRep = ((hiddenCIdx-1)*NFeatures + 1) :  ((hiddenCIdx-1)*NFeatures + NFeatures);

        for outputCIdx = 1:outputComponents

            % for every edge (task) in the bipartie graph, create a task
            % pathway
            if(A(hiddenCIdx, outputCIdx) > 0)

                % determine task ID
                task = (hiddenCIdx-1) * NPathways + outputCIdx;

                % adjust weights for task
                network.weights.W_TH(hiddenTaskRep, task) = w_TH;

            end

        end

    end

    % fix weights
    network.setFixedWeights([network.getFixedWeights() network.W_TASK_HIDDEN]);

end

end
