    LCA_settings.dt_tau = 0.1; %0.1
    LCA_settings.ITI = 4.0; % 0.5
    LCA_settings.maxTimeSteps = 50; %50;        
    LCA_settings.numSimulations = 100;%100   
    LCA_settings.lambda = 0.4;                      % leakage
    LCA_settings.alpha = 0.2;                         % recurrent excitation
    LCA_settings.beta = 0.2;                           % inhibition
    LCA_settings.responseThreshold = [0:0.1:1.5];         % tested threshold    [0:0.1:1.5]; 
    LCA_settings.T0 = 0.15;                                         % non-decision time
    if(exist('NFeatures'))
        LCA_settings.W_ext = eye(NFeatures) * 0.5;           % weight of external input
    end
    LCA_settings.c = 0.1;                                       % noise org: 0.1
