function [ws,sigGs,Gs] = gainNoiseNeuralModelParameterSweeps(simi,varargin)
%% gainNoiseNeuralModelParameterSweeps
%
%   [ws,sigGs,Gs] = gainNoiseNeuralModelParameterSweeps()
%
%   Sweeps through a default set of parameters and measures the resulting w
%   and sigG of circuit model output for comparison to behavior.
%
%%

%% Defaults
surround_weights_default = linspace(0,0.1,2);
thresholds_default = linspace(0,0.5,2);
exponentials_default = linspace(0.01,1,2);
gainNoises_default = linspace(0,1,2);

saveOpts_default.On = true;
saveOpts_default.location = ...
    ['/hpc/group/lisbergerlab/se138/Projects/MultiSizePursuit/Circuit/parameterSweeps/NeuralModel_v2_run_' datestr(now,'yyyymmdd') '_0'];
saveOpts_default.Figs = false;

thetas_default = 0;
speeds_default = 4:4:20;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'simi')
addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'surround_weights',surround_weights_default)
addParameter(Parser,'thresholds',thresholds_default)
addParameter(Parser,'exponentials',exponentials_default)
addParameter(Parser,'gainNoises',gainNoises_default)
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,simi,varargin{:})

simi = Parser.Results.simi;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
surround_weights = Parser.Results.surround_weights;
thresholds = Parser.Results.thresholds;
exponentials = Parser.Results.exponentials;
gainNoises = Parser.Results.gainNoises;
saveOpts = Parser.Results.saveOpts;

%% Set up parameter combinations
[SWs,Ts,Es,GNs] = ndgrid(surround_weights,thresholds,exponentials,gainNoises);

params = [SWs(:),Ts(:),Es(:),GNs(:)];

%% For each parameter combination, run simulation
sizeProps.minEccentricity = 1;
sizeProps.maxEccentricity = 30;
disp(['Sweep ' num2str(simi) ' of ' num2str(size(params,1))])
sizeProps.surround_weight = params(simi,1);
sizeProps.exponential = params(simi,3);
sizeProps.threshold = params(simi,2);

% With gain noise
[~, ~, ~, ~, ~, ws(simi), sigGs(simi), Gs(simi,:)] = NeuralModel_v2(...
    'thetas',thetas,'speeds',speeds,...
    'gainNoise',params(simi,4),'sizeProps',sizeProps,...
    'plotMT',false,'plotDecoding',false,'plotResults',false,...
    'saveOpts',saveOpts);
    
