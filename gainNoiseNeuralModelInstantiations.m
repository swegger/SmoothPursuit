%% gainNoiseNeuralModelInstantiations
%
%   Implements neural circuit model in different ways to explore the
%   effects of different parameterizations on expression of gain noise
%   effect.
%
%%

%% Defaults
thetas = 0;
speeds = 4:4:20;
mymakeaxisflg = false;

Cov.sigf = 0.36;
Cov.thetaLengthConstant = 0.4;
Cov.speedLengthConstant = 0.3;
Cov.separationLengthConstant = 0.3;
Cov.alpha = 0;
Cov.diffAlpha = 0;

sizeProps.minEccentricity = 1;
sizeProps.maxEccentricity = 30;

saveOpts.On = true;
saveOpts.Figs = true;
locationBase = '~/Projects/MultiSizePursuit/Circuit/Results/';

mymakeaxisflg = false;

%% Vanillia (e.g. Cosyne 2020 poster)
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

% Without gain noise
saveOpts.location = [locationBase '_vanilla_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_vanilla_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.4,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Nonlinear integration
sizeProps.surround_weight = 0;
sizeProps.exponential = 1/1000;
sizeProps.threshold = 0.5;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

% Without gain noise
saveOpts.location = [locationBase '_nonlinear_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'sizeProps',sizeProps,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_nonlinear_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.3,'sizeProps',sizeProps,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Surround suppression
sizeProps.surround_weight = 1;
sizeProps.exponential = 1;
sizeProps.threshold = 0.1;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

% Without gain noise
saveOpts.location = [locationBase '_sursup_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'sizeProps',sizeProps,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_sursup_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.3,'sizeProps',sizeProps,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Independent effect of threshold
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0.5;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

% Without gain noise
saveOpts.location = [locationBase '_threshold_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'sizeProps',sizeProps,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_threshold_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'sizeProps',sizeProps,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Differential correlations (diffAlpha = 0.1)
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0.1;
Cov.separationLengthConstant = 0.3;
% Without noise
saveOpts.location = [locationBase '_diffCorr010_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With noise
saveOpts.location = [locationBase '_diffCorr010_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Differential correlations (diffAlpha = 0.2)
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0.2;
Cov.separationLengthConstant = 0.3;
% Without noise
saveOpts.location = [locationBase '_diffCorr020_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With noise
saveOpts.location = [locationBase '_diffCorr020_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Differential correlations (diffAlpha = 0.4)
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0.4;
Cov.separationLengthConstant = 0.3;
% Without noise
saveOpts.location = [locationBase '_diffCorr040_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With noise
saveOpts.location = [locationBase '_diffCorr040_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Reduced correlation between neurons over RF location
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.03;
% Without noise
saveOpts.location = [locationBase '_sepCons003_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% Vanillia w/ motor noise
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;
motorNoise = 0.1;

% Without gain noise
saveOpts.location = [locationBase 'motorNoise_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts,'motorNoise',motorNoise)

% With gain noise
saveOpts.location = [locationBase 'motorNoise_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.4,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts,'motorNoise',motorNoise)


%% Vanlilla w/ different tuning widths
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-90,90,1800];
thetaTuning.amplitudeRange = [1,20,1000];
thetaTuning.widthRange = [10,60,1000];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [1,20,1000];
speedTuning.widthRange = [1,2,1000];
speedTuning.d = 0.1;

% Without gain noise
saveOpts.location = [locationBase '_heterogeneous_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'theta',thetaTuning,'speed',speedTuning,'N',10240,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts);

% With gain noise
saveOpts.location = [locationBase '_heterogeneous_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.4,...
    'theta',thetaTuning,'speed',speedTuning,'N',10240,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts);

%% Vanlilla w/ dir tuning from -180 to 180
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [10,10,1];
thetaTuning.widthRange = [45,45,1];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [10,10,1];
speedTuning.widthRange = [1.64,1.64,1];
speedTuning.d = 0.1;

% Without gain noise
saveOpts.location = [locationBase '_d360_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'theta',thetaTuning,'speed',speedTuning,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_d360_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.4,...
    'theta',thetaTuning,'speed',speedTuning,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% g*log2shat w/ dir tuning from -180 to 180
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [1,20,1000];
thetaTuning.widthRange = [10,60,1000];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [1,20,1000];
speedTuning.widthRange = [1,2,1000];
speedTuning.d = 0.1;

% Without gain noise
saveOpts.location = [locationBase '_g*log2shat_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','g*log2shat',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_g*log2shat_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.13,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','g*log2shat',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% g*2^shat w/ dir tuning from -180 to 180
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [10,10,1];
thetaTuning.widthRange = [45,45,1];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [10,10,1];
speedTuning.widthRange = [1.64,1.64,1];
speedTuning.d = 0.1;

% Without gain noise
saveOpts.location = [locationBase '_g*2^shat_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','g*2^shat',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_g*2^shat_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.13,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','g*2^shat',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% simpleSumVA^2 w/ dir tuning from -180 to 180
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [10,10,1];
thetaTuning.widthRange = [45,45,1];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [10,10,1];
speedTuning.widthRange = [1.64,1.64,1];
speedTuning.d = 0.1;

% Without gain noise
saveOpts.location = [locationBase '_simpleSumVA^2_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','simpleSumVA^2',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_simpleSumVA^2_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.4,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','simpleSumVA^2',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

%% simpleSumVA w/ dir tuning from -180 to 180
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0;
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [10,10,1];
thetaTuning.widthRange = [45,45,1];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [10,10,1];
speedTuning.widthRange = [1.64,1.64,1];
speedTuning.d = 0.1;

% Without gain noise
saveOpts.location = [locationBase '_simpleSumVA_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','simpleSumVA',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)

% With gain noise
saveOpts.location = [locationBase '_simpleSumVA_gainNoiseOn_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.4,...
    'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','simpleSumVA',...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)