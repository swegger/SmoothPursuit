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

saveOpts.On = false;
saveOpts.Figs = false;
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
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,...
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
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'sizeProps',sizeProps,...
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
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'sizeProps',sizeProps,...
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
Cov.separationLengthConstant = 0.3;
Cov.sigf = 0.0001;
% Without noise
saveOpts.location = [locationBase '_diffCorr040_gainNoiseOff_' datestr(now,'yyyymmdd')];
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts)
