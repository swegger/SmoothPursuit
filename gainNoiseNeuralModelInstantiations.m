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

%% Vanillia (e.g. Cosyne 2020 poster)
% Without gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,...
    'mymakeaxisflg',true)

% With gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,...
    'mymakeaxisflg',true)

%% Nonlinear integration
sizeProps.surround_weight = 0;
sizeProps.exponential = 1/1000;
sizeProps.threshold = 0.5;
% Without gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'sizeProps',sizeProps,...
    'mymakeaxisflg',true)

% With gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'sizeProps',sizeProps,...
    'mymakeaxisflg',true)

%% Surround suppression
sizeProps.surround_weight = 0.1;
sizeProps.exponential = 0;
sizeProps.threshold = 0.1;
% Without gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'sizeProps',sizeProps,...
    'mymakeaxisflg',true)

% With gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'sizeProps',sizeProps,...
    'mymakeaxisflg',true)

%% Independent effect of threshold
sizeProps.surround_weight = 0;
sizeProps.exponential = 1;
sizeProps.threshold = 0.5;
% Without gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'sizeProps',sizeProps,...
    'mymakeaxisflg',true)

% With gain noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'sizeProps',sizeProps,...
    'mymakeaxisflg',true)

%% Differential correlations (diffAlpha = 0.1)
Cov.diffAlpha = 0.1;
% Without noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',true)

% With noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'Cov',Cov,...
    'mymakeaxisflg',true)

%% Differential correlations (diffAlpha = 0.2)
Cov.diffAlpha = 0.2;
% Without noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',true)

% With noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'Cov',Cov,...
    'mymakeaxisflg',true)

%% Differential correlations (diffAlpha = 0.4)
Cov.diffAlpha = 0.4;
% Without noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
    'mymakeaxisflg',true)

% With noise
NeuralModel_v2('thetas',thetas,'speeds',speeds,'gainNoise',0.5,'Cov',Cov,...
    'mymakeaxisflg',true)