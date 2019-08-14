function [n, M, rNN, e] = NeuralModel_v1(varargin)
%% NeuralModel_v1
%
%
%
%%

%% Defaults
thetas_default = linspace(-90,90,4);
speeds_default = 2.^(linspace(-1,8,20)); %cumprod(2*ones(1,6));

[Ts, Ss] = meshgrid(thetas_default,speeds_default);

tuning_default.theta.Amp = 10;
tuning_default.theta.pref = Ts(:);
tuning_default.theta.sig = 45;
tuning_default.speed.Amp = 10;
tuning_default.speed.pref = Ss(:);
tuning_default.speed.sig = 1;
tuning_default.speed.d = 0.1;
tuning_default.n0 = 1;
tuning_default.Cov.sigf = 0.36;
tuning_default.Cov.thetaLengthConstant = 0.3;
tuning_default.Cov.speedLengthConstant = 4.4;
tuning_default.Cov.alpha = 0;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'tuning',tuning_default)

parse(Parser,varargin{:})

thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
tuning = Parser.Results.tuning;

%% Decoder properties
[Ds,Ss] = meshgrid(thetas,speeds);
s = cat(3,Ds',Ss');

%% Simulate MT and then decode

[n, M, rNN, ~, tuning] = SimpleMT(thetas,speeds,'trialN',100,'tuning',tuning);
e = DecodeMT(n,tuning,s,'gainNoise',0.1);