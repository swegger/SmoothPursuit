function [s, m, c] = SimulateGainNoiseModel(ws,G,sigma,wp,varargin)
%% SimulateGainNoiseModel
%
%   [s, m] = SimulateGainNoiseModel()
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'ws')
addRequired(Parser,'G')
addRequired(Parser,'sigma')
addRequired(Parser,'wp')
addParameter(Parser,'trials',1000)
addParameter(Parser,'svalues',linspace(4,20,5))

parse(Parser,ws,G,sigma,wp,varargin{:})

ws = Parser.Results.ws;
G = Parser.Results.G;
sigma = Parser.Results.sigma;
wp = Parser.Results.wp;
trials = Parser.Results.trials;
svalues = Parser.Results.svalues;

%% Simulate GainNoise model
s = randsample(svalues,trials,true);
s = s(:);
c = randsample(length(G),trials,true);
g = G(c);
g = g(:);
sm = s + ws*s.*randn(trials,1);
Gprime = g + sigma*randn(trials,1);
mp = Gprime.*sm;
m = mp + wp*mp.*randn(trials,1);