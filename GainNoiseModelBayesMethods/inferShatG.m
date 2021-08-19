function [shat, ghat] = inferShatG(e,g,sig,s,w,varargin)
%% inferShatG
%
%   [shat, ghat] = inferShatG(e,g,sig,s,w)
%       Infers the measured stimulus speed (shat) and the noisy gain (ghat)
%       based on the observed eye speed (e), mean gain (g), estimate of
%       gain noise (sig), actual stimulus speed (s), and the estimated 
%       Weber fraction for speed perception (w).
%
%       The algorithm assumes that the observed eye speed results from
%           e = (g + eta_g)(s + eta_s)
%       Where eta_g is noise drawn from a zero mean Gaussian with variance
%       sig^2 and eta_s is drawn from a zero mean Gaussian with variance
%       w^2s^2 (see Egger and Lisberger, 2020). Therefore,
%           ghat = g + eta_g,
%           shat = s + eta_s,
%       and
%           e = ghat*shat.
%
%       Algorithm makes use of the assumption that motor noise is
%       negligable to search the likelihhood function, p(ghat|g)*p(shat|s)
%       along only those values that are consistent with e = ghat*shat
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'e')
addRequired(Parser,'g')
addRequired(Parser,'sig')
addRequired(Parser,'s')
addRequired(Parser,'w')

parse(Parser,e,g,sig,s,w,varargin{:})

e = Parser.Results.e;
g = Parser.Results.g;
sig = Parser.Results.sig;
s = Parser.Results.s;
w = Parser.Results.w;

%% Infer G and shat from eye speed
minfun = @(p)( -logProbG(e./p,g,sig)-logProbS(p,s,w) );
shat = fminsearch(minfun,e);
ghat = e./shat;

%% Functions

%% probG
function p_g = logProbG(ghat,g,sig)
    p_g = -log(sqrt(2*pi*sig.^2)) - (ghat-g).^2./(2*sig^2);
    
%% probS
function p_s = logProbS(sm,s,w)
    p_s = -log(sqrt(2*pi*w^2*s^2)) - (sm-s).^2./(2*w^2*s^2);