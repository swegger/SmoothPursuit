function [C, vecs, vals] = EyeTemporalCovariance(v,varargin)
%% EyeTemporalCovariance
%
%   C = EyeTemporalCovariance(v)
%
%   Measures the covariance in eye movement velocity over time.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'v')
addParameter(Parser,'Windows',NaN)
addParameter(Parser,'dim',1)

parse(Parser,v,varargin{:})

v = Parser.Results.v;
Windows = Parser.Results.Windows;
dim = Parser.Results.dim;

%% Calculate sample covariance, full data matrix
dims = 1:2;
dimObs = dims(dims~=dim);
C = cov(permute(v,[dimObs dim]));



%% Perform PCA
if any(isnan(Windows(:)))
    [vecs,vals] = ...
        eig(C);
else
    for wi = 1:size(Windows,1)
        [vecs{wi},vals{wi}] = ...
            eig(C(Windows(wi,1):Windows(wi,2),Windows(wi,1):Windows(wi,2)));
    end
end