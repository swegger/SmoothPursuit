function [x, C, R] = simulateOrnsteinUhlenbeck(G,dt,input,eta,varargin)
%% simulateOrnsteinUhlenbeck
%
%
%%

%% Defaults

%% Parse inputs

Parser = inputParser;

addRequired(Parser,'G')
addRequired(Parser,'dt')
addRequired(Parser,'input')
addRequired(Parser,'eta')
addParameter(Parser,'x0',[0 0])

parse(Parser,G,dt,input,eta,varargin{:});

G = Parser.Results.G;
dt = Parser.Results.dt;
input = Parser.Results.input;
eta = Parser.Results.eta;
x0 = Parser.Results.x0;

if any(size(input) ~= size(eta))
    error('Size of input and eta must be identical')
end


%% Simulate process
T = size(input,2);
x = nan(size(input));
x(:,1) = x0(1) + x0(2)*randn(size(input,1),1);

for ti = 2:T
    dx(:,ti) = -(1-G(2))*x(:,ti-1) + G(1)*input(:,ti) + eta(:,ti);
    x(:,ti) = x(:,ti-1) + dx(:,ti)*dt;
end

%% Calculate covariance and correlation over time
C = cov(x);
R = corrcoef(x);