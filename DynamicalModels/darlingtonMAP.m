function r = darlingtonMAP(t,n,b,a,internal,varargin)
%% darlingtonMAP
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'t')
addRequired(Parser,'n')
addRequired(Parser,'b')
addRequired(Parser,'a')
addRequired(Parser,'internal')
addParameter(Parser,'r0',0)
addParameter(Parser,'dt',1/1000)
addParameter(Parser,'alpha',0.5)
addParameter(Parser,'beta',0.005)

parse(Parser,t,n,b,a,internal,varargin{:})

t = Parser.Results.t;
n = Parser.Results.n;
b = Parser.Results.b;
a = Parser.Results.a;
internal = Parser.Results.internal;
r0 = Parser.Results.r0;
dt = Parser.Results.dt;
alpha = Parser.Results.alpha;
beta = Parser.Results.beta;

%% Run dynamics
r = nan(length(t),1);
r(1) = r0;

for ti = 2:length(t)
    dr(ti) = alpha*(b(1:end-1)*permute(n(ti,:),[2,1]) + b(end)) - ...
        beta*(a(1:end-1)*permute(n(ti,:),[2,1]) + a(end) + internal(ti))*r(ti-1);
    r(ti) = r(ti-1) + dr(ti)*dt;
end