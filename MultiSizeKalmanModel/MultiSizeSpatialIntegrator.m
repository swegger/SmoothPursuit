function zstar = MultiSizeSpatialIntegrator(N,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'N')
addParameter(Parser,'trials',100)
addParameter(Parser,'z',10*ones(1,200))
addParameter(Parser,'H',NaN)
addParameter(Parser,'Wx',NaN)
addParameter(Parser,'zhat0',NaN)
addParameter(Parser,'wT0',NaN)
addParameter(Parser,'Wt',NaN)

parse(Parser,N,varargin{:})

N = Parser.Results.N;
trials = Parser.Results.trials;
z = Parser.Results.z;
H = Parser.Results.H;
Wx = Parser.Results.Wx;
zhat0 = Parser.Results.zhat0;
wT0 = Parser.Results.wT0;
Wt = Parser.Results.Wt;

if any(isnan(H))
    H = ones(N,1);
end

if any(isnan(Wx))
    Wx = ones(N,1);
end

if any(isnan(zhat0))
    zhat0 = zeros(N,1);
end

if any(isnan(wT0))
    wT0 = 1e6*ones(N,1);
end

if any(isnan(Wt))
    Wt = 0.5*ones(N,1);
end

%% Perform simulations

for triali = 1:trials
    x = MultiSizeSensorModel(z,H,Wx);
    zhat(:,:,triali) = MultiSizeTemporalIntegratorModel(x,zhat0,wT0,Wx,Wt);
end
zstar = mean(zhat,1);