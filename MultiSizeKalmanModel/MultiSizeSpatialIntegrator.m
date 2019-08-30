function zstar = MultiSizeSpatialIntegrator(varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'trials',100)
addParameter(Parser,'z',10*ones(1,200))
addParameter(Parser,'H',ones(100,1))
addParameter(Parser,'Wx',ones(100,1))
addParameter(Parser,'zhat0',zeros(100,1))
addParameter(Parser,'wT0',1e6*ones(100,1))
addParameter(Parser,'Wt',0.5*ones(100,1))

parse(Parser,varargin{:})

trials = Parser.Results.trials;
z = Parser.Results.z;
H = Parser.Results.H;
Wx = Parser.Results.Wx;
zhat0 = Parser.Results.zhat0;
wT0 = Parser.Results.wT0;
Wt = Parser.Results.Wt;

%% Perform simulations

for triali = 1:trials
    x = MultiSizeSensorModel(z,H,Wx);
    zhat(:,:,triali) = MultiSizeTemporalIntegratorModel(x,zhat0,wT0,Wx,Wt);
end
zstar = mean(zhat,1);