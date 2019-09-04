function [zstar, wT, zhat, K] = MultiSizeSpatialIntegrator(N,varargin)
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
addParameter(Parser,'t0',20)
addParameter(Parser,'Nmax',600)
addParameter(Parser,'gainNoise',0)

parse(Parser,N,varargin{:})

N = Parser.Results.N;
trials = Parser.Results.trials;
z = Parser.Results.z;
H = Parser.Results.H;
Wx = Parser.Results.Wx;
zhat0 = Parser.Results.zhat0;
wT0 = Parser.Results.wT0;
Wt = Parser.Results.Wt;
t0 = Parser.Results.t0;
Nmax = Parser.Results.Nmax;
gainNoise = Parser.Results.gainNoise;

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
    wT0 = 0.1*ones(N,1);
end

if any(isnan(Wt))
    Wt = 0*ones(N,1);
end

%% Perform simulations

for triali = 1:trials
    x = MultiSizeSensorModel(z,H,Wx);
    [zhat(:,:,triali), wT(:,:,triali), Ktemp] = ...
        MultiSizeTemporalIntegratorModel(...
        x,zhat0,wT0,Wx,Wt,'t0',t0);
    K(:,triali) = Ktemp(1,1,:);
    
    gain = sqrt(N)/sqrt(Nmax) + gainNoise*randn;
    zstar(:,:,triali) = gain*sum( ((1./wT(:,:,triali))./repmat(sum(1./wT(:,:,triali),1),[size(zhat,1),1,1])) .* zhat(:,:,triali),1);  % Mean, accounting for potential differences in reliabilities across temporal integrator units
end
% zstar = mean(zhat,1); % Stupid mean
wT = mean(wT,1);