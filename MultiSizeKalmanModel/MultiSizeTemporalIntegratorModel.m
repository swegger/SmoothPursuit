function zhat = MultiSizeTemporalIntegratorModel(x,zhat0,wT0,Wx,Wt)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'x')
addRequired(Parser,'zhat0')
addRequired(Parser,'wT0')
addRequired(Parser,'Wx')
addRequired(Parser,'Wt')

parse(Parser,x,zhat0,wT0,Wx,Wt)

x = Parser.Results.x;
zhat0 = Parser.Results.zhat0;
wT0 = Parser.Results.wT0;
Wx = Parser.Results.Wx;
Wt = Parser.Results.Wt;

%% Integrate noisy sensory inputs over time

% Preallocate
zhat = nan(size(x));
wT = nan(size(x));

% Initialize
zhat(:,1) = zhat0;
wT(:,1) = wT0;

% Estimate over time
for ti = 2:size(x,2)
    wT(:,ti) = ((wT(:,ti-1).^2+Wt.^2).*Wx.^2)./((wT(:,ti-1).^2+Wt.^2) + Wx.^2);
    K(:,:,ti) = diag( wT(:,ti).^2 ./ (wT(:,ti).^2 + Wx.^2) );
    zhat(:,ti) = zhat(:,ti-1) + K(:,:,ti)*(x(:,ti) - zhat(:,ti-1));
end