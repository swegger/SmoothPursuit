function MinimumInterventionPolicy(varargin)
%% MinimumInterventionPolicy
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

parse(Parser,varargin{:})

%%


%% Functions

%% latentDynamics
function latentDynamics(t,x0,A,B,H,Gamma,Sigma,R)

    x(:,ti) = x0;
    for ti = 2:length(t)
        dx(:,ti) = A*x(:,ti-1) + B*u(:,ti-1) + Gamma*randn(size(x,1),1);
        x(:,ti) = x(:,ti-1) + dx(:,ti);
        m(:,ti) = H*x(:,ti) + Sigma*randn(size(x(:,ti),1),1);
        xhat(:,ti) = kalmanUpdate(xhat(:,ti-1),u(:,ti-1),A,B,H,Gamma,Sigma);
        u(:,ti) = controlPolicy(xhat,R);
    end
    
%% kalmanUpdate
function kalmanUpdate(xhat,u,A,B,H,Gamma,Sigma)

%% controlPolicy
function controlPolicy(xhat,R)