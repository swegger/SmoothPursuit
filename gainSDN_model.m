function [Ve, C] = gainSDN_model(Vs,trials,varargin)
%% gainSDN_model
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'Vs')
addRequired(Parser,'trials')
addParameter(Parser,'t',-100:200)
addParameter(Parser,'gainParams',[1.5,0.02])
addParameter(Parser,'SDNparams',[1e6,1,0])
addParameter(Parser,'dt',1/1000);

parse(Parser,Vs,trials,varargin{:})

Vs = Parser.Results.Vs;
trials = Parser.Results.trials;
t = Parser.Results.t;
gainParams = Parser.Results.gainParams;
SDNparams = Parser.Results.SDNparams;
dt = Parser.Results.dt;

%% Generate sensory noise

% SDN
wt = wT(SDNparams(1),SDNparams(2),2,sum(t>=0));

[T1,T2] = meshgrid(t(t>=0));
if SDNparams(3) > 0
    Sigma = exp( -(T1-T2).^2/(2*SDNparams(3).^2) );
    temp = length(wt)*ones(length(wt),length(wt));
    temp2 = fliplr(ceil((T2-T1)/2+temp/2));
    W = wt(temp2);
    Sigma = Vs^2*nearestSPD(W.^2.*Sigma);
else
    Sigma = Vs.^2*diag(wt.^2);
end
etaVs = mvnrnd(zeros(1,length(wt)),Sigma,trials);

%% Generate gain noise
etaK = repmat(gainParams(2)*randn(trials,1),[1,sum(t>=0)]);

%% Generate visual speed measurements over time
VsM = [randn(trials,sum(t<0)) Vs*ones(size(etaVs)) + etaVs];

% reversal trial
VsMr = [randn(trials,sum(t<0)),...
    Vs*ones(trials,floor(size(etaVs,2)/2)) + etaVs(:,1:floor(size(etaVs,2)/2)),...
    -Vs*ones(trials,ceil(size(etaVs,2)/2)) + etaVs(:,floor(size(etaVs,2)/2)+1:end)];

%% Generate weight profile over time
Kbar = wt/SDNparams(2);
K = [etaK(:,1:sum(t<0)) repmat(Kbar(:)',[trials,1]) + etaK];

%% Generate eye speeds over time
% Ve = eyeSpeeds(t,VsM,K,dt);

K2 = repmat(exp(-t(t>0)/50),[trials,1]) + etaK(:,1:sum(t>0));
Ve = eyeSpeeds_filter(t,VsM,K2);
Ver = eyeSpeeds_filter(t,VsMr,K2);

%% Calculate covariance
C = cov(Ve);

Cr = cov(Ver);

%% Functions

%% wmT
function wt = wT(w0,w,t0,T)

    wt(1) = w0;
    for ti = 2:(t0+T-1)
        wt(ti,1) = sqrt((wt(ti-1,1)^2*w^2) / (wt(ti-1,1)^2 + w^2));
    end
    wt = wt(t0:end);
    
%% eyeSpeeds
function Ve = eyeSpeeds(t,VsM,K,dt)
    
    Ve = zeros(size(VsM));
    for ti = 2:length(t)
        Ve(:,ti) = Ve(:,ti-1) + K(:,ti).*VsM(:,ti)*dt;
    end

%% eyeSpeeds_filter
function Ve = eyeSpeeds_filter(t,VsM,K)
    
    for triali = 1:size(VsM,1)
        Ve(triali,:) = filter(K(triali,:),1,VsM(triali,:));
    end
    