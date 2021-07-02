function e = simpleClosedLoop(varargin)
%%
%
%
%
%%

%% Defaults
t_default = 0:1300;
speed_default = 10*ones(size(t_default));
coh_default = [60*ones(1,401), 100*ones(1,300), 20*ones(1,300), 60*ones(1,300)];

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'t',t_default)
addParameter(Parser,'speed',speed_default)
addParameter(Parser,'coh',coh_default)
addParameter(Parser,'tau',40)
addParameter(Parser,'tauG',20)
addParameter(Parser,'e0',0)
addParameter(Parser,'G0',[0; 0])
addParameter(Parser,'lag',80)
addParameter(Parser,'W',[-1 0; 0 -1])
addParameter(Parser,'betas',[0,0])
addParameter(Parser,'alphas',[0 0])
addParameter(Parser,'plotFlg',true)

parse(Parser,varargin{:})

t = Parser.Results.t;
speed = Parser.Results.speed;
coh = Parser.Results.coh;
tau = Parser.Results.tau;
tauG = Parser.Results.tauG;
e0 = Parser.Results.e0;
G0 = Parser.Results.G0;
lag = Parser.Results.lag;
W = Parser.Results.W;
betas = Parser.Results.betas;
alphas = Parser.Results.alphas;
plotFlg = Parser.Results.plotFlg;

%% Eye speed dynamics
speed_obs = speed;% + randn(size(speed))./coh*10;
s = nan(size(speed));
e = nan(size(speed));
de = nan(size(speed));
G = nan(2,size(speed,2));
dG = nan(2,size(speed,2));
% G1 = nan(size(speed));
% G2 = nan(size(speed));
% dG1 = nan(size(speed));
% dG2 = nan(size(speed));
for ti = 1:length(t)
    if ti == 1
        G(:,ti) = G0;
%         G1(ti) = 0; 
%         G2(ti) = 0;
        s(ti) = speed_obs(ti);
        e(ti) = e0;
    elseif ti <= lag
        dG(:,ti) = [0; 0];
%         dG1(ti) = 0; 
        G(:,ti) = G(:,ti-1) + dG(:,ti)/tauG;
%         G1(ti) = G1(ti-1) + dG1(ti)/tauG;
%         dG2(ti) = 0;
%         G2(ti) = G2(ti-1) + dG2(ti)/tauG;
        s(ti) = speed_obs(ti);
        de(ti) = 0;
        e(ti) = 0;
    else
%         dG1(ti) = g1function(speed(ti-lag),coh(ti-lag),G1(ti-1),betas(1),G2(ti-1),alphas(1)); 
%         G1(ti) = G1(ti-1) + dG1(ti)/tauG;
%         dG2(ti) = g2function(e(ti-1),coh(ti-lag),G2(ti-1),betas(2),G1(ti-1),alphas(2));
%         G2(ti) = G2(ti-1) + dG2(ti)/tauG;
        s(ti) = speed_obs(ti-lag) - e(ti-1);
%         de(ti) = G1(ti)*s(ti);
%         e(ti) = G1(ti)*e(ti-1) + de(ti)/tau;
        dG(:,ti) = W*G(:,ti-1) + [betas(1)./(1 + exp(-0.05*(coh(ti-lag)-50))); 1-betas(2)+betas(2)./(1 + exp(-0.05*(coh(ti-lag)+10)))];
        G(:,ti) = G(:,ti-1) + dG(:,ti)/tauG;
%         dG1(ti) = -alphas(1)*G1(ti-1) + 1./(1 + exp(-0.05*(coh(ti-lag)-50))) - alphas(2)*G2(ti-1);
%         G1(ti) = G1(ti-1) + dG1(ti)/tauG;
%         dG2(ti) = -alphas(1)*G2(ti-1) + 1./(1 + exp(-0.05*(coh(ti-lag)-0))) - alphas(2)*G1(ti-1);
%         G2(ti) = G2(ti-1) + dG2(ti)/tauG;
        de(ti) = (G(2,ti)-1)*e(ti-1) + G(1,ti)*s(ti);
        e(ti) = e(ti-1) + de(ti)/tau;
    end
end

%% Plotting
if plotFlg
    figure
    subplot(5,1,1)
    plot(t,speed_obs)
    
    subplot(5,1,2)
    plot(t,e)
    
    subplot(5,1,3)
    plot(t,s)
    
    subplot(5,1,4)
    plot(t,G(1,:))
    
    subplot(5,1,5)
    plot(t,G(2,:))
end
%% Functions

%%
function dG1 = g1function(speed,coh,G1,beta,G2,alpha)
dG1 = -G1 + ((1-beta)*log(coh)/log(100) + beta*speed/10) - alpha*G2;
% dG1 = -G1 + 1./(1 + exp(-0.05*(coh-50)));

%%
function dG2 = g2function(e,coh,G2,beta,G1,alpha)
dG2 = -G2 + (1-beta)*log(coh)/log(100) + beta*sqrt(abs(e))/sqrt(10) - alpha*G1;