function simpleClosedLoop_dynCoh_initCoh(varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'W',[-1 0; 0 -1])
addParameter(Parser,'betas',[1,1])
addParameter(Parser,'alphas',[0 0])
addParameter(Parser,'tauG',20)
addParameter(Parser,'G0',[0; 0])

parse(Parser,varargin{:})

W = Parser.Results.W;
betas = Parser.Results.betas;
alphas = Parser.Results.alphas;
tauG = Parser.Results.tauG;
G0 = Parser.Results.G0;

%% Dynamic Coherence
blockT = [450,300,300,300];
cT = cumsum(blockT);
dynCoh.t = 0:sum(blockT)-1;
dynCoh.coh = 60*ones(5,length(dynCoh.t));
dynCoh.coh(1,cT(1):cT(2)-1) = 100;
dynCoh.coh(1,cT(2):cT(3)-1) = 20;
dynCoh.coh(2,cT(1):cT(2)-1) = 20;
dynCoh.coh(2,cT(2):cT(3)-1) = 100;
dynCoh.coh(3,cT(1):cT(3)-1) = 100;
dynCoh.coh(4,cT(1):cT(3)-1) = 20;
dynCoh.speed = 10*ones(size(dynCoh.t));

for ci = 1:size(dynCoh.coh,1)
    dynCoh.e(ci,:) = simpleClosedLoop('t',dynCoh.t,'speed',dynCoh.speed,'G0',G0,...
        'coh',dynCoh.coh(ci,:),'W',W,'betas',betas,'alphas',alphas,'tauG',tauG,...
        'plotFlg',false);
end

figure('Name','Dynamic Coherence')
plot(dynCoh.t,dynCoh.e');
hold on
plotVertical(150);
plotVertical(cT);

%% Initiation coherence
initCoh.t = 0:300;
speeds = [5,10,20];
cohs = [20,60,100];
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        initCoh.speed(si,:,ci) = speeds(si)*ones(1,length(initCoh.t));
        initCoh.coh(si,:,ci) = cohs(ci)*ones(1,length(initCoh.t));
        
        initCoh.e(si,:,ci) = simpleClosedLoop('t',initCoh.t,'speed',initCoh.speed(si,:,ci),'G0',G0,...
        'coh',initCoh.coh(si,:,ci),'W',W,'betas',betas,'alphas',alphas,'tauG',tauG,...
        'plotFlg',false);
    end
end

figure('Name','Initial Coherence')
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        subplot(1,length(speeds),si)
        plot(initCoh.t,initCoh.e(si,:,ci));
        hold on
        plotHorizontal(speeds(si));
        ylim([0 1.1*max(speeds)])
    end
end

%% Behling init
binit.t = 0:1200;
binit.speed = 10*ones(1,length(binit.t));
cohs = [10,20,30,40,50,60,100];
for ci = 1:length(cohs)
    binit.coh(ci,:) = cohs(ci)*ones(1,length(binit.t));
    
    binit.e(ci,:) = simpleClosedLoop('t',binit.t,'speed',binit.speed,'G0',G0,...
        'coh',binit.coh(ci,:),'W',W,'betas',betas,'alphas',alphas,'tauG',tauG,...
        'plotFlg',false);
end

figure('Name','Behling Init')
plot(binit.t,binit.e')

%% Perturbations
blockT = [450,300,300,300];
onsets = [300, 450, 600, 750, 900 1050];
cT = cumsum(blockT);
perturb.t = 0:sum(blockT)-1;
perturb.coh = 60*ones(1,length(perturb.t));
perturb.speed = 10*ones(length(onsets),size(perturb.t,2));
period = 100;
duration = 100;
amplitude = 0.2*perturb.speed(1);

for oi = 1:length(onsets)
    perturb.speed(oi,onsets(oi):onsets(oi)+99) = perturb.speed(oi,onsets(oi):onsets(oi)+99) + ...
        amplitude*sin(2*pi/period*[0:duration-1]);
    
    perturb.e(oi,:) = simpleClosedLoop('t',perturb.t,'speed',perturb.speed(oi,:),'G0',G0,...
        'coh',perturb.coh,'W',W,'betas',betas,'alphas',alphas,'tauG',tauG,...
        'plotFlg',false);
end

figure('Name','Perturbation response')
plot(perturb.t,perturb.e')
