function gainNoiseNeuralModelParameterSweep_analysis(varargin)
%% gainNoiseNeuralModelParameterSweep_analysis
%
%   gainNoiseNeuralModelParameterSweep_analysis()
%   
%   Performs analysis of neural model output for default results directory
%
%%

%% Defaults
dataDirectory_default = '~/Projects/MultiSizePursuit/parameterSweeps';
dataDate_default = '20200326';


%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dataDirectory',dataDirectory_default)
addParameter(Parser,'dataDate',dataDate_default)
addParameter(Parser,'calcNeuronPursuitCorrelation',false)

parse(Parser,varargin{:})

dataDirectory = Parser.Results.dataDirectory;
dataDate = Parser.Results.dataDate;
calcNeuronPursuitCorrelation = Parser.Results.calcNeuronPursuitCorrelation;

%% Find files corresponding to desired date
potentialFiles = dir([dataDirectory '/*' dataDate '*.mat']);
files = potentialFiles(3:end);

%% For each file, load in parameter values and results
params = nan(length(files),4);
ws = nan(length(files),1);
sigGs = nan(length(files),1);
Gs = nan(length(files),3);

for filei = 1:length(files)
    results = load(files(filei).name,'sizeProps','w','sigG','G','gainNoise');
    params(filei,1) = results.sizeProps.exponential;
    params(filei,2) = results.sizeProps.threshold;
    params(filei,3) = results.sizeProps.surround_weight;
    params(filei,4) = results.gainNoise;
    ws(filei) = results.w;
    sigGs(filei) = results.sigG;
    Gs(filei,:) = results.G;
    
    if calcNeuronPursuitCorrelation
        results = load(files(filei).name,'e','n','s','tuning','epsilon','normalizer');
        n = results.n{1};
        e = results.e{1};
        s = results.s;
        
        % z-score
        nM = mean(n,3);
        pM = mean(e,3);
        nVar = var(n,[],3);
        pVar = var(e,[],3);
        
        nZ = (n - repmat(nM,[1,1,size(n,3),1]))./repmat(sqrt(nVar),[1,1,size(n,3),1]);
        pZ = (e - repmat(pM,[1,1,size(n,3),1]))./repmat(sqrt(pVar),[1,1,size(n,3),1]);
        
        % Correlate responses
        Rs = nan(size(s,1),size(s,2),size(n,4),2);
        PVals = nan(size(s,1),size(s,2),size(n,4),2);
        for thetai = 1:size(s,1)
            disp(thetai/size(s,1))
            for speedi = 1:size(s,2)
                for neuroni = 1:size(nZ,4)
                    [rtemp,ptemp] = corrcoef(...
                        [permute(nZ(thetai,speedi,:,neuroni),[3,4,1,2]),...
                        permute(pZ(thetai,speedi,:,:),[3,4,2,1])]);
                    
                    Rs(thetai,speedi,neuroni,:) = rtemp(1,2:3);
                    PVals(thetai,speedi,neuroni,:) = ptemp(1,2:3);
                end
            end
        end
        
        % Max corr
        maxCorr(filei,1) = max(Rs(:));
        
        % Now add noise
        [~, ~, ~, Rs, ~] = DecodeMT(n,results.tuning,s,...
            'gainNoise',0.4,'epsilon',results.epsilon,'b',results.normalizer,...
            'mymakeaxisflg',false,'plotflg',false,...
            'motorNoise',0);
        maxCorr_noise(filei,1) = max(Rs(:));
    end
end

%% Load in experimental fit
sigGobs = 0.1;

%%
figure('Name','Parameter effect on sigG')
subplot(2,3,1)
% scatter3(params(:,1),params(:,2),params(:,3),10,sigGs.^2*1e9)
% xlabel('Exp')
% ylabel('Thres')
% zlabel('SS')
% colorbar
% cax = caxis;
subplot(2,3,2)
exps = unique(params(:,1));
thres = unique(params(:,2));
sss = unique(params(:,3));
for expi = 1:length(exps)
    for thresi = 1:length(thres)
        msigG(expi,thresi) = mean(sigGs(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
        p1(expi,thresi) = exps(expi);
        p2(expi,thresi) = thres(thresi);
    end
end
scatter(p1(:),p2(:),10,msigG(:).^2*1e9)
caxis(cax)
%scatter(params(:,1),params(:,2),10,sigGs)
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
clear msigG p1 p2
for expi = 1:length(exps)
    for ssi = 1:length(sss)
        msigG(expi,ssi) = mean(sigGs(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
        p1(expi,ssi) = exps(expi);
        p2(expi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,msigG(:).^2*1e9)
caxis(cax)
%scatter(params(:,1),params(:,3),10,sigGs)
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
clear msigG p1 p2
for thresi = 1:length(thres)
    for ssi = 1:length(sss)
        msigG(thresi,ssi) = mean(sigGs(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
        p1(thresi,ssi) = thres(thresi);
        p2(thresi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,msigG(:).^2*1e9)
caxis(cax)
%scatter(params(:,2),params(:,3),10,sigGs)
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
histogram((sigGs.^2)*1e9,20,'Normalization','probability')
xlabel('\sigma^2')
ylabel('Relative frequency')
% plotVertical(mean(sigGs));
% hold on
% plotVertical(sigGobs);

%%
figure('Name','Parameter effect on w')
subplot(2,3,1)
scatter3(params(:,1),params(:,2),params(:,3),10,ws)
xlabel('Exp')
ylabel('Thres')
zlabel('SS')
subplot(2,3,2)
for expi = 1:length(exps)
    for thresi = 1:length(thres)
        mw = mean(ws(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
        scatter(exps(expi),thres(thresi),10,mw)
        hold on
    end
end
%scatter(params(:,1),params(:,2),10,ws)
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
for expi = 1:length(exps)
    for ssi = 1:length(sss)
        mw = mean(ws(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
        scatter(exps(expi),sss(ssi),10,mw)
        hold on
    end
end
%scatter(params(:,1),params(:,3),10,ws)
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
for thresi = 1:length(thres)
    for ssi = 1:length(sss)
        mw = mean(ws(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
        p1(thresi,ssi) = thres(thresi);
        p2(thresi,ssi) = sss(ssi);
    end
end
%scatter(params(:,2),params(:,3),10,ws)
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
histogram(ws,linspace(0,0.2,100))
% plotVertical(mean(ws));
hold on
% plotVertical(sigGobs);

%%
figure('Name','Parameter effect on Gain')
subplot(1,3,1)
plot3(params(:,1),params(:,2),Gs,'o')
xlabel('Exp')
ylabel('Thres')
subplot(1,3,2)
plot3(params(:,1),params(:,3),Gs,'o')
xlabel('Exp')
ylabel('SS')
subplot(1,3,3)
plot3(params(:,2),params(:,3),Gs,'o')
xlabel('Thres')
ylabel('SS')

%% Parameter effect on maximum neuron-behavior correlation
if calcNeuronPursuitCorrelation
    figure('Name','Parameter effect on maximum neuron-behavior correlation')
    subplot(2,3,1)
    % scatter3(params(:,1),params(:,2),params(:,3),10,sigGs.^2*1e9)
    % xlabel('Exp')
    % ylabel('Thres')
    % zlabel('SS')
    % colorbar
    % cax = caxis;
    subplot(2,3,2)
    exps = unique(params(:,1));
    thres = unique(params(:,2));
    sss = unique(params(:,3));
    for expi = 1:length(exps)
        for thresi = 1:length(thres)
            mR(expi,thresi) = mean(maxCorr(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
            p1(expi,thresi) = exps(expi);
            p2(expi,thresi) = thres(thresi);
        end
    end
    scatter(p1(:),p2(:),10,mR(:))
    cax = [0.53 0.55];
    caxis(cax)
    %scatter(params(:,1),params(:,2),10,sigGs)
    xlabel('Exp')
    ylabel('Thres')
    subplot(2,3,4)
    clear msigG p1 p2
    for expi = 1:length(exps)
        for ssi = 1:length(sss)
            mR(expi,ssi) = mean(maxCorr(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
            p1(expi,ssi) = exps(expi);
            p2(expi,ssi) = sss(ssi);
        end
    end
    scatter(p1(:),p2(:),10,mR(:))
    caxis(cax)
    %scatter(params(:,1),params(:,3),10,sigGs)
    xlabel('Exp')
    ylabel('SS')
    subplot(2,3,5)
    clear msigG p1 p2
    for thresi = 1:length(thres)
        for ssi = 1:length(sss)
            mR(thresi,ssi) = mean(maxCorr(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
            p1(thresi,ssi) = thres(thresi);
            p2(thresi,ssi) = sss(ssi);
        end
    end
    scatter(p1(:),p2(:),10,mR(:))
    caxis(cax)
    %scatter(params(:,2),params(:,3),10,sigGs)
    xlabel('Thres')
    ylabel('SS')
    
    subplot(2,3,3)
    histogram(maxCorr,linspace(0,0.6,50),'Normalization','probability')
    hold on
    histogram(maxCorr_noise,linspace(0,0.6,50),'Normalization','probability')
    xlabel('Maximum MT-behavior correalation')
    ylabel('Relative frequency')
    legend('\sigma_g^2 = 0','\sigma_g^2 = 0.16')
    
end
