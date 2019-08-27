function [n, M, rNN, e] = NeuralModel_v1(varargin)
%% NeuralModel_v1
%
%
%
%%

%% Defaults
thetas_default = linspace(-90,90,20);
speeds_default = 2.^(linspace(-1,8,20)); %cumprod(2*ones(1,6));

theta_default.range = [-90,90,1800];
theta_default.Amp = 10;
theta_default.sig = 100;

speed_default.range = [-1,8,1000];
speed_default.Amp = 10;
speed_default.sig = 1.45;
speed_default.d = 0.1;

Cov_default.sigf = 0.36;
Cov_default.thetaLengthConstant = 0.4;
Cov_default.speedLengthConstant = 0.3;
Cov_default.alpha = 0;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'sizes',[2,6,20])
addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'theta',theta_default)
addParameter(Parser,'speed',speed_default)
addParameter(Parser,'Cov',Cov_default)
addParameter(Parser,'n0',1)
addParameter(Parser,'epsilon',0.05)
addParameter(Parser,'gainNoise',0)

parse(Parser,varargin{:})

sizes = Parser.Results.sizes;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
theta = Parser.Results.theta;
speed = Parser.Results.speed;
Cov = Parser.Results.Cov;
n0 = Parser.Results.n0;
epsilon = Parser.Results.epsilon;
gainNoise = Parser.Results.gainNoise;

%% Generate population sizes
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for szi = 1:length(sizes)
    N(szi) = ceil(20*(integral(fun,0.5,sizes(szi))+9));
end

%% Initial simulation w/ maximum number of neurons
% Generate population tuning parameters
    tuning = tuningFunctions(max(N),theta,speed,Cov,n0);
    
    % Decoder properties
    [Ds,Ss] = meshgrid(thetas,speeds);
    s = cat(3,Ds',Ss');
    
    % Simulate MT and then decode    
    [n, M, rNN, ~, tuning] = SimpleMT(thetas,speeds(end),'trialN',400,'tuning',tuning,'plotflg',false);
%     normalizer = [mean(sqrt(sum(n(1,end,:,:),4)),3) 0];
    normalizer = [120*mean(sqrt(sum(n(1,end,:,:),4)),3) 0];
%     normalizer = (1+(1-tuning.Cov.sigf)*max(N));

%% Run simulations
for szi = 1:length(N)
    disp(N(szi))
    % Generate population tuning parameters
    tuning = tuningFunctions(N(szi),theta,speed,Cov,n0);
    
    % Decoder properties
    [Ds,Ss] = meshgrid(thetas,speeds);
    s = cat(3,Ds',Ss');
    
    % Simulate MT and then decode
    
    [n, M, rNN, ~, tuning] = SimpleMT(thetas,speeds,'trialN',400,'tuning',tuning,'plotflg',true);
    e{szi} = DecodeMT(n,tuning,s,'gainNoise',gainNoise,'epsilon',epsilon,'b',normalizer);%,'plotflg',false);
    
    eBar = mean(e{szi},3);
    eVar = var(e{szi},1,3);
    
    VeM(:,:,szi) = squeeze(eBar(:,:,:,2));
    VeVAR(:,:,szi) = squeeze(eVar(:,:,:,2));

%     figure(h)
%     subplot(1,length(N),szi)
%     plot(speeds,squeeze(e{szi}(1,:,:,2)),'ko')
%     hold on
%     plotUnity;
end

%% Fit gainSDN model across sizes
OPTIONS = optimset('Display','off');
[w,sigG] = fit_gainSDN(speeds,VeM,VeVAR,0.1,gainNoise,OPTIONS);

%% Estimate gain
for szi = 1:length(sizes)
    ys = e{szi}(1,:,:,2);
    xs = repmat(speeds,[size(ys,1),1,size(ys,3)]);
    betas(:,szi) = regress(ys(:),[xs(:) ones(size(xs(:)))]);
end

%% Plotting
h = figure('Name','Target v Eye speed','Position',[664 822 1530 387]);
for szi = 1:length(N)
    subplot(1,length(N),szi)
    plot(speeds,squeeze(e{szi}(1,:,:,2)),'ko')
    hold on
    plotUnity;
    xs = linspace(min(speeds),max(speeds),100);
    plot(xs,betas(1,szi)*xs+betas(2,szi),'r-')
    axis square
    xlabel('Target speed (deg/s)')
    ylabel('Eye speed (deg/s)')
    mymakeaxis(gca);
end

figure;
Ncolors = colormap('lines');

    szcolors = [0.8 0.8 0.8;...
                0.5 0.5 0.5;...
                  0   0   0];
              
Ncolors = [szcolors; Ncolors];

colors = projectColorMaps('speeds','samples',1:length(speeds),...
    'sampleDepth',length(speeds));

% figure
for di = 1:length(thetas)
    subplot(1,length(thetas),di)
    for szi = 1:length(N)
        for si = 1:length(speeds)
            %         for si = 1:length(speeds)
            plot(VeM(di,si,szi),VeVAR(di,si,szi),'s','Color',colors(si,:),...
                'MarkerFaceColor',Ncolors(szi,:),'MarkerSize',10)
            hold on
%             plot(VeM(di,si,szi),gainSDN(VeM(di,si,szi),speeds(si),w,sigG),'o',...
%                 'Color',colors(si,:),'MarkerFaceColor',Ncolors(szi,:),'MarkerSize',10)
            %         end
        end
        x = linspace(0,max(speeds),100);
        plot(betas(1,szi)*x+betas(2,szi),gainSDN(betas(1,szi)*x+betas(2,szi),x,w,sigG),'-','Color',Ncolors(szi,:))
        xlabel('Mean eye speed (deg/s)')
        ylabel('Eye speed variance (deg/s)^2')
    end
    mymakeaxis(gca);
end

%% Functions

%% tuningFunctions
function tuning = tuningFunctions(N,theta,speed,Cov,n0)

    thetas = linspace(theta.range(1),theta.range(2),theta.range(3));
    speeds = 2.^(linspace(speed.range(1),speed.range(2),speed.range(3))); 

    A(:,1) = randsample(thetas,N,true);
    A(:,2) = randsample(speeds,N,true);
    B = sortrows(A,[2,1]);
    
    tuning.theta= theta;
    tuning.theta.pref = B(:,1);
    
    tuning.speed = speed;
    tuning.speed.pref = B(:,2);
    
    tuning.n0 = n0;
    
    tuning.Cov = Cov;
    
%% fit_gainSDN
function [w, sigG] = fit_gainSDN(speeds,VeM,VeVAR,w0,sigG0,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeVAR - gainSDN(VeM,repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),p(2))).^2 ) ) );
    p = fmincon(minimizer,[w0 sigG0],[],[],[],[],[0,0],[Inf,Inf],[],OPTIONS);
    w = p(1);
    sigG = p(2);
    
%% gainSDN
function v = gainSDN(ve,vs,w,sigG)
    %%
    v = w.^2.*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;