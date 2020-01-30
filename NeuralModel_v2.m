function [n, M, rNN, e, tuning] = NeuralModel_v2(varargin)
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
Cov_default.separationLengthConstant = 0.3;
Cov_default.alpha = 0;

sizeProps_default.minEccentricity = 1;
sizeProps_default.maxEccentricity = 30;

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
addParameter(Parser,'sizeProps',sizeProps_default)
addParameter(Parser,'mymakeaxisflg',true)

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
sizeProps = Parser.Results.sizeProps;
mymakeaxisflg = Parser.Results.mymakeaxisflg;

%% Generate population size
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson cell density function (mm/deg)
%N = ceil(20*(integral(fun,sizeProps.minEccentricity,sizeProps.maxEccentricity)+9));  % By Pack Lab. Should probably multiply by 2!!!!!
N = ceil( 20* ( 6*( 10*sizeProps.maxEccentricity.^(1/10) - 10*sizeProps.minEccentricity.^(1/10) ) + 9 ) ); % Definite integral of cell density function; should probably multiply by 2!!!!!

%% Initial simulation w/ maximum number of neurons
% Generate population tuning parameters
tuning = tuningFunctions(N,theta,speed,Cov,n0,sizeProps);

% Simulate MT and then decode    
[nTemp, MTemp, rNNTemp, ~, ~] = DirSpeedSizeLocMT(thetas,speeds,max(sizes),'trialN',400,'tuning',tuning,'plotflg',false);

gainTemp = gainFunction(nTemp(1,:,:,:),tuning,0);
normalizer = [mean(gainTemp(:))/mean(log2(speeds)), 0]; % For -90 90 theta range

%% Run MT simulations
for szi = 1:length(sizes)
    
    [n{szi}, M{szi}, rNN{szi}, ~, ~] = DirSpeedSizeLocMT(thetas,speeds,sizes(szi),'trialN',400,'tuning',tuning,'plotflg',false,'mymakeaxisflg',mymakeaxisflg);
    
end

%% Run decoding simulations

% Decoder properties
[Ds,Ss] = meshgrid(thetas,speeds);
s = cat(3,Ds',Ss');

for szi = 1:length(sizes)
    e{szi} = DecodeMT(n{szi},tuning,s,'gainNoise',gainNoise,'epsilon',epsilon,'b',normalizer,'mymakeaxisflg',mymakeaxisflg);%,'plotflg',false);
    
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
for szi = 1:length(sizes)
    subplot(1,length(sizes),szi)
    plot(speeds,squeeze(e{szi}(1,:,:,2)),'ko')
    hold on
    plotUnity;
    xs = linspace(min(speeds),max(speeds),100);
    plot(xs,betas(1,szi)*xs+betas(2,szi),'r-')
    axis square
    xlabel('Target speed (deg/s)')
    ylabel('Eye speed (deg/s)')
    if mymakeaxisflg
        mymakeaxis(gca);
    end
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
    for szi = 1:length(sizes)
        plot(permute(VeM(di,:,szi),[2,3,1]),permute(VeVAR(di,:,szi),[2,3,1]),...
            'o','Color',Ncolors(szi,:),'MarkerFaceColor',Ncolors(szi,:),'MarkerSize',10)
        hold on
        for si = 1:length(speeds)
            %         for si = 1:length(speeds)
%             plot(VeM(di,si,szi),VeVAR(di,si,szi),'o','Color',colors(si,:),...
%                 'MarkerFaceColor',Ncolors(szi,:),'MarkerSize',10)
%             hold on
%             plot(VeM(di,si,szi),gainSDN(VeM(di,si,szi),speeds(si),w,sigG),'s',...
%                 'Color',colors(si,:),'MarkerFaceColor',Ncolors(szi,:),'MarkerSize',10)
            %         end
        end
        x = linspace(0,max(speeds),100);
        plot(betas(1,szi)*x+betas(2,szi),gainSDN(betas(1,szi)*x+betas(2,szi),x,w,sigG),'-','Color',Ncolors(szi,:))
        xlabel('Mean eye speed (deg/s)')
        ylabel('Eye speed variance (deg/s)^2')
    end
    if mymakeaxisflg
        mymakeaxis(gca);
    end
end

h = figure('Name','MT RF centers','Position',[713 316 1530 420]);
subplot(1,2,1)
plot(tuning.size.x,tuning.size.y,'k.')
hold on
for szi = 1:length(sizes)
    [x,y] = ellipse(sizes(szi)/2,sizes(szi)/2,0,0,pi/360);
    plot(x,y,'r')
end
axis([-1.2*sizeProps.maxEccentricity 1.2*sizeProps.maxEccentricity -1.2*sizeProps.maxEccentricity 1.2*sizeProps.maxEccentricity])
axis square
ind1 = find(sqrt(tuning.size.x.^2+tuning.size.y.^2) > 9 & sqrt(tuning.size.x.^2+tuning.size.y.^2) < 15,1);
[x,y] = ellipse(tuning.size.radius(ind1),tuning.size.radius(ind1),tuning.size.x(ind1),tuning.size.y(ind1),pi/360);
plot(x,y,'k')
ind2 = find(sqrt(tuning.size.x.^2+tuning.size.y.^2) < 6 & sqrt(tuning.size.x.^2+tuning.size.y.^2) > 3,1);
[x,y] = ellipse(tuning.size.radius(ind2),tuning.size.radius(ind2),tuning.size.x(ind2),tuning.size.y(ind2),pi/360);
plot(x,y,'k')
plotVertical(0);
plotHorizontal(0);
xlabel('Horizontal position (deg)')
ylabel('Vertical position (deg)')
if mymakeaxisflg
    mymakeaxis(gca);
end

subplot(1,2,2)
for szi = 1:length(sizes)
    plot(sqrt(tuning.size.x.^2+tuning.size.y.^2),squeeze(mean(M{szi},2)),'o')
    hold on
    plot(sqrt(tuning.size.x(ind1).^2+tuning.size.y(ind1).^2),squeeze(mean(M{szi}(:,:,ind1),2)),...
        'ko','MarkerFaceColor',[0 0 0])
    plot(sqrt(tuning.size.x(ind2).^2+tuning.size.y(ind2).^2),squeeze(mean(M{szi}(:,:,ind2),2)),...
        'ko','MarkerFaceColor',[0 0 0])
    plotVertical(sizes(szi)/2);
end

xlabel('Eccentricity (deg)')
ylabel('Mean response')
if mymakeaxisflg
    mymakeaxis(gca);
end

%% Functions

%% sampleEccentricities
function eccentricities = sampleEccentricities(N,minEccentricity,maxEccentricity)
    foveaN = 9*20; % 9*20 neurons from 0 to 1 deg
    u = rand(N-foveaN,1);
    eccentricities = ((maxEccentricity^0.1-1)*u+1).^10;
    foveaMin = 0.25;
    eccentricities = [(1-foveaMin)*rand(foveaN,1)+foveaMin; ...
        eccentricities];
    
    eccentricities = randsample(eccentricities,N);

%% tuningFunctions
function tuning = tuningFunctions(N,theta,speed,Cov,n0,sizeProps)

    thetas = linspace(theta.range(1),theta.range(2),theta.range(3));
    speeds = 2.^(linspace(speed.range(1),speed.range(2),speed.range(3))); 

    A(:,1) = randsample(thetas,N,true);
    A(:,2) = randsample(speeds,N,true);
    B = sortrows(A,[2,1]);
    
    tuning.theta = theta;
    tuning.theta.pref = B(:,1);
    
    tuning.speed = speed;
    tuning.speed.pref = B(:,2);
    
    eccentricities = sampleEccentricities(N,sizeProps.minEccentricity,sizeProps.maxEccentricity);
    rhos = 2*pi*rand(N,1);
    tuning.size.x = eccentricities.*cos(rhos);
    tuning.size.y = eccentricities.*sin(rhos);
    tuning.size.radius = (0.69*eccentricities+1)/sqrt(pi); % Raduis approixmation from Ferrera and Lisberger (1997)
    %tuning.size.radius = sqrt(eccentricities);     % Approximate radius is square root of eccentricity
    if isfield(sizeProps,'threshold')
        tuning.size.threshold = sizeProps.threshold.*ones(N,1);
    else
        tuning.size.threshold = zeros(N,1);
    end
    if isfield(sizeProps,'exponential')
        tuning.size.exponential = sizeProps.exponential.*ones(N,1);
    else
        tuning.size.exponential = ones(N,1);
    end
    if isfield(sizeProps,'surround_weight')
        tuning.size.surround_weight = sizeProps.surround_weight.*ones(N,1);
    else
        tuning.size.surround_weight = zeros(N,1);
    end
    
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
    
%% gainFunction
function out = gainFunction(n,tuning,gainNoise)
    %%
    for thetai = 1:size(n,1)
        for speedi = 1:size(n,2)
            out(thetai,speedi,:) = log2(tuning.speed.pref)' ...
                * permute(n(thetai,speedi,:,:),[4,3,1,2]) + ...
                gainNoise*randn(size(n,3),1)';

            % SDN
            %         out(thetai,speedi,:) = log2(tuning.speed.pref)' ...
            %             * permute(n(thetai,speedi,:,:),[4,3,1,2]);
            %         out(thetai,speedi,:) = out(thetai,speedi,:) + ...
            %             out(thetai,speedi,:)*gainNoise.*randn(1,1,size(n,3));
        end
    end