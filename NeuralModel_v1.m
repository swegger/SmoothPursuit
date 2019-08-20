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
theta_default.sig = 45;

speed_default.range = [-1,8,1000];
speed_default.Amp = 10;
speed_default.sig = 1;
speed_default.d = 0.1;

Cov_default.sigf = 0.36;
Cov_default.thetaLengthConstant = 0.3;
Cov_default.speedLengthConstant = 0.4;
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

parse(Parser,varargin{:})

sizes = Parser.Results.sizes;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
theta = Parser.Results.theta;
speed = Parser.Results.speed;
Cov = Parser.Results.Cov;
n0 = Parser.Results.n0;

%% Generate population sizes
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
for szi = 1:length(sizes)
    N(szi) = ceil(20*(integral(fun,0.5,sizes(szi))+9));
end

%% Run simulations
h = figure;
for szi = 1:length(N)
    disp(N(szi))
    % Generate population tuning parameters
    tuning = tuningFunctions(N(szi),theta,speed,Cov,n0);
    
    % Decoder properties
    [Ds,Ss] = meshgrid(thetas,speeds);
    s = cat(3,Ds',Ss');
    
    % Simulate MT and then decode
    
    [n, M, rNN, ~, tuning] = SimpleMT(thetas,speeds,'trialN',400,'tuning',tuning,'plotflg',true);
    e = DecodeMT(n,tuning,s,'gainNoise',0);%,'plotflg',false);
    
    eBar = mean(e,3);
    eVar = var(e,1,3);
    
    VeM(:,:,szi) = squeeze(eBar(:,:,:,2));
    VeVAR(:,:,szi) = squeeze(eVar(:,:,:,2));

    figure(h)
    subplot(1,length(N),szi)
    plot(speeds,squeeze(e(1,:,:,2)),'ko')
    hold on
    plotUnity;
end
%% Plotting
Ncolors = colormap('lines');

    szcolors = [0.8 0.8 0.8;...
                0.5 0.5 0.5;...
                  0   0   0];

Ncolors = [szcolors; Ncolors];
% figure
for di = 1:length(thetas)
    subplot(1,length(thetas),di)
    for szi = 1:length(N)
        %         for si = 1:length(speeds)
        plot(VeM(di,:,szi),VeVAR(di,:,szi),'o','Color',Ncolors(szi,:),...
            'MarkerFaceColor',Ncolors(szi,:))
        hold on
        %         end
        xlabel('Mean eye speed (deg/s)')
        ylabel('Eye speed variance (deg/s)^2')
    end
    mymakeaxis(gca);
end

%% Functions

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