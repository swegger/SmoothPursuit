function [nFEF, z, nMT, tuning] = dynamicalCohModel_v1(varargin)
%% dynamicalCohModel_v1
%
%
%
%%

%% Defaults
thetas_default = 0;
speeds_default = 10;

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

Sequences_default = [60*ones(1,45), 100*ones(1,30), 20*ones(1,30), 60*ones(1,30);...
    60*ones(1,45), 20*ones(1,30), 100*ones(1,30), 60*ones(1,30);
    60*ones(1,45), 100*ones(1,30), 100*ones(1,30), 60*ones(1,30);
    60*ones(1,45), 20*ones(1,30), 20*ones(1,30), 60*ones(1,30);
    60*ones(1,45), 60*ones(1,30), 60*ones(1,30), 60*ones(1,30)];
cohFunction_default = @(x)(x);

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'Sequences',Sequences_default)
addParameter(Parser,'cohFunction',cohFunction_default)
addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'theta',theta_default)
addParameter(Parser,'speed',speed_default)
addParameter(Parser,'Cov',Cov_default)
addParameter(Parser,'n0',1)
addParameter(Parser,'epsilon',500)
addParameter(Parser,'N',100)
addParameter(Parser,'trialN',400)
addParameter(Parser,'mymakeaxisflg',true)

parse(Parser,varargin{:})

Sequences = Parser.Results.Sequences;
cohFunction = Parser.Results.cohFunction;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
theta = Parser.Results.theta;
speed = Parser.Results.speed;
Cov = Parser.Results.Cov;
n0 = Parser.Results.n0;
epsilon = Parser.Results.epsilon;
N = Parser.Results.N;
trialN = Parser.Results.trialN;
mymakeaxisflg = Parser.Results.mymakeaxisflg;

%% Generate population sizes
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
% fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
% for szi = 1:length(sizes)
%     N = ceil(20*(integral(fun,0.5,sizes(szi))+9));
% end

%% Generate baseline tuning functions for MT neurons
tuning_orig = tuningFunctions(N,theta,speed,Cov,n0);

%% Run simulations
for seqi = 1:size(Sequences,1)
    disp(seqi)
    Coherence = Sequences(seqi,:);
    ampMultiplier = cohFunction(Coherence)/100;
    for ti = 1:length(Coherence)
        
        % Generate population tuning parameters
        tuning = tuning_orig;
        tuning.theta.Amp = theta.Amp*ampMultiplier(ti);
        
        % Decoder properties
        [Ds,Ss] = meshgrid(thetas,speeds);
        s = cat(3,Ds',Ss');
        
        % Simulate MT and then decode
        [n, M, rNN, ~, tuning] = SimpleMT(thetas,speeds,'trialN',trialN,...
            'tuning',tuning,'plotflg',false,'mymakeaxisflg',mymakeaxisflg);
        nMT(:,:,:,:,ti,seqi) = n;
    end
end

%% Model FEF response
alpha = 0.07;
a = 0.0003;
b = 0.0001;
c = 1;
tau = 1/alpha;
for seqi = 1:size(Sequences,1)
    for ti = 1:length(Coherence)
        if ti == 1
            z(ti,:,seqi) = epsilon*ones(1,trialN);%zeros(trialN,1);%1./(tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2]));
            nFEF(ti,:,seqi) = zeros(1,trialN);%(tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2]));
            
            dFEF(ti,:,seqi) = deltaFEF(nFEF(ti,:,seqi),nMT(:,:,:,:,ti,seqi),...
                z(ti,:,seqi),a,b,tuning);
            dz(ti,:,seqi) = deltaZ(nMT(:,:,:,:,ti,seqi),z(ti,:,seqi),c,tuning);
        else
            
            nFEF(ti,:,seqi) = nFEF(ti-1,:,seqi) + dFEF(ti-1,:,seqi)/tau;
            z(ti,:,seqi) = z(ti-1,:,seqi) + dz(ti-1,:,seqi)/tau;
            
            dFEF(ti,:,seqi) = deltaFEF(nFEF(ti,:,seqi),nMT(:,:,:,:,ti,seqi),...
                z(ti,:,seqi),a,b,tuning);
            dz(ti,:,seqi) = deltaZ(nMT(:,:,:,:,ti,seqi),z(ti,:,seqi),c,tuning);
        end
    end
end

%% Plotting

figure('Name','Responses','Position',[265 360 2099 660]);
for seqi = 1:size(Sequences,1)
    subplot(2,size(Sequences,1),seqi)
    imagesc(squeeze(mean(nMT(:,:,:,:,:,seqi),3)))
    hold on
    xlabel('time')
    ylabel('Perf dir')
end

for seqi = 1:size(Sequences,1)
    subplot(2,size(Sequences,1),size(Sequences,1)+seqi)
    plotyy(1:size(Sequences,2),mean(nFEF(:,:,seqi),2),...
        1:size(Sequences,2),mean(z(:,:,seqi),2))
    xlabel('time')
end

figure('Name','Sequence response comparison','Position',[1000 855 1179 474])
controlSeq = 5;
for seqi = 1:size(Sequences,1)
    subplot(2,2,1)
    plot(1:size(Sequences,2),mean(nFEF(:,:,seqi),2))
    hold on
    xlabel('time')
    
    subplot(2,2,3)
    plot(1:size(Sequences,2),mean(nFEF(:,:,seqi),2)-mean(nFEF(:,:,controlSeq),2))
    hold on
    xlabel('time')
    ylabel('Response minus cotnrol sequence')
    
    subplot(2,2,[2,4])
    plot(mean(nFEF(:,:,seqi),2),mean(z(:,:,seqi),2),'.-')
    hold on
    xlabel('nFEF')
    ylabel('z')
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

%% Dynamical systems
% dFEF
function out = deltaFEF(nFEF,nMT,z,a,b,tuning)
    out = -nFEF +...
        a*( (tuning.speed.pref' * permute(nMT(:,:,:,:),[4,3,1,2])) ) - ...
        b*( (z + (tuning.speed.pref' * permute(nMT(:,:,:,:),[4,3,1,2])) ) );
    
% dz
function out = deltaZ(nMT,z,c,tuning)
    out = -z + c*( (tuning.speed.pref' * permute(nMT(:,:,:,:),[4,3,1,2])) );
    
    
% dFEF2
function out = deltaFEF2(nFEF,nMT,z,a,b,tuning)
    out = -nFEF +...
        a*( (tuning.speed.pref' * permute(nMT(:,:,:,:),[4,3,1,2])) ) - ...
        b*( (z + (tuning.speed.pref' * permute(nMT(:,:,:,:),[4,3,1,2])) ) );