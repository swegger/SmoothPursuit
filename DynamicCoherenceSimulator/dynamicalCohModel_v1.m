function [n, M, rNN, e, tuning] = dynamicalCohModel_v1(varargin)
%% dynamicalCohModel_v1
%
%
%
%%

%% Defaults
thetas_default = 0;%linspace(-90,90,20);
speeds_default = 10;%2.^(linspace(-1,8,20)); %cumprod(2*ones(1,6));

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

Coherence_default = [60*ones(1,20), 100*ones(1,20)];
cohFunction_default = @(x)(x);

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'Coherence',Coherence_default)
addParameter(Parser,'cohFunction',cohFunction_default)
addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'theta',theta_default)
addParameter(Parser,'speed',speed_default)
addParameter(Parser,'Cov',Cov_default)
addParameter(Parser,'n0',1)
addParameter(Parser,'epsilon',500)
addParameter(Parser,'gainNoise',0)
addParameter(Parser,'N',100)
addParameter(Parser,'mymakeaxisflg',true)

parse(Parser,varargin{:})

Coherence = Parser.Results.Coherence;
cohFunction = Parser.Results.cohFunction;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
theta = Parser.Results.theta;
speed = Parser.Results.speed;
Cov = Parser.Results.Cov;
n0 = Parser.Results.n0;
epsilon = Parser.Results.epsilon;
gainNoise = Parser.Results.gainNoise;
N = Parser.Results.N;
mymakeaxisflg = Parser.Results.mymakeaxisflg;

%% Generate population sizes
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
% fun=@(x1) (6*x1.^-0.9); % Erickson (mm/deg)
% for szi = 1:length(sizes)
%     N = ceil(20*(integral(fun,0.5,sizes(szi))+9));
% end

%% Generate theta amplidue multipliers over time
ampMultiplier = cohFunction(Coherence)/100;
tuning_orig = tuningFunctions(N,theta,speed,Cov,n0);

%% Run simulations
theta_temp = theta;
for ti = 1:length(Coherence)
    disp(ti)
    % Generate population tuning parameters
    tuning = tuning_orig;
    tuning.theta.Amp = theta.Amp*ampMultiplier(ti);
    
    % Decoder properties
    [Ds,Ss] = meshgrid(thetas,speeds);
    s = cat(3,Ds',Ss');
    
    % Simulate MT and then decode
    [n, M, rNN, ~, tuning] = SimpleMT(thetas,speeds,'trialN',400,'tuning',tuning,'plotflg',false,'mymakeaxisflg',mymakeaxisflg);    
    nMT(:,:,:,:,ti) = n;
end

%% Model FEF response
alpha = 0.3;
for ti = 1:length(Coherence)
    if ti == 1
        z(ti,:) = 1./(tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2]));
    else
        z(ti,:) = (1-alpha)*z(ti-1,:) + alpha*...
            (1./(tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2])) + 1./z(ti-1,:)) ./ ...
            (1./(tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2])).*(1./z(ti-1,:)));
    end
    nFEF(ti,:) = mean((tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2])) ./ ...
        (epsilon + sum(permute(nMT(:,:,:,:,ti),[4,3,1,2]),1)),2);
    
    if ti == 1
        nFEF2(ti,:) = nFEF(ti,:);
    else
        nFEF2(ti,:) = (1-alpha)*nFEF2(ti-1,:) + alpha*nFEF(ti,:);
    end
%     nFEF2(ti,:) = mean((tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2])) ./ ...
%         (z(ti,:)/1e5 + (tuning.speed.pref' * permute(nMT(:,:,:,:,ti),[4,3,1,2]))),2);
end

%% Plotting

figure;
subplot(2,2,1)
imagesc(squeeze(mean(nMT,3)))
hold on
xlabel('time')
ylabel('Perf dir')

subplot(2,2,3)
% plot(1:length(Coherence),squeeze(tuning.speed.pref' * mean(permute(nMT,[4,5,3,1,2]),3)) ./ ...
%     (epsilon+squeeze(sum(mean(nMT,3),4)))')
plot(1:length(Coherence),epsilon+squeeze(sum(mean(nMT,3),4)))
hold on
plot(1:length(Coherence),squeeze(tuning.speed.pref' * mean(permute(nMT,[4,5,3,1,2]),3)))
xlabel('time')
ylabel('numerator and denominator')

subplot(2,2,[2 4])
plotyy(1:length(Coherence),nFEF,1:length(Coherence),nFEF2)
xlabel('time')


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
    