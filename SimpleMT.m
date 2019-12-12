function [n, M, rNN, pNN, tuning] = SimpleMT(thetas,speeds,varargin)
%% SimpleMT
%
%
%%

%% Defaults
tuning_default = defaultTuning();

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dirs')
addRequired(Parser,'speeds')
addParameter(Parser,'trialN',1)
addParameter(Parser,'tuning',tuning_default)
addParameter(Parser,'plotflg',true)
addParameter(Parser,'mymakeaxisflg',true)

parse(Parser,thetas,speeds,varargin{:})

thetas = Parser.Results.dirs;
speeds = Parser.Results.speeds;
trialN = Parser.Results.trialN;
tuning = Parser.Results.tuning;
plotflg = Parser.Results.plotflg;
mymakeaxisflg = Parser.Results.mymakeaxisflg;

%% Generate noisy bump for each trial
% Preallocate
f = nan(length(thetas),length(speeds),length(tuning.theta.pref));
n = nan(length(thetas),length(speeds),trialN,length(tuning.theta.pref));

% Generate covariance matrix
[C,Q] = neuralCov(tuning);

% Cycle through input thetas and speeds
for thetai = 1:length(thetas)
    for speedi = 1:length(speeds)
        f(thetai,speedi,:) = neuralTuning(thetas(thetai),speeds(speedi),tuning);
        
%         n(thetai,speedi,:,:) = mvnrnd(permute(f(thetai,speedi,:),[1,3,2]),C,trialN);
        for triali = 1:trialN
            y = permute(Q*randn(length(tuning.theta.pref),1),[2,3,4,1]);
            n(thetai,speedi,triali,:) = permute(f(thetai,speedi,:),[1,2,4,3]) ...
                + y.*sqrt(1.5*permute(f(thetai,speedi,:),[1,2,4,3]));
        end
    end
end

%% Summary stats
% Means
M = permute(nanmean(n,3),[1,2,4,3]);
Mtheta = nansum(nansum(n,2),3)/(size(n,2)+size(n,3));
Mspeed = nansum(nansum(n,1),3)/(size(n,1)+size(n,3));

% Signal correlations

% Noise correlations
noise = n - repmat(permute(M,[1 2 4 3]),[1,1,size(n,3),1]);
for neuroni = 1:size(n,4)
    temp1 = reshape(noise(:,:,:,neuroni),[numel(noise(:,:,:,neuroni)),1]);
    temp1 = temp1/std(temp1);
    for neuronj = 1:size(n,4)
        temp2 = reshape(noise(:,:,:,neuronj),[numel(noise(:,:,:,neuronj)),1]);
        temp2 = temp2/std(temp2);
        [tempR, tempP] = corrcoef(temp1,temp2);
        rNN(neuroni,neuronj) = tempR(1,2);
        pNN(neuroni,neuronj) = tempP(1,2);
    end
end

%% Plotting
if plotflg
    figure('Name','Noise correlations','Position',[150 157 1089 641])
    h = subplot(2,2,[1 3]);
    imagesc(1:length(tuning.theta.pref),1:length(tuning.theta.pref),rNN - diag(diag(rNN)))
    colormap gray
    axis square
    h.YDir = 'normal';
    xlabel('Neuron i')
    ylabel('Neuron j')
    colorbar
    if mymakeaxisflg
        mymakeaxis(gca)
    end
    
    subplot(2,2,2)
    [S1,S2] = meshgrid(tuning.speed.pref);
    dS = S1-S2;
    mask = logical(tril(ones(size(rNN)),-1));
    hS = scatter(abs(dS(mask)),rNN(mask),20,'k');
    hS.MarkerFaceColor = [0 0 0];
    hS.MarkerFaceAlpha = 0.01;
    hS.MarkerEdgeAlpha = 0.01;
    hold on
    axis([0 150 -0.25 1])
    plotHorizontal(0);
    xlabel('\Delta Speed perference')
    ylabel('Noise Correlation')
    if mymakeaxisflg
        mymakeaxis(gca)
    end
    
    subplot(2,2,4)
    [D1,D2] = meshgrid(tuning.theta.pref);
    dD = D1-D2;
    mask = logical(tril(ones(size(rNN)),-1));
    hS = scatter(abs(dD(mask)),rNN(mask),20,'k');
    hS.MarkerFaceColor = [0 0 0];
    hS.MarkerFaceAlpha = 0.01;
    hS.MarkerEdgeAlpha = 0.01;
    hold on
    axis([0 180 -0.25 1])
    plotHorizontal(0);
    plotHorizontal(0.2);
    xlabel('\Delta Direction perference')
    ylabel('Noise Correlation')
    if mymakeaxisflg
        mymakeaxis(gca)
    end
end

%% Functions

%% neuralTuning
function f = neuralTuning(theta,speed,tuning)
f = tuning.n0 + ...
    tuning.theta.Amp ...
    .*exp( -(theta - tuning.theta.pref).^2./tuning.theta.sig.^2/2 ) .* ...
    tuning.speed.Amp ...
    .*exp( -( log2(speed./tuning.speed.pref) ).^2 ./ ...
    (tuning.speed.sig).^2/2 );

%% neuralCov
function [C,Q] = neuralCov(tuning)
%% Model the covariance as a Gaussian process
N = length(tuning.theta.pref);
C = zeros(length(tuning.theta.pref),length(tuning.theta.pref));
Q = zeros(size(C));
for i = 1:length(tuning.theta.pref)
    for j = 1:length(tuning.theta.pref)
        
        if i > j
            if isfield(tuning.Cov,'structured') && ~tuning.Cov.structured
                C(i,j) = tuning.Cov.sigf;
            else
                C(i,j) = tuning.Cov.sigf * exp( ...
                    -( ...
                    ((1-tuning.Cov.alpha)*(tuning.theta.pref(i)-tuning.theta.pref(j)) +...
                    tuning.Cov.alpha*(360*rand-180)).^2 ./ ...
                    (180^2*tuning.Cov.thetaLengthConstant.^2) + ...
                    ((1-tuning.Cov.alpha)*(tuning.speed.pref(i)-tuning.speed.pref(j)) + ...
                    tuning.Cov.alpha*256*rand).^2 ./ ...
                    (255.5^2*tuning.Cov.speedLengthConstant.^2) ...
                    ) );
                
            end
            Q(i,j) = 1/sqrt(N) * ...
                sqrt( 2/N + C(i,j) - 2*C(i,j)/N - 2/N * ...
                sqrt( (1-C(i,j))*(1-C(i,j)+C(i,j)*N) ) );
        end
    end
end
Q = real(Q);
for i = 1:length(tuning.theta.pref)
    for j = 1:length(tuning.theta.pref)
        if i < j
            C(i,j) = C(j,i);
            Q(i,j) = Q(j,i);
        elseif i ==j
            C(i,j) = tuning.Cov.sigf;%C(i,j);
            
            Q(i,j) = 1/C(i,j)/sqrt(N) * ...
                sqrt( 2/N + C(i,j) - 2*C(i,j)/N -2/N * ...
                sqrt( (1-C(i,j))*(1-C(i,j)+C(i,j)*N) ) ) * ...
                (1 + sqrt( (1-C(i,j))*(1 - C(i,j) + C(i,j)*N) ) );
        end
        
    end
end
Q = real(Q);
QQ = Q*Q';
% QQ = QQ - diag(diag(QQ)) + diag(ones(size(QQ,1)));
Q = real(sqrtm(QQ));
% Q = real(sqrtm(C));

%% defaultTuning
function tuning = defaultTuning()
%%
thetas = linspace(-180,180,36);
speeds = 2.^(linspace(-1,8,20));
[Ts, Ss] = meshgrid(thetas,speeds);
tuning.theta.Amp = 10;
tuning.theta.pref = Ts(:);
tuning.theta.sig = 45;

tuning.speed.Amp = 10;
tuning.speed.pref = Ss(:);
tuning.speed.sig = 1;
tuning.speed.d = 0.1;

tuning.n0 = 1;

tuning.Cov.sigf = 0.36;
tuning.Cov.thetaLengthConstant = 0.3;
tuning.Cov.speedLengthConstant = 0.4;
tuning.Cov.alpha = 0;