function [n, M, rNN, pNN, tuning] = DirSpeedSizeLocMT(thetas,speeds,sz,varargin)
%% DirSpeedSizeLocMT
%
%
%%

%% Defaults
tuning_default = defaultTuning();

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dirs')
addRequired(Parser,'speeds')
addRequired(Parser,'sz')
addParameter(Parser,'trialN',1)
addParameter(Parser,'tuning',tuning_default)
addParameter(Parser,'plotflg',true)
addParameter(Parser,'mymakeaxisflg',true)

parse(Parser,thetas,speeds,sz,varargin{:})

thetas = Parser.Results.dirs;
speeds = Parser.Results.speeds;
sz = Parser.Results.sz;
trialN = Parser.Results.trialN;
tuning = Parser.Results.tuning;
plotflg = Parser.Results.plotflg;
mymakeaxisflg = Parser.Results.mymakeaxisflg;

%% Generate noisy bump for each trial
% Preallocate
f = nan(length(thetas),length(speeds),length(tuning.theta.pref));
n = nan(length(thetas),length(speeds),trialN,length(tuning.theta.pref));

% Cycle through input thetas and speeds
for thetai = 1:length(thetas)
    for speedi = 1:length(speeds)
        f(thetai,speedi,:) = neuralTuning(thetas(thetai),speeds(speedi),sz,tuning);
        [C,Q] = neuralCov(tuning,thetas(thetai),sz,speeds(speedi));
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
    h = subplot(3,2,[1 3 5]);
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
    
    subplot(3,2,2)
    [S1,S2] = meshgrid(tuning.speed.pref);
    dS = S1-S2;
    mask = logical(tril(ones(size(rNN)),-1));
    alldS = abs(dS(mask));
    allrNN = rNN(mask);
    exInds = randsample(length(alldS),200);
    hS = scatter(alldS(exInds),allrNN(exInds),20,'k');
    hS.MarkerFaceColor = [0 0 0];
%     hS.MarkerFaceAlpha = 0.01;
%     hS.MarkerEdgeAlpha = 0.01;
    hold on
    axis([0 150 -0.25 1])
    plotHorizontal(0);
    xlabel('\Delta Speed perference')
    ylabel('Noise Correlation')
    if mymakeaxisflg
        mymakeaxis(gca,'yticks',[-0.25,0,0.25,0.5,0.75,1])
    end
    
    subplot(3,2,4)
    [D1,D2] = meshgrid(tuning.theta.pref);
    dD = D1-D2;
    mask = logical(tril(ones(size(rNN)),-1));
    alldD = abs(dD(mask));
    allrNN = rNN(mask);
    exInds = randsample(length(alldD),200);
    hS = scatter(alldD(exInds),allrNN(exInds),20,'k');
    hS.MarkerFaceColor = [0 0 0];
%     hS.MarkerFaceAlpha = 0.01;
%     hS.MarkerEdgeAlpha = 0.01;
    hold on
    axis([0 180 -0.25 1])
    plotHorizontal(0);
    plotHorizontal(0.2);
    xlabel('\Delta Direction perference')
    ylabel('Noise Correlation')
    if mymakeaxisflg
        mymakeaxis(gca,'yticks',[-0.25,0,0.25,0.5,0.75,1])
    end
    
    subplot(3,2,6)
    [X1,X2] = meshgrid(tuning.size.x);
    [Y1,Y2] = meshgrid(tuning.size.y);
    dX = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    mask = logical(tril(ones(size(rNN)),-1));
    alldX = abs(dX(mask));
    allrNN = rNN(mask);
    exInds = randsample(length(alldX),200);
    hS = scatter(alldX(exInds),allrNN(exInds),20,'k');
    hS.MarkerFaceColor = [0 0 0];
%     hS.MarkerFaceAlpha = 0.01;
%     hS.MarkerEdgeAlpha = 0.01;
    hold on
    axis([0 60 -0.25 1])
    plotHorizontal(0);
    plotHorizontal(0.2);
    xlabel('\Delta Position')
    ylabel('Noise Correlation')
    if mymakeaxisflg
        mymakeaxis(gca,'yticks',[-0.25,0,0.25,0.5,0.75,1])
    end
end

%% Functions

%% neuralTuning
function f = neuralTuning(theta,speed,sz,tuning)
f = tuning.n0 + ...
    dirTuning(theta,tuning) .* ...
    speedTuning(speed,tuning) .* ...
    sizeTuning(sz,tuning);

%% dirTuning
function f = dirTuning(theta,tuning)
    f = tuning.theta.Amp ...
        .*exp( -(theta - tuning.theta.pref).^2./tuning.theta.sig.^2/2 );
    
%% speedTuning
function f = speedTuning(speed,tuning)
    f = tuning.speed.Amp ...
        .*exp( -( log2(speed./tuning.speed.pref) ).^2 ./ ...
        (tuning.speed.sig).^2/2 );

%% sizeTuning
function f = sizeTuning(sz,tuning)
    sz = sz/2;  % Convert diameter to radius
    
    % Determine response due to classical receptive field
    d = sqrt(tuning.size.x.^2 + tuning.size.y.^2);
    a = tuning.size.radius.^2 .* acos( (d.^2 + tuning.size.radius.^2 - sz.^2)./(2*d.*tuning.size.radius) );
    b = sz.^2 .* acos( (d.^2 + sz.^2 - tuning.size.radius.^2)./(2*d.*sz) );
    c = sqrt( (-d+tuning.size.radius+sz).*(d+tuning.size.radius-sz).*(d-tuning.size.radius+sz).*(d+tuning.size.radius+sz) )/2;
    z = real((a+b-c)./(pi*tuning.size.radius.^2));
%     f = z;
    
    y = (1./(1-tuning.size.threshold)).*( sqrt(z) - tuning.size.threshold);
    y(y<0) = 0;
    r_crf = y.^tuning.size.exponential;
    
    % Determine surround response
    surround_radius = 3*tuning.size.radius;
    d = sqrt(tuning.size.x.^2 + tuning.size.y.^2);
    a = surround_radius.^2 .* acos( (d.^2 + surround_radius.^2 - sz.^2)./(2*d.*surround_radius) );
    b = sz.^2 .* acos( (d.^2 + sz.^2 - surround_radius.^2)./(2*d.*sz) );
    c = sqrt( (-d+surround_radius+sz).*(d+surround_radius-sz).*(d-surround_radius+sz).*(d+surround_radius+sz) )/2;
    zsur = real((a+b-c)./(pi*surround_radius.^2));
    
    y = (1./(1-tuning.size.threshold)).*( sqrt(zsur) - tuning.size.threshold);
    y(y<0) = 0;
    r_sur = y.^tuning.size.exponential;
    
    % Set epsilon such that maximal response occurs at units' radii
    d = sqrt(tuning.size.x.^2 + tuning.size.y.^2);
    a = surround_radius.^2 .* acos( (d.^2 + surround_radius.^2 - tuning.size.radius.^2)./(2*d.*surround_radius) );
    b = tuning.size.radius.^2 .* acos( (d.^2 + tuning.size.radius.^2 - surround_radius.^2)./(2*d.*tuning.size.radius) );
    c = sqrt( (-d+surround_radius+tuning.size.radius).*(d+surround_radius-tuning.size.radius).*(d-surround_radius+tuning.size.radius).*(d+surround_radius+tuning.size.radius) )/2;
    zepsi = real((a+b-c)./(pi*surround_radius.^2));
    
    y = (1./(1-tuning.size.threshold)).*( sqrt(zepsi) - tuning.size.threshold);
    y(y<0) = 0;
    epsilon = (1-y.^tuning.size.exponential)/2;
    
    % Implement divisive normalization
    if (1-tuning.size.surround_weight) + tuning.size.surround_weight.*(epsilon + (r_crf + r_sur)/2) == 0
        f = r_crf;
    else
        f = r_crf ./ ((1-tuning.size.surround_weight) + tuning.size.surround_weight.*(epsilon + (r_crf + r_sur)/2));
    end
    
    
    
%% neuralCov
function [C,Q] = neuralCov(tuning,theta,sz,s)
%% Model the covariance as a Gaussian process
N = length(tuning.theta.pref);
C = zeros(length(tuning.theta.pref),length(tuning.theta.pref));
Q = zeros(size(C));
diffC = tuning.Cov.sigf*derivativeWRTspeed(tuning,theta,sz,s)*derivativeWRTspeed(tuning,theta,sz,s)';
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
                    (255.5^2*tuning.Cov.speedLengthConstant.^2) + ...
                    ((1-tuning.Cov.alpha)*(sqrt( (tuning.size.x(i)-tuning.size.x(j)).^2 + (tuning.size.y(i)-tuning.size.y(j)).^2 )) + ...
                    tuning.Cov.alpha*60*rand).^2 ./ ...
                    (60^2*tuning.Cov.separationLengthConstant.^2) ...
                    ) );
                
%                 C(i,j) = tuning.Cov.sigf * exp( ...
%                     -( ...
%                     ((1-tuning.Cov.alpha)*(tuning.theta.pref(i)-tuning.theta.pref(j)) +...
%                     tuning.Cov.alpha*(360*rand-180)).^2 ./ ...
%                     (180^2*tuning.Cov.thetaLengthConstant.^2) ...
%                     ) );
                
%                 C(i,j) = tuning.Cov.sigf * exp( ...
%                     -( ...
%                     ((1-tuning.Cov.alpha)*(tuning.speed.pref(i)-tuning.speed.pref(j)) + ...
%                     tuning.Cov.alpha*256*rand).^2 ./ ...
%                     (255.5^2*tuning.Cov.speedLengthConstant.^2) ...
%                     ) );

%                 C(i,j) = tuning.Cov.sigf * exp( ...
%                     -( ...
%                     ((1-tuning.Cov.alpha)*(sqrt( (tuning.size.x(i)-tuning.size.x(j)).^2 + (tuning.size.y(i)-tuning.size.y(j)).^2 )) + ...
%                     tuning.Cov.alpha*60*rand).^2 ./ ...
%                     (60^2*tuning.Cov.separationLengthConstant.^2) ...
%                     ) );
                
                C(i,j) = (1-tuning.Cov.diffAlpha)*C(i,j) + ...
                    tuning.Cov.diffAlpha*diffC(i,j);
                
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

%% Derivative wrt speed
function df = derivativeWRTspeed(tuning,theta,sz,s)
%%
f_theta = dirTuning(theta,tuning);
f_sz = sizeTuning(sz,tuning);

df_s = (-tuning.speed.Amp.*log2(s./tuning.speed.pref)./(s*tuning.speed.sig.^2*log(2))) .* ...
    exp( -(log2(s./tuning.speed.pref)).^2./(2*tuning.speed.sig.^2) );

df = f_theta.*f_sz.*df_s;

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

tuning.size.x = 0;
tuning.size.y = 0;
tuning.size.radius = 7.9/sqrt(pi);
tuning.size.threshold = 0;
tuning.size.exponential = 1;
tuning.size.surround_weight = 0;
%tuning.size.epsilon = 0.1;
% tuning.size.Ae = 1;
% tuning.size.Ai = 0.1;
% tuning.size.se = sqrt(tuning.loc.x.^2+tuning.loc.y.^2);
% tuning.size.si = sqrt(tuning.loc.x.^2+tuning.loc.y.^2);

tuning.n0 = 1;

tuning.Cov.sigf = 0.36;
tuning.Cov.thetaLengthConstant = 0.3;
tuning.Cov.speedLengthConstant = 0.4;
tuning.Cov.alpha = 0;
tuning.Cov.diffAlpha = 0;