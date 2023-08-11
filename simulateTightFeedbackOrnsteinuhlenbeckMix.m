function [x, C, R] = simulateTightFeedbackOrnsteinuhlenbeckMix(G,dt,input,eta,varargin)
%% simulateOrnsteinUhlenbeck
%
%
%%

%% Defaults
plotOpts_default.On = true;

%% Parse inputs

Parser = inputParser;

addRequired(Parser,'G')
addRequired(Parser,'dt')
addRequired(Parser,'input')
addRequired(Parser,'eta')
addParameter(Parser,'x0',[0 0])
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,G,dt,input,eta,varargin{:});

G = Parser.Results.G;
dt = Parser.Results.dt;
input = Parser.Results.input;
eta = Parser.Results.eta;
x0 = Parser.Results.x0;
plotOpts = Parser.Results.plotOpts;

if any(size(input) ~= size(eta))
    error('Size of input and eta must be identical')
end


%% Simulate process
T = size(input,2);
x = nan(size(input));
x(:,1) = x0(1) + x0(2)*randn(size(input,1),1);

for ti = 2:T
    slip(:,ti) = input(:,ti) - x(:,ti-1);
    slip(isnan(slip(:,ti)),ti) = 0;
    dx(:,ti) = -(1-G(2))*x(:,ti-1) + G(1)*slip(:,ti) + eta(:,ti);
    x(:,ti) = x(:,ti-1) + dx(:,ti)*dt;
end

%% Calculate covariance and correlation over time
C = cov(x);
R = corrcoef(x);

%% Plot results
if plotOpts.On
    figure('Name','First and second order statistics','Position',[1387 873 1134 449])
    subplot(1,2,1)
    samps = randsample(size(x,1),20);
    plot(linspace(0,T*dt,T),x(samps,:)','Color',[0.6 0.6 0.6])
    hold on
    plot(linspace(0,T*dt,T),mean(x,1),'k','LineWidth',2)
    plot(linspace(0,T*dt,T),mean(input,1),'k--')
    axis tight
    ylim([0 max(input(:))*1.2])
    axis square
    xlabel('Time')
    ylabel('Simulated output')
    title(['G(1) = ' num2str(G(1)) ', G(2) = ' num2str(G(2))])
    
    subplot(1,2,2)
    imagesc(linspace(0,T*dt,T),linspace(0,T*dt,T),C)
    quants = quantile(C(:),[1,99]/100);
    caxis(quants)
    axis square
    xlabel('Time')
    ylabel('Time')
end