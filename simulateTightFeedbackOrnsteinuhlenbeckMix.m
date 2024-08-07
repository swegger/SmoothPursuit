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
addParameter(Parser,'latency',0)
addParameter(Parser,'x0',[NaN,NaN])
addParameter(Parser,'t',NaN)
addParameter(Parser,'assumeSteadyState',true)
addParameter(Parser,'openLoop',false)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,G,dt,input,eta,varargin{:});

G = Parser.Results.G;
dt = Parser.Results.dt;
input = Parser.Results.input;
eta = Parser.Results.eta;
latency = Parser.Results.latency;
x0 = Parser.Results.x0;
t = Parser.Results.t;
assumeSteadyState = Parser.Results.assumeSteadyState;
openLoop = Parser.Results.openLoop;
plotOpts = Parser.Results.plotOpts;

if any(size(input) ~= size(eta))
    error('Size of input and eta must be identical')
end

if any(isnan(x0))
    x0(1) = G(1)*input(1,1)/(1-G(2)+G(1));
    x0(2) = 0;
end

T = size(input,2);
if any(isnan(t))
    t = linspace(0,T*dt,T);
end


%% Simulate process
x = nan(size(input));
x(:,1) = x0(1) + x0(2)*randn(size(input,1),1);

for ti = 2:T
    if ti > latency
        if openLoop
            slip(:,ti) = input(:,ti-latency);
        else
            slip(:,ti) = input(:,ti-latency) - x(:,ti-1);
        end
        slip(isnan(slip(:,ti)),ti) = 0;
        slip(:,ti) = slip(:,ti) + eta(:,ti);
        dx(:,ti) = -(1-G(2))*x(:,ti-1) + G(1)*slip(:,ti);
        x(:,ti) = x(:,ti-1) + dx(:,ti)*dt;
    else
        if openLoop
            slip(:,ti) = x0(1)-x(:,ti-1);
        else
            if assumeSteadyState
                tempState = input(1,1);
            else
                tempState = 0;
            end
            slip(:,ti) = tempState*ones(size(input,1),1)-x(:,ti-1);
        end
        slip(isnan(slip(:,ti)),ti) = 0;
        slip(:,ti) = slip(:,ti) + eta(:,ti);
        dx(:,ti) = -(1-G(2))*x(:,ti-1) + G(1)*slip(:,ti);
        x(:,ti) = x(:,ti-1) + dx(:,ti)*dt;
    end
end

%% Calculate covariance and correlation over time
C = cov(x);
R = corrcoef(x);

%% Plot results
if plotOpts.On
    figure('Name','First and second order statistics','Position',[1387 873 1134 449])
    subplot(1,2,1)
    samps = randsample(size(x,1),20);
    plot(t,x(samps,:)','Color',[0.6 0.6 0.6])
    hold on
    plot(t,mean(x,1),'k','LineWidth',2)
    plot(t,mean(input,1),'k--')
    axis tight
    ylim([0 max(input(:))*1.2])
    axis square
    xlabel('Time')
    ylabel('Simulated output')
    title(['G(1) = ' num2str(G(1)) ', G(2) = ' num2str(G(2))])
    
    subplot(1,2,2)
    imagesc(t,t,C)
    quants = quantile(C(:),[1,99]/100);
    caxis(quants)
    axis square
    xlabel('Time')
    ylabel('Time')
end