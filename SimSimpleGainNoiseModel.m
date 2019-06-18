function SimSimpleGainNoiseModel(varargin)
%% SimSimpleGainNoiseModel
%
%   SimSimpleGainNoiseModel()
%
%   Simulates simple sensorimotor gain noise model of the form:
%           e = (g + eta_g)(s + eta_s);
%
%   with the covariance between eta_g and eta_s set to zero. Simulates a
%   number of values of g and s for the default noise set to Gaussian,
%   independent, with mean 0 and standard deviation set to 0.1 and 1,
%   respectively.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'sigma_g',0.1)
addParameter(Parser,'sigma_s',1)
addParameter(Parser,'sigma_gs',0)
addParameter(Parser,'gs',linspace(0,1,10))
addParameter(Parser,'ss',linspace(0,10,10))
addParameter(Parser,'sampleN',1000)

parse(Parser,varargin{:})

sigma_g = Parser.Results.sigma_g;
sigma_s = Parser.Results.sigma_s;
sigma_gs = Parser.Results.sigma_gs;
gs = Parser.Results.gs;
ss = Parser.Results.ss;
sampleN = Parser.Results.sampleN;

%% For each combination of g and s generate samples and eye outputs

% All combinations
[Gs,Ss] = meshgrid(gs,ss);

% Generate noise
COV = [sigma_g^2, sigma_gs^2; sigma_gs^2, sigma_s^2];
Noise = mvnrnd([0,0],COV,sampleN);

% Calculate e
Es = (repmat(Gs,[1,1,sampleN]) + ...
    repmat(permute(Noise(:,1),[3,2,1]),[size(Gs,1),size(Gs,2),1]))...
    .*...
    (repmat(Ss,[1,1,sampleN]) + ...
    repmat(permute(Noise(:,2),[3,2,1]),[size(Ss,1),size(Ss,2),1]));

% Calculate stats
mE = mean(Es,3);
vE = var(Es - repmat(mE,[1,1,sampleN]),[],3);

%% Theoretical results
M = Gs.*Ss;
VAR = sigma_s^2*Gs.^2 + sigma_s^2*sigma_g^2 + sigma_g^2*Ss.^2;

%% Analyze
% figure
% subplot(1,2,1)
% surf(Gs,Ss,mE,'EdgeColor','none')
% xlabel('Gain')
% ylabel('Speed')
% 
% subplot(1,2,2)
% surf(Gs,Ss,M,'EdgeColor','none')
% xlabel('Gain')
% ylabel('Speed')

%%
figure
subplot(1,2,1)
surf(Gs,Ss,vE,'EdgeColor','none')
xlabel('Gain')
ylabel('Speed')

subplot(1,2,2)
surf(Gs,Ss,VAR,'EdgeColor','none')
xlabel('Gain')
ylabel('Speed')
