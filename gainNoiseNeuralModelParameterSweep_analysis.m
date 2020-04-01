function gainNoiseNeuralModelParameterSweep_analysis(varargin)
%% gainNoiseNeuralModelParameterSweep_analysis
%
%   gainNoiseNeuralModelParameterSweep_analysis()
%   
%   Performs analysis of neural model output for default results directory
%
%%

%% Defaults
dataDirectory_default = '~/Projects/MultiSizePursuit/parameterSweeps';
dataDate_default = '20200326';


%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dataDirectory',dataDirectory_default)
addParameter(Parser,'dataDate',dataDate_default)

parse(Parser,varargin{:})

dataDirectory = Parser.Results.dataDirectory;
dataDate = Parser.Results.dataDate;

%% Find files corresponding to desired date
potentialFiles = dir([dataDirectory '/*' dataDate '*.mat']);
files = potentialFiles(3:end);

%% For each file, load in parameter values and results
params = nan(length(files),4);
ws = nan(length(files),1);
sigGs = nan(length(files),1);
Gs = nan(length(files),3);

for filei = 1:length(files)
    results = load(files(filei).name,'sizeProps','w','sigG','G','gainNoise');
    params(filei,1) = results.sizeProps.exponential;
    params(filei,2) = results.sizeProps.threshold;
    params(filei,3) = results.sizeProps.surround_weight;
    params(filei,4) = results.gainNoise;
    ws(filei) = results.w;
    sigGs(filei) = results.sigG;
    Gs(filei,:) = results.G;
end

%% Load in experimental fit
sigGobs = 0.1;

%%
figure('Name','Parameter effect on sigG')
subplot(2,3,1)
scatter3(params(:,1),params(:,2),params(:,3),10,sigGs)
xlabel('Exp')
ylabel('Thres')
zlabel('SS')
subplot(2,3,2)
scatter(params(:,1),params(:,2),10,sigGs)
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
scatter(params(:,1),params(:,3),10,sigGs)
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
scatter(params(:,2),params(:,3),10,sigGs)
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
histogram(sigGs,linspace(0,0.2,100))
plotVertical(mean(sigGs));
hold on
plotVertical(sigGobs);

%%
figure('Name','Parameter effect on w')
subplot(2,3,1)
scatter3(params(:,1),params(:,2),params(:,3),10,ws)
xlabel('Exp')
ylabel('Thres')
zlabel('SS')
subplot(2,3,2)
scatter(params(:,1),params(:,2),10,ws)
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
scatter(params(:,1),params(:,3),10,ws)
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
scatter(params(:,2),params(:,3),10,ws)
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
histogram(ws,linspace(0,0.2,100))
% plotVertical(mean(ws));
hold on
% plotVertical(sigGobs);

%%
figure('Name','Parameter effect on Gain')
subplot(1,3,1)
plot3(params(:,1),params(:,2),Gs,'o')
xlabel('Exp')
ylabel('Thres')
subplot(1,3,2)
plot3(params(:,1),params(:,3),Gs,'o')
xlabel('Exp')
ylabel('SS')
subplot(1,3,3)
plot3(params(:,2),params(:,3),Gs,'o')
xlabel('Thres')
ylabel('SS')