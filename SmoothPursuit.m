function [P, pVar, P2, p2Var] = SmoothPursuit(t,e,varargin)
%% DecodeMT
%
%
%
%%

%% Defaults
% motorNoise_default.timeConstant = [0, 0];
% motorNoise_default.signalDependent = [0, 0.2];
% motorNoise_default.signalIndependent = [0.1, 0];

plantModel_default.f = @defaultPlant;
plantModel_default.parameters(1) = -1;  % velocity resistence
plantModel_default.parameters(2) = 0;    % acceleration resistence
plantModel_default.parameters(3) = 60;    % time constant
plantModel_default.noiseFraction = 0.5;
plantModel_default.noiseInd = 0;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'t')
addRequired(Parser,'e')
% addParameter(Parser,'motorNoise',motorNoise_default)
addParameter(Parser,'plantModel',plantModel_default)
addParameter(Parser,'plotflg',true)

parse(Parser,t,e,varargin{:})

t = Parser.Results.t;
e = Parser.Results.e;
% motorNoise = Parser.Results.motorNoise;
plantModel = Parser.Results.plantModel;
plotflg = Parser.Results.plotflg;


%% Direction and speed
Input = permute(repmat(e,[1,1,1,1,length(t)]),[1,2,3,5,4]);
sz = size(Input);
P = executePlant(Input,zeros(sz(1:3)),plantModel);

plantModelNoNoise = plantModel;
plantModelNoNoise.noiseFraction = 0;
plantModelNoNoise.noiseInd = 0;
P2 = executePlant(Input,zeros(sz(1:3)),plantModelNoNoise);

%% Statistics
pVar = var(P,[],3);
p2Var = var(P2,[],3);

%% Velocity and position


%% Plot
if plotflg
    figure('Name','Example speed/dir')
    subplot(1,2,1)
    plot(t,squeeze(P2(1,4,:,:,2)),'Color',[0.6 0.6 0.6]);
    hold on
    plot(t,squeeze(mean(P(1,4,:,:,2),3)),'Color',[0 0 0],'LineWidth',2);
    plotHorizontal(mean(e(1,4,:,2),3));
    axis([min(t) max(t) -2 mean(e(1,4,:,2),3)*1.5])
    
    subplot(1,2,2)
    plot(t,squeeze(P(1,4,:,:,2)),'Color',[0.6 0.6 0.6]);
    hold on
    plot(t,squeeze(mean(P(1,4,:,:,2),3)),'Color',[0 0 0],'LineWidth',2);
    plotHorizontal(mean(e(1,4,:,2),3));
    axis([min(t) max(t) -2 mean(e(1,4,:,2),3)*1.5])
    
    figure('Name','Example speeds')
    plot(t,squeeze(P(1,4,:,:,2)),'Color',[1 0.6 0.6]);
    hold on
    plot(t,squeeze(P(1,5,:,:,2)),'Color',[0.6 0.6 1]);
    plot(t,squeeze(mean(P(1,4,:,:,2),3)),'Color',[1 0 0],'LineWidth',2);
    plot(t,squeeze(mean(P(1,5,:,:,2),3)),'Color',[0 0 1],'LineWidth',2);
end

%% Functions

% %% motorNoiseDefault
% function COV = assignCOV(t,motorNoise)
% if motorNoise.timeConstant(1) == 0
%     COV(:,:,1) = diag(ones(length(t),1));
% else
%     COV(:,:,1) = nan(length(t),length(t),1);
%     for ti1 = 1:length(t)
%         for ti2 = 1:length(t)
%             COV(ti1,ti2) = (1/sqrt(2*pi*motorNoise.timeConstant(1))) * ...
%                 exp( -(ti1 - ti2)^2 / (2*motorNoise.timeConstant(1)) );
%         end
%     end
% end
% 
% if motorNoise.timeConstant(2) == 0
%     COV(:,:,2) = diag(ones(length(t),1));
% else
%     COV(:,:,2) = nan(length(t),length(t),1);
%     for ti1 = 1:length(t)
%         for ti2 = 1:length(t)
%             COV(ti1,ti2) = (1/sqrt(2*pi*motorNoise.timeConstant(2))) * ...
%                 exp( -(ti1 - ti2)^2 / (2*motorNoise.timeConstant(2)) );
%         end
%     end
% end
% 
% %% assignNoise
% function noise = assignNoise(t,e,COV,motorNoise)
% %%
% noise = nan([size(e), length(t)]);
% noise = permute(noise,[1 2 3 5 4]);
% for i = 1:size(e,1)
%     for j = 1:size(e,2)
%         noise(i,j,:,:,1) = (motorNoise.signalDependent(1) * ...
%             repmat(permute(e(i,j,:,1),[3,1,2,4]),[1,length(t)]) + ...
%             motorNoise.signalIndependent(1)) .* ...
%             mvnrnd(zeros(size(e,3),length(t)),COV(:,:,1));
%         noise(i,j,:,:,2) = (motorNoise.signalDependent(2) * ...
%             repmat(permute(e(i,j,:,1),[3,1,2,4]),[1,length(t)]) + ...
%             motorNoise.signalIndependent(1)).* ...
%             mvnrnd(zeros(size(e,3),length(t)),COV(:,:,2));
%     end
% end

%% defaultPlant
function out = defaultPlant(in,p)
out = (p(1)*in(:,:,:,2) + p(2)*(in(:,:,:,2)-in(:,:,:,1)));

%% executePlant
function out = executePlant(Input,out0,plantModel)
%%
sz = size(Input);
sz = sz(1:3);
out = nan(size(Input));
out(:,:,:,:,1) = Input(:,:,:,:,1);
out(:,:,:,1:2,2) = repmat(out0,[1,1,1,2]);
for ti = 3:size(Input,4)
    dOut = plantModel.f(out(:,:,:,ti-2:ti-1,2),plantModel.parameters);
    out(:,:,:,ti,2) = out(:,:,:,ti-1,2) + ( dOut + ...
        Input(:,:,:,ti-1,2) + ...
        (plantModel.noiseFraction.*Input(:,:,:,ti-1,2) + ...
        plantModel.noiseInd).*randn(sz) )/plantModel.parameters(3);
end