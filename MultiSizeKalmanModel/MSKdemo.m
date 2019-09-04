function MSKdemo(varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'T',200)
addParameter(Parser,'Ns',[2,4,8])
addParameter(Parser,'zs',4:4:20)
addParameter(Parser,'gainNoise',0)
addParameter(Parser,'t0',1)
addParameter(Parser,'Wx',1)
addParameter(Parser,'wT0',0.1)
addParameter(Parser,'Wt',0)


parse(Parser,varargin{:})

T = Parser.Results.T;
Ns = Parser.Results.Ns;
zs = Parser.Results.zs;
gainNoise = Parser.Results.gainNoise;
t0 = Parser.Results.t0;
Wx = Parser.Results.Wx;
wT0 = Parser.Results.wT0;
Wt = Parser.Results.Wt;

Zs = repmat(zs(:),[1,T]);

colors = projectColorMaps('speeds','samples',1:length(zs),'sampleDepth',5);

%% Run simulation
for ni = 1:length(Ns)
    for zi = 1:length(zs)
        [zstar(:,:,:,zi,ni), wT(:,:,:,zi,ni), ~, K(:,:,zi,ni)] = ...
            MultiSizeSpatialIntegrator(Ns(ni),'z',Zs(zi,:),'Nmax',max(Ns),...
            'gainNoise',gainNoise,'t0',t0,...
            'Wx',Wx*ones(Ns(ni),1),'Wt',Wt*ones(Ns(ni),1),'wT0',wT0*ones(Ns(ni),1));
    end
end

%% Summary statistics
mZstar = mean(zstar,3);
varZstar = var(zstar,1,3);

for ni = 1:length(Ns)
    for zi = 1:length(zs)
        zstarCOV(:,:,zi,ni) = cov(permute(zstar(:,:,:,zi,ni),[3,2,1,4,5]));
    end
end

%% Plotting
figure('Name','Estimate over time','Position',[405 281 1885 420])
for ni = 1:length(Ns)
    subplot(1,length(Ns),ni)
    
    for zi = 1:length(zs)
        patchProps.FaceColor = colors(zi,:);
        patchProps.FaceAlpha = 0.5;
        h = myPatch([0:T-1]',squeeze(mean(zstar(:,:,:,zi,ni),3))',...
            squeeze(std(zstar(:,:,:,zi,ni),1,3))',...
            'patchProperties',patchProps);
        hold on
    end
    
    for zi = 1:length(zs)
        plot(0:T-1,squeeze(mean(zstar(:,:,:,zi,ni),3)),'-',...
            'Color',colors(zi,:),'LineWidth',2)
        hold on
    end
    axs(ni,:) = axis;
end
for ni = 1:length(Ns)
    subplot(1,length(Ns),ni)
    axis([min(axs(:,1)) max(axs(:,2)) min(axs(:,3)) max(axs(:,4))])
    xlabel('Time')
    ylabel('Estimate')
    mymakeaxis(gca,'xytitle',['N = ' num2str(Ns(ni))])
end

%%
figure('Name','Mean vs variance in estimates','Position',[405 281 1885 420])
for ni = 1:length(Ns)
    subplot(1,length(Ns),ni)
    for zi = 1:length(zs)
        plot(squeeze(mean(zstar(:,:,:,zi,ni),3)),squeeze(var(zstar(:,:,:,zi,ni),1,3)),'o',...
            'Color',colors(zi,:),'MarkerFaceColor',colors(zi,:))
        hold on
    end
    axs(ni,:) = axis;
end
for ni = 1:length(Ns)
    subplot(1,length(Ns),ni)
    axis([min(axs(:,1)) max(axs(:,2)) min(axs(:,3)) max(axs(:,4))])
    if ni == length(Ns)
        text(0.8*max(axs(:,2)),0.1*max(axs(:,4)),['\sigma_g = ' num2str(gainNoise)])
        text(0.8*max(axs(:,2)),0.2*max(axs(:,4)),['w_x = ' num2str(Wx)])
        text(0.8*max(axs(:,2)),0.3*max(axs(:,4)),['w_T(0) = ' num2str(wT0)])
    end
    xlabel('Estimate mean')
    ylabel('Estimate variance')
    mymakeaxis(gca,'xytitle',['N = ' num2str(Ns(ni))])
end

%% Relative estimates over time
figure('Name','Relative mean estimates')
for ni = 2:3
    plot(squeeze(mean(zstar(:,:,:,:,1),3)),squeeze(mean(zstar(:,:,:,:,ni),3)),...
        'o','Color',(1-(ni-1)/length(Ns))*ones(1,3),'MarkerFaceColor',(1-(ni-1)/length(Ns))*ones(1,3))
    hold on
end
plotUnity;
axis square

%% Temporal integration gain over time
figure('Name','Temporal integration gain')
for ni = 1:length(Ns)
    plot(0:T-1,log(K(:,1,1,1)),'k')
    hold on
end
xlabel('Time')
ylabel('$\ln (K)$')
mymakeaxis(gca,'interpreter','latex')

%% Covariance
figure('Name','Estimate temporal covariance','Position',[341 76 1214 1183])
ind = 0;
for zi = 1:length(zs)
    for ni = 1:length(Ns)
        ind = ind+1;
        subplot(length(zs),length(Ns),ind)
        imagesc(0:T-1,0:T-1,zstarCOV(:,:,zi,ni))
        colorbar
        axis square
        title(['N = ' num2str(Ns(ni)) ', z = ' num2str(zs(zi))])
    end
end