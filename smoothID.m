function smoothID(d,varargin)
%% smoothID
%
%
%
%%

%% Defaults
t0 = 500;
twin = 100;
tvec = -twin:twin;

%% Unpack variables
hvel = [d.eye.hvel{:}];
saccades = [d.eye.saccades{:}];

x = [d.target.x{:}];
y = [d.target.y{:}];
s = [d.target.s{:}];
dir = [d.target.d{:}];
sz = [d.target.sz{:}];

xs = unique(x);
ys = unique(y);
ss = unique(s);
ds = unique(dir);
szs = unique(sz);

frac_test = 0.25;

%% Conditional means
v = hvel(t0-twin:t0+twin,:);
sac = saccades(t0-twin:t0+twin,:);
v(sac) = NaN;
for xi = 1:length(xs)
    for yi = 1:length(ys)
        for si = 1:length(ss)
            for di = 1:length(ds)
                for szi = 1:length(szs)
                    logicvec = x == xs(xi) & ...
                               y == ys(yi) & ...
                               s == ss(si) & ...
                               dir == ds(di) & ...
                               sz == szs(szi);
                    logic2 = ~any(sac,1);
                    vs{xi,yi,si,di,szi} = v(:,logicvec & logic2);
                    vs_sac{xi,yi,si,di,szi} = v(:,logicvec & ~logic2);
                    mu(:,xi,yi,si,di,szi) = nanmean(v(:,logicvec),2);
                    %Sig(:,:,xi,yi,si,di,szi) = cov(v(:,logicvec & logic2)');
                    res{xi,yi,si,di,szi} = vs{xi,yi,si,di,szi} - mu(:,xi,yi,si,di,szi);
                    res_sac{xi,yi,si,di,szi} = vs_sac{xi,yi,si,di,szi} - mu(:,xi,yi,si,di,szi);
                    test{xi,yi,si,di,szi} = false(1,size(res{xi,yi,si,di,szi},2));
                    test{xi,yi,si,di,szi}(randsample(size(test{xi,yi,si,di,szi},2),ceil(size(test{xi,yi,si,di,szi},2)*frac_test))) = true;
                end
            end
        end
    end
end

RES = [res{:}];
TEST = [test{:}];
SigAll = cov(RES(:,~TEST)');

%% Plot mean and covariance of each condition
% for xi = 1:length(xs)
%     for yi = 1:length(ys)
%         for si = 1:length(ss)
%             for di = 1:length(ds)
%                 for szi = 1:length(szs)
%                     figure
%                     subplot(1,2,1)
%                     plot(mu(:,xi,yi,si,di,szi))
%                     subplot(1,2,2)
%                     imagesc(Sig(:,:,xi,yi,si,di,szi))
%                 end
%             end
%         end
%     end
% end

%% Fit model to covariance matrix
% [T1,T2] = meshgrid(tvec);
% SigModel = @(p)( exp(-(T1-T2).^2/2/p(1)) .* cos( (T1-T2)*p(2) ) + p(3) );
% FitModel = @(p)(sum( (reshape(SigModel(p),[numel(SigAll),1]) - reshape(SigAll,[numel(SigAll),1])).^2 ));
% 
% p = fminsearch(FitModel,[4000,1/40,0]);

%% Now predict smooth pursuit in saccade windows for an example
z = res_sac{1,1,1,1,1}(:,79);
obsInds = ~isnan(z);
fObs = z(obsInds);
tObs = tvec(obsInds);
tUn = tvec(~obsInds);
sign = 0;
% [mu_, Sig_] = conditionalGaussian(zeros(size(tvec))',SigModel(p)+sign*eye(size(SigAll)),fObs(:),'x_indices',obsInds);
[mu_, Sig_] = conditionalGaussian(zeros(size(tvec))',SigAll,fObs(:),'x_indices',obsInds);

% Predict smooth pursuit w/out observation noise
N = length(z);
q = sum(obsInds);
% SigTemp = SigModel(p);
SigTemp = SigAll;
    indsOrig = 1:N;
    indsNew = [indsOrig(~obsInds) indsOrig(obsInds)];
    [Inew,Jnew] = meshgrid(indsNew);
    SigNew = nan(size(SigTemp));
    for i = 1:N
        for j = 1:N
            SigNew(i,j) = SigTemp(Inew(i,j),Jnew(i,j));     % Now covariance matrix is of form SigNew = [SigUnUn, SigUnObs; SigObsUn, SigObsObs];
        end
    end
    SigObs = SigNew((N-q+1):end,(N-q+1):end);
mu_z = SigObs/(SigObs + sign*eye(size(SigObs)))*fObs;

%% Plot predictions
figure
plot(tvec,z)
hold on
plot(tUn,mu_,'o')
plot(tObs,mu_z,'o')
%plot(tUn,mu_2)
xlabel('Time (ms)')
ylabel('Residual speed (deg/s)')

%% Now quantify the performance using trials without saccades
med_saclength = median(sum(sac(:,sum(sac,1)>0),1));
obsInds = true(size(tvec));
obsInds(ceil(length(tvec)/2)-ceil(med_saclength/2):ceil(length(tvec)/2)+ceil(med_saclength/2)) = false;
tStart = -ceil(med_saclength/2)-1;
tEnd = ceil(med_saclength/2)+1;
indStart = ceil(length(tvec)/2)-ceil(med_saclength/2)-1;
indEnd = ceil(length(tvec)/2)+ceil(med_saclength/2)+1;
tObs = tvec(obsInds);
tUn = tvec(~obsInds);

extrapModel = @(tStart,zStart,tEnd,zEnd,tUn)(zStart + ...
    (zEnd-zStart)/(tEnd-tStart)*(tUn-tStart));

ind = 1;
for xi = 1:length(xs)
    for yi = 1:length(ys)
        for si = 1:length(ss)
            for di = 1:length(ds)
                for szi = 1:length(szs)
                    for triali = 1:sum(test{xi,yi,si,di,szi})
                        ind = ind+1;
                    end
                end
            end
        end
    end
end
errs = nan(length(tUn),ind);
errs2 = nan(length(tUn),ind);
condition = nan(5,ind);

ind = 1;
for xi = 1:length(xs)
    for yi = 1:length(ys)
        for si = 1:length(ss)
            for di = 1:length(ds)
                for szi = 1:length(szs)
                    for triali = 1:size(res{xi,yi,si,di,szi},2)
                        if test{xi,yi,si,di,szi}(triali)
                            z = res{xi,yi,si,di,szi}(:,triali);
                            fObs = z(obsInds);
                            fUn = z(~obsInds);
                            %                         [mu_, Sig_] = conditionalGaussian(zeros(size(tvec))',...
                            %                             SigModel(p)+sign*eye(size(SigAll)),fObs(:),...
                            %                             'x_indices',obsInds);
                            [mu_, Sig_] = conditionalGaussian(zeros(size(tvec))',...
                                SigAll,fObs(:),...
                                'x_indices',obsInds);
                            errs(:,ind) = fUn(:) - mu_(:);
                            condition(:,ind) = [xi,yi,si,di,szi];
                            
                            Extrap = extrapModel(tStart,z(indStart),tEnd,z(indEnd),tUn);
                            errs2(:,ind) = fUn(:) - Extrap(:);
                            ind = ind+1;
                        end
                    end
                end
            end
        end
    end
end

%% Ex trial
figure
plot(tvec,z)
hold on
plotVertical(tStart);
plotVertical(tEnd);
plot(tUn,mu_,'ro')
plot(tUn,Extrap,'rx')
xlabel('Time in window (ms)')
ylabel('Residual (deg/s)')

%% Plot improvement of fit
figure('Name','Validation','Position',[119 257 1221 429])
subplot(3,3,[1,2,4,5,7,8])
plot(sqrt(mean(errs.^2,1)),sqrt(mean(errs2.^2,1)),'o')
hold on
plotUnity;
axis square
xlabel('RMSE GP model (deg/s)')
ylabel('RMSE EM (deg/s)')

subplot(3,3,3)
histogram(sqrt(mean(errs.^2,1)),linspace(0,6,100));
xlabel('RMSE GP model (deg/s)')
ax = axis;

subplot(3,3,6)
histogram(sqrt(mean(errs2.^2,1)),linspace(0,6,100));
xlabel('RMSE EM (deg/s)')
axis(ax)

subplot(3,3,9)
histogram(sqrt(mean(errs2.^2,1))-sqrt(mean(errs.^2,1)))%,linspace(-30,30,100));
hold on
plotVertical(0);
xlabel('RMSE EM - RMSE GP (deg/s)')

figure('Name','sqrt(VAR)')
for xi = 1:length(xs)
    for yi = 1:length(ys)
        for si = 1:length(ss)
            for di = 1:length(ds)
                for szi = 1:length(szs)
                    for triali = 1:sum(test{xi,yi,si,di,szi})
                        logicvec = condition(1,:) == xi &...
                                    condition(2,:) == yi &...
                                    condition(3,:) == si &...
                                    condition(4,:) == di &...
                                    condition(5,:) == szi;
                         es = errs(:,logicvec);
                         es2 = errs2(:,logicvec);
%                          plot(nanmean(es(:)),sqrt(nanvar(es(:))),'o','Color',[0.5 0.5 0.5])
%                          hold on
%                          plot(nanmean(es2(:)),sqrt(nanvar(es2(:))),'x','Color',[0.5 0.5 0.5])
                         plot(sqrt(nanvar(es(:))),sqrt(nanvar(es2(:))),'o')%,'Color',[0.5 0.5 0.5])
                         hold on
                    end
                end
            end
        end
    end
end
% plot(abs(nanmean(errs(:))),sqrt(nanvar(errs(:))),'ko','MarkerSize',10,'MarkerFaceColor',[0 0 0])
% hold on
% plot(abs(nanmean(errs2(:))),sqrt(nanvar(errs2(:))),'kx','MarkerSize',10)
plot(sqrt(nanvar(errs(:))),sqrt(nanvar(errs2(:))),'ko','MarkerSize',10,'MarkerFaceColor',[0 0 0])
% xlabel('BIAS (deg/s)')
xlabel('\sqrt{VAR} GP (deg/s)')
ylabel('\sqrt{VAR} EM (deg/s)')
axis tight
plotUnity;
axis square
% ax = axis;
% axis([0 max(ax) 0 max(ax)])