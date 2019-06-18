function EvalBIAS_VAR()
%%
%
%
%%

%% Defaults

%% Variables

t = 0:250;
profileN = 10000;
s = 12;
e0 = 0;
ws = 0.1;
we = 0.05;

% Gain prior
a = 0.0006;
b = 5.33;
c = 123.88;
d = 1;

sigf = 2;
sign = 0;
tau = 50;

%% Select temporal gain profiles
COV = setCOV(t,sigf,sign,tau);
MU = real(( d*(t-a).^b ./ ( (t-a).^b + c^b ) ));
MU = [0 diff(MU)];
MU = 5*MU;

Gs = mvnrnd(MU,COV,profileN);
% Gs = mvnrnd(10*ones(size(t)),COV,profileN);
% Gs(:,t<50) = 0;

%% Run each profile and evaluate RMSE
for i = 1:profileN
    [e(:,i) BIAS(:,i), VAR(:,i)] = eyeDynamics(s,e0,Gs(i,:),t,ws,we);
    RMSE(:,i) = sqrt( BIAS(:,i).^2 + VAR(:,i).^2 );
end


%% Plot something interesting
% figure
% h = plot(t,Gs);
sRMSE = log(nanmean(RMSE,1));
for i = 1:profileN
    if isnan(sRMSE(i)) || isinf(sRMSE(i))
        wts(i) = 0;
    else
        wts(i) = min(sRMSE)/sRMSE(i);
    end
end
% for i = 1:profileN
%     h(i).Color = [wts(i) 0.5 0.5];
% end

figure
plot(t,nansum( repmat(wts(:)/sum(wts),[1,length(t)]).*Gs,1 ))
hold on
bestInd = find(wts == 1);
plot(t,Gs(bestInd,:),'r')
plot(t,MU,'k')

%%
mG = nansum( repmat(wts(:)/sum(wts),[1,length(t)]).*Gs,1);
[me mBIAS, mVAR] = eyeDynamics(s,e0,mG,t,ws,we);

%% Functions
function [e, BIAS, VAR] = eyeDynamics(s,e0,G,t,ws,we)
    e(1) = e0;
    BIAS(1) = s;
    VAR(1) = 0;
    for ti = 2:length(t)
        if t(ti) < 50
            e(ti) = 0;
            BIAS(ti) = s;
            VAR(ti) = 0;
        else
            e(ti) = e(ti-1) + G(ti)*(s - e(ti-1));
            BIAS(ti) = evalBIAS(s,e(ti),G(ti));
            VAR(ti) = evalVAR(s,e(ti),G(ti),ws,we);
        end
    end
    
function BIAS = evalBIAS(s,e,G)
    BIAS = abs((1-G).*(s-e));
    
function VAR = evalVAR(s,e,G,ws,we)
    VAR = we^2*e.^2 + ws^2*G.^2.*s.^2 + we^2*G.^2.*e.^2 - 2*we^2*G.*e.^2;
    
function COV = setCOV(t,sigf,sign,tau)
    [T1,T2] = meshgrid(t);
    COV = reshape( sigf^2 * exp( -(T1(:) - T2(:)).^2/(2*tau^2) ), ...
        size(T1) );
    COV = COV + diag(sign^2*ones(size(t)));
    