function [wms, sig, e_varHat, sumSqErr, wm0, wm] = ...
    FitGainSDN_vel_var_wm(wm0,t0,sig0,s,e_vel,e_var,varargin)
%%
%
%
%%

%% Defaults
fmin_options_default = optimset('Display','iter');


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'wm0')
addRequired(Parser,'t0')
addRequired(Parser,'sig0')
addRequired(Parser,'s')
addRequired(Parser,'e_vel')
addRequired(Parser,'e_var')
addParameter(Parser,'fmin_options',fmin_options_default)

parse(Parser,wm0,t0,sig0,s,e_vel,e_var,varargin{:})

wm0 = Parser.Results.wm0;
t0 = Parser.Results.t0;
sig0 = Parser.Results.sig0;
s = Parser.Results.s;
e_vel = Parser.Results.e_vel;
e_var = Parser.Results.e_var;
fmin_options = Parser.Results.fmin_options;

%% Set up minimization
p0 = [wm0, sig0];

minimizant = @(p)(varErr(p,s,e_vel,e_var,t0));
p = fmincon(minimizant,p0,[],[],[],[],zeros(size(p0)),inf(size(p0)),[],fmin_options);

%% Best fitting estimate of variance
wms = wmT(1e100,p(1),t0,size(e_vel,1));
sig = p(end);
wm0 = 1e100;
wm = p(1);

% Ks = wms.^2./(wms.^2 + repmat(wm.^2,[size(wms,1),1]));
% Ss(1,:) = zeros(size(s(1,:)));
% for ti = 2:length(e_vel)
%     Ss(ti,:) = Ss(ti-1,:) + Ks(ti)*s(ti,:);
% end
% e_varHat = varHat(wms,sig,Ss,e_vel);

e_varHat = varHat(wms,sig,s,e_vel);
sumSqErr = sum((e_var(:) - e_varHat(:)).^2);

%% Functions

%% varErr
function sumSqErr = varErr(p,s,e_vel,e_var,t0)
    
    wms = wmT(1e100,p(1),t0,size(e_vel,1));
    sig = p(end);
    
%     Ks = wms.^2./(wms.^2 + repmat(p(1).^2,[size(wms,1),1]));
%     e_varHat = varHat(wms,sig,Ss,e_vel);
%     Ss(1,:) = zeros(size(s(1,:)));
%     for ti = 2:length(e_vel)
%         Ss(ti,:) = Ss(ti-1,:) + Ks(ti)*s(ti,:);
%     end
%     e_varHat = varHat(wms,sig,Ss,e_vel);
    
    e_varHat = varHat(wms,sig,s,e_vel);
    
    sumSqErr = sum((e_var(:) - e_varHat(:)).^2);

%% varHat
function e_varHat = varHat(wms,sig,s,e_vel)
    
    e_varHat = wms.^2.*e_vel.^2 + (sig^2 + wms.^2*sig^2).*s.^2;
    
%% wmT
function wms = wmT(wm0,wm,t0,T)

    wms(1) = wm0;
    for ti = 2:(t0+T-1)
        wms(ti,1) = sqrt((wms(ti-1,1)^2*wm^2) / (wms(ti-1,1)^2 + wm^2));
    end
    wms = wms(t0:end);