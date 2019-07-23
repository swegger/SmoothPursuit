function [wms, sig, e_varHat, sumSqErr] = ...
    FitGainSDN_vel_var(wm0,sig0,s,e_vel,e_var,varargin)
%%
%
%
%%

%% Defaults
fmin_options_default = optimset('Display','iter');


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'wm0')
addRequired(Parser,'sig0')
addRequired(Parser,'s')
addRequired(Parser,'e_vel')
addRequired(Parser,'e_var')
addParameter(Parser,'fmin_options',fmin_options_default)

parse(Parser,wm0,sig0,s,e_vel,e_var,varargin{:})

wm0 = Parser.Results.wm0;
sig0 = Parser.Results.sig0;
s = Parser.Results.s;
e_vel = Parser.Results.e_vel;
e_var = Parser.Results.e_var;
fmin_options = Parser.Results.fmin_options;

%% Set up minimization
p0 = [wm0(:)', sig0];

minimizant = @(p)(varErr(p,s,e_vel,e_var));
p = fmincon(minimizant,p0,[],[],[],[],zeros(size(p0)),inf(size(p0)),[],fmin_options);

%% Best fitting estimate of variance
wms = reshape(p(1:end-1),[size(e_vel,1), size(e_vel,3)]);
sig = p(end);

e_varHat = varHat(wms,sig,s,e_vel);
sumSqErr = sum((e_var(:) - e_varHat(:)).^2);

%% Functions

%% varErr
function sumSqErr = varErr(p,s,e_vel,e_var)
    
    wms = reshape(p(1:end-1),[size(e_vel,1), size(e_vel,3)]);
    sig = p(end);
    
    e_varHat = varHat(wms,sig,s,e_vel);
    sumSqErr = sum((e_var(:) - e_varHat(:)).^2);

%% varHat
function e_varHat = varHat(wms,sig,s,e_vel)
    
    e_varHat = wms.^2.*e_vel.^2 + (sig^2 + wms.^2*sig^2).*s.^2;