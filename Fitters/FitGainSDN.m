function [g, sigma_g, sigma_s, w, wp] = FitGainSDN(e,s,varargin)
%% FitGainSDN
%
%   [g, sigma_g, sigma_s, w] = FitGainSDN(e,s)
%
%   Fits a model of the form:
%       e = (g + eta_g)m
%
%   where
%       eta_g ~ N(0,sigma_g^2)
%       m = s + eta_s
%       eta_s ~ N(0,sigma_s^2 + w^2s^2)
%
%   to data using maximium likelihood.
%
%%

%% Defaults
fmin_options_default = optimset('Display','iter');
options_default.dx = 0.05;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'e')
addRequired(Parser,'s')
addParameter(Parser,'MinMax',NaN)
addParameter(Parser,'p0',[1 0 0 0 0])
addParameter(Parser,'fmin_options',fmin_options_default)
addParameter(Parser,'options',options_default)

parse(Parser,e,s,varargin{:})

e = Parser.Results.e;
s = Parser.Results.s;
MinMax = Parser.Results.MinMax;
p0 = Parser.Results.p0;
fmin_options = Parser.Results.fmin_options;
options = Parser.Results.options;

%% Set measurements
if any(isnan(MinMax(:)))
    MinMax = [0.1*min(s), 1.5*max(s);...
              0.05,       1.5];
end

%% Maximize log likelihood
minimizant = @(p)(-loglikelihood(p,e,s,MinMax,options));
p = fmincon(minimizant,p0,[],[],[],[],[0 0 0 0 0],[Inf Inf Inf Inf Inf],[],fmin_options);

g = p(1);
sigma_g = p(2);
sigma_s = p(3);
w = p(4);
wp = p(5);

%% Functions

%% loglikelihood
function ll = loglikelihood(p,e,s,MinMax,options)
%%

    functionHandle = @(x)( likelihoodf(x,p,e,s) );
    like = ndintegrate(functionHandle,MinMax,'method','quad',...
        'options',options);
    
    like(like == 0) = realmin;
    ll = nansum( log(like), 2 );
    
 %% likelihood
 function out = likelihoodf(x,p,e,s)
     %%
     % x(:,1) - measurement; x(:,2) - gain
    X = repmat(permute(x,[1,3,2]),[1,numel(e),1]);
    
    S = repmat(s(:)',[size(X,1), 1, 1]);
    E = repmat(e(:)',[size(X,1), 1, 1]);
    Ehat = X(:,:,2).*X(:,:,1);
        
    g = p(1);
    sigma_g = p(2);
    sigma_s = p(3);
    wm = p(4);
    wp = p(5);
    
    p_e_take_ehat = (1./sqrt(2*pi*wp^2*Ehat)) .* exp( -(E-Ehat).^2./(2*wp^2*Ehat) );;
    
    p_g = (1./sqrt(2*pi*sigma_g.^2)) .* ...
        exp( -(g - X(:,:,2)).^2./(2*sigma_g.^2) );
    
    p_m_take_s = (1./sqrt(2*pi*(sigma_s.^2 + wm.^2.*s.^2))) .* ...
        exp( -(X(:,:,1) - S).^2 ./(2*(sigma_s.^2 + wm.^2.*s.^2)) );
    
    out = p_g.*p_e_take_ehat.*p_m_take_s;
%     out = p_g.*p_m_take_s;