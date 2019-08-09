function [e, likelihood] = ScalarBayesSpeedEstimator(m,wm,sigm,xmin,xmax,varargin)
%% ScalarBayesEstimators
%
%   e = ScalarBayesSpeedEstimator(m,wm,sigm,xmin,xmax)
%       Computes the BLS estimates (e) of x, given a set of measurements
%       (m). Assumes uniform prior over x between xmin and xmax
%
%   e = ScalarBayesSpeedEstimator(m,wm,sigm,xmin,xmax,'method',method_opts)
%       Computes the BLS estimate using the the method specifiied by
%       method_opts.
%           method_opts.type = 'integral'
%               Use integral functions, default
%           method_opts.type = 'trapz'
%               Approximates the integral using trapz using the step
%               size method_opts.dx
%           method_opts.type = 'quad'
%               Approximates the integral using Simpson's quadrature using
%               the step size method_opts.dx
%
%   e = ScalarBayesSpeedEstimator(...,'estimator',estimator_opts)
%       Computes the estimate specified by estimator_opts.type
%           estimator_opts.type = 'BLS'
%               Bayes least-squares; default.
%           estimator_opts.type = 'ObsAct'
%               Observer-Actor model. User must supply the weber fraction
%               on production in the estimator_opts.wp
%
%   e = ScalarBayesSpeedEstimator(m,wm,sigm,xmin,xmax,'prior',prior_opts)
%       Computes the BLS estimate using a prior specified by the structure
%       prior_opts.
%           TODO: prior types
%
%   e = ScalarBayesSpeedEstimator(...,'estimator',estimator_opts)
%       with
%           estimator_opts.type = 'weightedMean'
%       and
%           estimator_opts.weights = weights
%       Computs the BLS estimate using a weighted average of the
%       measurements according a set of weights.
%           TODO: integral and trapz support
%
%   e = ScalarBayesSpeedEstimator(...,'estimator',estimator_opts)
%       with
%           estimator_opts.type = 'sequential'
%       and
%           estimator_opts.transitionFunction = function handle for
%           transition probabilities
%       and
%           estimator_opts.p = vector of parameters for transition function
%       Computs the estimate using a sequential message passing
%       algorithm. In brief, updates the joint probability distribution of
%       the measurement and sample according to
%           p(x(t),m(1),m(2),...,m(t)) = p(m(t)|x(t))...
%               *int/(p(x(t-1),m(1),m(2),...,m(t-1))*p(x(t)|x(t-1))dx(t-1))
%       Algorithm is equivalent to the foward-backward message passing
%       algorithm (see Bishop, 2006; chapter 13), without the backward
%       passes. Estimate is then calculated as the expected value of x(t). 
%           TODO: integral and trapz support
%
%   NOTE: Methods requiring mmx.mex are mostly obsolete.
%
%%

%% Defaults
method_opts.type = 'integral';
estimator_opts.type = 'BLS';
prior_opts.type = 'uniform';

%% Parse inputs
inputs = inputParser;
addRequired(inputs,'m');
addRequired(inputs,'wm');
addRequired(inputs,'sigm');
addRequired(inputs,'xmin');
addRequired(inputs,'xmax');
addParameter(inputs,'method',method_opts,@isstruct);
addParameter(inputs,'estimator',estimator_opts,@isstruct);
addParameter(inputs,'prior',prior_opts,@isstruct);

parse(inputs,m,wm,sigm,xmin,xmax,varargin{:})

m = inputs.Results.m;
wm = inputs.Results.wm;
sigm = inputs.Results.sigm;
xmin = inputs.Results.xmin;
xmax = inputs.Results.xmax;
method = inputs.Results.method;
estimator = inputs.Results.estimator;
prior = inputs.Results.prior;       %% TODO generalize to aribitrary prior


if isfield(estimator,'wy')
    wy = estimator.wy;
else
    wy = 0;
end

if ~isfield(estimator,'ObsAct')
    estimator.ObsAct = 0;
end

%% Compute the estimate
switch estimator.type
    case {'BLS','BLSbiasedLapse','BLS_wm_wp_sigp'}
        switch method.type
            case 'integral'
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4)/(1+wy.^2)^estimator.ObsAct;
                e = permute(e,[2 1]);
                
            case 'quad'
                switch prior.type
                    case 'uniform'
                        % Number of measurements
                        N = size(m,2);
                        
                        % Create x-vector
                        dx = method.dx;
                        x = xmin:dx:xmax;
                        
                        % Create Simpson's nodes
                        l = length(x);
                        h = (xmax - xmin)/l;
                        w = ones(1,l);
                        w(2:2:l-1) = 4;
                        w(3:2:l-1) = 2;
                        w = w*h/3;
                        
                        % Reshape measurements for processing
                        M = permute(m,[2 3 1]);
                        M = repmat(M,[1,1,1,l]);
                        x = reshape(x,[1 1 1 l]);
                        X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                        
                        % Generate estimate
                        w = reshape(w,[1 1 1 l]);
                        w = repmat(w,[1 1 size(m,1) 1]);
                        likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(sum((X-M).^2,1))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
                        %                likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
                        e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                        e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                        
                    case 'Gaussian'
                        % Parameterize prior
                        mu = prior.mu;
                        sig = prior.sig;
                        
                        % Number of measurements
                        N = size(m,2);
                        
                        % Create x-vector
                        dx = method.dx;
                        x = mu-3*sig:dx:mu+3*sig;
                        
                        % Create Simpson's nodes
                        l = length(x);
                        h = (xmax - xmin)/l;
                        w = ones(1,l);
                        w(2:2:l-1) = 4;
                        w(3:2:l-1) = 2;
                        w = w*h/3;
                        
                        % Reshape measurements for processing
                        M = permute(m,[2 3 1]);
                        M = repmat(M,[1,1,1,l]);
                        x = reshape(x,[1 1 1 l]);
                        X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                        P = 1/sqrt(2*pi*sig^2) * exp( -(X(1,:,:,:)-mu).^2/(2*sig^2));
                        
                        % Generate estimate
                        w = reshape(w,[1 1 1 l]);
                        w = repmat(w,[1 1 size(m,1) 1]);
                        likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(sum((X-M).^2,1))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
                        %                likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
                        e = sum(w.*X(1,:,:,:).*likelihood.*P,4)./sum(w.*likelihood.*P,4);
                        e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                       
                    otherwise
                        error(['Prior ' prior.type ' not recognized!'])
                end
                
            case 'MonteCarlo'
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo_batch'
                % Set up integration variables
                options.N = method.N;
                options.batch_sz = method.batch_sz;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
        end
        
    
        
    
end

%% Functions

function out = MonteCarloIntegrand_numerator(x,m,wm)

% Reshape measurements for processing
N = size(m,2);
l = length(x);
M = permute(m,[2 3 1]);
M = repmat(M,[1,1,1,l]);
x = reshape(x,[1 1 1 l]);
X = repmat(x,[size(M,1) 1 size(M,3) 1]);

% compute likelihood
likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(sum((X-M).^2,1))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
%likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
out = squeeze(permute(X(1,:,:,:).*likelihood,[4 3 2 1]));


function out = MonteCarloIntegrand_denominator(x,m,wm)

% Reshape measurements for processing
N = size(m,2);
l = length(x);
M = permute(m,[2 3 1]);
M = repmat(M,[1,1,1,l]);
x = reshape(x,[1 1 1 l]);
X = repmat(x,[size(M,1) 1 size(M,3) 1]);

% compute likelihood
likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(sum((X-M).^2,1))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
%likelihood = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
out = squeeze(permute(likelihood,[4 3 2 1]));

function out = posteriorMAP(x,m,wm,xmin,xmax)

N = size(m,2);
out = -(1./sqrt(2*pi)./wm./x).^N .* exp( -(sum((x-m).^2))./(2*wm.^2.*x.^2) );
out(x < xmin | x > xmax) = 0;

function out = posteriorMLE(x,m,wm,xmin,xmax)

N = size(m,2);
out = -(1./sqrt(2*pi)./wm./x).^N .* exp( -(sum((x-m).^2))./(2*wm.^2.*x.^2) );


function l = likelihoods(M,X,wms)

ls = nan(size(X));
for wmi = 1:length(wms)
    wm = wms(wmi);
    ls(wmi,:,:,:) = ( (1./sqrt(2*pi*(sigm^2+wm^2*X(1,:,:,:).^2))) .* ...
        exp( -(X(1,:,:,:)-M(wmi,:,:,:)).^2 ./...
        (2*(sigm^2+wm^2*X(1,:,:,:).^2)) ) );
end
l = prod(ls,1);
