function TestFitGainSDN(varargin)
%% TestFitGainSDN
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'sigma_g',0.1)
addParameter(Parser,'sigma_s',0.1)
addParameter(Parser,'w',0.1)
addParameter(Parser,'wp',0.01)
addParameter(Parser,'sigma_gs',0)
addParameter(Parser,'gs',0.5)
addParameter(Parser,'ss',linspace(4,20,5))
addParameter(Parser,'sampleN',1000)
addParameter(Parser,'dm',0.05)

parse(Parser,varargin{:})

sigma_g = Parser.Results.sigma_g;
sigma_s = Parser.Results.sigma_s;
w = Parser.Results.w;
wp = Parser.Results.wp;
sigma_gs = Parser.Results.sigma_gs;
gs = Parser.Results.gs;
ss = Parser.Results.ss;
sampleN = Parser.Results.sampleN;
dm = Parser.Results.dm;

%% Simulate gainSDN model

% Generate noise
COV = [sigma_g^2, sigma_gs^2; sigma_gs^2, sigma_s^2];
Noise = mvnrnd([0,0],COV,sampleN);

% Calculate e
[Ss,Gs] = meshgrid(ss,gs);
Es = (repmat(Gs,[1,1,sampleN]) + ...
    repmat(permute(Noise(:,1),[3,2,1]),[size(Gs,1),size(Gs,2),1]))...
    .*...
    (repmat(Ss,[1,1,sampleN]) + ...
    repmat(permute(Noise(:,2),[3,2,1]),[size(Ss,1),size(Ss,2),1]) + ...
    + w*repmat(Ss,[1,1,sampleN]).*randn(size(Ss,1),size(Ss,2),sampleN));

S = repmat(Ss,[1,1,sampleN]);
s = S(:)';
e = Es(:)' + wp*Es(:)'.*randn(size(Es(:)'));

% Calculate stats
mE = mean(Es,3);
vE = var(Es - repmat(mE,[1,1,sampleN]),[],3);

%% Fit gainSDN model
options.dx = 0.01;
p0 = [gs,sigma_g,sigma_s,w,wp];
[phat(1), phat(2), phat(3), phat(4), phat(5)] = FitGainSDN(e,s,'p0',p0,'options',options);

%% Model predictions
VARa = (sigma_s^2+w^2*ss.^2).*gs.^2 + sigma_s^2*sigma_g^2 + sigma_g^2*ss.^2 + wp^2*ss.^2;
VAR = (phat(3)^2+phat(4)^2*ss.^2).*phat(1).^2 + phat(3)^2*phat(2)^2 + phat(2)^2*ss.^2 + phat(5)^2*ss.^2;

%%
sig_gs = linspace(0.05,5*sigma_g,10);
sig_ss = linspace(0.05,2*sigma_s,10);
logLike = llSurf(sig_gs,sig_ss,gs,w,e,s,options);

%% Functions
%% llSurf
function logLike = llSurf(sig_gs,sig_ss,g,w,e,s,options)
%%
    MinMax = [0.1*min(s) 1.5*max(s);...
              0.01       1.1        ];
    
    for i = 1:length(sig_gs)
        for j = 1:length(sig_ss)
            disp([i,j])
            logLike(i,j) = loglikelihood([g,sig_gs(i),sig_ss(j),w],e,s,...
                MinMax,options);
        end
    end

%% loglikelihood
function ll = loglikelihood(p,e,s,MinMax,options)
%%

    functionHandle = @(x)( likelihoodf(x,p,e,s) );
    like = ndintegrate(functionHandle,MinMax,'method','quad',...
        'options',options);
    
    
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