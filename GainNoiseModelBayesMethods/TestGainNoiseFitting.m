function TestGainNoiseFitting(varargin)
%% TestGainNoiseFitting
%
%   TestGainNoiseFitting()
%
%%

%% Defaults
modelparams_default.variant = 'simple';
fitparams_default.variant = 'simple';

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'ws',0.1)
addParameter(Parser,'gains',[0.4,0.6,0.8])
addParameter(Parser,'sigma',0.1)
addParameter(Parser,'wp',0)
addParameter(Parser,'trialN',2500)
addParameter(Parser,'modelparams',modelparams_default)
addParameter(Parser,'fitparams',fitparams_default)
addParameter(Parser,'testparameter','none')
parse(Parser,varargin{:})

ws = Parser.Results.ws;
gains = Parser.Results.gains;
sigma = Parser.Results.sigma;
wp = Parser.Results.wp;
trialN = Parser.Results.trialN;
modelparams = Parser.Results.modelparams;
fitparams = Parser.Results.fitparams;
testparameter = Parser.Results.testparameter;

%% Run simulation
switch modelparams.variant
    case 'simple'
        [s, m, c] = SimulateGainNoiseModel(ws,gains,sigma,0,'trials',trialN);
        g = gains(c)';
    
    case 'flexible_ws_constrained_g'
        s = [];
        m = [];
        c = [];
        for ci = 1:length(modelparams.ws_list)
            ws = modelparams.ws_list(ci);
            gains(ci) = modelparams.w0^2./(modelparams.w0^2 + ws^2);
            [stemp, mtemp, ~] = SimulateGainNoiseModel(...
                ws,gains(ci),sigma,wp,'trials',trialN);
            s = [s; stemp];
            m = [m; mtemp];
            c = [c; ci*ones(size(mtemp))];
        end      
        
        g = gains(c)';
        
    case 'full'
        [s, m, c] = SimulateGainNoiseModel(ws,gains,sigma,wp,'trials',trialN);
        g = gains(c)';
        
    otherwise
        error(['Model variant ' modelparams.variant ' not recognized!'])
end

% Calculate basic statistics
ss = unique(s);
cs = unique(c);
for si = 1:length(ss)
    for ci = 1:length(cs)
        mean_m(si,ci) = mean(m(s == ss(si) & c == cs(ci)));
        std_m(si,ci) = std(m(s == ss(si) & c == cs(ci)));
        
        z(s == ss(si) & c == cs(ci)) = (m(s == ss(si) & c == cs(ci)) - ...
            mean_m(si,ci))/std_m(si,ci);
    end
end

%% Fit data
sfit = s(m>0);
mfit = m(m>0);
cfit = c(m>0);
switch fitparams.variant
    case 'full'
        if isfield(fitparams,'params0')
            [ws_fit, sigma_fit, wp_fit, G_fit, ll] = GainNoise_fitter(sfit',mfit',cfit',...
                'params0',fitparams.params0,...
                'lb',fitparams.lb,...
                'ub',fitparams.ub);
        else
            [ws_fit, sigma_fit, wp_fit, G_fit, ll] = GainNoise_fitter(sfit',mfit',cfit',...
                'params0',[0.1,0.05,0.1,ones(1,length(gains))],...
                'lb',[0.03,0.03,0.03,zeros(1,length(gains))],...
                'ub',inf(1,3+length(gains)));
            
        end
        
        % Evaluate likelihood as a function of parameter space
        switch testparameter
            case {'none',false}
                
            case 'ws'
                % Test over ws
                wss = linspace(0.05,0.25,20);
                for i = 1:length(wss)
                    ll(i) = LogLikelihood_data_take_gainNoise(sfit',mfit',G_fit(cfit)',wss(i),wp_fit,sigma_fit);
                end
                figure;
                plot(wss,ll,'o')
                hold on
                plotVertical(ws_fit);
                
            otherwise
                error(['Likelihood evaluation parameter ' testparameter ' not recognized!'])
        end
        
    case 'simple'
        if isfield(fitparams,'params0')            
            [ws_fit, sigma_fit, G_fit, ll] = GainNoiseSimple_fitter(sfit',mfit',cfit',...
                'params0',fitparams.params0,...
                'lb',fitparams.lb,...
                'ub',fitparams.ub);
            wp_fit = 0;
            
        else
            [ws_fit, sigma_fit, G_fit, ll] = GainNoiseSimple_fitter(sfit',mfit',cfit',...
                'params0',[0.1,0.05,ones(1,length(gains))],...
                'lb',[0.03,0.03,zeros(1,length(gains))],...
                'ub',inf(1,2+length(gains)));
            wp_fit = 0;
        end
        
        % Evaluate likelihood as a function of parameter space
        switch testparameter
            case {'none',false}
                
            case 'ws'
                % Test over ws
                wss = linspace(0.05,0.25,20);
                g_fit = G_fit(c);
                for i = 1:length(wss)
                    ll(i) = LogLikelihood_data_take_gainNoiseSimple(sfit',mfit',g_fit,wss(i),sigma_fit);
                end
                figure;
                plot(wss,ll,'o')
                hold on
                plotVertical(ws_fit);
                
            otherwise
                error(['Likelihood evaluation parameter ' testparameter ' not recognized!'])
        end
        
        
    case 'flexible_ws_constrained_g'
        if isfield(fitparams,'params0')
            [ws_fit, w0_fit, wp_fit, ll] = FlexWsConG_fitter(sfit',mfit',cfit',...
                'params0',fitparams.params0,...
                'lb',fitparams.lb,...
                'ub',fitparams.ub);
            G_fit = w0_fit^2./(w0_fit^2 + ws_fit.^2);
        else
            [ws_fit, w0_fit, wp_fit, ll] = FlexWsConG_fitter(sfit',mfit',cfit',...
                'params0',[0.05,0.15,0.1*ones(1,length(gains))],...
                'lb',[0.03,0.03,0.03*ones(1,length(gains))],...
                'ub',inf(1,2+length(gains)));
            G_fit = w0_fit^2./(w0_fit^2 + ws_fit.^2);
        end
        
        % Evaluate likelihood as a function of parameter space
        switch testparameter
            case {'none',false}
                
            case 'ws'
                % Test over w0
                w0s = linspace(0.05,0.25,20);
                for i = 1:length(w0s)
                    ll(i) = LogLikelihood_data_take_FlexWsConG(sfit',mfit',cfit',w0s(i),ws_fit,wp_fit);
                end
                figure;
                plot(w0s,ll,'o')
                hold on
                plotVertical(w0_fit);
                
            otherwise
                error(['Likelihood evaluation parameter ' testparameter ' not recognized!'])
        end
        
    otherwise
        error(['Fit model variant ' fitparams.variant ' not recognized!'])
end

%% Plotting
figure
subplot(1,2,1)
for ci = 1:length(gains)
    plot(s(c==ci),m(c==ci),'o')
    hold on
end
plotUnity;
plot([min(s),max(s)],G_fit'.*repmat([min(s),max(s)],length(gains),1),'k')

xlabel('s')
ylabel('m')
mymakeaxis(gca,'interpreter','latex')

subplot(1,2,2)
for ci = 1:length(gains)
    switch fitparams.variant
        case 'full'
            
        case 'simple'
            ss = linspace(min(s),max(s),1000);
            plot(G_fit(ci)*ss,gainSDN(G_fit(ci)*ss,ss,ws_fit,sigma_fit),'k-')
        case 'flexible_ws_constrained_g'
            ss = linspace(min(s),max(s),1000);
            plot(G_fit(ci)*ss,(ws_fit(ci).^2+wp_fit(ci).^2)*G_fit(ci).^2*ss.^2,'k-')            
    end
    hold on
end
plot(mean_m,std_m.^2,'o')

xlabel('$\langle m \rangle$')
ylabel('$\langle (m - \langle m\rangle)^2 \rangle$')
mymakeaxis(gca,'interpreter','latex')

%% Functions
%% gainSDN
function v = gainSDN(ve,vs,w,sigG)
    %%
    v = w.^2.*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;
