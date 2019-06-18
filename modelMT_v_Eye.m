function [Rs, PVals] = modelMT_v_Eye(n,P,s,tuning,varargin)
%% modelMT_v_Eye
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'n')
addRequired(Parser,'P')
addRequired(Parser,'s')
addRequired(Parser,'tuning')

parse(Parser,n,P,s,tuning,varargin{:})

n = Parser.Results.n;
P = Parser.Results.P;
s = Parser.Results.s;
tuning = Parser.Results.tuning;

%% Z score
nM = mean(n,3);
pM = mean(P,3);
nVar = var(n,[],3);
pVar = var(P,[],3);

nZ = (n - repmat(nM,[1,1,size(n,3),1]))./repmat(sqrt(nVar),[1,1,size(n,3),1]);
pZ = (P - repmat(pM,[1,1,size(n,3),1]))./repmat(sqrt(pVar),[1,1,size(n,3),1]);

%% Correlate responses
Rs = nan(size(s,1),size(s,2),size(P,4),size(n,4),size(P,5));
RVals = nan(size(s,1),size(s,2),size(P,4),size(n,4),size(P,5));
for thetai = 1:size(s,1)
    disp(thetai/size(s,1))
    for speedi = 1:size(s,2)
        for ti = 1:size(P,4)
            for neuroni = 1:size(nZ,4)
                [rtemp,ptemp] = corrcoef(cat(2,...
                    permute(nZ(thetai,speedi,:,neuroni),[3,4,1,2]),...
                    permute(pZ(thetai,speedi,:,ti,:),[3,5,1,2,4])));
                
                Rs(thetai,speedi,ti,neuroni,:) = rtemp(1,2:3);
                PVals(thetai,speedi,ti,neuroni,:) = ptemp(1,2:3);
            end
        end
    end
end

%% Sort by % preferred speed
speedPercent = 100*repmat(s(:,:,2),[1,1,size(P,4),size(Rs,4),size(Rs,5)]) ./ ...
    repmat(permute(tuning.speed.pref,[4,2,3,1]),[size(s,1),size(s,2),size(P,4),1,2]);

dirDiff = repmat(s(:,:,1),[1,1,size(P,4),size(Rs,4),size(Rs,5)]) - ...
    repmat(permute(tuning.theta.pref,[4,2,3,1]),[size(s,1),size(s,2),size(P,4),1,2]);