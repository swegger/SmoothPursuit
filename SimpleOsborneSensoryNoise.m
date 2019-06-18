%% Parameters
t = 1:100;
A = -1;
B = -0.1;
I = 1;
trialN = 10000;
x(1,1,1:trialN) = 0;
x(1,2,1:trialN) = 0;
tau = 200;
W = 0.1;
t0 = 80;
WI = 2;

%% Simulate
for triali = 1:trialN
    Itemp = I + randn/WI;
    t0temp = t0 + floor(5*randn);
    for ti = 2:length(t)
        if ti > t0temp
            x(ti,1,triali) = B*x(ti-1,1,triali) + A*x(ti-1,2,triali) + Itemp + W*randn;
            x(ti,2,triali) = x(ti-1,2,triali) + x(ti-1,1,triali)/tau;
        else
            x(ti,1,triali) = B*x(ti-1,1,triali) + A*x(ti-1,2,triali) + W*randn;
            x(ti,2,triali) = x(ti-1,2,triali) + x(ti-1,1,triali)/tau;
        end
    end
end

%%
M = mean(x,3);
C = cov(permute(x(:,1,:)-M(:,1),[3,1,2]));
C2 = C - diag(diag(C));

% figure
subplot(1,2,1)
plot(M)
subplot(1,2,2)
imagesc(C)