function performanceMonitor(date,h2o,trials,weight)
%%
%
%
%
%%

d = datetime(date,'ConvertFrom','yyyymmdd');
days = -caldays(between(d,d(1),'days'));
sh2o = smooth(days,h2o,5);
strials = smooth(days,trials,5);

figure('Name','Performance','Position',[71 124 1258 700])
subplot(2,1,1)
% ax = plotyy(days,h2o/weight,d,trials,'stem','stem');
ax = plot(days,h2o/weight,'ro');
hold on
plot(days,sh2o/weight,'r-')
hold on
plotHorizontal(20);
plotVertical(3:7:max(days));
% xlabel('Date')
xlabel('Day')
% ylabel(ax(1),'ml/kg')
ylabel('ml/kg')
% ylabel(ax(2),'trials')
subplot(2,1,2)
% plot(days,log(h2o./trials),'o')
plot(days,trials,'bo')
hold on
plot(days,strials,'b-')
plotVertical(3:7:max(days));
% xlabel('Date')
xlabel('Days')
% ylabel('ln(ml/trial)')