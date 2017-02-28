%% test_identify
var = param.V1;
est = struct2array(estimate); est = est(15+1:15+24);

figure;hold on;
scatter(var,est);
plot([0 100],[0 100],'k');
lsline;
ax = gca;
ax.YLim = [0 max([est,var])];
ax.YLim = [0 max([est,var])];
xlabel('generative parameter');
ylabel('estimated parameter');