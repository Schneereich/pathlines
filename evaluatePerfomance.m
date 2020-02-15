function evaluatePerfomance
close all
%*** start with first, end with last = N-1
first = 100;
N = 110;

methods = ["explicitEuler","implicitEuler","betterEuler","eulerHeun", ...
           "crankNicolson","rungeKutta4Newton", "rungeKutta4Classic"];
nMethods = length(methods);
lineStyle = {'--o','-d','-.*','--s',':p','-.x','-^'};

%*** solve and time pathlines for multiple substeps
subSteps = 50:50:500; % [50, 100, 150];
nSubSteps = length(subSteps);
solveTimes = zeros(nSubSteps, nMethods);
for m = 1:nMethods
    for s = 1:nSubSteps
        solveTimes(s,m) = mean( ...
            postprocessing_mini('method',methods(m),...
                                'first',first, ...
                                'N',N, ...
                                'subSteps',subSteps(s)));
    end
    %*** compare runtime of solvers for increasing substeps
    figure(3)
    plot(subSteps,solveTimes(:,m),lineStyle{m},'LineWidth',1)
    hold on
end

legend(methods,'Location','northwest')
title('Runtime Comparison')
xlabel('substeps in each timestep')
ylabel('mean runtime in seconds')
ax = gca;
ax.YAxis.Exponent = 0;

%*** compare runtime for fixed substeps, sort by descendent runtime
[sortedTimes, sortedIdx] = sort(solveTimes(end,:),'descend');
figure
stem(sortedTimes,'linewidth',2)
title(['Mean Runtime of Pathline Iteration for ', ...
       num2str(subSteps(end)),' substeps'])
ylabel('runtime in seconds')
xticklabels(methods(sortedIdx))
xtickangle(45)
ax = gca;
ax.YAxis.Exponent = 0;