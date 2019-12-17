function integratedStreakComparison(filenameConv,filenameUnc)
%% integratedStreakComparison 
% plots the streak data of converted an unconverted light over each other and
% also the difference.  
%
%   Input Arguments:
%       filenames: The files should be in .mat format and contain time and
%       integrated intensity. 

cd('integrated-data');
load(filenameConv);
convTime = time;
convInt = Int;
load(filenameUnc);
unTime = time;
unInt = Int; 

% normalises the data and determines position of maximum 
[convImax,convNormInt,convP] = normshift(convTime,convInt);
[unImax,unNormInt,unP] = normshift(unTime,unInt);

% shifts the data so that their maxima lie at the same time
timeDiff = convImax - unImax;
if timeDiff >0
    unPshifted = cat(1,zeros(timeDiff,1), unP);
    convPshifted = cat(1,convP,zeros(timeDiff,1));  
else
    unPshifted = cat(1,unP,zeros(timeDiff,1));
    convPshifted = cat(1,zeros(timeDiff,1),convP);  
end
 newTime = 1:mean(diff(convTime)):length(unPshifted)*mean(diff(convTime));  

%% plot both functions 
cd('..');
plot(newTime, unPshifted, 'linewidth',2);
hold on;
plot(newTime, convPshifted, 'linewidth',2);
set(gca, 'Ylim',[-0.1 1.1]);
hold off;
graphicsSettings;
fontName = 'Times New Roman';
fontSize = 22;
set(gca,'DefaultTextInterpreter','latex');
legend('unconverted','converted');
ylabel('Integrated Intensity (a.u.)','FontSize',fontSize,'FontName',fontName);
xlabel('time (ps)','FontSize',fontSize,'FontName',fontName);
print([filenameConv '-IntegratedPlot-both.png'],'-dpng','-r300');
savefig([filenameConv '-IntegratedPlot-both.fig']);
clf();

%% plot difference
plot(newTime, convPshifted-unPshifted, 'linewidth',2);
%     set(gca, 'Ylim',[0 1.1]);
hold off;
set(gca,'DefaultTextInterpreter','latex');
graphicsSettings;
ylim([-1 1]);
legend('difference');
ylabel('Difference Intensity (a.u.)','FontSize',fontSize,'FontName',fontName);
xlabel('time (ps)','FontSize',fontSize,'FontName',fontName);
title('converted - unconverted');
print([filenameConv '-IntegratedPlot-diff.png'],'-dpng','-r300');
savefig([filenameConv '-IntegratedPlot-diff.fig']);
    
    
end
