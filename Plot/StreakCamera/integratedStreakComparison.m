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

%use not smoothed data
convP = convNormInt;
unP = unNormInt;

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
plot(newTime/1000, unPshifted, 'linewidth',2);
hold on;
plot(newTime/1000, convPshifted, 'linewidth',2);
plot(newTime/1000, convPshifted-unPshifted, 'linewidth',2);
set(gca, 'Ylim',[min(convPshifted-unPshifted) 1.1]);
hold off;
legend('Unconverted','Converted','Difference');
ylabel('Integrated intensity (a.u.)');
xlabel('Time (ns)');
graphicsSettings;
grid;
print([filenameConv '-IntegratedPlot-both.png'],'-dpng','-r300');
savefig([filenameConv '-IntegratedPlot-both.fig']);
clf();

%% plot difference
% plot(newTime/1000, convPshifted-unPshifted, 'linewidth',2);
% %     set(gca, 'Ylim',[0 1.1]);
% hold off;
% graphicsSettings;
% ylim([-1 1]);
% ylabel('Difference intensity (a.u.)');
% xlabel('Time (ns)');
% print([filenameConv '-IntegratedPlot-diff.png'],'-dpng','-r300');
% savefig([filenameConv '-IntegratedPlot-diff.fig']);
    
    
end
