function saveA4(filename)
%SAVEA4 Save current figure as A4 pdf

fig = gcf;
fig.Units = 'centimeters';
fig.Position = [1,1,21,29.7];
fig.PaperPositionMode = 'auto';
print(filename,'-dpdf');

end

