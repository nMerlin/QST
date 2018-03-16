function saveA5Landscape(filename)
%SAVEA5Landscape Save current figure as A5 landscape pdf

fig = gcf;
fig.Units = 'centimeters';
fig.Position = [1,1,21,14.8];
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
print(filename,'-dpdf');

end

