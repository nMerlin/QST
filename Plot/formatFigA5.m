function formatFigA5(fig)
%FORMATFIGA5 Format the given figure to A5 format. If no handle is given,
%the current figure will be formatted.

if nargin == 0
    fig = gcf;
end
fig.Units = 'centimeters';
fig.Position = [1,1,21,14.8];
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';

end

