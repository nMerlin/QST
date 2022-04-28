function plotHusimi(O1,O2,iSelect,varargin)
%PLOTHUSIMI Plot histogram of Husimi-Q function
%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename);
defaultLimits = [-7,7];
addParameter(p,'Limits',defaultLimits);
defaultStyle = 'auto'; % 'simple'
addParameter(p,'Style',defaultStyle);
defaultTitleStr = 'H(q,p)';
addParameter(p,'Title',defaultTitleStr,@isstr);
defaultnBinsA = 1000;
addParameter(p,'nBinsA',defaultnBinsA,@isnumeric);
defaultnBinsB = 1000;
addParameter(p,'nBinsB',defaultnBinsB,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,limits,nBinsA,nBinsB,style,titlestr] = c{:};

%% Plotting
[H, binsO1, binsO2] = histogram2D(O1,O2,'nBinsA',nBinsA,'nBinsB',nBinsB);
imagesc(binsO1,binsO2,H); axis on; colormap hot; hold on;
scatter(O1(iSelect),O2(iSelect),'.g'); hold off;
set(gca,'XLim',limits,'YLim',limits);
if strcmp(style,'simple')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    title(titlestr,'FontSize',48);
else
    xlabel('q');
    ylabel('p');
    title(titlestr);
end
pbaspect([1 1 1]); % ensure that a circle doesn't look elliptic

saveFig(filename);

end

