function [Max, integratedInt] = plotFullChipSpectrumAndFit(filenameSIG,filenameBG,varargin)
%% PLOTSPECTRUM plots an optical spectrum with full chip ROI (data has form x, y, Intensity), subtracts background, 
% fits a Gaussian to the peak. 
%
%   Input Arguments:
%       filenameSIG: file with the data. The file should be located in folder
%       'raw-data'. The data should consist of two columns; 
%       first column contains wavelength and second column contains
%       intensities
%       filenameBG: file with background data.
%       'XLim': limits for x-axis around the peak.
%       'XUnit': x unit can be nm (wavelength) or Hz (frequency).

%% Validate and parse input arguments
parser = inputParser;
defaultFit = 'yes'; %
addParameter(parser,'Fit',defaultFit);
defaultSave = 'yes'; %
addParameter(parser,'Save',defaultSave);
defaultInterpolate = 'yes';
addParameter(parser,'Interpolate',defaultInterpolate);
defaultSubtract = 'yes'; %
addParameter(parser,'Subtract',defaultSubtract);
defaultXLim = []; %
addParameter(parser,'XLim',defaultXLim,@isnumeric);
defaultXUnit = 'nm';
addParameter(parser,'XUnit',defaultXUnit);
defaultXrange = [];
addParameter(parser,'Xrange',defaultXrange);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[fitoption,intp,save,subtract,xLim,xRange,xUnit] = c{:};

%%
fontsize = 24;
lightVelocity = 299792458;

%% load data
    cd('raw-data')
    %data = textread(filenameSIG);
    data = textread(filenameSIG,'','delimiter',',','headerlines',1); 
    
    X = data(:,1); % wavelength
    Y = data(:,2); % pixel position
    Int = data(:,3); %intensity
    
    %% subtract background
    if strcmp(subtract, 'yes')
        dataBG = textread(filenameBG);
        IntBG = dataBG(:,3);
        Int = Int - IntBG;
    else
        %bg = min(M);
        %M = M-bg;
    end
    
    %% set negative values to zero
    Int(Int<0) = 0;
    
%% sort data
    if Y(1) == 0
        Y=Y+1;
    end
     n = length(Y)/max(Y);
     X = X(1:n);
     Y = 1:max(Y);
     Y = Y';
     Int = reshape(Int, [n max(Y)]);
 
% %% confine data in pixel range if wished
% if not(isempty(xRange))
%     Int = Int((w>= min(xRange)) & (w<= max(xRange))); 
%     w = w((w>= min(xRange)) & (w<= max(xRange)));
% end
         
%% get maximum and peak position
    [~,IXs] = max(Int); %indices of the rows with the maxima of each column
    [Max,IY] = max(max(Int)); %IY: column of the absolute maximum
    IX = IXs(IY); %IX: row of the absolute maximum
    
%% confine data in pixel range if wished 
if not(isempty(xLim))
%     Int = Int(IX-xLim:IX+xLim, IY-xLim:IY+xLim); 
%     X = X(IX-xLim:IX+xLim);
%     Y = Y(IY-xLim:IY+xLim);
    midX = round(length(X)/2);
    midY = round(length(Y)/2);
    Int = Int(midX-xLim:midX+xLim, midY-xLim:midY+xLim); 
    X = X(midX-xLim:midX+xLim);
    Y = Y(midY-xLim:midY+xLim);
end
 
%% integrated data
    integratedInt = sum(sum(Int));

         
 %% make 2D surface plot

surf(X, Y, Int');
colorbar;
view(0,90);
shading flat;
axis tight;

fontsize = 24;
graphicsSettings;
fontName = 'Times New Roman';
set(gca,'XGrid','on','YGrid','on');
xlabel('X (pixel)','FontSize',fontsize,'FontName',fontName);
ylabel('Y (pixel)','FontSize',fontsize,'FontName',fontName);

set(gca,'DefaultTextInterpreter','latex');
text(max(X)/4, max(Y)/4,['max Int ' num2str(Max,'%.0f') ' counts' char(10) ...
    'integr. Int ' num2str(integratedInt,'%.0f') ' counts' char(10)],...
    'FontSize',fontsize-4);

cd('..');
print([filenameSIG '-2Dsurfplot.png'],'-dpng','-r300');
savefig([filenameSIG '-2Dsurfplot.fig']);
clf();

end
