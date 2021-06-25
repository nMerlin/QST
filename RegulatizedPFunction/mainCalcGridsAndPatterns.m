function [] = mainCalcGridsAndPatterns(varargin)
%calculates all coordinate grids and pattern functions for given
%parameters and saves them.
%% This is for a normalization such that the fluctuation of a vacuum field Delta X = 1. 

%% Validate and parse input arguments
p = inputParser;
defaultDirectory = 'C:\Users\Public\Documents\archived-data';
% Directory in which to save the grids and patterns
%directory = 'C:\Users\lab\Documents\@archived-data';
addParameter(p,'Directory',defaultDirectory,@isstr);
defaultMaxQuad = 20; %max abs value of the new quadrature coordinates 
addParameter(p,'MaxQuad',defaultMaxQuad,@isnumeric);
defaultMaxX= 20; %max abs value of the old quadrature coordinates 
addParameter(p,'MaxX',defaultMaxX,@isnumeric);
defaultPhiStep= 0.1; %step size of the phase grid
addParameter(p,'PhiStep',defaultPhiStep,@isnumeric);
defaultR= 0.7; %filter Parameter R
addParameter(p,'R',defaultR,@isnumeric);
defaultXStep= 1; %step size of the old quadrature coordinates grid
addParameter(p,'XStep',defaultXStep,@isnumeric);
defaultResolution= 1; %step size of the new quadrature coordinates grid
addParameter(p,'Resolution',defaultResolution,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[directory,maxQuad,maxX,phiStep,R,Resolution,xStep] = c{:};

directory = [directory,'\PFunction-R-',num2str(R),'-maxQuad-',num2str(maxQuad),'-Resolution-',num2str(Resolution)];
if ~exist(directory,'dir')
    mkdir(directory)
end

tic;
[xGrid,phiGrid] = makeGridxAndPhi(maxX,xStep,phiStep);
makeXGridAndPattern(xGrid,phiGrid,R,maxQuad,Resolution,directory);
toc;

end