function [] = testCalcGrids()

maxX = 10;
xStep = 0.5;
phiStep = 0.1;
MaxQuad = 10;
Resolution = 0.5;
R = 0.5;
directory = ['C:\Users\Public\Documents\archived-data\PFunction-R-' num2str(R)];
%directory = 'C:\Users\lab\Documents\@archived-data\PFunction-XParameter';
if ~exist(directory,'dir')
    mkdir(directory)
end

tic;
[xGrid,phiGrid] = makeGridxAndPhi(maxX,xStep,phiStep);
makeGridX(xGrid,phiGrid,R,MaxQuad,Resolution,directory);
calcPatternFunctions(directory);
toc;

end