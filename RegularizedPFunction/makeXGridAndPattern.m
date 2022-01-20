function [] = makeXGridAndPattern(xGrid,phiGrid,R,MaxQuad,Resolution,directory)
%for each position in xGrid, phiGrid, this makes a matrix in q,p
%coordinates and saves it (XGrid). It also computes the corresponding pattern function and
%saves it. 

%% This is for a normalization such that the fluctuation of a vacuum field Delta X = 1. 

    QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);
    [QAxis,PAxis]=meshgrid(QuadVals,QuadVals);
    PhaseMatrix = atan2(PAxis,QAxis) - pi/2; %-pi/2 for the correct axis 
    AmplMatrix = sqrt(QAxis.^2 + PAxis.^2);
    alphaAmplMatrix = AmplMatrix/2; % relationship between X,P and alpha has factor 2
    
    for x = xGrid
        for phi = phiGrid
            % compute the X parameter, see Jans notes, below equ. 14
            XGrid = 2*R*(x - 2*alphaAmplMatrix.*cos(phi + PhaseMatrix)); 
            pattern = patternFunction(XGrid,R);
            filename = strcat(directory,'\x-',num2str(x),'-phi-',num2str(phi),'.mat');
            save(filename, 'XGrid','pattern');
        end
    end
   


end