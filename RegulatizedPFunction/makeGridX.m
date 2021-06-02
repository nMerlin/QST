function [] = makeGridX(xGrid,phiGrid,R,MaxQuad,Resolution,directory)
%for each position in xGrid, phiGrid, this makes a matrix in q,p
%coordinates and saves it. 

%% Normierung???

    QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);
    [QAxis,PAxis]=meshgrid(QuadVals,QuadVals);
    PhaseMatrix = atan2(PAxis,QAxis); %returns values in the closed interval [-pi,pi???
    AmplMatrix = sqrt(QAxis.^2 + PAxis.^2);
    
    for x = xGrid
        for phi = phiGrid
            % compute the X parameter, see Jans notes, below equ. 14
            XGrid = 2*R*(x - 2*AmplMatrix.*cos(phi + PhaseMatrix));
            filename = strcat(directory,'\x-',num2str(x),'-phi-',num2str(phi),'-XParameter.mat');
            save(filename, 'XGrid');
        end
    end
   


end