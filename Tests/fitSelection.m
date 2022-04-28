phaseShiftF = fitCoeff.phaseShift;
amplitudeF = fitCoeff.amplitude;
positionF = fitCoeff.position;
for i= length(fitCoeff.position): -1 : 1
    if ~selection(i)
        phaseShiftF(i) = [];
        amplitudeF(i) = [];
        positionF(i) = [];
    end
end


  
nDs_norm=nDs./nDs_amplSel;   
T=table(phaseselection', nDs', nDsStd', nDs_norm', 'Variablenames',{'phaseselection', 'nDs', 'nDsStd', 'nDs_norm'});
writetable(T,[num2str(positions(j)) 'mm-photonNumberVsPhase.txt']);

A = [phaseselection';nDs' ;nDsStd'; nDs_norm'];  
fileID = fopen([num2str(positions(1)) 'mm-photonNumberVsPhase.txt'],'w');
fprintf(fileID,'%5d %5d %5d %5d\n',phaseselection',nDs',nDsStd', nDs_amplSel);
%f  % 5.2f
fclose(fileID);


