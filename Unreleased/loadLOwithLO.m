filenameLO = '07-5mW-LOonly.raw';
filenameSIG = '08-5mW-LOwithLO.raw';
[X1,piezoSign1] = prepare1ChData(filenameLO,filenameSIG,'Channel',1);
[X2,piezoSign2] = prepare1ChData(filenameLO,filenameSIG,'Channel',2);
[X3,piezoSign3] = prepare1ChData(filenameLO,filenameSIG,'Channel',3);
[theta1,ys1] = computePhase(X1,1,piezoSign1);
[theta2,ys2] = computePhase(X2,1,piezoSign2);
[theta3,ys3] = computePhase(X3,1,piezoSign3);