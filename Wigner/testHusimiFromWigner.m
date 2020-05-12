function [cohWF,thermWF,fockWF,cohHF,thermHF,fockHF,cohHS,thermHS,fockHS,cohDiff,thermDiff,fockDiff] = testHusimiFromWigner()

q = -20:0.125:20;
p = q;
nPhotons = 4; 
q0 = 5;
p0 = 3;

%% Theory functions
cohWF = cohWigner(q, p, nPhotons,'Q0',q0,'P0',p0);
thermWF = thermWigner( q, p, nPhotons); 
fockWF = fockWigner(q,p,nPhotons);

cohHF = cohHusimi(q, p, nPhotons,'Q0',q0,'P0',p0);
thermHF = thermHusimi( q, p, nPhotons ); 
fockHF = fockHusimi(q,p,nPhotons);

%% compute Husimi from Wigner functions
cohHS = HusimiFromWigner(cohWF,'PQ',q);
thermHS = HusimiFromWigner(thermWF,'PQ',q);
fockHS = HusimiFromWigner(fockWF,'PQ',q);

%% compare theory Husimi with Husimi from Wigner
cohDiff = cohHS - cohHF;
thermDiff = thermHS - thermHF;
fockDiff = fockHS - fockHF;

plotWigner(cohDiff,'Filename','coherentHusimiDifference.fig');
plotWigner(thermDiff,'Filename','thermalHusimiDifference.fig');
plotWigner(fockDiff,'Filename','fockHusimiDifference.fig');

plotWigner(cohHF,'Filename','coherentHusimiFromTheory.fig');
plotWigner(thermHF,'Filename','thermalHusimiFromTheory.fig');
plotWigner(fockHF,'Filename','fockHusimiFromTheory.fig');

plotWigner(cohHS,'Filename','coherentHusimiFromWigner.fig');
plotWigner(thermHS,'Filename','thermalHusimiFromWigner.fig');
plotWigner(fockHS,'Filename','fockHusimiFromWigner.fig');

end