%for LO:
folder C:\Users\lab\Documents\@data\2019-08-09-LO
filename = 'SpotsWithLens2.bmp';
A = imread(filename);
A = A(:,500:end); %cut the other spot off that was from fiber 
fcam = 0.3; %focal lenght of detector camera lens in m
%correct for telescopes in the way:
%Signal is after telescopes on the camera in k space only 375 px and
%without telescopes 500 px. So telescopes reduce k space of signal on
%camera by 0.75. So vice versa, if we project LO to signal k space before
%MO, this should be 1/0.75 bigger. 
telFactor = 1/0.75;

%for zwischenabbildung signal pinhole:
folder C:\Users\lab\Documents\@data\2020-10-08-M3396-92\camera
filename = '258mW-signalspot-zwischenabbildung-unpolarized-pinhole.bmp';
A = imread(filename);
fcam = 0.15; %focal lenght of zwischenabbildung camera lens in m
telFactor = 75/50;


%% get real space width
A = double(A);
pixToMum = 4.5; %pixels to micrometer of the camera
[Xlength, Ylength] = size(A);
x = pixToMum * (1:Xlength);
[~,I] = max(A);
%for LO
Acut = A(:,316); %better... 
%for phinhole
 Acut = A(:,1033);
FWHMx = fwhm(x,Acut); % x width in micrometer 
FWHMxOnSample =  FWHMx*fMO/fcam/telFactor;

%% make fourier transform 
Y = fft2(A);
% shift center
Y = abs(fftshift(Y));

%k space vectors in inverse micrometers 
kxCut = 2*pi*(1:Xlength)/pixToMum /Xlength;
kyCut = 2*pi*(1:Ylength)/pixToMum /Ylength;
[kx,ky] = meshgrid(kxCut,kyCut);
kx = kx'; ky = ky';
%surf(kx,ky,Y); shading flat;

% make cut 
Ycut = Y(:,round(Ylength/2));
%have zero in center
kxCut = kxCut-max(kxCut)/2;
plot(kxCut,Ycut);
% for LO
Ycut(Ycut>1.659e5) = 1.659e5;%remove peak at kx = 0 %make this better

plot(kxCut,Ycut);
xlabel('k_{x} (\mum ^{-1})');
savefig('kSpace.fig');
clf();
% get k space width
FWHMk = fwhm(kxCut,Ycut);

%We want to get the k space as big as it was in the sample emission, that
%is with the microscope objective, because the dispersions k space is
%calibrated to the NA of the MO. 
%correct for focal length of camera / focal length of MO
fMO = 0.01; %focal lenght of MO in m
kxCutCor = kxCut*fcam/fMO;
plot(kxCutCor,Ycut);
xlabel('k_{x} (\mum ^{-1})');
savefig('kSpaceCorrectedForMO.fig');
clf();
FWHMkCor = fwhm(kxCutCor,Ycut);

FWHMkCorTel = FWHMkCor /telFactor;

%% test other formula; I think, somewehre factor are missing to convert between FWHM and sigma 
%get D = FWHM in collimated space bebvfore the camera lens
lambda = 770e-9; % light wavelength in m 
B = FWHMx*1e-6 /sqrt(2*log(2)); %picture width in m; devide by /sqrt(2*log(2)) to go from FWHM to radius
D = 4*lambda *fcam/pi/B; % in m
% from https://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=2940
 % convert this to the angle we would get when ficusing with MO
 theta = atan(D/2/fMO);
% get half k width frim this in micrometer 
k = 2*pi/lambda *sin(theta) /1e6;
FWHMkTest = 2*k;  %What about factor /sqrt(2*log(2))??

% Alternative: 
theta = 2*lambda/pi/B;
% from https://en.wikipedia.org/wiki/Gaussian_beam
k = 2*pi/lambda *sin(theta) /1e6;
FWHMkTest = 2*k *fcam/ fMO;%What about factor /sqrt(2*log(2))??

%ave results
save('results-fft.mat','A','Y','Ycut','kxCut','FWHMx','FWHMk','FWHMkCor',...
    'FWHMkCorTel','FWHMkTest');



