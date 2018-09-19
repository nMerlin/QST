function [g1,t,t_coh] = g1FromPSD(lambda,psd)
%G1FROMPSD computes the unnormalized first order coherence function from
%the power spectral density of light via the Wiener-Khintchine theorem. For
%mathematical details see, for example, book "Optical coherence and quantum
%optics" by Leonard Mandel and Emil Wolf.
%
% Inputs:
%   psd: unnormalized power spectral density of investigated light
%   lambda: wavelength axis of psd in meters
%   nRep: multiples of the 
%
% Outputs:
%   g1: unnormalized first order coherence function
%   t: time axis for g1
%   t_coh: coherence time as normalized mean root-mean square width

%% Constants and parameters
c = 299792458; % speed of light
nMult = 50; % adapts resolution of g1 in time, should be even

%% Compute inverse fourier transform
% Calculate angular frequencies
om1 = 2*pi*c/max(lambda);
om2 = 2*pi*c/min(lambda);
dOm = om2 - om1;

% add zeros to the end of the psd for increased resolution
n = nMult*length(lambda);
om = linspace(om1,om1+nMult*dOm,n);
if ~isrow(psd)
    psd = psd';
end
psd = [psd,zeros(1,(nMult-1)/nMult*n)];

% fourier transform and calculate time axis
g1 = abs(fftshift(ifft(psd)));
dt = 1/(mean(diff(om))*n)*2*pi;
t = linspace(-n/2*dt,n/2*dt,n);

%% Compute coherence time t_coh as normalized root-mean square width of g1
t_coh = sqrt(sum(t.^2.*g1.^2)/sum(g1.^2));

end

