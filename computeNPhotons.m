function [nPhotons] = computeNPhotons(powerLO, lambda, frequencyLO)
%COMPUTENPHOTONS Computes the number of photons per pulse of a laser beam
%with power POWERLO, wavelength LAMBDA and repetition rate FREQUENCYLO

SPEED_OF_LIGHT = 299792459; % m / s
PLANCK = 6.626070040e-34; % J * s

nPhotons = powerLO / ...
    (PLANCK * SPEED_OF_LIGHT / lambda * frequencyLO); % per pulse

end

