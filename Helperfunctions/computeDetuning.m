function [detuningMeV, detuning2g0,Ecav,lambdaCav] = computeDetuning(ELP,R,EX)
% computes the detuning of a polariton cavity. Input parameters are: 
% ELP : lower polariton energy at k = 0 (in eV)
% R: Rabisplitting 2g0 of the sample (in meV)
% EX : Exciton energy (at k = 0) (in eV)
% Example: [detuningMeV, detuning2g0,Ecav,lambdaCav] = computeDetuning(1.6093,9.5,1.6195)

%% Sample char. for M3396
% R = 9.5; %Rabi splitting in meV
% EX = 1.6195; % Exciton energy in eV
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;

u = (EX - ELP)/(R/1000);

% detuning in units of 2g0
detuning2g0 = 1./(4*u) - u;

% detuning in meV
detuningMeV = detuning2g0 * R;

% Energy (eV) and wavelength (nm) of the cavity
Ecav = EX + detuningMeV/1000;
lambdaCav = h *c0/(Ecav*e0)/10^-9;

end