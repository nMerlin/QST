function [detuningMeV, detuning2g0] = computeDetuning(ELP,R,EX)
% computes the detuning of a polariton cavity. Input parameters are: 
% ELP : lower polariton energy at k = 0 (in eV)
% R: Rabisplitting 2g0 of the sample (in meV)
% EX : Exciton energy (at k = 0) (in eV)

u = (EX - ELP)/(R/1000);

% detuning in units of 2g0
detuning2g0 = 1./(4*u) - u;

% detuning in meV
detuningMeV = detuning2g0 * R;

end