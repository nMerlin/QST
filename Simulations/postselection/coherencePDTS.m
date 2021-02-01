function Coherence = coherencePDTS(nTherm,nCoherent)
% computes coherence of a phase-averaged displaced thermal state.
% nTherm: thermal photon number
% nCoherent: coherent photon number
%     a = 2*nTherm + 1;
    a = (nTherm + 1).^2 - nTherm.^2;
    Coherence = (1 - exp(-2*nCoherent./a) .* besseli(0,2*nCoherent./a ) )./a;
end

