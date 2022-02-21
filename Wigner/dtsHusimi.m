function [ HF ] = dtsHusimi(q,p,nTh,nCoh,varargin)
%THERMHUSIMI Returns Q(q,p) for a thermal state with NPHOTONS.
%   Detailed explanation goes here

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultP0 = sqrt(2*nCoh); % Change that for another normalization!
addParameter(parser,'P0',defaultP0,@isnumeric);
defaultQ0 = 0;
addParameter(parser,'Q0',defaultQ0,@isnumeric);
defaultPhaseAveraged = false;
addParameter(parser,'PhaseAveraged',defaultPhaseAveraged,@islogical);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,p0,phaseAveraged,q0] = c{:};

if ~phaseAveraged
    disc = mean(diff(q));
    HF = zeros(length(q),length(p));
    for iP = 1:length(p)
        HF(:,iP)=disc^2*norm^2*1/(pi*(nTh+1))* ...
            exp(-((q - q0).^2 + (p(iP) - p0).^2)/((2*norm)^2*(nTh+1)));
    end
    HF = HF./sum(sum(HF)); % Usually not necessary

else
    p = q';
    alphaSquare = (q.^2 + p.^2)/(2*norm)^2;
    alpha0Square = (q0.^2 + p0.^2)/(2*norm)^2;
    HF = exp(-(alphaSquare + alpha0Square)./(nTh+1))/(pi*(nTh+1)) .* besseli(0,2*sqrt(alphaSquare.*alpha0Square)/(nTh+1));
    HF = HF./sum(sum(HF));
end


end

