function WF = thermWigner( q, p, nPhotons, varargin )
%THERMWIGNER Returns W(q,p) for a thermal state with NPHOTONS photons
%
%   Important: Normalization A = 1/sqrt(2) for quadrature operator q =
%   A*(a^dagger + a) is used.

%% Validate and parse input arguments
parser = inputParser;
defaultTheta = 0; % Rotation angle in degrees
addParameter(parser,'Theta',defaultTheta,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[theta] = c{:};

WF = zeros(length(q),length(p));
for iP = 1:length(p)
    qrot = q*cosd(theta)-p(iP)*sind(theta);
    prot = q*sind(theta)+p(iP)*cosd(theta);
    WF(:,iP) = 1/pi*1/(2*nPhotons+1)* ...
        exp(-(qrot.^2+prot.^2)/(2*nPhotons+1));
end
WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end

