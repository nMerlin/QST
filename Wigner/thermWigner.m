function WF = thermWigner( q, p, nPhotons, varargin )
%THERMWIGNER Returns W(q,p) for a thermal state with NPHOTONS photons
%
% Formulation is in alpha-space but choice of normalization has to be taken
% into account.

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultTheta = 0; % Rotation angle in degrees
addParameter(parser,'Theta',defaultTheta,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,theta] = c{:};

disc = mean(diff(q));
WF = zeros(length(q),length(p));
for iP = 1:length(p)
    qrot = q*cosd(theta)-p(iP)*sind(theta);
    prot = q*sind(theta)+p(iP)*cosd(theta);
    WF(:,iP) = disc^2*norm^2*2/pi*1/(2*nPhotons+1)* ...
        exp(-2*((qrot*norm).^2+(prot*norm).^2)/(2*nPhotons+1));
end
%WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end

