function [HS,pq] = HusimiFromWigner(WF,varargin)
% computes the Husimi Q function (HS) from a given Wigner function WF. WF must
% be a quadratic matrix, e.g. computed with mainWignerFromRho. 
% Formula from equation (6) from V. A. Andreev: A TRANSFORMATIONAL
% PROPERTY..., Theoretical and Mathematical Physics, 166(3): 356–368 (2011)
% https://link.springer.com/content/pdf/10.1007/s11232-011-0028-8.pdf
% Optional Input Arguments:
%   'PQ': Specify p axis. Default is -20:0.125:20, as this is used in
%   mainWigner. 

%% Validate and parse input arguments
parser = inputParser;
defaultPQ = -20:0.125:20;
addParameter(parser,'PQ',defaultPQ,@isvector);
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,pq] = c{:};

%% prepare p, q axis
[p,q] = deal(pq);
nP = length(p);
nWF = length(WF);
shift = (nWF-nP)/2;
WFs = WF(shift+1:end-shift,shift+1:end-shift);

%% make p and q matrix
qMatrix = zeros(nP,nP);
for x = 1:nP
    qMatrix(x,:) = q;
end
pMatrix = qMatrix';

%% prepare Husimi matrix 
HS = zeros(size(WFs));
[maxX,maxY]=size(WFs);

%% evaluate each value of the Husimi matrix
for x = 1:maxX   
    for y = 1:maxY
        gaussian = exp(-2/(2*norm)^2 * (( q(y) - qMatrix).^2 + (p(x) - pMatrix).^2) );
        Conv = WFs.*gaussian;
        HS(x,y) = sum(sum(Conv));        
    end;
end;

HS = HS * 2/(2*norm)^2/pi;
HS = HS/sum(sum(HS));

end