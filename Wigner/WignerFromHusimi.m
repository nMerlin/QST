function [WF,pq] = WignerFromHusimi(HS,varargin)
% computes the Wigner function (WF) from a given Husimi function HS. HS must
% be a quadratic matrix
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
nHS = length(HS);
shift = (nHS-nP)/2;
HSs = HS(shift+1:end-shift,shift+1:end-shift);

%% make p and q matrix
qMatrix = zeros(nP,nP);
for x = 1:nP
    qMatrix(x,:) = q;
end
pMatrix = qMatrix';

%% prepare Wigner matrix 
WF = zeros(size(HSs));
[maxX,maxY]=size(HSs);

%% evaluate each value of the Wigner matrix
for  y = 1:maxY 
    for  x = 1:maxX 
        gaussian = exp(1/norm^2*( 1/8 *(qMatrix.^2 + pMatrix.^2) + 1i*(q(y)* qMatrix + p(x) * pMatrix ))) ;
        Conv = HSs.*gaussian;
        WF(x,y) = sum(sum(Conv));        
    end;
end;

WF = WF/sum(sum(WF));

end