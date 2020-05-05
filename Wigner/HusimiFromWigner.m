function [HS,pq] = HusimiFromWigner(WF,varargin)
% computes the Husimi Q function (HS) from a given Wigner function WF. WF must
% be a matrix, e.g. computed with mainWignerFromRho. 
% Optional Input Arguments:
%   'PQ': Specify p axis. Default is -20:0.125:20, as this is used in
%   mainWigner. 

%% Validate and parse input arguments
p = inputParser;
defaultPQ = -20:0.125:20;
addParameter(p,'PQ',defaultPQ,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[pq] = c{:};

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
        gaussian = exp(-( q(y) - qMatrix).^2 - (p(x) - pMatrix).^2 );
        Conv = WFs.*gaussian;
        HS(x,y) = sum(sum(Conv));        
    end;
end;

HS = HS/pi;
HS = HS/sum(sum(HS));

end