function WF = mainWignerFromRho(rho,varargin)
% calculates the Wigner function from the input density matrix rho.
% Optional Input: Directory: The directory, where the WignerTables WigTab
% are located that should be used.  

%% Validate and parse input arguments
p = inputParser;
defaultDirectory = 'C:\Users\lab\Documents\@archived-data\Wigner';
%directory = 'Z:\freefilesync-lab\matlab\QST\Wigner';
%directory = 'D:\@archive\2016-08-30-wigner-test';
%directory = 'C:\Users\lab\Documents\@archived-data\Wigner-Resolution-0.25';
addParameter(p,'Directory',defaultDirectory,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[directory] = c{:};

maxn=size(rho,1);
load(strcat(directory,'\n0m0'));
WF=zeros(size(WigTab,1),size(WigTab,2));
for nVal=0:maxn-1
    for mVal=0:maxn-1
        loadstring=strcat(directory,'\n',int2str(nVal),'m',int2str(mVal));
        load (loadstring);
        WF=WF+(rho(nVal+1,mVal+1)*WigTab);
    end;
end
%WF = real(WF);
WF = WF/sum(sum(WF));

end