function WF = mainWignerFromRho(rho)

%directory = 'Z:\freefilesync-lab\matlab\QST\Wigner';
%directory = 'D:\@archive\2016-08-30-wigner-test';
directory = 'C:\Users\lab\Documents\@data\Wigner';

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