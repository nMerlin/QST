function WF= WignerFromDensityMatrix(rho)

maxn=size(rho,1);
load('WignerTable\n0m0');
WF=zeros(size(WigTab,1),size(WigTab,2));
for nVal=0:maxn-1
    for mVal=0:maxn-1
        loadstring=strcat('WignerTable\n',int2str(nVal),'m',int2str(mVal));
        load (loadstring);
        WF=WF+(rho(nVal+1,mVal+1)*WigTab);
    end;
end