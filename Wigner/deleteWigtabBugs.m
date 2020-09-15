directory = 'C:\Users\lab\Documents\@archived-data\Wigner-Resolution-0.25';
maxn=101;
load(strcat(directory,'\n0m0'));
for nVal=0:maxn-1
    for mVal=0:maxn-1
        loadstring=strcat(directory,'\n',int2str(nVal),'m',int2str(mVal),'.mat');
        if (nVal > 85 || mVal >85)
            delete(loadstring);
        end
      
    end;
end