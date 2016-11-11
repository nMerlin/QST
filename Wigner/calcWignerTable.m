function calcWignerTable(nVal,mVal,directory,minq,maxq,qintstep,overwrite)
%CALCWIGNERTABLE calculates FT(<q+1/2*q'|nVal>*<mVal|q-1/2*q'>). The
%resulting matrix will be saved in the given DIRECTORY. The function will
%be discretized in phase space from p,q = MINQ to MAXQ in steps of
%QINTSTEP. FT is the fourier transform.

assert((exist(directory,'dir')>0),'The given path is not a directory');
filename = strcat(directory,'\n',int2str(nVal),'m',int2str(mVal));
isFile = (exist(filename,'file')==2);

if (~isFile || overwrite)
    WigTab=zeros(size(minq:qintstep:maxq,2),size(minq:qintstep:maxq,2));
    nWigTab = size(WigTab,1);

    %disp(size(WigTab);
    for q=minq:qintstep:maxq
        for p=minq:qintstep:maxq
            WigTab(int16((q-minq)/qintstep)+1,int16((p-minq)/qintstep)+1) = ...
                wignernm(nVal,mVal,q,p,qintstep,2*minq,2*maxq);
        end;
        if parallel==0
            dispstat(strcat(num2str(100*((q-minq)/qintstep+1)/nWigTab),'%'));
        end
        %disp(q);
    end;

    save(filename, 'WigTab');
end

end