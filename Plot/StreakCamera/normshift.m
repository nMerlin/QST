%%normalize and shift
    function [Imax,normInt,normP] = normshift(timeFunc,IntFunc)
        %make spline interpolation
       % p = csaps(timeFunc, IntFunc,1e-4,timeFunc);
        p = IntFunc;
        %get max and startpoint 
        [maxP,Imax]=max(p);
        %startindex = min(find(p>0));
        % normalize
        normInt = IntFunc/maxP;
        normP = p/maxP;
        %shift time 
    %     shiftedTime = time(startindex-10:end);
    %     shiftedInt = time(startindex-10:end);
    end