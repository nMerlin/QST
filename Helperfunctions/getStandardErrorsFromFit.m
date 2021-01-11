function [se]= getStandardErrorsFromFit(f,gof,method)
% returns the standard deviations on fitting coefficients.
% f, gof are results from fit 
% both methods apparently give the same results. 

    if strcmp(method,'method1')
        %source:
        %https://de.mathworks.com/matlabcentral/answers/34234-how-to-obtain-std-of-coefficients-from-curve-fitting,
        % comment from Tom Lane
%         The 1 comes from wanting 1 standard error. The negative sign is to 
%         get the level associated with 1 standard error below zero. 
%         The multiplication by 2 is to include the values beyond 2 standard 
%         error above the mean, by symmetry.
        level = 2*tcdf(-1,gof.dfe);
        m = confint(f,1-level);    
        se = (m(2,:)-m(1,:))/2;    
    else
        % source: https://de.mathworks.com/matlabcentral/answers/153547-how-can-i-compute-the-standard-error-for-coefficients-returned-from-curve-fitting-functions-in-the-c
        alpha = 0.95;
        ci = confint(f, alpha);
        t = tinv((1+alpha)/2, gof.dfe); 
        se = (ci(2,:)-ci(1,:)) ./ (2*t);
    end

end