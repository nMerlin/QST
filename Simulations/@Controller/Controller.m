classdef Controller < handle
    %CONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hnext;
        errold;
        reject;
    end
    
    methods
        function obj = Controller()
            obj.reject = false;
            obj.errold = 1.0e-4;
        end
        function [status, h] = success(obj, err, h)
            beta = 0.4/8.0;
            alpha = 1.0/8.0-beta*0.75;
            safe = 0.9;
            minscale = 0.333;
            maxscale = 6.0;
            if (err <= 1.0)
                if (err == 0.0)
                    scale = maxscale;
                else
                    scale = safe*err^(-alpha)*obj.errold^beta;
                    if (scale<minscale); scale=minscale; end;
                    if (scale>maxscale); scale=maxscale; end;
                end
                if (obj.reject)
                    obj.hnext = h*min(scale,1.0);
                else
                    obj.hnext = h*scale;
                end
                obj.errold = max(err,1.0e-4);
                obj.reject=false;
                status = true;
            else
                scale = max(safe*err^(-alpha),minscale);
                h = h * scale;
                obj.reject = true;
                status = false;
            end
        end
    end
    
end

