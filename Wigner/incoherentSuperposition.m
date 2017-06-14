function WF = incoherentSuperposition( nPhotons, alpha, q, varargin )
%INCOHERENTSUPERPOSITION WF = alpha*WF_thermal+(1-alpha)*WF_coherent
%
%   Paramters:
%       nPhotons - average number of photons
%       alpha - mixing paramter
%       q - WF will be evaluated in both directions at the given points
%
%   Optional Parameters:
%       gif - creates a gif-file with a 3D-animation for different alpha

% Optional input arguments
verbose = 0;
gif = 0;
quiet = 'notquiet';
if nargin > 3
    for i = 4:nargin
        eval([varargin{i-3} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end


WF_thermal = thermWigner(q,q,nPhotons);
WF_coherent = cohWigner(q,q,nPhotons);
WF = alpha*WF_thermal + (1-alpha)*WF_coherent;

if gif == 1
    h = figure;
    axis tight manual
    filename = 'testAnimated.gif';
    for alpha = 0:0.01:1
        WF = alpha*WF_thermal + (1-alpha)*WF_coherent;
        plotWigner(WF,'surface','narrow');
        xlim([-5 5]);
        ylim([-5 5]);
        zlim([0 max(max(WF_coherent))]);
        drawnow
        
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Write to GIF
        if alpha == 0
            imwrite(imind,cm,filename,'gif','Loopcount',inf, ...
                'DelayTime',0.1);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', ...
                'DelayTime',0.03);
        end
    end
end

end

