function plot3ChMovie(O1,O2,O3,theta,varargin)
%PLOT3CHMOVIE Creates a MP4-file while scanning different regions in a
%3Ch-Dataset
%
% Arguments:
%   (O1,O2,O3,theta) - 3-Channel dataset where O1 and O2 are orthogonal
% Optional Arguments:
%   plot3ChMovie(~,~,~,~,filename) - 'filename' is the name of the output
%       file, the default is '3ChGif.gif'
%   plot3ChMovie(~,~,~,~,~,scanMode) - scanMode selects the mode of
%       scanning. The default 'square' scans a small square in lines along
%       the O1 and O2 directions.
%   plot3ChMovie(~,~,~,~,~,~,nBins) - number of bins for the 1D histogram
%       (default: 100)

%% Handle optional input arguments and default values
nVarargin = length(varargin);
optArgs = {'3ChMovie.mp4' 'square' 100 5};
optArgs(1:nVarargin) = varargin;
[filename,scanMode,nBins,delays] = optArgs{:};

switch scanMode
    case 'square'
        region = {'rectangle' 0.5 0.5};
        x = repmat(-5:0.5:4.5,1,length(-5:0.5:4.5)); x = x';
        y = repmat(4.5:-0.5:-5,length(4.5:-0.5:-5),1); y = y(:);
end

h = figure;
axis tight  % set axis limit to the range of the data
axis manual  % keep the current or manually chosen axis limits
[H, binsO1, binsO2] = histogram2D(O1,O2);

% Create movie object
movie = VideoWriter(filename, 'MPEG-4');
movie.Quality = 100;
open(movie);

for k = 1:length(x)
    clf(h);
    
    % Main plot
    XSel = selectRegion(O1,O2,O3,theta,region{1},x(k),y(k), ...
        region{2},region{3});
    edges = linspace(-10,10,nBins);
    histogram(XSel,edges,'Normalization','probability');
    set(gca,'YLim',[0 0.05],'XLim',[-10 10]);
    xlabel('X3');
    
    % Inset
    insetAx = axes('Parent',gcf,'Position',[0.2 0.6 0.25 0.25]);
    imagesc(binsO1,binsO2,H); axis on; colormap hot;
    set(insetAx,'FontSize',8,'XLim',[-5 5],'YLim',[-5 5],'XTickLabel','');
    title('Selected Region');
    xlabel('X1');
    ylabel('X2');
    hold on;
    switch scanMode
        case 'square'
            fill([x(k) x(k)+region{2} x(k)+region{2} x(k)], ...
                [y(k) y(k) y(k)+region{3} y(k)+region{3}],'b');
    end
    hold off;
    
    %% Write movie
    % Capture the plot as an image
    %frame = getframe(h); (no control over resolution)
    frame = print(h,'-r150','-RGBImage');
    for i = 1:delays
        writeVideo(movie, frame);
    end
    
%     % Write to GIF
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1
%         imwrite(imind,cm,filename,'gif','Loopcount',inf, ...
%             'DelayTime',0.1);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append', ...
%             'DelayTime',0.1);
%     end
end

close(movie);

end

