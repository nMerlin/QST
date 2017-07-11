function XSel = threeChannelPhase(X1, X2, X3, timestamps)
%THREECHANNELPHASE Summary of this function goes here
%   Detailed explanation goes here

%% Constants
ORTH_WIDTH = 0.05;
PEAK_OPTS.MinPeakDistance = 200000;

%% Remove offsets, cut data into piezo segments and remove NaNs
X1 = X1 - mean(mean(X1));
X2 = X2 - mean(mean(X2));
X3 = X3 - mean(mean(X3));
X1 = piezoSegments(timestamps,X1);
X2 = piezoSegments(timestamps,X2);
X3 = piezoSegments(timestamps,X3);
X1 = X1(:,1:end-4,:);
X2 = X2(:,1:end-4,:);
X3 = X3(:,1:end-4,:);

%% Choose product for orthogonality post-selection
XProd = X1.*X2;

%% Approximate each segment with a cubic smoothing spline
[nPulses,nPieces,nSegments] = size(XProd);
x = 1:(nPulses*nPieces);
y = reshape(XProd,[nPulses*nPieces nSegments]);
ys = transpose(csaps(x,y',0.0000000001,x));
ys = ys(:);

%% Identify the intervals where the two channels are orthogonal
[~,maxlocs]=findpeaks(ys,PEAK_OPTS);
[~,minlocs]=findpeaks(-ys,PEAK_OPTS);
maxlocs = maxlocs(2:end-1);
minlocs = minlocs(2:end-1);
meanMax = mean(ys(maxlocs));
meanMin = mean(ys(minlocs));
dev = (meanMax-meanMin)*(1-ORTH_WIDTH)/2;
lowerBnd = meanMin + dev;
upperBnd = meanMax - dev;
iOrth = find(ys<upperBnd & ys>lowerBnd);
X1 = X1(iOrth); X2 = X2(iOrth); X3 = X3(iOrth);
% x = 1:length(ys);
% plot(x,ys,x,repmat(upperBnd,1,length(x)),x,repmat(lowerBnd,1,length(x)),x,repmat(mean(ys),1,length(x)))
% legend('signal','up','low','mean');
% waitforbuttonpress

%% 2D histogram of X1 and X2 with selected region
uniq1 = unique(X1(:));
hDisc = min(diff(uniq1)); % discretization
binsX1 = min(uniq1):hDisc:max(uniq1); % bin centers
uniq2 = unique(X2(:));
hDisc = min(diff(uniq2));
binsX2 = min(uniq2):hDisc:max(uniq2); % bin centers
numBins1 = numel(binsX1);
numBins2 = numel(binsX2);
% map X1 and X2 to bin indices
X1i = round(interp1(binsX1,1:numBins1,X1,'linear','extrap'));
X2i = round(interp1(binsX2,1:numBins2,X2,'linear','extrap'));
% limit indices to [1,numBins]
X1i = max(min(X1i,numBins1),1);
X2i = max(min(X2i,numBins2),1);
% Count number of elements in each bin
H = accumarray([X1i(:) X2i(:)], 1, [numBins1 numBins2]);
% Plot
% imagesc(binsX1,binsX2,H), axis on
% colormap hot; colorbar;

%% Phase & Amplitude Selection
h = figure;
axis tight manual
filename = 'postSelections.gif';
% for X2min = -5:0.5:4.5
% for X1min = -5:0.1:4.5
    clf(h);
    %bndX1 = [X1min,X1min+0.5]; bndX2 = [X2min,X2min+0.5];
    %iSelect = find(X1<bndX1(2) & X1>bndX1(1) & X2<bndX2(2) & X2>bndX2(1));
    iSelect = find(sqrt((X1.^2+X2.^2))<6 & sqrt((X1.^2+X2.^2))>5);
    XSel = X3(iSelect);
    uniq = unique(XSel(:));
    maxValue = max(-min(uniq),max(uniq));
    hDisc = min(diff(uniq)); % discretization
    histEdges = (-maxValue-hDisc/2):hDisc:(maxValue+hDisc/2);
    histogram(XSel,histEdges,'Normalization','probability');
    set(gca,'YLim',[0 0.02],'XLim',[-10 10]);
    xlabel('X3');
    % Inset
    insetAx = axes('Parent',gcf,'Position',[0.2 0.6 0.25 0.25]);
    imagesc(binsX1,binsX2,H), axis on
    colormap hot;
    set(insetAx,'FontSize',8,'XLim',[-5 5],'YLim',[-5 5],'XTickLabel','');
    title('Selected Region');
    xlabel('X1');
    ylabel('X2');
%     hold on;
%     fill([bndX1(1) bndX1(1) bndX1(2) bndX1(2)], ...
%         [bndX2(1) bndX2(2) bndX2(2) bndX2(1)],'b');
%     hold off;
    
%     % Capture the plot as an image
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     
%     % Write to GIF
%     if X1min == -5 && X2min == -5
%         imwrite(imind,cm,filename,'gif','Loopcount',inf, ...
%             'DelayTime',0.1);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append', ...
%             'DelayTime',0.1);
%     end
% end
% end

end
