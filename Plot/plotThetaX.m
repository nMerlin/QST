function [varX,delX] = plotThetaX(theta,X, ys, filename, varargin)

%theta, X, ys should contain only one segment

%% Validate and parse input arguments
%method for computing the variance and stuff
p = inputParser;
defaultMethod = 'Fit';
addParameter(p,'Method',defaultMethod,@isstr);
defaultPlotDeviation = 'No';
addParameter(p,'PlotDeviation',defaultPlotDeviation,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[method,plotDeviation] = c{:};

% remove NaN values
X = X(~isnan(X));
theta = theta(~isnan(theta));

% calculate stuff

if strcmp(method,'Intervals')

    nIntervals = 200;
    [ XD, thetaD ] = discretizeTheta( X, theta,nIntervals );
    [Xline, varX] = deal(zeros(nIntervals, 1));
    ThetaLine = zeros(nIntervals, 1);
    for iInterval = 1 : nIntervals
            Xline(iInterval) = mean(XD(~isnan(XD(:,iInterval)), iInterval));
            ThetaLine(iInterval) = mean(thetaD(~isnan(thetaD(:, iInterval)),iInterval));
            varX(iInterval) = var(XD(~isnan(XD(:,iInterval)),iInterval));
    end 

    varX = mean(varX);
    delX = mean(sqrt(varX));
    
elseif strcmp(method,'Fit')
    
    [ fitParams, fitFunction,~,~] = fitSinus( theta,X);
    Xline=fitFunction(fitParams,theta);
    [ThetaLine,I] = sort(theta);
    deviation = abs(X - Xline);
    varX = sum(deviation.^2)/(length(deviation)-1) ;
    delX = sqrt(varX);
    Xline = Xline(I);  

elseif strcmp(method,'Spline')
    
    Xline=ys;
    [ThetaLine,I] = sort(theta);
    deviation = abs(X - Xline);
    varX = sum(deviation.^2)/(length(deviation)-1) ;
    delX = sqrt(varX);
    Xline = Xline(I);  
end

[maxX,maxI]=max(Xline);

%plot
  
plot(theta,X,'r.');
xlim([0 max(theta)]);
fontname = 'Times New Roman';
fontsize1 =22;
fontsize2 =20;
xlabel('$\theta$','FontSize',fontsize1,'Interpreter','latex');
ylabel('$X_{\theta}$','FontSize',fontsize1,'Interpreter','latex');
set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', 'k','FontName',fontname);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi]);
set(gca,'XTickLabel',{'$0$','$\frac{1}{2} \pi$','$\pi$','$\frac{3}{2} \pi$','$2 \pi$'});
grid on;
hold on;

plot(ThetaLine,Xline,'k-','LineWidth',1.5);

if strcmp(plotDeviation,'yes')
    h = plot(theta,deviation.^2*5,'Color',[132 184 24]/255);
    uistack(h,'bottom');
    l=legend('$(X_{i} - X_{fit,i})^{2}$','$X$ gemessen','Fit');
    l.Interpreter = 'latex';
else
    
%plot Annotations
    h = annotation('doublearrow');
    set(h,'parent',gca,'Position',[ThetaLine(maxI) maxX-delX 0 2*delX],'LineWidth',1.5);
    set(gca,'DefaultTextInterpreter','latex');
    text(ThetaLine(maxI)+0.2,maxX + 0.4,'$\pm \Delta X$','FontSize',fontsize2);
end

print(['plotThetaX-' filename '-method-' method],'-dpng');

end