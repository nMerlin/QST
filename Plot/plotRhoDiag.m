function plotRho(rho, filename, varargin)

%% Validate and parse input arguments
parser = inputParser;
defaultPart = 'real'; 
addParameter(parser,'Part',defaultPart);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[part] = c{:};


maxN = length(rho) -1;
n = 0:maxN;


if strcmp(part,'real')
    pn = real(diag(rho));
elseif strcmp(part,'abs')
    pn = abs(diag(rho));
elseif strcmp(part,'imag')
    pn = imag(diag(rho));
end

meanN = n*pn;
varN = (n.^2)*pn - meanN^2;

bar(n,pn,'r');

xlim([-1 maxN+1]);

fontname = 'Times New Roman';
fontsize1 =22;
fontsize2 =20;
xlabel('$n$','FontSize',fontsize1,'Interpreter','latex');

if strcmp(part,'real')
    ylabel('$Re(\rho_{nn})$','FontSize',fontsize1,'Interpreter','latex');
elseif strcmp(part,'abs')
    ylabel('$|\rho_{nn}|$','FontSize',fontsize1,'Interpreter','latex');
elseif strcmp(part,'imag')
    ylabel('$Im(\rho_{nn})$','FontSize',fontsize1,'Interpreter','latex');
end

set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', 'k','FontName',fontname);

set(gca,'DefaultTextInterpreter','latex');
text(10,0.1,['$\bar{n} =$ ' num2str(meanN,'%.2f') char(10) 'Var$(n) =$ ' ...
    num2str(varN,'%.2f') ],'FontSize',fontsize2);

% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi]);
% set(gca,'XTickLabel',{'$0$','$\frac{1}{2} \pi$','$\pi$','$\frac{3}{2} \pi$','$2 \pi$'});

%grid on; 

print(['PlotRho-' filename '-' num2str(part) '.png'],'-dpng');

end