function plotG2vsTime(times, g2vec, ada, filename)

%   FILENAME - this string will be included in the filename of 
%       the png-file
%
% Optional Parameters:
%   'limBins': Limit the number of bins in the histogram to 1000

%PLOTG2 Plot g2 and ada values against recording times (expects SI units)
%
% Input Arguments:
%   times: vector with time axis
%   g2vec: vector with g2 values
%   FILENAME - this string will be included in the filename of 
%       the png-file

%% Create plot of g2 and ada against time
[times,units] = convenientUnits(times,'s');

% plot g2
subplot(2,1,1)
plot(times,g2vec,'r');
axis([0 max(times) 0.5 2.5]);
%axis([0 max(times) 0.5 1.5]);
fontname = 'Times New Roman';
fontsize1 =22;
fontsize2 =20;
%xlabel(['$t$(',units,')'],'FontSize',fontsize1,'Interpreter','latex');
ylabel('$g^{(2)}(0,t)$','FontSize',fontsize1,'Interpreter','latex');
set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', 'k','FontName',fontname);
grid on; 

%%
% plot ada
subplot(2,1,2)
plot(times,ada,'k');
axis([0 max(times) 0 max(ada)+1]);
fontname = 'Times New Roman';
fontsize1 =22;
fontsize2 =20;
xlabel(['$t$ (',units,')'],'FontSize',fontsize1,'Interpreter','latex');
ylabel('$\bar{n} (t)$','FontSize',fontsize1,'Interpreter','latex');
set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', 'k','FontName',fontname);
grid on; 
%%

%Save figure 
print(strcat(filename,'G2-', '.png'), '-dpng');
savefig(strcat(filename,'G2-', '.fig'));
% bw = imread(strcat('G2-',filename,'-nResolution-',num2str(nResolution)...
%     ,'-n-',num2str(nAv),'.png'));
% bw2 = imcomplement(bw);
% imwrite(bw2,strcat('G2-',filename,'-nResolution-',num2str(nResolution)...
%     ,'-n-',num2str(nAv),'-inverted.png'));
%

clf;
%plot histogram of g2
histogram(g2vec,150,'Normalization','probability','FaceColor','r','EdgeColor','r');
xlabel('$g^{(2)}(\tau = 0)$','FontSize',fontsize1,'Interpreter','latex');
ylabel('Probability','FontSize',fontsize1,'Interpreter','latex');
xlim([0 3]);
set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', [0 0 0],'FontName',fontname);
grid on; 
print(strcat(filename,'G2-histogram','.png'), '-dpng');
savefig(strcat(filename,'G2-histogram','.fig'));
% bw = imread(strcat('G2-histogram',filename,'-nResolution-',num2str(nResolution),'.png'));
% bw2 = imcomplement(bw);
% imwrite(bw2,strcat('G2-histogram',filename,'-nResolution-',num2str(nResolution),'-inverted.png'));
clf;


end

