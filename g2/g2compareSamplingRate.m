function [nRes, g2Av, g2Std, g2weighted, slowg2] = g2compareSamplingRate(X, filename)

maxRes = length(X(:));
k = 1:0.5:log(maxRes);
nRes =  exp(k);
nRes = floor(nRes);
fVector = 75.4e6 ./ nRes;

%nRes = 100:100:1000;
X = X(:);

[g2Av, g2Std,g2weighted,slowg2] = deal(zeros(length(nRes),1));

for i = 1:length(nRes)
%     nSegments = floor(length(X)/nRes(i));
%     Xres = X(1:nSegments*nRes(i));
%     Xres = reshape(Xres,[nRes(i) nSegments]); % Consecutive values in columns
    [g2vec,ada,~] = g2(X, nRes(i));
    g2Std(i) = sqrt(var(g2vec));
    g2Av(i) = mean(g2vec,'omitnan');
    g2weighted(i) =sum(ada.^2.*g2vec')/sum(ada.^2);
    slowg2(i)=mean(ada.^2)/mean(ada).^2;
end

fontname = 'Times New Roman';
fontsize1 =22;
fontsize2 =20;

%% plot av
loglog(nRes,g2Av,'o','markerSize',9,'lineWidth',1.5);
hold on;
loglog(nRes,g2weighted,'o','markerSize',9,'lineWidth',1.5);
loglog(nRes,slowg2,'o','markerSize',9,'lineWidth',1.5);
l = legend('simple mean','weighted with $n^{2}$','slow $g^{(2)}(0)$','Location','best');
l.Interpreter = 'latex';
l.FontSize = 15;
hold off;
set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', 'k','FontName',fontname);
ylabel('$< g^{(2)}(0)>$','FontSize',fontsize1,'Interpreter','latex');
xlabel('Anzahl Quadraturen $ m $ in Statistik','FontSize',fontsize1,'Interpreter','latex');

print([filename '.-g2Av.png'],'-dpng');
savefig([filename '-g2Av.fig']);
clf();

%% plot std
loglog(nRes,g2Std,'ro','markerSize',9,'lineWidth',1.5);
xlim([min(nRes) max(nRes)]);
set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
    'YColor', 'k','FontName',fontname);

ylabel('$\Delta g^{(2)}(0)$','FontSize',fontsize1,'Interpreter','latex');
xlabel('Anzahl Quadraturen $ m $ in Statistik','FontSize',fontsize1,'Interpreter','latex');

print([filename '-g2std.png'],'-dpng');
savefig([filename '-g2std.fig']);
clf();

%% plot frequency std
loglog(fVector,g2Std,'ro','markerSize',9,'lineWidth',1.5);
graphicsSettings;
ylabel('$\Delta g^{(2)}(0)$','FontSize',fontsize1,'Interpreter','latex');
xlabel('Frequency (Hz)','FontSize',fontsize1,'Interpreter','latex');
print([filename '-g2std-freq.png'],'-dpng');
savefig([filename '-g2std-freq.fig']);
clf();
%% plot frequency Av
loglog(fVector,g2Av,'o','markerSize',9,'lineWidth',1.5);
hold on;
loglog(fVector,g2weighted,'o','markerSize',9,'lineWidth',1.5);
loglog(fVector,slowg2,'o','markerSize',9,'lineWidth',1.5);
l = legend('simple mean','weighted with $n^{2}$','slow $g^{(2)}(0)$','Location','best');
l.Interpreter = 'latex';
l.FontSize = 15;
hold off;
graphicsSettings;
ylabel('$< g^{(2)}(0) >$','FontSize',fontsize1,'Interpreter','latex');
xlabel('Frequency (Hz)','FontSize',fontsize1,'Interpreter','latex');
print([filename '-g2Av-freq.png'],'-dpng');
savefig([filename '-g2Av-freq.fig']);
clf();
end