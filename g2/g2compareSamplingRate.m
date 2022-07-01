function [nRes, g2Av, g2Std, g2weighted, slowg2] = g2compareSamplingRate(X, filename)

% maxRes = length(X(:));
% k = 1:0.5:log(maxRes);
% nRes =  exp(k);
% nRes = floor(nRes);
% fVector = 75.4e6 ./ nRes;

LOsamplerate = 75.39e6;

minF = LOsamplerate/1e7;
maxF = LOsamplerate/100;
numberOfF = 50;
step = (log(maxF)-log(minF))/numberOfF;
logs = log(minF):step:log(maxF);
fVector = exp(logs);
nRes = round(75.39e6 ./ fVector);

X = X(:);

[g2Av, g2Std,g2weighted,slowg2] = deal(zeros(length(nRes),1));

for i = 1:length(nRes)
    [g2vec,ada,~] = g2(X, nRes(i));
    g2Std(i) = sqrt(var(g2vec));
    g2Av(i) = mean(g2vec,'omitnan');
    g2weighted(i) =sum(ada.^2.*g2vec')/sum(ada.^2);
    slowg2(i)=mean(ada.^2)/mean(ada).^2;
end


%% plot av
semilogx(nRes,g2Av,'o');
hold on;
semilogx(nRes,g2weighted,'o');
semilogx(nRes,slowg2,'o');
l = legend('simple mean','weighted with n^{2}','slow g^{(2)}(0)','Location','best');
l.FontSize = 15;
hold off;
ylabel('\langle g^{(2)}(0) \rangle');
xlabel('Anzahl Quadraturen m in Statistik');
graphicsSettings;
ylim([1 3]);
print([filename '-g2Av.png'],'-dpng');
savefig([filename '-g2Av.fig']);
clf();

%% plot std
semilogx(nRes,g2Std,'ro');
xlim([min(nRes) max(nRes)]);
ylabel('\Delta g^{(2)}(0)');
xlabel('Anzahl Quadraturen  m  in Statistik');
graphicsSettings;
ylim([1 3]);
print([filename '-g2std.png'],'-dpng');
savefig([filename '-g2std.fig']);
clf();

%% plot frequency std
semilogx(fVector,g2Std,'ro');
ylabel('\Delta g^{(2)}(0)');
xlabel('f_{av} (Hz)');
graphicsSettings;
ylim([1 3]);
print([filename '-g2std-freq.png'],'-dpng');
savefig([filename '-g2std-freq.fig']);
clf();
%% plot frequency Av
semilogx(fVector,g2Av,'o');
hold on;
semilogx(fVector,g2weighted,'o');
semilogx(fVector,slowg2,'o');
l = legend('simple mean','weighted with n^{2}','slow g^{(2)}(0)','Location','best');
hold off;
ylabel('\langle g^{(2)}(0) \rangle');
xlabel('f_{av} (Hz)');
ylim([1 3]);
graphicsSettings;
print([filename '-g2Av-freq.png'],'-dpng');
savefig([filename '-g2Av-freq.fig']);
clf();
end