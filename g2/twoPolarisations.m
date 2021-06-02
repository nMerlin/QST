function [] = twoPolarisations()
%% Validate and parse input arguments
folder1 = 'g2data-X1';
folder2 = 'g2data-X2';
%%
targetfolder = 'TwoPolarisationsPlots';
if ~exist(targetfolder,'dir')
    mkdir(targetfolder)
end

[filenames,~,~]= getParametersFromFilenames('Folder',folder1);
%[filenames,numbers,Is]= getParametersFromFilenames('Folder2',folder);
%[period_mean_ons,period_mean_offs,fliprates] = deal(zeros(length(filenames),1));

for fileI = 1:length(filenames)
    load([folder1 '\' cell2mat(filenames(fileI))], 'times', 'ada', 'g2vec');
    times1 = times;
    ada1 = ada;
    g2vec1 = g2vec;
    load([folder2 '\' cell2mat(filenames(fileI))], 'times', 'ada', 'g2vec');
    times2 = times;
    ada2 = ada;
    g2vec2 = g2vec;
    plotFilename = [targetfolder '\' cell2mat(filenames(fileI)) ];
    %% Create plot of g2 and ada against time
    [times1, ~] = convenientUnits(times1,'s');
    [times2,units] = convenientUnits(times2,'s');

    % plot g2
    subplot(2,1,1)
    plot(times1,g2vec1,'b', times2,g2vec2,'r');
    axis([0 max(times1) 0.5 2.5]);
    %axis([0 max(times) 0.5 1.5]);
    fontname = 'Times New Roman';
    fontsize1 =22;
    fontsize2 =20;
    %xlabel(['$t$(',units,')'],'FontSize',fontsize1,'Interpreter','latex');
    ylabel('$g^{(2)}(0,t)$','FontSize',fontsize1,'Interpreter','latex');
    set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
        'YColor', 'k','FontName',fontname);
    legend('g2Polarisation1','g2Polarisation2');
    grid on;
    %%
    % plot ada
    subplot(2,1,2)
    plot(times1,ada1,'k', times2, ada2, 'r');
    max1=max(ada1);
    max2=max(ada2);
    axis([0 max(times2) 0 max(max1, max2)+1]);
    fontname = 'Times New Roman';
    fontsize1 =22;
    fontsize2 =20;
    xlabel(['$t$ (',units,')'],'FontSize',fontsize1,'Interpreter','latex');
    ylabel('$\bar{n} (t)$','FontSize',fontsize1,'Interpreter','latex');
    set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
        'YColor', 'k','FontName',fontname);
    legend('Polarisation1', 'Polarisation2');
    grid on;

    %%

    %Save figure 
    print(strcat(plotFilename,'G2-', '.png'), '-dpng');
    savefig(strcat(plotFilename,'G2-', '.fig'));

    %%
    % plot ada
    subplot(2,1,2)
    plot(times1,ada1,'k', times2, ada2, 'r');
    max1=max(ada1);
    max2=max(ada2);
    axis([0 max(times2) 0 max(max1, max2)+1]);
    fontname = 'Times New Roman';
    fontsize1 =22;
    fontsize2 =20;
    xlabel(['$t$ (',units,')'],'FontSize',fontsize1,'Interpreter','latex');
    ylabel('$\bar{n} (t)$','FontSize',fontsize1,'Interpreter','latex');
    set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
        'YColor', 'k','FontName',fontname);
    legend('Polarisation1', 'Polarisation2');
    grid on;

end



