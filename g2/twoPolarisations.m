function [] = twoPolarisations()
%% Validate and parse input arguments
folder1 = 'g2data-X1';
folder2 = 'g2data-X2';
%%
targetfolder = 'TwoPolarisationsPlotsWithG2';
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
    plot(times1,g2vec1,'k', times2,g2vec2,'r');
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
%      %scale ada
%      max1=max(ada1);
%      max2=max(ada2);
%      s1 = 1/max1;
%      s2 = 1/max2;
%      ada1 = s1* ada1;
%      ada2 = s2* ada2;
     %added adas
     %ada1 = ada1 +ada2;
    % plot ada
    subplot(2,1,2)
   % plot(times1, ada1, 'k');
    plot(times1, ada1, 'k', times2, ada2, 'r');
    %title('Photon numbers added')
    max1=max(ada1);
    max2=max(ada2);
    %axis([0 max(times2) 0 1.1]);
    axis([0 max(times2) 0 max(max1, max2)+1]);
    fontname = 'Times New Roman';
    fontsize1 =22;
    fontsize2 =20;
    xlabel(['$t$ (',units,')'],'FontSize',fontsize1,'Interpreter','latex');
    ylabel('$\bar{n} (t)$','FontSize',fontsize1,'Interpreter','latex');
    set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
        'YColor', 'k','FontName',fontname);
   % legend('Polarisation1', 'Polarisation2');
   %legend('Location', 'northoutside');
    grid on;
    grid minor;

    %%

    %Save figure 
    print(strcat(plotFilename,'G2-', '.png'), '-dpng');
    savefig(strcat(plotFilename,'G2-', '.fig'));

    %%
    %scale ada
%     max1=max(ada1);
%     max2=max(ada2);
%     s1 = 1/max1;
%     s2 = 1/max2;
%     ada1 = s1* ada1;
%     ada2 = s2* ada2;
    % plot ada
%     subplot(2,1,2)
%     plot(times1,ada1,'k', times2, ada2, 'r');
%     axis([0 max(times2) 0 2]);
%     fontname = 'Times New Roman';
%     fontsize1 =22;
%     fontsize2 =20;
%     xlabel(['$t$ (',units,')'],'FontSize',fontsize1,'Interpreter','latex');
%     ylabel('$\bar{n} (t)$','FontSize',fontsize1,'Interpreter','latex');
%     set(gca, 'LineWidth', 1.5,'FontSize',fontsize2, 'XColor',[0 0 0],...
%         'YColor', 'k','FontName',fontname);
%     legend('Polarisation1', 'Polarisation2');
%     grid on;

end



