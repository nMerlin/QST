function QSeries(O1,O2,O3,oTheta)

Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 

%make befor some steps e g select orthogonal
%make here switch case for several region selections eg fullcircle

Q = linspace(0,6,13);
[ expQmax,  expQ2max, delQmean, Unc, N]=deal(zeros(length(Q),1));

for i = 1:length(Q)
    [XQ,thetaQ] = selectRegion(O1,O2,O3,oTheta,'Type','Pline','Position',[Q(i) 0.5],'Plot','show');
    [XQdis,thetaQdis]=discretizeTheta(XQ,thetaQ,180);
    [ expQ, ~, expQ2, ~, delQ, ~, meanUnc,~ , meanN, ~] =...
        computeExpectations2( XQdis, thetaQdis, ['11-Qline-' strrep(num2str(Q(i)),'.',',')] );
    expQmax(i) = max(expQ);
    expQ2max(i) = max(expQ2);
    delQmean(i) = mean(delQ);
    Unc(i) = meanUnc;
    N(i) = meanN;
end

close all;
hold off;
plot(Q, expQmax, Q, expQ2max, Q, delQmean, Q, Unc, Q, N, 'linewidth', 2);

xlabel('Q_{PS}');
%axis([0 max(Q) 0  ]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<Q>_{max}$', '$<Q^{2}>_{max}$', '$<\Delta Q>$',...
     '$<\Delta Q \cdot \Delta P>$ ', '$<n>$ ', 'location', 'best');
title(strcat('11-Qline'));

    % Save plot
print('11-Qline-Qdependence', '-dpng');

end