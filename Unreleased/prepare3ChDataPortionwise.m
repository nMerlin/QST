function [ X1, X2, X3 ] = prepare3ChDataPortionwise( filenameLO, filenameSIG )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

maxRecordings = 1581; %maximum number of segments in one portion

[data8bitLO,configLO,~]=load8BitBinary(filenameLO,'dontsave');

[data8bitSIG,configSIG,~,portions]=load8BitBinaryPortion(filenameSIG,1,'dontsave');

if portions == 1
    for iCh = 1:3
        % Compute number of LO photons
        dispstat(['Computing number of LO photons for Channel ' num2str(iCh)],'timestamp','keepthis',0);
        XLO = computeQuadratures(data8bitLO(:,:,iCh), configLO, CALIBRATION_CH1);

        % Calculate the variance piece-wise to compensate slow drifts (e.g. piezos)
        NLO = mean(var(XLO));

        % Compute quadratures for target quantum state
        dispstat(['Computing quadratures for target quantum state for Channel ' num2str(iCh)],...
          'timestamp','keepthis',0);
        X = computeQuadratures(data8bitSIG(:,:,iCh), configSIG, CALIBRATION_CH1);

        % Calibration of quadratures to vacuum state
        X = Norm * X / sqrt(NLO);

        switch iCh
            case 1
                X1 = X;
            case 2
                X2 = X;
            case 3
                X3 = X;
        end
    end

else
    
    X = zeros(999, portions*maxRecordings,3);
    % Compute quadratures for target quantum state
    dispstat('Computing quadratures for target quantum state ...',...
          'timestamp','keepthis',0);
    
    for iCh = 1:3   %first portion
        X(:,1:maxRecordings,iCh) = computeQuadratures(data8bitSIG(:,:,iCh),...
            configSIG, CALIBRATION_CH1);
    end
        for portion = 2:portions  %load and process remaining portions
            dispstat(['Portion ' num2str(portion)],...
          'timestamp','keepthis',0);
            [data8bitSIGnext,~,~,~]=load8BitBinaryPortion(filenameSIG,portion,'dontsave');
            for iCh = 1:3
                Xnext = computeQuadratures(data8bitSIGnext(:,:,iCh), configSIG, CALIBRATION_CH1);
                X(:,1+(portion -1)*maxRecordings : portion*maxRecordings, iCh) = Xnext;
            end
        end
        
     for iCh = 1:3
        % Compute number of LO photons
        dispstat(['Computing number of LO photons for Channel ' num2str(iCh)],'timestamp','keepthis',0);
        XLO = computeQuadratures(data8bitLO(:,:,iCh), configLO, CALIBRATION_CH1);
        % Calculate the variance piece-wise to compensate slow drifts (e.g. piezos)
        NLO = mean(var(XLO));

        % Calibration of quadratures to vacuum state
        X(:,:,iCh) = Norm * X(:,:,iCh) / sqrt(NLO);
     end
     X1 =  X(:,1:find(mean(X(:,:,1)),1,'last'),1); %cut off trailing zeros
     X2 =  X(:,1:find(mean(X(:,:,2)),1,'last'),2); %cut off trailing zeros
     X3 =  X(:,1:find(mean(X(:,:,3)),1,'last'),3); %cut off trailing zeros
        
end

end

