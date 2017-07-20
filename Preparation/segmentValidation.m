function [X1,X2,X3,theta] = segmentValidation(X1,X2,X3,theta,varargin)
%SEGMENTVALIDATION manual selection of correctly reconstructed phases
%
% Input Parameters:
%   X1,X2,X3 - quadrature data in the format [nPoints,nRecords,nSegments]
%   theta - reconstructed phase with cross correlation X1.*X3
%
% Output Paramters:
%   X1,X2,X3,theta - only the manually accepted parts of the corresponding
%       input arrays are returned

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
parse(p,varargin{:});

%% Let the user choose the segments to reject
% First, we need to calculate the smoothed cross correlation that was used
% for phase reconstruction.
ys = smoothCrossCorr(X1,X3);

% Second, we loop over all piezo segments and ask the user whether he wants
% to reject the segment or not
[~,~,nSegments] = size(X1);
segSel = zeros(nSegments,1);
for iSeg = 1:nSegments
    plot(ys(:,iSeg)); hold on;
    plot(theta(:,iSeg)); hold off;
    ylabel(['Segment ',num2str(iSeg)]);
    title({'Mousebutton: Reject','Keyboard: Accept'});
    key = waitforbuttonpress;
    if key == 1 % keyboard
        segSel(iSeg) = 1;
    else % e.g. mouse
        segSel(iSeg) = 0;
    end
end % iSeg

% Third, apply the selection to the output parameters
X1 = X1(:,:,segSel==1);
X2 = X2(:,:,segSel==1);
X3 = X3(:,:,segSel==1);
theta = theta(:,segSel==1);

%% Optional: Save the output paramters
if ~isempty(p.Results.Filename)
    save(['validatedSegments-',p.Results.Filename,'.mat'],'X1','X2','X3','theta');
end

end

