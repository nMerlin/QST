function [X,theta,iSelect] = selectRegionAroundZeroChangeHusimi(O1,O2,O3,oTheta,oThetaMira,varargin)
%function [X,theta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,varargin)
%SELECTREGION Return all values of O3 in the specified (O1,O2)-region
%
%   Rectangle:
%   selectRegion(O1,O2,O3,theta,'rectangle',[x y w h]) - lower left corner
%   at (x,y) and (width,height) = (w,h)
%
% Options:
%   Output - 'print' saves plotted figure to file

%% Validate and parse input arguments
p = inputParser;
defaultType = 'dot';
defaultPosition = [2 2 0.25];
defaultPlotOpt = 'hide';
defaultOutput = 'figure';
defaultFilename = '';
addParameter(p,'Type',defaultType,@isstr);
addParameter(p,'Position',defaultPosition,@isvector);
addParameter(p,'Plot',defaultPlotOpt,@isstr);
addParameter(p,'Output',defaultOutput,@isstr);
addParameter(p,'Filename',defaultFilename,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,output,plotopt,position,type] = c{:};



%% Selection process for husimiPhase only
% switch type
%     case 'phase'
%        % For the case of a phase region around zero or
%         %2*pi problems are solved with a 2*pi-shift of the boundaries of thetaMin/thetaMax
%         mTheta = position(1);
%         dTheta = position(2);
%         thetaMin = mTheta-dTheta/2;
%         thetaMax = mTheta+dTheta/2;
% 
%         husimiPhase = pi+atan2(O2,O1);
%         if thetaMin<0
%             if thetaMax <0
%                 iSelect = find(husimiPhase>(thetaMin+2*pi) & ...
%             husimiPhase<(thetaMax+2*pi) & (O1.^2+O2.^2)>0.5); 
%             else
%              iSelect1 = find(husimiPhase>(thetaMin+2*pi) & (O1.^2+O2.^2)>0.5);
%             iSelect2 = find(husimiPhase<thetaMax & (O1.^2+O2.^2)>0.5);
%             iSelect = cat(1,iSelect1,iSelect2);
%             end
%         elseif thetaMax > 2*pi
%             if thetaMin > 2*pi
%                 iSelect = find(husimiPhase>(thetaMin-2*pi) & ...
%             husimiPhase<(thetaMax-2*pi) & (O1.^2+O2.^2)>0.5);
%             else
%              iSelect1 = find(husimiPhase>thetaMin & (O1.^2+O2.^2)>0.5); 
%              iSelect2 = find(husimiPhase <(thetaMax-2*pi)&(O1.^2+O2.^2)>0.5);
%              iSelect = cat(1,iSelect1,iSelect2);
%             end
%         else
%              iSelect = find(husimiPhase>thetaMin & ...
%             husimiPhase<thetaMax & (O1.^2+O2.^2)>0.5);
%         end
%         
%     case 'phaseAndAmplitude'
%         mTheta = position(1);
%         dTheta = position(2);
%         thetaMin = mTheta-dTheta/2;
%         thetaMax = mTheta+dTheta/2;
%         r = position(3);
%         w = position(4);
%         husimiPhase = pi+atan2(O2,O1);
%         
%          if thetaMin<0
%             if thetaMax <0
%                 iSelect = find(husimiPhase>(thetaMin+2*pi) & ...
%             husimiPhase<(thetaMax+2*pi) & (O1.^2+O2.^2)>0.5 & ...
%             sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2); 
%             else
%             iSelect1 = find(husimiPhase>(thetaMin+2*pi) & (O1.^2+O2.^2)>0.5 & ...
%                 sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
%             iSelect2 = find(husimiPhase<thetaMax & (O1.^2+O2.^2)>0.5 & ...
%                 sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
%             iSelect = cat(1,iSelect1,iSelect2);
%             end
%         elseif thetaMax > 2*pi
%             if thetaMin > 2*pi
%                 iSelect = find(husimiPhase>(thetaMin-2*pi) & ...
%             husimiPhase<(thetaMax-2*pi) & (O1.^2+O2.^2)>0.5 & ...
%             sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
%             else
%              iSelect1 = find(husimiPhase>thetaMin & (O1.^2+O2.^2)>0.5 & ...
%                 sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2); 
%              iSelect2 = find(husimiPhase <(thetaMax-2*pi)&(O1.^2+O2.^2)>0.5 & ...
%                 sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
%              iSelect = cat(1,iSelect1,iSelect2);
%             end
%         else
%              iSelect = find(husimiPhase>thetaMin & ...
%                 husimiPhase<thetaMax & (O1.^2+O2.^2)>0.5 & ...
%                 sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
%          end
% 
% end
% X = O3(iSelect);
%%
%Selection process for totalPhase
switch type
    case 'phase'
        %For the case of a phase region around zero or
        %2*pi problems are solved with a 2*pi-shift of the boundaries of thetaMin/thetaMax
        mTheta = position(1);
        dTheta = position(2);
        thetaMin = mTheta-dTheta/2;
        thetaMax = mTheta+dTheta/2;

        husimiPhase = atan2(O2,O1);
        thetaTotal = mod(husimiPhase - oTheta + oThetaMira,2*pi);
        if thetaMin<0
            if thetaMax <0
                iSelect = find(thetaTotal>(thetaMin+2*pi) & ...
            thetaTotal<(thetaMax+2*pi) & (O1.^2+O2.^2)>0.5); 
            else
             iSelect1 = find(thetaTotal>(thetaMin+2*pi) & (O1.^2+O2.^2)>0.5);
            iSelect2 = find(thetaTotal<thetaMax & (O1.^2+O2.^2)>0.5);
            iSelect = cat(1,iSelect1,iSelect2);
            end
        elseif thetaMax > 2*pi
            if thetaMin > 2*pi
                iSelect = find(thetaTotal>(thetaMin-2*pi) & ...
            thetaTotal<(thetaMax-2*pi) & (O1.^2+O2.^2)>0.5);
            else
             iSelect1 = find(thetaTotal>thetaMin & (O1.^2+O2.^2)>0.5); 
             iSelect2 = find(thetaTotal <(thetaMax-2*pi)&(O1.^2+O2.^2)>0.5);
             iSelect = cat(1,iSelect1,iSelect2);
            end
        else
             iSelect = find(thetaTotal>thetaMin & ...
            thetaTotal<thetaMax & (O1.^2+O2.^2)>0.5);
        end
        
    case 'phaseAndAmplitude'
        mTheta = position(1);
        dTheta = position(2);
        thetaMin = mTheta-dTheta/2;
        thetaMax = mTheta+dTheta/2;
        r = position(3);
        w = position(4);
        husimiPhase = atan2(O2,O1);
        thetaTotal = mod(husimiPhase - oTheta + oThetaMira,2*pi);
        
         if thetaMin<0
            if thetaMax <0
                iSelect = find(thetaTotal>(thetaMin+2*pi) & ...
            thetaTotal<(thetaMax+2*pi) & (O1.^2+O2.^2)>0.5 & ...
            sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2); 
            else
            iSelect1 = find(thetaTotal>(thetaMin+2*pi) & (O1.^2+O2.^2)>0.5 & ...
                sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
            iSelect2 = find(thetaTotal<thetaMax & (O1.^2+O2.^2)>0.5 & ...
                sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
            iSelect = cat(1,iSelect1,iSelect2);
            end
        elseif thetaMax > 2*pi
            if thetaMin > 2*pi
                iSelect = find(thetaTotal>(thetaMin-2*pi) & ...
            thetaTotal<(thetaMax-2*pi) & (O1.^2+O2.^2)>0.5 & ...
            sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
            else
             iSelect1 = find(thetaTotal>thetaMin & (O1.^2+O2.^2)>0.5 & ...
                sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2); 
             iSelect2 = find(thetaTotal <(thetaMax-2*pi)&(O1.^2+O2.^2)>0.5 & ...
                sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
             iSelect = cat(1,iSelect1,iSelect2);
            end
        else
             iSelect = find(thetaTotal>thetaMin & ...
                thetaTotal<thetaMax & (O1.^2+O2.^2)>0.5 & ...
                sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
         end

end
X = O3(iSelect);
%  
%% Phase calculation
theta = mod(oTheta + atan2(O2,O1),2*pi);
theta = theta(iSelect);

%% Show selection
if strcmp(plotopt,'show')
    assessTheta(theta,X,'Husimi',{O1,O2,iSelect},'VarBins', ...
        200,'PhaseBins',200,'Output',output,'Filename',filename);
end

end
