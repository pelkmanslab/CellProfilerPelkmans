function handles = ImageBySegmentation(handles)

% ImageBySegmentation
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Generates a masked image based upon the segmentation of objects. Values
% outside of this mask will be set to 0.
% *************************************************************************
%
% See also MeasureImageIntensity.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   TS
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4526 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = What do you want to call the Image derived from the Segmentation?
%defaultVAR01 = SegmenationImage
%infotypeVAR01 = imagegroup indep
iOutputName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Which objects should be used?
%infotypeVAR02 = objectgroup
iObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = Which Spot amplification do you want to use?
%choiceVAR03 = Complete Object
%choiceVAR03 = Centroid of Object
iObjectDefinition = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Which intensities should be used?
%choiceVAR04 = Zero And One
%choiceVAR04 = Zero And From Image
iScaleMethod = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = (Optional: Template for Intensities)?
%infotypeVAR05 = imagegroup
iInputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu



%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Retrieves the label matrix image that contains the segmented objects which
%%% will be measured with this module.
LabelMatrixImage = CPretrieveimage(handles,['Segmented', iObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Reads (opens) the image you want to analyze and assigns it to a variable,
%%% "OrigImage".
switch iScaleMethod
    case 'Zero And From Image'
        OrigImage = CPretrieveimage(handles,iInputImageName,ModuleName);
        if any(size(OrigImage) ~= size(LabelMatrixImage))
            error(['Image processing was canceled in the ', ModuleName, ' module. The size of the image you want to measure is not the same as the size of the image from which the ',ObjectName,' objects were identified.'])
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAKE MEASUREMENTS & SAVE TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ImageMask = zeros(size(LabelMatrixImage));

switch iObjectDefinition
    case 'Complete Object'
        f = LabelMatrixImage>1;
        ImageMask(f) = 1;
    case 'Centroid of Object'
        if max(LabelMatrixImage(:)>0)
            tmp = regionprops(LabelMatrixImage,'Centroid');
            Centroid = cat(1,tmp.Centroid); clear tmp:
            Centroid = round(Centroid);
            Centroid = [Centroid(:,2), Centroid(:,1)]; %reorganize that row is first so that corresponds with subindices of image
            ImageMask(Centroid) = 1;
        end
end

switch iScaleMethod
    case 'Zero And One'
        OutputImage = ImageMask;
    case 'Zero And From Image'
        OutputImage = OrigImage.*ImageMask;
end

%%% Save measurements
handles.Pipeline.(iOutputName) = OutputImage;
% 
% %%% Report measurements
% if any(findobj == ThisModuleFigureNumber);
%     CPimagesc(OutputImage,handles);
%     title(sprintf('Intensity Image created from %s Segmentation', iObjectName));
% end
end