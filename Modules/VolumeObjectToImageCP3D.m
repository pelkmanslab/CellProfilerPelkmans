function handles = VolumeObjectToImageCP3D(handles)

% Help for VolumeObjectToImageCP3D:
% Category: Other
%
% SHORT DESCRIPTION:
% Transforms a 3D object segmentation into a 2D intensity image, where
% intensities show the number of voxels occupied with the object. E.g.:
% simple quantification of 3D properties by 2D modules 
%
% Example: Measure Volume Of Nuclei grown in monolayer
% AddSegmentationVolumeCP3D -> VolumeObjectToImageCP3D ->
% MeasureObjectIntensity
% *************************************************************************
%
% Does an image projection
%
%   Authors:
%   Thomas Stoeger
%   Nico Battich
%   Lucas Pelkmans
%
% Website: http://www.pelkmanslab.org
% ***********************************************************************
%
% $Revision: 3524 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What is the name of the 3D object?
%choiceVAR01 = Nuclei3D
%infotypeVAR01 = objectgroup
strSegmentation3D = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the image?
%defaultVAR02 = VolumeImage
%infotypeVAR02 = imagegroup indep
strOccupancyImage = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% Check / Load input segmetnation
fieldname = ['Segmented', strSegmentation3D];
if isfield(handles.Pipeline.(fieldname), 'Format')
    switch handles.Pipeline.(fieldname).Format
        case 'SegmentationCC'
            SegmentationCC = handles.Pipeline.(fieldname).Label;
        otherwise
            error('Current 3D format not supported for this module');
    end
else
    error('This module only works with 3D segmentations');
end
    
% Compute occupancy image
ViewDirectionOfInterest = 'Z';
ViewFromZ = createOccupancyImage2(SegmentationCC,ViewDirectionOfInterest);
ViewFromZ = ViewFromZ./65535; % Adhere to CP-2D convention of intensities of images

%%%%%%%%%%%%%%%%%%%%%%%
%%% STORE RESULTS   %%%
%%%%%%%%%%%%%%%%%%%%%%%

handles.Pipeline.(strOccupancyImage) = ViewFromZ;

%%%%%%%%%%%%%%%
%%% DISPLAY  %%
%%%%%%%%%%%%%%%

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if ~CPisHeadless()
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure(ViewFromZ,'TwoByTwo',ThisModuleFigureNumber);
        end
        
        CPimagesc(ViewFromZ,handles);
        colormap('jet')
        title(sprintf('Layers containing objects within %s , cycle #%d', strSegmentation3D ,handles.Current.SetBeingAnalyzed));
    end
end
end