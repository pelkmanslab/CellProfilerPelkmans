function handles = SaveSegmentedCellsOverlay(handles)

% Help for the SaveSegmentedCellsOverlay module:
% Category: Other
%
% SHORT DESCRIPTION:
% Identifies objects (e.g. cell edges) using "seed" objects identified by
% an Identify Primary module (e.g. nuclei) and overlays their outline on
% the original image (e.g. DAPI staining).
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1725 $
% 
% [Markus Herrmann, 06/2013]


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%defaultVAR01 = Cells
%infotypeVAR01 = objectgroup indep
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = How did you call the image you want to use for overlay?
%defaultVAR02 = OrigBlue
%infotypeVAR02 = objectgroup indep
strOverlayImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%%%VariableRevisionNumber = 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    drawnow

    fieldname = ['Segmented',strObjectName];
    
    fieldname2 = ['OverlaySegmented',strObjectName];
    

    % create new image name
    % compatibility with both old and new CellProfiler measurement storing
    % schemas.
    if isfield(handles.Measurements.Image, 'FileNames')
        strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
    else
        cellstrImageFieldNames = fieldnames(handles.Measurements.Image);
        strImageFieldName = cellstrImageFieldNames{find(strncmp(cellstrImageFieldNames,'FileName_',9),1,'first')};
        strOrigImageName = char(handles.Measurements.Image.(strImageFieldName){handles.Current.SetBeingAnalyzed});
    end    
    
    matDotIndices=strfind(strOrigImageName,'.');
    % new CP apparently removes file extensions from image names
    if ~isempty(matDotIndices)
        strOrigImageName = strOrigImageName(1,1:matDotIndices(end)-1);
    end
    strTmpFileName = [strOrigImageName,'_',fieldname2,'.png'];
    
    % if analyzing subdirectories, do not take along the subdirectory
    % structure... but replaces file separators with underscores to
    % preserver the additional information of the path.
    if ~isempty(strfind(strOrigImageName,filesep))
        strTmpFileName = strrep(strTmpFileName,filesep,'_');
    end

    % determine output directory, if BATCH, create SEGMENTATIONOVERLAY
    % directory, otherwise, use the default output directory to dump the
    % segmentation images.
    strOutputDir = handles.Current.DefaultOutputDirectory;
    if strcmp(getlastdir(strOutputDir),'BATCH')
        strOutputDir = strrep(strOutputDir, [filesep,'BATCH'],[filesep,'SEGMENTATIONOVERLAY']);
        if ~fileattrib(strOutputDir)
            disp(sprintf('%s: creating default output directory SEGMENTATIONOVERLAY in %s',mfilename,getbasedir(strOutputDir)))
            mkdir(strOutputDir)
        end
    end

    strFullFileName = fullfile(strOutputDir, strTmpFileName);
    
    % retrieve segmentation image
    LabelMatrixImage = CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale');
    % retrieve image used for overlay
    GrayImage = CPretrieveimage(handles,strOverlayImageName,ModuleName,'MustBeGray','DontCheckScale');
    % create segmentation outline image
    if ~max(max(LabelMatrixImage))==0
        cellB = bwboundaries(LabelMatrixImage,'holes');
        matB = cell2mat(cellB);
        linInd = sub2ind(size(GrayImage),matB(:,1),matB(:,2));
        imSegOver = zeros(size(GrayImage));
        imSegOver(linInd) = 1;
        imSegOver = logical(imSegOver);
    else
        imSegOver = LabelMatrixImage;
    end
    % create and write overlay image
    OverlayImage = imoverlay(GrayImage,imSegOver,[1,1,0]);
    imwrite(OverlayImage,strFullFileName,'png');
    
end
