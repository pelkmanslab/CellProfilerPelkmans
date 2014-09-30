function handles = SaveSegmentedCells_MPcycle(handles)

% Help for the SaveSegmentedCells module:
% Category: Other
%
% SHORT DESCRIPTION:
% Identifies objects (e.g. cell edges) using "seed" objects identified by
% an Identify Primary module (e.g. nuclei).
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object?
%defaultVAR01 = Cells
%infotypeVAR01 = objectgroup indep
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%%%VariableRevisionNumber = 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    drawnow
    
    %%% retrieve shift descriptors from handles
    shift = handles.shiftDescriptor;
    
    %%% get index of current image
    strOrigImageName = char(handles.Measurements.Image.FileNames{handles.Current.SetBeingAnalyzed}{1,1});
    
    fieldname = ['Segmented',strObjectName];
    
    matDotIndices=strfind(strOrigImageName,'.');
    % new CP apparently removes file extensions from image names
    if ~isempty(matDotIndices)
        strOrigImageName = strOrigImageName(1,1:matDotIndices(end)-1);
    end
    
    % if analyzing subdirectories, do not take along the subdirectory
    % structure... but replaces file separators with underscores to
    % preserver the additional information of the path.
    if ~isempty(strfind(strOrigImageName,filesep))
        strOrigImageName = strrep(strOrigImageName,filesep,'_');
    end
    
    % built absolute SEGMENTATION path from relative path stored in handles
    strSegmentationDir = [strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], filesep), shift.SegmentationDirectory];
    
    % built full segmentation image filename from filename trunk stored in handles
    strSegmentationFileName = [regexprep(strOrigImageName,'.+(_\w{1}\d{2}_)',sprintf('%s$1',shift.SegmentationFileNameTrunk)),'_Segmented',strObjectName,'.png'];
    
    % determine output directory, if BATCH, create SEGMENTATION directory,
    % otherwise, use the default output directory to dump the segmentation
    % images.
    strOutputDir = strSegmentationDir;
    if strcmp(getlastdir(strOutputDir),'SEGMENTATION')
        strOutputDir = strrep(strOutputDir, [filesep,'SEGMENTATION'],[filesep,'SEGMENTATION_ALIGNED']);
        if ~fileattrib(strOutputDir)
            fprintf('%s: creating default output directory SEGMENTATION in %s',mfilename,getbasedir(strOutputDir))
            mkdir(strOutputDir)
        end
    end

    strFullFileName = fullfile(strOutputDir, strSegmentationFileName);
    LabelMatrixImage = CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale');
    imwrite(uint16(LabelMatrixImage),strFullFileName,'png');
    
end
