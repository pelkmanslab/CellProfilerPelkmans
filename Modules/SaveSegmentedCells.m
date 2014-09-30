function handles = SaveSegmentedCells(handles)

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

    fieldname = ['Segmented',strObjectName];

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
    strTmpFileName = [strOrigImageName,'_',fieldname,'.png'];
    
    % if analyzing subdirectories, do not take along the subdirectory
    % structure... but replaces file separators with underscores to
    % preserver the additional information of the path.
    if ~isempty(strfind(strOrigImageName,filesep))
        strTmpFileName = strrep(strTmpFileName,filesep,'_');
    end

    % determine output directory, if BATCH, create SEGMENTATION directory,
    % otherwise, use the default output directory to dump the segmentation
    % images.
    strOutputDir = handles.Current.DefaultOutputDirectory;
    if strcmp(getlastdir(strOutputDir),'BATCH')
        strOutputDir = strrep(strOutputDir, [filesep,'BATCH'],[filesep,'SEGMENTATION']);
        if ~fileattrib(strOutputDir)
            disp(sprintf('%s: creating default output directory SEGMENTATION in %s',mfilename,getbasedir(strOutputDir)))
            mkdir(strOutputDir)
        end
    end

    strFullFileName = fullfile(strOutputDir, strTmpFileName);
    LabelMatrixImage = CPretrieveimage(handles,fieldname,ModuleName,'MustBeGray','DontCheckScale');
    imwrite(uint16(LabelMatrixImage),strFullFileName,'png');
    
end
