function handles = MeasureGPperSingleCell(handles)

%%%MF on 120111 in pelkmans lab
%%%This module works only on laurdan stained images, the GP formula has
%%%been taken from 'Visualizing lipid structure and raft domains in living
%%%cells with two-photon microscopy, Gaus & al 2003 PNAS

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the greyscale image for the 400 to 460 window?
%infotypeVAR01 = imagegroup
ImageName1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the greyscale image for the 470 to 530 window?
%infotypeVAR02 = imagegroup
ImageName2 = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the objects that you want to define the GP values?
%choiceVAR03 = Do not use
%infotypeVAR03 = objectgroup
ObjectNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow

    %%% Reads (opens) the images you want to analyze and assigns it to a variable,
    %%% "OrigImage".
    OrigImage1 = CPretrieveimage(handles,ImageName1,ModuleName);
    OrigImage2 = CPretrieveimage(handles,ImageName2,ModuleName);
    ObjectName = ObjectNameList{1};
    GPimage = zeros(size(OrigImage1));
    GPimageDetailed = zeros(size(OrigImage1));
    
    %%% Retrieves the label matrix image that contains the segmented objects which
    %%% will be measured with this module.
    LabelMatrixImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');

    %%% For the cases where the label matrix was produced from a cropped
    %%% image, the sizes of the images will not be equal. So, we crop the
    %%% LabelMatrix and try again to see if the matrices are then the
    %%% proper size. Removes Rows and Columns that are completely blank.
    if any(size(OrigImage1) < size(LabelMatrixImage))
        ColumnTotals = sum(LabelMatrixImage,1);
        RowTotals = sum(LabelMatrixImage,2)';
        warning off all
        ColumnsToDelete = ~logical(ColumnTotals);
        RowsToDelete = ~logical(RowTotals);
        warning on all
        drawnow
        CroppedLabelMatrix = LabelMatrixImage;
        CroppedLabelMatrix(:,ColumnsToDelete,:) = [];
        CroppedLabelMatrix(RowsToDelete,:,:) = [];
        clear LabelMatrixImage
        LabelMatrixImage = CroppedLabelMatrix;
        %%% In case the entire image has been cropped away, we store a single
        %%% zero pixel for the variable.
        if isempty(LabelMatrixImage)
            LabelMatrixImage = 0;
        end
    end

    if any(size(OrigImage1) ~= size(LabelMatrixImage))
        error(['Image processing was canceled in the ', ModuleName, ' module. The size of the image you want to measure is not the same as the size of the image from which the ',ObjectName,' objects were identified.'])
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MAKE MEASUREMENTS & SAVE TO HANDLES STRUCTURE %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow

    %%% Initialize measurement structure
    GPvaluePerObject = [];
    BasicFeatures    = {'MeanGPvaluePerSingleCell'};

    %%% Get pixel indexes (fastest way), and count objects
    props = regionprops(LabelMatrixImage,'PixelIdxList');
    ObjectCount = length(props);
    
    if ObjectCount > 0
        TotalGPvaluePerPixel = [];
        GPvaluePerObject = zeros(ObjectCount,12);
        
        for Object = 1:ObjectCount
            
            if isempty(props(Object).PixelIdxList)
                disp(sprintf('%s: BEREND BUGFIX: DETECTED EMPTY OBJECT, SKIPPING MEASUREMENT.',mfilename))
                continue
            end
            
           
           PixelCount = length(props(Object).PixelIdxList);
           GPvaluePerPixel = zeros(PixelCount,1);
           for Pixel = 1:PixelCount
               
               Pix400460 = OrigImage1(props(Object).PixelIdxList(Pixel));
               Pix470530 = OrigImage2(props(Object).PixelIdxList(Pixel));
               
               
               GPvaluePerPixel(Pixel)= (Pix400460-Pix470530)/(Pix400460+Pix470530);
               GPimageDetailed(props(Object).PixelIdxList(Pixel)) = GPvaluePerPixel(Pixel);
               
           end    
           %%%%%%%This one is for plotting the GP distribution%%%%%%%
           TotalGPvaluePerPixel = [TotalGPvaluePerPixel;GPvaluePerPixel];
           %%%%%%%This one is for the usual 'against LCD' plotting%%%    
           GPvaluePerObject(Object) = nanmean(GPvaluePerPixel);
           for Pixel = 1:PixelCount
           GPimage(props(Object).PixelIdxList(Pixel)) = GPvaluePerObject(Object);     
           end 
        end
    else
        TotalGPvaluePerPixel = 0;
        GPvaluePerObject(1) = 0;
    end
     
    ImageName = 'OrigMixed';
    
    %%% Save measurements %%%GPImage creation and storing
    handles.Measurements.(ObjectName).(['SingleCellGPanalysis_',ImageName,'Features']) = BasicFeatures;
    handles.Measurements.(ObjectName).(['SingleCellGPanalysis_',ImageName])(handles.Current.SetBeingAnalyzed) = {GPvaluePerObject};
    
    handles.Measurements.(ObjectName).(['PixelGP_',ImageName,'Features']) = 'GPforEachPixelForDistribAnalysis';
    handles.Measurements.(ObjectName).(['PixelGP_',ImageName])(handles.Current.SetBeingAnalyzed) = {TotalGPvaluePerPixel};
    
    fieldname = ['GP_in_',ObjectName];

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

    % determine output directory, if BATCH, create GP directory,
    % otherwise, use the default output directory to dump the segmentation
    % images.
    strOutputDir = handles.Current.DefaultOutputDirectory;
    if strcmp(getlastdir(strOutputDir),'BATCH')
        strOutputDir = strrep(strOutputDir, [filesep,'BATCH'],[filesep,'GP']);
        if ~fileattrib(strOutputDir)
            disp(sprintf('%s: creating default output directory GP in %s',mfilename,getbasedir(strOutputDir)))
            mkdir(strOutputDir)
        end
    end

    strFullFileName = fullfile(strOutputDir, strTmpFileName);
    matCmap = redgreencmap(255);
    matCmap(1,:)=[1,1,1];
    GPimage = (GPimage + 1);
    GPimage = (GPimage / 2);
    GPimage = uint8(GPimage * 255);
    imwrite(GPimage,matCmap,strFullFileName,'png');
    
    strTmpFileName2 = ['Pix_' strTmpFileName];
    strFullFileName2 = fullfile(strOutputDir, strTmpFileName2);
    matCmap = redgreencmap(255);
    matCmap(1,:)=[1,1,1];
    GPimageDetailed = (GPimageDetailed + 1);
    GPimageDetailed = (GPimageDetailed / 2);
    GPimageDetailed = uint8(GPimageDetailed * 255);
    imwrite(GPimageDetailed,matCmap,strFullFileName2,'png');
    
end
   
    
