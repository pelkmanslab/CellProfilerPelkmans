function handles = MeasureObjectEnvironment(handles)

% Help for the Measure Object Environment module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Calculates intensity and texture features in a defined environment 
% of a cell.
% *************************************************************************
%
% Given an image with objects identified (e.g. nuclei or cells), this
% module determines how many neighbors each object has. The user selects
% the radius/radii of a circle(s) around the centroid of an object (cell).
% The module then measures the intensity and texture within the defined 
% circles. 

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

k = handles.Current.SetBeingAnalyzed;

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the greyscale images you want to measure?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What second image do you want to overlay in the visualization?
%infotypeVAR02 = imagegroup
SecondImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the objects whose environment you want to measure?
%infotypeVAR03 = objectgroup
%defaultVAR03 = cells
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = How many pixels do you expand away from the cells centroid for each environment? (in ascending ranking) Examples 
%defaultVAR04 = 300 600 1200
strEnvironmentRadii = char(handles.Settings.VariableValues{CurrentModuleNum,4});


%%%VariableRevisionNumber = 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a variable,
%%% "OrigImage".
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','DontCheckScale');

%%% Reads (opens) the image which is used in the visualization of the defined environments and assigns it to a variable,
%%% "SecondImage".
SecondImage = CPretrieveimage(handles,SecondImageName,ModuleName,'MustBeGray','DontCheckScale');

% %%% Retrieves the label matrix image that contains the segmented objects
% %%% whose environment will be measured with this module.
% ObjectLabelMatrixImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');

% Creates vector containing the radii of environments to be analysed.
EnvironmentRadii = double(cell2mat(textscan(strEnvironmentRadii,'%d','delimiter',' ')));
nEnvironments = numel(EnvironmentRadii);

intObjectCountColumn = find(strcmpi(handles.Measurements.Image.ObjectCountFeatures,ObjectName));
intNrOfObjects = handles.Measurements.Image.ObjectCount{k}(intObjectCountColumn);

% initialize measurement output for current image
EnvironmentFeatures = cell(1,nEnvironments);
matMeasurement = NaN(intNrOfObjects,nEnvironments);

% creating matrix holding a disk-shaped template with specified radius 
for j = 1: nEnvironments
    matDisk = fspecial('disk',EnvironmentRadii(j)) > 0;
    % zero-centered pixel id list for the first ball (disk, whateva)
    matPixelIDList = [];
    [matPixelIDList(:,1),matPixelIDList(:,2)] = find(matDisk);
    matPixelIDList = matPixelIDList - EnvironmentRadii(j);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MAKE MEASUREMENTS & SAVE TO HANDLES STRUCTURE %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    EnvironmentFeatures{1,j} = sprintf('MeanIntensity with R=%d',EnvironmentRadii(j));

    if k == 1
        % initialize output
        handles.Measurements.(ObjectName).(sprintf('Environment_%s',ImageName)) = cell(1,handles.Current.NumberOfImageSets);
    end
    handles.Measurements.(ObjectName).(sprintf('Environment_%s',ImageName)){k} = matMeasurement;

    matObjectCentroids = handles.Measurements.(ObjectName).Location{k};

    if intNrOfObjects<1
        return
    end

    matImSize = size(OrigImage);% this is ROW & COLUMN

    matMeasurementCounterImage = zeros(size(OrigImage));

    %loop over all objects of an image
    for i = 1:intNrOfObjects
        matCurObjCentroids = matObjectCentroids(i,:);

        % shift disk to be centered on current object (nucleus centroid)
        matCurObjDiskPixels = matPixelIDList;

        % double check, X and Y coordinates are not always the same as row and
        % column... 
        matCurObjDiskPixels(:,1) = matCurObjDiskPixels(:,1) + round(matCurObjCentroids(1,2)); % again, check X & Y vs ROWS & COLUMNS... :)
        matCurObjDiskPixels(:,2) = matCurObjDiskPixels(:,2) + round(matCurObjCentroids(1,1));

        % remove pixel values outside of image dimensions

        matBadPixels = any(matCurObjDiskPixels<1,2) | matCurObjDiskPixels(:,2)>matImSize(2) | matCurObjDiskPixels(:,1)>matImSize(1);
        matCurObjDiskPixels(matBadPixels,:) = [];

        matLinearPixelIX = sub2ind2(matImSize,matCurObjDiskPixels);
        % measure stuff
        matMeasurement(i,j) = mean(OrigImage(matLinearPixelIX));

        % keep track for visualization
        matMeasurementCounterImage(matLinearPixelIX) = matMeasurementCounterImage(matLinearPixelIX) + 1;
    end
    
    % show figure
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);

        plotsize = (nEnvironments+1)/2;
        subplot(plotsize,2,j);
        % merge pictures
        matRGB = cat(3,OrigImage,matMeasurementCounterImage,SecondImage);
        % rescale each channel to go from 0 to 1
        for i = 1:3
            matRGB(:,:,i) = matRGB(:,:,i) - min(lin(matRGB(:,:,i)));
            matRGB(:,:,i) = (matRGB(:,:,i) / quantile(lin(matRGB(:,:,i)),.99));
        end
        matRGB(matRGB<0)=0;
        matRGB(matRGB>1)=1;
        %figure;
        CPimagesc(matRGB,handles)
        title(sprintf('overlayed with Environment of size %d',EnvironmentRadii(j)));
        drawnow
     end
end

% store measurements

handles.Measurements.(ObjectName).EnvironmentFeatures = EnvironmentFeatures;
handles.Measurements.(ObjectName).(sprintf('Environment_%s',ImageName)){k} = matMeasurement;


