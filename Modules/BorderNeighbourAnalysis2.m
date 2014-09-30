function handles = BorderNeighbourAnalysis2(handles)

% Help for the Neighbours module:
%
% Category: Measurement
%
% SHORT DESCRIPTION:
% Analyses the neighbours of objects in a segmented grayscale image.
% Neighbours are defined by presence of at least one pixel within an
% extension surrounding object of interest. Returns
%   - basic features of the neigbhourhood to objects and background
%   - objectIDs of objects in neighbourhood
% Note: This module has several alterations in input-parameters, output,
% image analysis and robustness, when compared to James initial module.


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%


drawnow

[~, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to analyse the borders of?
%infotypeVAR01 = objectgroup
%defaultVAR01 = Cells
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = How many pixels do you maximally expand away from object?
%defaultVAR02 = 2
%infotypeVAR02 = objectgroup indep
distanceToObjectMax = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = How many pixels do you minimally expand away from object?
%defaultVAR03 = 2
%infotypeVAR03 = objectgroup indep
distanceToObjectMin = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%%VariableRevisionNumber = 5


%% Checking Inputs
distanceToObjectMax=str2double(distanceToObjectMax);
distanceToObjectMin=str2double(distanceToObjectMin);

if isnan(distanceToObjectMax)
    distanceToObjectMax=2;
    fprintf('Maximal extension must be numerical value. Used default value of 2.');
end
if isnan(distanceToObjectMin)
    distanceToObjectMin=2;
    fprintf('Minimal extension must be numerical value. Used default value of 2.');
end

if distanceToObjectMin < 1
    distanceToObjectMin = 2;
    fprintf('Minimal extension must be positive value. Used default value of 2.');
end

if distanceToObjectMax < 1
    distanceToObjectMax = 2;
    fprintf('Maxmial extension must be positive value. Used default value of 2.');
end

if distanceToObjectMin > distanceToObjectMax
    distanceToObjectMax = distanceToObjectMin;
    fprintf(['Minimal extension must not be larger than maximal extension. Set both values to larger value:' num2str(distanceToObjectMax)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing of the Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBTAIN SEGMENTATION IMAGE
SegmentedObjectImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'MustBeGray','DontCheckScale');
loadedImage=SegmentedObjectImage;

% IDENTIFY IDs OF OBJECT OF INTEREST
% note: As opposed to James' version IDs are not derived from segmentation
% file, but instead from object count

% Obtain (first) Index, where Object name is found for Object Count (note:
% look for name because size of object count influenced by different
% classes of primary objects, such as spots and nuclei)
IXObject = find(cell2mat(cellfun(@(x) strcmp(x, ObjectName), handles.Measurements.Image.ObjectCountFeatures, 'uniformoutput',false)),1,'first');
ObjCounts = handles.Measurements.Image.ObjectCount{1,handles.Current.SetBeingAnalyzed}(1,IXObject); 

if ObjCounts > 0
    ImageEmpty = false;
    listUniq = 1:ObjCounts;
else
    ImageEmpty = true;
end


% ENSURE BACKGROUND VALUE OF 0
% note: added by T since James mentioned images with non-0-background
if ImageEmpty == false
    listUniqInImage = unique(loadedImage);
    backgroundID = listUniqInImage(~ismember(listUniqInImage',listUniq)');
    isBackground=false(size(loadedImage));
    for k=1:length(backgroundID)
        isBackground(loadedImage==backgroundID(k))=true;
    end
    loadedImage(isBackground)=0;
end



% DEFINE FEATURES
% In cell (where varying amount of measurements per object)
FeaObjIDs={'Neighbour_IDs'};            % excludes background: IDs of adjacent objects
FeaObjAbs={'Overlap_objects_pixels'};   % Overlap with objects, in structure
FeaObjFra={'Overlap_objects_fraction'}; % Overlap with objects, in structure
FeaRelations=[FeaObjIDs FeaObjAbs FeaObjFra];  %position within each cell identifies same adjacent object
% In different columns of one matrix (where fixed number of
% measurements per object)
FeaBackNum={'Border_to_Background'}; % says, whether indvidual object is next to background
FeaBackAbs={'Overlap_background_pixels'}; % absolute pixel amount next to background
FeaBackFra={'Overlap_background_fraction'}; % absolute fraction next to background
FeaLength={'Area_extension_pixels'};  % absolute number of pixels covered by extension
FeaObjNum={'Number_neighbouring_objects'};  % Amount of neigbouring objects (note: this module is indepenent of SVM cleanup)
FeaExtOutImg={'Extension_out_of_image'};
FeaGlobal = [FeaBackNum FeaBackAbs FeaBackFra FeaLength FeaObjNum FeaExtOutImg];

% INITIALIZE FEATURES
if handles.Current.SetBeingAnalyzed==1
    % Labels
    handles.Measurements.(ObjectName).NeigbourGeneralFeatures = FeaGlobal;
    handles.Measurements.(ObjectName).NeigbourRelationsFeatures = FeaRelations;
    % Numerical Measurements
    handles.Measurements.(ObjectName).NeigbourGeneral = cell(1,handles.Current.NumberOfImageSets);
    handles.Measurements.(ObjectName).NeigbourRelations = cell(1,handles.Current.NumberOfImageSets);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OBTAIN MEASUREMENTS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ImageEmpty==1
    % If there is no object, set measurment to empty cell (as recommended
    % by Berend, contrasting James' pecursor version with NaN)
    NeigbourGeneral = [];
    NeigbourRelations = {[] [] []};
else
    % Obtain measuremnts per object; in contrast to James' version:
    % distinguish between background and real objects as different classes
    % of neighbours
    [neighbours, borderLength, totSum, totFra] = arrayfun(@(x) getParameters(x,loadedImage,distanceToObjectMax,distanceToObjectMin), listUniq,'uniformoutput',false);
    
    % Identify background and non backround neigbouring object
    fIXBackground=cellfun(@(x) ismember(x,0)', neighbours,'uniformoutput',false);        % determine indices of Background (should either be first or all zero)
    fIXNoBackground=cellfun(@(x) ~(x), fIXBackground,'uniformoutput',false);            % determine indices of non -background
    
    % Select Measurements
    totalBorder=borderLength';
    numNeighboursBack=cellfun(@sum, fIXBackground)';     % Will return 1 or 0 dependent upon whether cell touches background
    numNeighboursNoBack=cellfun(@sum, fIXNoBackground)';        % will return number of adjacent objects
    
    NeigbourID=(cellfun(@(x,y) (x(y))',neighbours,fIXNoBackground,'uniformoutput',false))';      %Obtain ID from neighbouring objects (filter for non-background)
    
    borderBackAbs=(cellfun(@(x,y) sum((x(y))),totSum,fIXBackground,'uniformoutput',false))';         % sum over individual fraction -> case empty: overlap will be 0 insead of []
    borderNoBackAbs=cellfun(@(x,y) x(y),totSum,fIXNoBackground,'uniformoutput',false)';
    borderBackRel=(cellfun(@(x,y) sum(x(y)),totFra,fIXBackground,'uniformoutput',false))';      % sum over individual fraction -> case empty: overlap will be 0 insead of []
    borderNoBackRel=cellfun(@(x,y) x(y),totFra,fIXNoBackground,'uniformoutput',false)';
    
    % Convert individual measurments to matrix
    if iscell(numNeighboursBack); numNeighboursBack = cell2mat(numNeighboursBack); end;
    if iscell(borderBackAbs); borderBackAbs = cell2mat(borderBackAbs); end;
    if iscell(borderBackRel); borderBackRel = cell2mat(borderBackRel); end;
    if iscell(totalBorder); totalBorder = cell2mat(totalBorder); end
    if iscell(numNeighboursNoBack); numNeighboursNoBack = cell2mat(numNeighboursNoBack); end;
    
    
    % Check, which objects would expand out of the image (equals ID within
    % border close to image with width of defined maximum distance)   
    
    cA = [1 size(loadedImage,1)-(distanceToObjectMax-1) size(loadedImage,2)-(distanceToObjectMax-1)];
    cB = [distanceToObjectMax size(loadedImage,1) size(loadedImage,2)];
    
    ImBorder = unique([reshape(loadedImage(cA(1):cB(1),:),1,[]),...
        reshape(loadedImage(cA(2):cB(2),:),1,[]),...
        reshape(loadedImage(:,cA(1):cB(1)),1,[]),...
        reshape(loadedImage(:,cA(3):cB(3)),1,[])]);
    ImBorder=ImBorder(ImBorder~=0)';         % remove background and transpose
    outImage = zeros(size(listUniq,2),1);
    outImage(ImBorder)=1;        %Set indices, of objects with extenstion out of image to 1
    
    % Concatanate measurements with one value per object of interest (note this contasts James' version, where they were saved individually)
    NeigbourGeneral=[numNeighboursBack borderBackAbs borderBackRel totalBorder numNeighboursNoBack outImage];
    % Concatanate measurements with corresponding number of elements per per object of interest (note this contasts James' version, where they were saved individually)
    NeigbourRelations = [NeigbourID borderNoBackAbs borderNoBackRel];
    
end


    
%%%%%%%%%%%%%%%%%%%%
%%% Save Results %%%
%%%%%%%%%%%%%%%%%%%%

handles.Measurements.(ObjectName).NeigbourGeneral{handles.Current.SetBeingAnalyzed} = NeigbourGeneral;
handles.Measurements.(ObjectName).NeigbourRelations{handles.Current.SetBeingAnalyzed} = NeigbourRelations;

end

function [neighbours, bordLength, totSum, totFra] = getParameters(listUniq,loadedImage,distanceToObjectMax,distanceToObjectMin)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions Defined %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% in contrast to James' version: one single larger function instead of
% separate small functions; introduce userd-defined size of extension &
% filtering done exclusively with logicals

Cell = loadedImage==listUniq;       % filter for object of interest
dilateCellMax= bwmorph(Cell,'dilate',distanceToObjectMax);    % expand from object by given pixel-size
dilForMin=(distanceToObjectMin-1);      %since later the difference of the masks is permitted, expand by one pixel less as than specified minimum for inclusion
dilateCellMin= bwmorph(Cell,'dilate',dilForMin);

maskedCell= dilateCellMax & ~dilateCellMin;     % filter marking expansion from object
bordList=loadedImage(maskedCell);   % Turns pixels within expansion into a linearised matrix
neighbours = unique(bordList);
numN = size(neighbours,1);          % obtain number of neighbours
if isempty(neighbours)              % mark single cells, which would fill the complete image by nan
    neighbours=nan(1,1);
    numN = nan(1,1);
end
bordLength=size(bordList,1);        % Obtain complete number of pixels incorporated by expansion

totFra=nan(1,numN);
totSum=nan(1,numN);

for k=1:numN
    totSum(k)=sum(bordList==neighbours(k));
    totFra(k)=totSum(k)./bordLength;
end


end

