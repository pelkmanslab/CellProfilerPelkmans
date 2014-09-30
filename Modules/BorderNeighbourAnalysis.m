function handles = BorderNeighbourAnalysis(handles)

% Help for the Neighbours module:
%
% Category: Measurement
%
% SHORT DESCRIPTION:
% Analyses the neighbours of objects in a segmented grayscale image.
% Neighbours are defined by presence of at least one pixel within an
% extension surrounding object of interest. Returns
%   - A) basic features of the neigbhourhood to objects and background
%   - B) relational features about objects in neighbourhood
% Note: This module has several alterations in input-parameters, output,
% image analysis and robustness, when compared to James initial module,
% which is attached below as a comment.
%
% A) Outputs are a MATRIX. They can be loaded with the standard
%    getrawprobmodeldata2
%    1] Border_to_Background            1 if cell (partially) borders the background
%    2] Overlap_background_pixels       Absolute pixels bordering background
%    3] Overlap_background_fraction     Fraction of extension bodering background
%    4] Area_extension_pixels           Amount of pixels covered by extension
%    5] Number_neighbouring_objects     Amount of neigbouring objects (note: this module is indepenent of SVM cleanup)
%    6] Extension_out_of_image          1 if extension reaches outside of the image
%
% B) Outputs are a CELL. They require the special getrawprobmodeldata2cell
%    They allow to construct arbitraray (weighted) relational features.
%    1] Neighbour_IDs                   IDs of adjacent objects (excludes background)
%    2] Overlap_objects_pixels          Amount of pixels of a given Neighbour within extension (same order as Neighbour_IDs)
%    3] Overlap_objects_fraction        Fraction of extension with a given Neighbour (same order as Neighbour_IDs)
%
%
% $Revision: 4533 $


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

%%%VariableRevisionNumber = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Checking Inputs
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
FeaExtOutImg={'Extension_out_of_image'};    % says, wheter the extension would reach out of the image
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
    % Get Cropping coordinates for each object   
    props = regionprops(loadedImage,'BoundingBox');
    BoxPerObj = cat(1,props.BoundingBox);

    N = floor(BoxPerObj(:,2)-distanceToObjectMax-1);                    f = N < 1;                      N(f) = 1;  
    S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+distanceToObjectMax+1);   	f = S > size(loadedImage,1);    S(f) = size(loadedImage,1);
    W = floor(BoxPerObj(:,1)-distanceToObjectMax-1);                    f = W < 1;                      W(f) = 1;
    E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+distanceToObjectMax+1);      f = E > size(loadedImage,2);    E(f) = size(loadedImage,2);
    
    % Obtain measuremnts per object;
    
    [neighbours, borderLength, totSum, totFra,qcIsPresent] = ...
        arrayfun(@(a,b,c,d,e) getParameters(a,b,c,d,e,loadedImage,distanceToObjectMax,distanceToObjectMin), listUniq,N',S',W',E','uniformoutput',false);
    
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
    
    % Set data, where nucleating object seems suspicious (eg. hinting at error in secondary object module on brutus) , to nan
    qcIsPresent = cell2mat(qcIsPresent');
    
    NeigbourGeneral(~qcIsPresent,:) = nan(sum(~qcIsPresent),size(NeigbourGeneral,2));
    NeigbourRelations(~qcIsPresent,:) = repmat({nan(1)},[sum(~qcIsPresent) size(NeigbourRelations,2)]);
    
end



%%%%%%%%%%%%%%%%%%%%
%%% Save Results %%%
%%%%%%%%%%%%%%%%%%%%

handles.Measurements.(ObjectName).NeigbourGeneral{handles.Current.SetBeingAnalyzed} = NeigbourGeneral;
handles.Measurements.(ObjectName).NeigbourRelations{handles.Current.SetBeingAnalyzed} = NeigbourRelations;

end

function [neighbours, bordLength, totSum, totFra,qcIsPresent] = getParameters(listUniq,N,S,W,E,loadedImage,distanceToObjectMax,distanceToObjectMin)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions Defined %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% in contrast to James' version: one single larger function instead of
% separate small functions; introduce user-defined size of extension &
% filtering done exclusively with logicals; also work on cropped image to
% reduce computing time 


% Reduce computationally heavy image morphing, by reducing area of
% operation - Crop images

loadedImage = loadedImage(N:S,W:E);


Cell = loadedImage==listUniq;       % filter for object of interest

if any(Cell(:))
    % Initialize quality control: in the map object distances module it has
% occured that matlab / CP on brutus would allocate cytoplasm to the wrong
% cells / fuse it with others. This error would not be reproducible
% locally. Thus include safety check for presence of nucleating object;
% Strategy is to temporarily create harmless measuremeents and use
% qcIsPresent towards end of calling module to set output to nan.
    qcIsPresent = true;    % ground assumption is that object is present
else
    qcIsPresent = false;
end

switch qcIsPresent
    case true     % If cell is really present
        dilateCellMax= bwmorph(Cell,'dilate',distanceToObjectMax);    % expand from object by given pixel-size
        dilForMin=(distanceToObjectMin-1);      %since later the difference of the masks is permitted, expand by one pixel less as than specified minimum for inclusion
        dilateCellMin= bwmorph(Cell,'dilate',dilForMin);
        
        maskedCell= dilateCellMax;
        maskedCell(dilateCellMin) = false;     % filter marking expansion from object
        bordList=loadedImage(maskedCell);   % Turns pixels within expansion into a linearised matrix
        neighbours = unique(bordList);
        numN = size(neighbours,1);          % obtain number of neighbours
        if isempty(neighbours)              % mark single cells, which would fill the complete image by nan
            neighbours=nan(1,1);
            numN = nan(1,1);
        end
        bordLength=size(bordList,1);        % Obtain complete number of pixels within expansion
        
        totFra=nan(1,numN);
        totSum=nan(1,numN);
        
        for k=1:numN
            totSum(k)=sum(bordList==neighbours(k));
            totFra(k)=totSum(k)./bordLength;
        end
        
    case false   % If cell is not present create harmless outputs.
        neighbours=0;
        bordLength=1;
        totSum=1; 
        totFra=1;       
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  JAMES VERSION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [TS] The code below is the original code for the neigbhour analysis as
% written by [JP]. It was completely rewritten to increase speed,
% robustness and amount of measuremetns.
%
% Ultimately James module was only used in the beginning to pioneer
% neighbour analysis.
%
% To prevent confusion among multiple modules, it was replaced by the
% module formerly known as BorderNeigbhourAnalysis2 on 2012/07/12.


% function handles = BorderNeighbourAnalysis(handles)
%
% % Help for the Neighbours module:
% %
% % Category: Measurement
% %
% % SHORT DESCRIPTION:
% % Analyses the neighbours of cells in a segmented grayscale image. Produces the names
% % (i.e. intensities) of cells neighbouring a cell X, the number of cells neighbouring X,
% % the length in pixels occupied by a particular cell on cell X's border and the fraction
% % of X's total border which this length takes up.
% % Output into 4 cells, first column of each cell is for the background's border analysis.
%
% %%%%%%%%%%%%%%%%%
% %%% VARIABLES %%%
% %%%%%%%%%%%%%%%%%
%
% drawnow
%
% [CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);
%
% %textVAR01 = What do you want to analyse the borders of?
% %infotypeVAR01 = objectgroup
% %defaultVAR01 = Cells
% NeighboursName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
% %inputtypeVAR01 = popupmenu
%
%
% %%VariableRevisionNumber = 5
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Processing of the Images %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %Retrieve images of segmented cells specified by user:
% SegmentedObjectImage = CPretrieveimage(handles,['Segmented', NeighboursName],ModuleName,'MustBeGray','DontCheckScale');
% loadedImage=SegmentedObjectImage;
% %find number of cells:
% %1. List unique elements in loadedImage
% listUniq = unique(loadedImage);
%
% %Define Features:
% ObjectIDFeatures={'Object_IDs'};
% neighboursFeatures={'Neighbours_of_Cell'};
% numNeighboursFeatures={'Number_of_Neighbours_of_Cell'};
% totSumFeatures={'Absolute_Pixel_Length_Border_Covered'};
% totFracFeatures={'Fraction_Length_Border_Covered'};
%
% %ensure there are cells within the image:
% if sum(listUniq)==0;
%     if handles.Current.SetBeingAnalyzed==1;
%         handles.Measurements.(NeighboursName).ObjectID = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).ObjectIDFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).BorderNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).BorderNeighbourFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).NumBorderNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).NumBorderNeighbourFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).SumNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).SumNeighbourFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).FracNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).FracNeighbourFeatures = cell(1,1);
%     end;
%
%     handles.Measurements.(NeighboursName).ObjectID{handles.Current.SetBeingAnalyzed} = NaN;
%     handles.Measurements.(NeighboursName).ObjectIDFeatures = ObjectIDFeatures;
%
%     handles.Measurements.(NeighboursName).BorderNeighbour{handles.Current.SetBeingAnalyzed} = NaN;
%     handles.Measurements.(NeighboursName).BorderNeighbourFeatures = neighboursFeatures;
%
%     handles.Measurements.(NeighboursName).NumBorderNeighbour{handles.Current.SetBeingAnalyzed} = NaN;
%     handles.Measurements.(NeighboursName).NumBorderNeighbourFeatures = numNeighboursFeatures;
%
%     handles.Measurements.(NeighboursName).SumNeighbour{handles.Current.SetBeingAnalyzed} = NaN;
%     handles.Measurements.(NeighboursName).SumNeighbourFeatures = totSumFeatures;
%
%     handles.Measurements.(NeighboursName).FracNeighbour{handles.Current.SetBeingAnalyzed} = NaN;
%     handles.Measurements.(NeighboursName).FracNeighbourFeatures = totFracFeatures;
% else
%
%     %Calculate Variables
%     [ObjectID,uniqCell] = arrayfun(@(x) cellUniq(loadedImage,x),listUniq,'uniformoutput',false);%output logical cell regions
%
%     [neighbours,maskedCell] = cellfun(@(x) bordercells(loadedImage,x),uniqCell,'uniformoutput',false);%output which cells are neighbours of a cell, and a logical cell border
%
%     numNeighbours = cellfun(@(x) numbordercell(x),neighbours,'uniformoutput',false);%output the number of neighbours a particular cell has
%
%     [totSum, totFrac] = cellfun(@(x) valEdge(loadedImage,x),maskedCell,'uniformoutput',false);%save fraction of cell bordered and absolute part of cell bordered by its other cells
%
%     ObjectID=cell2mat(ObjectID);
%     numNeighbours=cell2mat(numNeighbours);
%
%     %%%%%%%%%%%%%%%%%%%%
%     %%% Save Results %%%
%     %%%%%%%%%%%%%%%%%%%%
%
%     % let's initialize if it doesnt exist
%     if handles.Current.SetBeingAnalyzed==1;
%         handles.Measurements.(NeighboursName).ObjectID = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).ObjectIDFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).BorderNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).BorderNeighbourFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).NumBorderNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).NumBorderNeighbourFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).SumNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).SumNeighbourFeatures = cell(1,1);
%         handles.Measurements.(NeighboursName).FracNeighbour = cell(1,handles.Current.NumberOfImageSets);
%         handles.Measurements.(NeighboursName).FracNeighbourFeatures = cell(1,1);
%     end
%
%     handles.Measurements.(NeighboursName).ObjectID{handles.Current.SetBeingAnalyzed} = ObjectID;
%     handles.Measurements.(NeighboursName).ObjectIDFeatures = ObjectIDFeatures;
%
%     handles.Measurements.(NeighboursName).BorderNeighbour{handles.Current.SetBeingAnalyzed} = neighbours;
%     handles.Measurements.(NeighboursName).BorderNeighbourFeatures = neighboursFeatures;
%
%     handles.Measurements.(NeighboursName).NumBorderNeighbour{handles.Current.SetBeingAnalyzed} = numNeighbours;
%     handles.Measurements.(NeighboursName).NumBorderNeighbourFeatures = numNeighboursFeatures;
%
%     handles.Measurements.(NeighboursName).SumNeighbour{handles.Current.SetBeingAnalyzed} = totSum;
%     handles.Measurements.(NeighboursName).SumNeighbourFeatures = totSumFeatures;
%
%     handles.Measurements.(NeighboursName).FracNeighbour{handles.Current.SetBeingAnalyzed} = totFrac;
%     handles.Measurements.(NeighboursName).FracNeighbourFeatures = totFracFeatures;
%
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Fucntions Defined %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%
%
% %Generate image of only an individual cell, represented by 1s
% function [ObjectID,uniqCell] = cellUniq(loadedImage,listUniq);
% uniqCell = loadedImage==listUniq;
% ObjectID=listUniq;
%
% %Determine which cells lie on border of some Cell X
% function [neighbours,maskedCell] = bordercells(loadedImage,matObjectPosition);
% dilateCell= bwmorph(matObjectPosition,'dilate',1);
% maskedCell= dilateCell-matObjectPosition;
% neighbours = unique(loadedImage(maskedCell==1));%find the unique cells in the masked border
%
% %Determine number of neighbours of some cell X
% function [numN] = numbordercell(neighbours);
% [numN,numNjunk]=size(neighbours);
%
% %Determine absolute and fraction of the border of some cell X taken up by a
% %cell Y
% function [totSum, totFrac] = valEdge(loadedImage,maskedCell);
% maskedCell=logical(maskedCell);
% bordList=loadedImage(maskedCell);%Turns cells making up border into a linearised matrix
% bordLength=size(bordList,1);%Length of border
% neighInt=unique(bordList);%Cells making up border
% for j=1:size(neighInt,1);
%     sumOG=size(find(bordList==neighInt(j)),1);%Find length of border take up by one of the cells neighbours
%     totSum(1,j)=sumOG;
%     totFrac(1,j)=sumOG/bordLength;
% end
%
% end
%




