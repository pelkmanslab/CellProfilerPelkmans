function handles = SeparatePrimary(handles)

% Help for the Separate Primary module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Separate large obeject identified by the Identify Prim Automatic module
% *************************************************************************
% How it works>
% 0) Insert description here
%
% $Revision: 1808 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the primary objects you want to separate?
%infotypeVAR01 = objectgroup
PrimaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = SepNuclei
%infotypeVAR02 = objectgroup indep
SeparatedObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the images to be used to find the separation edges of the objects?
%infotypeVAR03 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Size threshold in pixels to try to separate objects. Objects smaller than the threshodl will not be considered for separation.
%defaultVAR04 = 1500
ThresholdCorrection = str2num(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Threshold for the angle acutness factor of the poits at which the objects is separated [0 Inf]. Point with a higher angle acutness will not be considered for separation.
%defaultVAR05 = 0.4
ThreAngle = str2num(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Strigtness factor of the separation line [1 Inf]. If objects are splited by a line with higher threshold than the specified the object will not be separated.
%defaultVAR06 = 4
ThreStraigth = str2num(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Separation intesity threshold. Fraction of the initial object intesiy that the final objects must at least acount for
%defaultVAR07 = 0.2
ThreIntensisty = str2num(handles.Settings.VariableValues{CurrentModuleNum,7});

%%%VariableRevisionNumber = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
PrimaryLabelImage = CPretrieveimage(handles,['Segmented', PrimaryObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));
SeparatedPrimaryLabelImage=PrimaryLabelImage>0;

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%filter the original image
%f=fspecial('gaussian',[5,5],1);
%IntensityImage=imfilter(OrigImage,f);
IntensityImage=OrigImage;
%figure;imagesc(IntensityImage)

%get the water sed transform
L = watershed(imcomplement(IntensityImage));
L(~PrimaryLabelImage)=0;
%figure;imagesc(L)

%get the Id of the primary objects
ObjectsId=unique(PrimaryLabelImage);
ObjectsId=ObjectsId(2:end);

%get the size of the objects
ObjectSize=regionprops(PrimaryLabelImage,'area');
ObjectSize=cat(1,ObjectSize(:).Area);

%Obects to separate
BigObjects=find(ObjectSize>ThresholdCorrection);

%Fraction of image covered with primary objects
fractionCoveredWithObjects = sum(SeparatedPrimaryLabelImage(:))./numel(SeparatedPrimaryLabelImage(:));

%Only process Image, if less than 90 percent of the image consists of
%primary objects. If the majority of the image would consist of false
%positive nuclei, this leads to a memory requirement of more than 30GB
%resulting in timeouts on Brutus

fractionLimit = 0.90;
if fractionCoveredWithObjects <= fractionLimit
    
    %loop through objects to be separated
    for iObj=BigObjects'
        %try % In case of any Problem, just do not do any separation
        
        %get the current object segmanetation
        ImCurrentObject=PrimaryLabelImage==iObj;
        [Ix,Jx]=find(ImCurrentObject);
        
        %get the current object watershed, intensity and segmentation image
        CurrentIntensity=IntensityImage(min(Ix):max(Ix),min(Jx):max(Jx));
        CurrentSegmentation=ImCurrentObject(min(Ix):max(Ix),min(Jx):max(Jx));
        CurrentWatershed=L(min(Ix):max(Ix),min(Jx):max(Jx));
        %figure;imagesc(CurrentSegmentation)
        
        %get the lines
        CurrentPreLines=zeros(size(CurrentIntensity));
        CurrentPreLines(~CurrentWatershed)=CurrentSegmentation(~CurrentWatershed);
        CurrentPreLines(~CurrentSegmentation)=0;
        
        %define lines and crossing points
        f=[0 1 0; 1 0 1; 0 1 0;];
        
        CurrentPreLines2=CurrentPreLines;
        CurrentPreLines2(~CurrentSegmentation)=5;
        CurrentLinesAndNodes=imfilter(CurrentPreLines2,f);
        CurrentLinesAndNodes(~CurrentPreLines2)=0;
        CurrentLinesAndNodes(~CurrentSegmentation)=0;
        %figure;imagesc(CurrentLinesAndNodes)
        
        
        %define lines and measure the area
        CurrentLines=bwlabel(CurrentLinesAndNodes<3 & CurrentLinesAndNodes>0,4);
        %figure;imagesc(CurrentLines)
        
        LineAreas=regionprops(CurrentLines,'area');
        LineAreas=cat(1,LineAreas(:).Area);
        LineIds=unique(CurrentLines(:));
        LineIds(LineIds==0)=[];
        %LineMeanInt=arrayfun(@(a,b) sum(sum(IntTemp(matLines==a)))/b,LineIds,LineAreas);
        
        %define nodes and measure their centroids
        CurrentNodes=bwlabel(CurrentLinesAndNodes>2,4);
        NodesCentroids=regionprops(CurrentNodes,'centroid');
        NodesCentroids=cat(1,NodesCentroids(:).Centroid);
        NodesIds=unique(CurrentNodes(:));
        NodesIds=NodesIds(2:end)';
        
        
        %build the connection matrix used to measure paths
        f={[0 1 0; 0 0 0; 0 0 0;] [0 0 0; 1 0 0; 0 0 0;] [0 0 0; 0 0 1; 0 0 0;] [0 0 0; 0 0 0; 0 1 0;]};
        DisplacedLines=cellfun(@(x) imfilter(CurrentLines,x),f,'uniformoutput',false);
        DisplacedLines=cat(3,DisplacedLines{:});
        
        
        NodeType=zeros(size(NodesIds));
        matNodesLines=zeros(length(NodesIds),length(LineAreas));
        for iNode=1:length(NodesIds)
            tmpid=NodesIds(iNode);
            
            tmpix=CurrentNodes==tmpid;
            temptype=unique(CurrentLinesAndNodes(tmpix));
            NodeType(iNode)=max(temptype)>5;
            
            tmpix=repmat(tmpix(:),4,1);
            templineids=unique(DisplacedLines(tmpix));
            templineids(templineids==0)=[];
            matNodesLines(iNode,templineids')=templineids';
        end
        
        
        matNodesNodes=zeros(length(NodesIds));
        matNodesNodesLabel=zeros(length(NodesIds));
        for iNode=1:length(NodesIds)
            tmplines=unique(matNodesLines(iNode,:));
            tmplines(tmplines==0)=[];
            for l=tmplines
                tmpnodes=find(matNodesLines(:,l));
                matNodesNodes(iNode,tmpnodes)=LineAreas(l);
                %matNodesNodesInt=LineMeanInt(l);
                matNodesNodesLabel(iNode,tmpnodes)=l;
            end
            
        end
        matNodesNodes(sub2ind([length(NodesIds) length(NodesIds)],NodesIds,NodesIds))=0;
        matNodesNodesLabel(sub2ind([length(NodesIds) length(NodesIds)],NodesIds,NodesIds))=0;
        
        %define the border nodes, Sourse Nodes and Target Nodes
        NodeToTest=NodesIds(logical(NodeType));
        NodeixS=repmat(NodeToTest,length(NodeToTest),1);
        NodeixT=repmat(NodeToTest,1,length(NodeToTest));
        NodeixS=NodeixS(:);
        NodeixT=NodeixT(:);
        
        %Measure the direct distance
        directdist=pdist2(NodesCentroids(NodeToTest,:),NodesCentroids(NodeToTest,:));
        directdist=directdist(:);
        
        %get the shortest paths
        %[dist, path]=arrayfun(@(x) function(sparse(matNodesNodes),x,NodeToTest),NodeToTest,'uniformoutput',false);->graphshor_testpath
        %dist=cat(2,dist{:});
        
        if isempty(NodeToTest) %[GG] checks whether the objects as a watershed line, if not it exits the loop
            continue
        end

        [dist,path] = dijkstraCP(matNodesNodes>0,matNodesNodes,NodeToTest,NodeToTest);
        dist=dist(:)';
        
        %measure the angle factor for the nodes
        AngleFactor=[];
        for i =1:length(NodeToTest)
            tmpcentroid=round(NodesCentroids(NodeToTest(i),:));
            tmpi=tmpcentroid(2)-5:tmpcentroid(2)+5;
            tmpj=tmpcentroid(1)-5:tmpcentroid(1)+5;
            tmpim=padarray(CurrentSegmentation,[10 10]);
            tmpim=tmpim(tmpi+10,tmpj+10);
            cyto=sum(tmpim(:));
            backgr=(length(tmpim(:))-cyto)/cyto;
            AngleFactor(i)=backgr;
        end
        
        QuantileDistance=quantile(dist(dist~=0 & dist~=Inf),1);
        %ThreAngle
        %ThreStraigth
        %ThreIntensisty
        
        %filter acording to sise of distance threshold
        thrix=dist~=0 & dist<QuantileDistance;
        
        NodeixS2=NodeixS(thrix);
        NodeixT2=NodeixT(thrix);
        
        SelectedLines={};
        InitialInt=sum(CurrentIntensity(CurrentSegmentation(:)).^2);
        for i=1:length(NodeixT2)
            
            tmpimage=zeros(size(CurrentNodes));
            
            % find the path
            tmppath=path{NodeToTest==NodeixS2(i),NodeToTest==NodeixT2(i)};
            tmpimage(ismember(CurrentNodes(:),tmppath(:)))=1;
            %[AngleFactor(NodeToTest==NodeixS2(i)) AngleFactor(NodeToTest==NodeixT2(i))]
            
            
            if AngleFactor(NodeToTest==NodeixS2(i))<ThreAngle & AngleFactor(NodeToTest==NodeixT2(i))<ThreAngle
                % builthe path image
                for j=1:length(tmppath)-1
                    tmpimage(find(CurrentLines==matNodesNodesLabel(tmppath(j),tmppath(j+1))))=1;
                end
                
                %get the straigt line: this can b put in a subfunction perhaps
                tmpcentroid1=round(NodesCentroids(NodeToTest(NodeToTest==NodeixS2(i)),:));
                tmpcentroid2=round(NodesCentroids(NodeToTest(NodeToTest==NodeixT2(i)),:));
                m=(tmpcentroid1(2)-tmpcentroid2(2))/(tmpcentroid1(1)-tmpcentroid2(1));
                
                if m ~= -Inf & m ~= Inf & ~isnan(m)
                    
                    x=(min([tmpcentroid1(1),tmpcentroid2(1)]):max([ tmpcentroid1(1),tmpcentroid2(1)]));
                    y=(min([tmpcentroid1(2),tmpcentroid2(2)]):max([ tmpcentroid1(2),tmpcentroid2(2)]));
                    c=tmpcentroid1(2)-m*tmpcentroid1(1);
                    py=round(m.*x+c);
                    px=round((y-c)/m);
                    StraigtLineix=sub2ind(size(tmpimage),[y py],[px x]);
                    
                    %get the straigness ratio
                    tmprim=tmpimage;
                    tmprim(StraigtLineix(~isnan(StraigtLineix)))=1;
                    tmprim=imfill(tmprim);
                    tmpratio=sum(tmprim(:))/length(unique(StraigtLineix));
                    
                    tmpnewseg=CurrentSegmentation;
                    tmpnewseg(logical(tmpimage))=0;
                    tmpnewseg=bwlabel(tmpnewseg,4);
                    
                    k=unique(tmpnewseg(:));
                    k(k==0)=[];
                    
                    tmpobject=arrayfun(@(x) tmpnewseg==x,k,'uniformoutput',false);
                    tmpint=cellfun(@(x) sum(CurrentIntensity(x(:)).^2),tmpobject,'uniformoutput',false);
                    tmpint=cat(1,tmpint{:});
                    
                    tmpint=sum(tmpint>InitialInt*ThreIntensisty)==length(tmpint);
                    
                    if tmpratio<ThreStraigth & tmpint
                        SelectedLines{end+1}=tmpimage;
                        %CurrentSegmentation=tmpnewseg;
                        %CurrentSegmentation=CurrentSegmentation>0;
                    end
                    
                end
            end
            
        end
        
        %filter selected lines
        
        ToBeDrawn=ones(1,length(SelectedLines));
        if ~isempty(SelectedLines)
            
            for iLine=1:length(SelectedLines)
                tmpline=SelectedLines{iLine};
                cellcomparison=cellfun(@(x) any((x(:)+tmpline(:))==2),SelectedLines);
                for iLine2=1:length(SelectedLines);
                    if cellcomparison(iLine2)
                        if sum(SelectedLines{iLine}(:))>sum(SelectedLines{iLine2}(:))
                            ToBeDrawn(iLine)=0;
                        end
                    end
                end
            end
            
            % CREATE SEGMENTATION, WHICH HAS BEEN CUT
            % make prospective Segmentation, which will be passed to additional
            % checks
            prospectiveSegmentation=CurrentSegmentation;
            for iLine=1:length(SelectedLines)
                if ToBeDrawn(iLine)
                    prospectiveSegmentation(logical(SelectedLines{iLine}))=0;
                end
            end
            
            % Check if new subobject resulting from two cuts, is too small to
            % qualify as individual object. If not, try to make educated guess,
            % which of the possible cuts should be done
            try   % just a measure of precaution
                if sum(ToBeDrawn) == 4       % only trigger if 2 pairs of line
                    % get areas of newly generated objects
                    SubSegmentation = bwconncomp(prospectiveSegmentation>0);
                    SubAreas = cell2mat(cellfun(@numel, SubSegmentation.PixelIdxList,'UniformOutput',false));
                    MinFraction = 0.2;          % Fraction, of area of parent, which should be within an child; Consider to be selectable as input
                    TooSmallArea = SubAreas < (sum(SubAreas).*MinFraction);
                    if any(TooSmallArea)      % Only correct for cuts, if any object is too small.
                        % Preferentially use smallest line  to cut
                        LineIdx = find(ToBeDrawn);
                        LineLengths = NaN(size(LineIdx));
                        for k = 1:4 % note that only triggered for 4
                            LineLengths(k) = sum(SelectedLines{LineIdx(k)}(:));
                        end
                        [LineMin, IdMin] = min(LineLengths);
                        if sum(LineLengths == LineMin) == 2;     % if only one pair has minimal length
                            LineIdx(IdMin) = [];                % set every candidate lines exept the used one to not be drawn
                            ToBeDrawn(LineIdx) = 0;
                        else
                            % Or the dimmest line of two equally long cuts
                            MeanIntensityPerLine = NaN(size(LineIdx));
                            for k=1:4
                                f = logical(SelectedLines{LineIdx(k)}(:));
                                MeanIntensityPerLine(k) = sum(CurrentIntensity(f));  % not that this corresponds to the absolute and mean intensity if the lines have same length
                            end
                            [~, IdMin] = min(MeanIntensityPerLine);
                            % make cut at minimum or in case of multiple minima at
                            % first minimum (which is the 2nd output of min)
                            LineIdx(IdMin) = [];       % set every candidate lines exept the used one to not be drawn
                            ToBeDrawn(LineIdx) = 0;
                        end
                        % Make NewObjectSegmenation, in case of too small object
                        NewObjectSegmentation=CurrentSegmentation;
                        for iLine=1:length(SelectedLines)
                            if ToBeDrawn(iLine)
                                NewObjectSegmentation(logical(SelectedLines{iLine}))=0;
                            end
                        end
                    else
                        % Make NewObjectSegmenation, in case of NO too small object
                        NewObjectSegmentation = prospectiveSegmentation;
                    end
                else
                    % Make NewObjectSegmenation, in case of not 4 lines
                    NewObjectSegmentation = prospectiveSegmentation;
                end
            catch someCutDidNotWork
                % In case of unforseen error, catch by doing segmenation with all
                % cuts instead of provoking a crash
                NewObjectSegmentation = prospectiveSegmentation;
                fprintf('A cut did not work \n')
            end
            
            if length(unique(bwlabel(NewObjectSegmentation,4)))>2
                %save the separated object
                ImCurrentObject2=ImCurrentObject;
                ImCurrentObject2(min(Ix):max(Ix),min(Jx):max(Jx))=NewObjectSegmentation;
                ImFinalFilter=xor(ImCurrentObject2,ImCurrentObject);
                SeparatedPrimaryLabelImage(ImFinalFilter)=0;
            end
        end
        
        %  catch CouldNotPerformSeparation % catch any error by skipping separation (might not best performance, but no error, at least as good as before)
        %      fprintf(['Error in SeparatePrimay for object ' num2str(iObj) '. Skipping separation of this object. \n']);
        %  end
    end
    
else
    fprintf(['Did not perform seapration of Primary objects as more than  ' num2str(fractionLimit) ' of the image consists of objects to prevent out of memory \n']);
end

%FinalLabelMatrixImage=bwlabel(SeparatedPrimaryLabelImage,4);

%%% Calculates the 2 sobel filters.  The sobel filter is directional, so it
%%% is used in both the horizontal & vertical directions and then the
%%% results are combined.
filter1 = [1 1 1;0 0 0;-1 -1 -1];%fspecial('sobel');
filter2 = filter1';
%%% Applies each of the sobel filters to the original image.
I1 = imfilter(SeparatedPrimaryLabelImage, filter1);
I2 = imfilter(SeparatedPrimaryLabelImage, filter2);
%%% Adds the two images.
%%% The Sobel operator results in negative values, so the absolute values
%%% are calculated to prevent errors in future steps.
AbsSobeledImage = (abs(I1) + abs(I2))>0;
clear I1; clear I2;                  %%% [NB] hack. save memory
SeparatedPrimaryLabelImage(AbsSobeledImage)=0;

%label the final image
FinalLabelMatrixImage=bwlabel(SeparatedPrimaryLabelImage,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final, segmented label matrix image of secondary objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented',SeparatedObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

fieldname = ['UneditedSegmented',SeparatedObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

fieldname = ['SmallRemovedSegmented',SeparatedObjectName];
handles.Pipeline.(fieldname) = FinalLabelMatrixImage;

%%% Saves the ObjectCount, i.e. the number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = strcmpi(handles.Measurements.Image.ObjectCountFeatures,SeparatedObjectName);%[GG] changed strfind with strcmpi
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SeparatedObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(FinalLabelMatrixImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(SeparatedObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(FinalLabelMatrixImage,'Centroid');
%%% [NB] hack. save memory.
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)
    Centroid = [0 0];
end
handles.Measurements.(SeparatedObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

if any(findobj == ThisModuleFigureNumber)
    
    %%% Calculates the ColoredLabelMatrixImage for displaying in the figure
    %%% window in subplot(2,2,2).
    ColoredLabelMatrixImage = CPlabel2rgb(handles,FinalLabelMatrixImage);
    %%% Calculates OutlinesOnOrigImage for displaying in the figure
    %%% window in subplot(2,2,3).
    %%% Note: these outlines are not perfectly accurate; for some reason it
    %%% produces more objects than in the original image.  But it is OK for
    %%% display purposes.
    %%% Maximum filters the image with a 3x3 neighborhood.
    MaxFilteredImage = ordfilt2(FinalLabelMatrixImage,9,ones(3,3),'symmetric');
    %%% Determines the outlines.
    IntensityOutlines = FinalLabelMatrixImage - MaxFilteredImage;
    %%% [NB] hack.s ave memory.
    clear MaxFilteredImage;
    %%% Converts to logical.
    warning off MATLAB:conversionToLogical
    LogicalOutlines = logical(IntensityOutlines);
    %%% [NB] hack.s ave memory.
    clear IntensityOutlines;
    warning on MATLAB:conversionToLogical
    %%% Determines the grayscale intensity to use for the cell outlines.
    %[NB-HACK] so that images are not so dim!!!!
    ObjectOutlinesOnOrigImage = OrigImage;
    ObjectOutlinesOnOrigImage=ObjectOutlinesOnOrigImage-quantile(OrigImage(:), 0.025);
    ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage<0)=0;
    ObjectOutlinesOnOrigImage(ObjectOutlinesOnOrigImage>quantile(ObjectOutlinesOnOrigImage(:), 0.95))=quantile(ObjectOutlinesOnOrigImage(:), 0.95);
    LineIntensity = quantile(ObjectOutlinesOnOrigImage(:), 0.99);
    
    %%% Overlays the outlines on the original image.
    %ObjectOutlinesOnOrigImage = OrigImage;
    
    ObjectOutlinesOnOrigImage(LogicalOutlines) = LineIntensity;
    %%% Calculates BothOutlinesOnOrigImage for displaying in the figure
    %%% window in subplot(2,2,4).
    %%% Creates the structuring element that will be used for dilation.
    StructuringElement = strel('square',3);
    %%% Dilates the Primary Binary Image by one pixel (8 neighborhood).
    DilatedPrimaryBinaryImage = imdilate(FinalLabelMatrixImage, StructuringElement);
    %%% Subtracts the PrelimPrimaryBinaryImage from the DilatedPrimaryBinaryImage,
    %%% which leaves the PrimaryObjectOutlines.
    PrimaryObjectOutlines = DilatedPrimaryBinaryImage - FinalLabelMatrixImage;
    %%% [NB] hack. save memory.
    clear DilatedPrimaryBinaryImage
    BothOutlinesOnOrigImage = ObjectOutlinesOnOrigImage;
    BothOutlinesOnOrigImage(PrimaryObjectOutlines == 1) = LineIntensity;
    %%% [NB] hack. save memory.
    clear PrimaryObjectOutlines LineIntensity;
    
    %%%%%%%%%%%%%%%%%%%%%%%% END OF INITIATION OF VISUALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%% (rearrangement)
    
    
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    ObjectCoverage = 100*sum(sum(FinalLabelMatrixImage > 0))/numel(FinalLabelMatrixImage);
    
    
    %%% A subplot of the figure window is set to display the original image.
    subplot(2,2,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the colored label
    %%% matrix image.
    subplot(2,2,2);
    CPimagesc(ColoredLabelMatrixImage,handles);
    clear ColoredLabelMatrixImage
    title(['Outlined ',SeparatedObjectName]);
    %%% A subplot of the figure window is set to display the original image
    %%% with secondary object outlines drawn on top.
    subplot(2,2,3);
    CPimagesc(ObjectOutlinesOnOrigImage,handles);
    clear ObjectOutlinesOnOrigImage
    title([SeparatedObjectName, ' Outlines on Input Image']);
    %%% A subplot of the figure window is set to display the original
    %%% image with outlines drawn for both the primary and secondary
    %%% objects.
    subplot(2,2,4);
    CPimagesc(BothOutlinesOnOrigImage,handles);
    clear BothOutlinesOnOrigImage;
    title(['Outlines of ', PrimaryObjectName, ' and ', SeparatedObjectName, ' on Input Image']);
end

