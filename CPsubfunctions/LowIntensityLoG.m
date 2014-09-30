function varargout=LowIntensityLoG(OrigImage,LoGSigma,LoGSize,LoGCutoff,MarkerCutoff,MarkerFraction,MinArea)
%LOWINTENSITYLOG simple but powerfull Laplacion of Gaussian object detection
    %
    %OBJECTMASK=LOWINTENSITYLOG(ORIGIMAGE,LOGSIGMA,LOGSIZE,LOGCUTOFF,MARKERCUTOFF,MARKERFRACTION,MINAREA)
    %[OBJECTMASK,LOGIMAGE,MARKERIMAGE]=LOWINTENSITYLOG(ORIGIMAGE,LOGSIGMA,LOGSIZE,LOGCUTOFF,MARKERCUTOFF,MARKERFRACTION,MINAREA)
    %
    %INTRODUCTION:
    %
    %The result of the LoG filter shows a sign transition at the boundary of
    %objects, but also for background structures. LOGIMAGE can visualize
    %this finding. As the gradient is measured, object boundaries will have 
    %(absolute) high values, the inner part of objects lower values.
    %Background structures in general have a higher value than objects.
    %
    %Only positive parts of the LOGIMAGE are further analyzed.
    %Setting LOGCUTOFF below the lowest (positive) background value can 
    %separate objects from background but will inavertably also remove
    %parts of the objects or even complete low intensity objects. If no
    %further segmentation is performed, clumped objects have to be 
    %separated using a sufficiently low LOGCUTOFF. In order to detect even
    %the lowest intensity objects, while taking some clumping into account
    %the cutoff should be disabled(LOGCUTOFF=0). Further segmentation may
    %have to be performed in this case.
    %
    %For either setting of LOGCUTOFF(=0 or >0) the background can
    %effectively removed without risking removing object parts!
    %A MARKERIMAGE is first created via a MARKERCUTOFF, containing smaller
    %regions within objects. In the next step structures (positive valued 
    %in the LoG result) are selected as objects based on wether a marker lies 
    %inside. To prevent background structures from being selected by very few
    %marker pixels, a minimum marker to object area ratio(MARKERFRACTION) is
    %used. This procedure allows a very stringent cutoff while still recovering
    %the original object shape. MARKERCUTOFF is usually chosen smaller than
    %LOGCUTOFF (if enabled).
    %
    %STRATEGY FOR CHOSING FILTER VALUES:
    %
    %This function is primarily designed to provide a coarse (under)segmentation 
    %for extremely dim images (living cells, hoechst etc.). The aim is to
    %recover all objects contained in the image and later sub-segment via shape
    %analysis (PerimeterAnalysis.m, PerimeterSegmentation.m). In this case
    %proceed as follows:
    %
    %SIGMA should be chosen as low as possible, but such that no
    %sub-segmentation/erosion of object structures occurs! Lower values will
    %provide a more detailed object outline and better segmentation. Higher
    %values give smoother outlines, better detection of low intensity objects
    %but a tendency to merge objects.(Yokogawa 20x: try 10)
    %
    %FILTERSIZE should be square and at least 8xSIGMA (eg [80 80]) large
    %filters can demand a lot of computational time.
    %
    %MARKERCUTOFF should be chosen by looking at imshow(LOGIMAGE) and picked
    %slightly above the lowest value of the object with the overall highest
    %values (to still place a marker in this object). The value has to be lower
    %than the lowest value of background structures. A too high cutoff will
    %lead to detecting backround as objects, a too low value will lead to
    %objects not being detected. (eg 16)
    %
    %MARKERFRACTION should be chosen by looking at
    %imagesc(OBJECTMASK+MARKERIMAGE). A very agressive cutoff may lead
    %to small markers and requires a smaller MARKERFRACTION. A too small value
    %may lead to detection of background structures as objects. Values may
    %range from 0 to 1. (Try: 0.1)
    %
    %LOGCUTOFF is disabled by setting it to zero, or a rather high value is
    %chosen if a high detection fidelity for low intensity objects is
    %still required.
    %
    %MINAREA allows to discard too small objects.
    %
    %Example:
    %TestMask=LowIntensityLoG(OrigImage,10,[80,80],16.5,15.8,0.1,80);
    %
    %[Anatol Schwab 21.11.12]


    LoGSpecial=fspecial('log',LoGSize,LoGSigma);
    LoG=imfilter(double(OrigImage),LoGSpecial,'replicate');
    LoGScaled=log(abs(LoG)).*sign(LoG);
    if LoGCutoff~=0
        LoGPositive=(LoGScaled>0)&(LoGScaled<LoGCutoff);%positive valued parts of image,apply cutoff to segemt objects
    else
        LoGPositive=(LoGScaled>0);%do not apply cutoff for objects
    end
    LoGLabeled=bwlabel(LoGPositive);%all structures (background and objects), labeled
    LoGMarkers=(LoGScaled<MarkerCutoff) & LoGPositive;%generate markers
    %get rid of background  structures via Marker vs Object size, also discard objects that do not contain a marker
    LoGMarkersLabeled=LoGMarkers.*LoGLabeled;%Project structure labels on markers - this allows structure/marker comparisson via ID
    ObjectsSize=regionprops(LoGLabeled,'Area');
    MarkerSize=regionprops(LoGMarkersLabeled,'Area');%list may be shorter!
    TrueLabels=[];
    for i=1:length(MarkerSize)
      if MarkerSize(i).Area>0%object i has overlap with the marker
          if MarkerSize(i).Area/ObjectsSize(i).Area > MarkerFraction%check how much area the marker covers
              if ObjectsSize(i).Area>=MinArea
                TrueLabels=[TrueLabels,i];
              end
          end
      end
    end

    %Labels=setdiff(unique(FirstLoGLabeled.*FirstLoGMarkers),0);
    ObjectMask=ismember(LoGLabeled,TrueLabels);%clever trick to select objects via their ID from the labeled image
    ObjectMask=imfill(ObjectMask,'holes');
    
    %Output
    if nargout==1
        varargout{1}=ObjectMask;
    elseif nargout==3
        varargout{1}=ObjectMask;
        varargout{2}=LoGScaled;
        varargout{3}=LoGMarkers;
    end
end


