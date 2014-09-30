function handles = SetImageBackgroundToZero(handles)

% Help for the Set Image Background To Zero module:
% Category: Image Processing
%
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4076 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the image to be rescaled?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the objects you want to identify background around?
%infotypeVAR02 = objectgroup
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the rescaled image?
%defaultVAR03 = RescaledBlue
%infotypeVAR03 = imagegroup indep
RescaledImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Filter size (segmentation expansion of X pixels, as well as median filtering with X pixels)
%defaultVAR04 = 5
intFilterSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%%%VariableRevisionNumber = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image to be analyzed and assigns it to a variable,
%%% "OrigImage".
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'DontCheckColor','CheckScale');

% Retreives the segmentation of the selected object. Converts to logical
% immediately.
matObjectSegmentation = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage))>0;

% must be a round number.
intFilterSize = round(intFilterSize);

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow


% calculate bin-edges
matHistogramEdges = linspace(min(OrigImage(:)),max(OrigImage(:)),200);

% slightly expand current object segmentation, and filter out single pixel
% values...
matSE = strel('disk',intFilterSize,0);
matObjectSegmentation = imdilate(matObjectSegmentation,matSE);

% let's do median filtering, to get rid of smaller areas between objects
matObjectSegmentation = medfilt2(matObjectSegmentation, [intFilterSize,intFilterSize]);


% if there is enough empty space to determine background...
if mean(matObjectSegmentation(:)==0) > 0.1

    % do histigram binning of background (non-segmented) values
    matHistCount = histc(OrigImage(matObjectSegmentation==0),matHistogramEdges);
    matHistCountSignal = histc(OrigImage(matObjectSegmentation>0),matHistogramEdges);
    
    % find peak of background distribution, set this to 0.
    intBackgroundValue = matHistogramEdges(findmax(matHistCount));
    
    % let's also calculate a rue signal value
    intForegroundvalue = matHistogramEdges(findmax(matHistCountSignal));
    

else
    
    % set this to NaN, it's not that important.
    intForegroundvalue = NaN;
    
    % if we already have other background estimates, average out those
    if isfield(handles.Measurements.Image,['Background_',ImageName])
        matPreviousBackgroundValues = cell2mat(handles.Measurements.Image.(['Background_',ImageName])');
        intBackgroundValue = nanmean(matPreviousBackgroundValues(:,1));
        fprintf('%s: warning, not enough empty space available, using mean background value so far\n',mfilename)
    else
        % otherwise, use a very low (conservative) quantile as backrgound
        % estimate
        intBackgroundValue = quantile(OrigImage(:),0.0001);
        fprintf('%s: warning, not enough empty space available, using minimum (0.0001 quantile) image intensity\n',mfilename)
    end
    
end
    
fprintf('%s: background value for image ''%s'' (new name ''%s'') object ''%s'' = %g.\n',mfilename,ImageName,RescaledImageName,ObjectName,intBackgroundValue)

% do the rescale
RescaledImage = OrigImage - intBackgroundValue;
% set pixels below 0 to 0...
RescaledImage(RescaledImage<0)=0;

% right now we just do a shift, but we could also do a true rescale...


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow




ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
%%% Check whether that figure is open. This checks all the figure handles
%%% for one whose handle is equal to the figure number for this module.
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
%     if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
%         CPresizefigure(OrigImage,'TwoByOne',ThisModuleFigureNumber)
%     end
    %%% A subplot of the figure window is set to display the original image.
    subplot(2,2,1); 
    imagesc(OrigImage,[0 quantile(OrigImage(:),0.999)]);
    colormap(gca,hot)
    title([ImageName,', cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% A subplot of the figure window is set to display the Rescaled
    %%% Image.
    subplot(2,2,2); 
    imagesc(RescaledImage,[0 quantile(RescaledImage(:),0.999)]);
    colormap(gca,hot)
    title(RescaledImageName);

    % also histogram true signal, for plotting purposes only.
    matHistCountSignal = histc(OrigImage(matObjectSegmentation>0),matHistogramEdges);
    
    subplot(2,2,3:4); 
    hold on
    cla reset
    plot(matHistogramEdges,matHistCount,matHistogramEdges,matHistCountSignal); 
    legend({'background','signal'})
    vline(intBackgroundValue,'-r')
    if ~isnan(intForegroundvalue)
        vline(intForegroundvalue,'-y');
    end
    % also draw background distribution
    if handles.Current.SetBeingAnalyzed>1
        matPreviousBackgroundValues = cell2mat(handles.Measurements.Image.(['Background_',ImageName])');
        vline(nanmean(matPreviousBackgroundValues(:,1)),'-k');
    end
    hold off
    title('Background histogram');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The Rescaled image is saved to the handles structure so it can be
%%% used by subsequent modules.
handles.Pipeline.(RescaledImageName) = RescaledImage;

% store background value...
handles.Measurements.Image.(['Background_',ImageName]){handles.Current.SetBeingAnalyzed} = [intBackgroundValue,intForegroundvalue];
handles.Measurements.Image.(['Background_',ImageName,'Features']) = {'Background','Foreground'};

%%%
% let's do a little check here...
RescaledImage2 = CPretrieveimage(handles,RescaledImageName,ModuleName);
isequal(RescaledImage,RescaledImage2)
%%%


function [intMaxXindex, intMaxX] = findmax(X)
%
% [Berend Snijder]. Returns the absolute maximum value
%
% [intMaxX, intMaxXindex] = absmax(X)

[intMaxX, intMaxXindex]=nanmax(abs(X));

