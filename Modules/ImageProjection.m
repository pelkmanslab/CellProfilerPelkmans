function handles = ImageProjection(handles)

% Help for the Combine module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Takes 1 to 7 color or grayscale images and combines them into 1. Each image's
% intensity can be adjusted independently.
% *************************************************************************
%
% Does an image projection
%
% Authors:
%   Berend Snijder
%
% Website: http://www.cellprofiler.org
%
% $Revision: 3524 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to call the projection image?
%defaultVAR01 = ProjectionImage
%infotypeVAR01 = imagegroup indep
ProjectionImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Select what projection method shold be applied
%choiceVAR02 = Mean
%choiceVAR02 = Median
%choiceVAR02 = Minimum
%choiceVAR02 = Maximum
%choiceVAR02 = Quantile_01
%choiceVAR02 = Quantile_05
%choiceVAR02 = Quantile_10
%choiceVAR02 = Quantile_20
%choiceVAR02 = Quantile_80
%choiceVAR02 = Quantile_90
%choiceVAR02 = Quantile_95
%choiceVAR02 = Quantile_99
%choiceVAR02 = Sum
ProjectionMethod = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu custom

%textVAR03 = Which image?
%choiceVAR03 = Leave this blank
%infotypeVAR03 = imagegroup
ImageNames{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Which image?
%choiceVAR04 = Leave this blank
%infotypeVAR04 = imagegroup
ImageNames{2} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Which image?
%choiceVAR05 = Leave this blank
%infotypeVAR05 = imagegroup
ImageNames{3} = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 = Which image?
%choiceVAR06 = Leave this blank
%infotypeVAR06 = imagegroup
ImageNames{4} = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 = Which image?
%choiceVAR07 = Leave this blank
%infotypeVAR07 = imagegroup
ImageNames{5} = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%textVAR08 = Which image?
%choiceVAR08 = Leave this blank
%infotypeVAR08 = imagegroup
ImageNames{6} = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 = Which image?
%choiceVAR09 = Leave this blank
%infotypeVAR09 = imagegroup
ImageNames{7} = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = Which image?
%choiceVAR10 = Leave this blank
%infotypeVAR10 = imagegroup
ImageNames{8} = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%textVAR11 = Which image?
%choiceVAR11 = Leave this blank
%infotypeVAR11 = imagegroup
ImageNames{9} = char(handles.Settings.VariableValues{CurrentModuleNum,11});
%inputtypeVAR11 = popupmenu

%textVAR12 = Which image?
%choiceVAR12 = Leave this blank
%infotypeVAR12 = imagegroup
ImageNames{10} = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%textVAR13 = Which image?
%choiceVAR13 = Leave this blank
%infotypeVAR13 = imagegroup
ImageNames{11} = char(handles.Settings.VariableValues{CurrentModuleNum,13});
%inputtypeVAR13 = popupmenu

%textVAR14 = Which image?
%choiceVAR14 = Leave this blank
%infotypeVAR14 = imagegroup
ImageNames{12} = char(handles.Settings.VariableValues{CurrentModuleNum,14});
%inputtypeVAR14 = popupmenu

%textVAR15 = Which image?
%choiceVAR15 = Leave this blank
%infotypeVAR15 = imagegroup
ImageNames{13} = char(handles.Settings.VariableValues{CurrentModuleNum,15});
%inputtypeVAR15 = popupmenu

%textVAR16 = Which image?
%choiceVAR16 = Leave this blank
%infotypeVAR16 = imagegroup
ImageNames{14} = char(handles.Settings.VariableValues{CurrentModuleNum,16});
%inputtypeVAR16 = popupmenu

%textVAR17 = Which image?
%choiceVAR17 = Leave this blank
%infotypeVAR17 = imagegroup
ImageNames{15} = char(handles.Settings.VariableValues{CurrentModuleNum,17});
%inputtypeVAR17 = popupmenu

%textVAR18 = Which image?
%choiceVAR18 = Leave this blank
%infotypeVAR18 = imagegroup
ImageNames{16} = char(handles.Settings.VariableValues{CurrentModuleNum,18});
%inputtypeVAR18 = popupmenu

%textVAR19 = Which image?
%choiceVAR19 = Leave this blank
%infotypeVAR19 = imagegroup
ImageNames{17} = char(handles.Settings.VariableValues{CurrentModuleNum,19});
%inputtypeVAR19 = popupmenu

%textVAR20 = Which image?
%choiceVAR20 = Leave this blank
%infotypeVAR20 = imagegroup
ImageNames{18} = char(handles.Settings.VariableValues{CurrentModuleNum,20});
%inputtypeVAR20 = popupmenu

%textVAR21 = Which image?
%choiceVAR21 = Leave this blank
%infotypeVAR21 = imagegroup
ImageNames{19} = char(handles.Settings.VariableValues{CurrentModuleNum,21});
%inputtypeVAR21 = popupmenu

%textVAR22 = Which image?
%choiceVAR22 = Leave this blank
%infotypeVAR22 = imagegroup
ImageNames{20} = char(handles.Settings.VariableValues{CurrentModuleNum,22});
%inputtypeVAR22 = popupmenu




%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% If selected, load images and create existance flag matrix for them
Images = {};
for i = 1:20
    if ~strcmpi(ImageNames{i},'Leave this blank')
        try
            Images{end+1} = CPretrieveimage(handles,ImageNames{i},ModuleName,'DontCheckColor','CheckScale'); %#ok Ignore MLint
        catch objFoo
            error(['Image processing was canceled in the ' ModuleName 'module because an error occurred while trying to load the image ' ImageNames{i} '. Please make sure it is a valid input image. Perhaps you chose an image that will be created later in the pipeline, in which case you should relocate the Combine module or the other one.']);
        end
    end
end
drawnow


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% If any of the images are binary/logical format, they must be converted
%%% to a double first before immultiply
Images = cellfun(@double,Images,'UniformOutput',0);


%%% Do projection over all images in the third dimension
switch ProjectionMethod
    case 'Mean'
        CombinedImage = mean(cat(3,Images{:}),3);
    case 'Median'
        CombinedImage = median(cat(3,Images{:}),3);
    case 'Minimum'
        CombinedImage = min(cat(3,Images{:}),[],3);
    case 'Maximum'
        CombinedImage = max(cat(3,Images{:}),[],3);
    case 'Quantile_01'
        CombinedImage = quantile(cat(3,Images{:}),0.01,3);
    case 'Quantile_05'
        CombinedImage = quantile(cat(3,Images{:}),0.05,3);
    case 'Quantile_10'
        CombinedImage = quantile(cat(3,Images{:}),0.10,3);
    case 'Quantile_20'
        CombinedImage = quantile(cat(3,Images{:}),0.20,3);
    case 'Quantile_80'
        CombinedImage = quantile(cat(3,Images{:}),0.80,3);
    case 'Quantile_90'
        CombinedImage = quantile(cat(3,Images{:}),0.90,3);
    case 'Quantile_95'
        CombinedImage = quantile(cat(3,Images{:}),0.95,3);
    case 'Quantile_99'
        CombinedImage = quantile(cat(3,Images{:}),0.99,3);
    case 'Sum'    
        CombinedImage = sum(cat(3,Images{:}),3);
        
    otherwise
        error('unknown projection method %s\n',ProjectionMethod)
end




%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

matSubPlotIX = [4,8,9,10,11,12];

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(CombinedImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    %%% A subplot of the figure window is set to display the Combined Image
    %%% image.  Using CPimagesc or image instead of imshow doesn't work when
    %%% some of the pixels are saturated.
    subplot(3,4,[1,2,3,5,6,7]);
    CPimagesc(CombinedImage,handles);
    title(sprintf('''%s'', %s projection, cycle #%d',ProjectionImageName,ProjectionMethod,num2str(handles.Current.SetBeingAnalyzed)));
    %%% A subplot of the figure window is set to display Image 1.
    for i = 1:(min(length(Images),length(matSubPlotIX)))
        subplot(3,4,matSubPlotIX(i))
        CPimagesc(Images{i},handles);
        title(sprintf('Image %d: %s',i,ImageNames{i}));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the adjusted image to the handles structure so it can be used by
%%% subsequent modules.
handles.Pipeline.(ProjectionImageName) = CombinedImage;
