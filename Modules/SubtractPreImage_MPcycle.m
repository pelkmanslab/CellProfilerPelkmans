function handles = SubtractPreImage_MPcycle(handles)

% Help for the SUBTRACTPREIMAGE_MPCYCLE module:
% Category: Other
%
% SHORT DESCRIPTION: 
% Obtains intensity images of two subsequent multiplexing cycles from the
% handles and subtracts the image of the previous cycle from the currently
% analysed cycle.
%
% DETAILS:
% This module requires both images to be loaded from different TIFF folders
% (using an additional LoadImages.m module and providing an alternative
% pathname to TIFF folder). 
%
% Note that the filenames of both images must only differ in one character
% (the cycle number). This ensures that the correct images (from the same
% well, site, and channel) are used for subtraction. Make sure that you
% have named your image files accordingly.
% Example:
%   Current Cycle Image Filename: 
%           '20140207MPcycle2_B02_T0001F001L01A02Z01C02.png'
%   Previous Cycle Image Filename:
%           '20140207MPcycle1_B02_T0001F001L01A02Z01C02.png'
%   
% *************************************************************************
%
% Author:
%    Markus Herrmann <markus.herrmann@imls.uzh.ch>
%

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the aligned intensity images of the current cycle?
%infotypeVAR01 = imagegroup
CurrImNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = 
%choiceVAR02 = Do not use
%infotypeVAR02 = imagegroup
CurrImNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What did you call the aligned intensity images of the previous cycle?
%infotypeVAR03 = imagegroup
PreImNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = 
%choiceVAR04 = Do not use
%infotypeVAR04 = imagegroup
PreImNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = How do you want to call the subtracted intensity images?
%defaultVAR05 = SubtrGreen
%infotypeVAR05 = imagegroup indep
OutputNameList{1} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = 
%defaultVAR06 = /
%infotypeVAR06 = imagegroup indep
OutputNameList{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});



%%%%%%%%%%%%%%%%%%
%%% PROCESSING %%%
%%%%%%%%%%%%%%%%%%

CurrImages = cell(1,length(CurrImNameList));
PreImages = cell(1,length(PreImNameList));
OutputImages = cell(1,length(CurrImNameList));
for j = 1:length(CurrImNameList)
    CurrImageName = CurrImNameList{j};
    if strcmpi(CurrImageName,'Do not use')
        continue
    end
    PreImageName = PreImNameList{j};
    if strcmpi(PreImageName,'Do not use')
        continue
    end
    
    %%% make sure the correct images are loaded for later subtraction
    CurrFilename = char(handles.Pipeline.(sprintf('Filename%s',CurrImageName)));
    PreFilename = char(handles.Pipeline.(sprintf('Filename%s',PreImageName)));
    % the two images from the two different multiplexing cycles should be
    % from identical sites/channels, hence their filenames should only
    % differ in one character, namely the cycle number
    if sum(~(CurrFilename == PreFilename))~=1 %difference of more than one character between cycles
        error('%s: images from two different sites/channels were selected',mfilename)
    end
    
    %%% retrieve intensity images
    CurrImages{j} = CPretrieveimage(handles,CurrImageName,ModuleName,'MustBeGray','CheckScale');
    PreImages{j} = CPretrieveimage(handles,PreImageName,ModuleName,'MustBeGray','CheckScale');
    
    %%% subtract intensity images
    OutputImages{j} = imsubtract(CurrImages{j},PreImages{j});
    OutputImages{j}(OutputImages{j}<0) = 0;
   
end



%%%%%%%%%%%%%%%
%%% DISPLAY %%%
%%%%%%%%%%%%%%%


drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
       
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    subplot(2,3,1); 
    CPimagesc(CurrImages{1},handles);
    title(sprintf('Current cycle image ''%s''', CurrImNameList{1}));
    
    subplot(2,3,2); 
    CPimagesc(PreImages{1},handles);
    title(sprintf('Previous cycle image ''%s''', PreImNameList{1}));
    
    subplot(2,3,3); 
    CPimagesc(OutputImages{1},handles);
    title(sprintf('Subtracted image ''%s''', OutputNameList{1}));
    
    if sum(cellfun(@isempty,OutputImages))==0
        
        subplot(2,3,4);
        CPimagesc(CurrImages{2},handles);
        title(sprintf('Current cycle image ''%s''', CurrImNameList{2}));
        
        subplot(2,3,5);
        CPimagesc(PreImages{2},handles);
        title(sprintf('Previous cycle image ''%s''', PreImNameList{2}));
        
        subplot(2,3,6);
        CPimagesc(OutputImages{2},handles);
        title(sprintf('Subtracted image ''%s''', OutputNameList{2}));
        
    end
    
    drawnow
end



%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%


%%% save shifted images to handles
for j = 1:length(CurrImNameList)
    CurrImageName = CurrImNameList{j};
    if strcmpi(CurrImageName,'Do not use')
        continue
    end
    PreImageName = PreImNameList{j};
    if strcmpi(PreImageName,'Do not use')
        continue
    end
    
    if strcmpi(OutputNameList{j},'/')
        error('name for output intensity image not specified')
    else
        handles.Pipeline.(OutputNameList{j}) = OutputImages{j};
    end
end
    
end
