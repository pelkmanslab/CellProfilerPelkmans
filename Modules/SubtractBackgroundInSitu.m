function handles = SubtractBackgroundInSitu(handles)

% Help for the Substract Background InSitu module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Calculates a background thresholds by iteration and subtarct that value
% per image
% *************************************************************************
% Developed by Nico Battich
%
% $Revision: 4129 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the image to be corrected?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the corrected image?
%defaultVAR02 = SubsBlue
%infotypeVAR02 = imagegroup indep
CorrectedImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Select an a method for calculating the background value.
%choiceVAR03 = Otsu Global
%choiceVAR03 = Iterative Threshold
%choiceVAR03 = Single Value
ThresholdMethod = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu custom

%textVAR04 = For the Iterative Threshold method select the number of iterations.
%defaultVAR04 = 500
Iterations = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = For the Otsu Global and Iterative Threshold methods select a correction factor.
%defaultVAR05 = 1
CorrectionFactor = str2double(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = For the Otsu Global and Iterative Threshold methods select a minimum allowed value.
%defaultVAR06 = 0
Minimum = str2double(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = For the Otsu Global and Iterative Threshold methods select a maximum allowed value.
%defaultVAR07 = 1
Maximum = str2double(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = For the Single Value method enter the theshhold value.
%defaultVAR08 = 0.0018
SingleValue = str2double(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Do you want to correct bleed through from another channel?
%choiceVAR09 = No
%choiceVAR09 = Yes
Debleed = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = If yes. Enter the image you want to use for correction.
%infotypeVAR10 = imagegroup
CorrImageName = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%textVAR11 = If yes. Enter the bleeding ration, which the image above will be multiplied to before subtracting.
%defaultVAR11 = 0.065
BleedFactor = str2double(handles.Settings.VariableValues{CurrentModuleNum,11});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%check values
if Minimum>Maximum
    error('The minimum allowed value must be smaller than the maximum allowed value')
end

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

%do bleed through correction
if strcmpi(Debleed,'yes')
    
    %get the image to use for bleed correction
    BleedImage = CPretrieveimage(handles,CorrImageName,ModuleName,'MustBeGray','CheckScale');
    
    CorrectedImage=OrigImage;
    
    %subtract the bleed through
    CorrectedImage=CorrectedImage-(BleedImage.*BleedFactor);
    CorrectedImage(CorrectedImage<0)=0;
    
else
    CorrectedImage=OrigImage;
end


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

switch ThresholdMethod
    
    case 'Otsu Global'
        % calculate the threshold
        im=CorrectedImage(:);
        %[handles,OrigThreshold] = CPthreshold(handles,ThresholdMethod,[],0,1,CorrectionFactor,OrigImage,ImageName,ModuleName,);
        if max(im) == min(im)
            level = im(1);
        elseif isempty(im)
            level = 1;
        else
            %%% We want to limit the dynamic range of the image to 256. Otherwise,
            %%% an image with almost all values near zero can give a bad result.
            minval = max(im)/256;
            im(im < minval) = minval;
            im = log(im);
            minval = min (im);
            maxval = max (im);
            im = (im - minval) / (maxval - minval);
            level = exp(minval + (maxval - minval) * graythresh(im));
        end
        clear im
        
        level=level*CorrectionFactor;
        
        %check if level is within allowed values
        if level>Maximum
            level=Maximum;
        elseif level<Minimum
            level=Minimum;
        end
        
        
    case 'Iterative Threshold'
        %calculate threshold
        im=CorrectedImage(:);
        level=get_threshold_iterative(im,Iterations);
        level=level*CorrectionFactor;
        clear im
        
        %check if level is within allowed values
        if level>Maximum
            level=Maximum;
        elseif level<Minimum
            level=Minimum;
        end
        
     case 'Single Value'
        level=SingleValue*CorrectionFactor;
        
    otherwise
        error('enter thresholding method or value')
        
end

%do the background subtraction
CorrectedImage=CorrectedImage-level;
CorrectedImage(CorrectedImage<0)=0;

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByOne',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original
    %%% image, some intermediate images, and the final corrected image.
    %define the colormap
    matColorMap=colormap('jet');
    matColorMap(1,:)=[1 1 1];
    
    subplot(2,1,1);
    matCLim = [0 quantile(OrigImage(:),0.999)];
    CPimagesc(OrigImage,handles,matCLim);
    colormap(matColorMap)
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    %%% The mean image does not absolutely have to be present in order to
    %%% carry out the calculations if the illumination image is provided,
    %%% so the following subplot is only shown if MeanImage exists in the
    %%% workspace.
    subplot(2,1,2);
    matCLim = [0 quantile(CorrectedImage(:),0.999)];
    CPimagesc(CorrectedImage,handles,matCLim);
    title('Substracted Image');
    colormap(matColorMap)
    %%% Displays the text.
    if isempty(findobj('Parent',ThisModuleFigureNumber,'tag','DisplayText'))
        displaytexthandle = uicontrol(ThisModuleFigureNumber,'tag','DisplayText','style','text', 'position', [0 0 200 20],'fontname','helvetica','backgroundcolor',[0.7 0.7 0.9],'FontSize',handles.Preferences.FontSize);
    else
        displaytexthandle = findobj('Parent',ThisModuleFigureNumber,'tag','DisplayText');
    end
    displaytext = sprintf('Background threshold used:%.7f', level);
    set(displaytexthandle,'string',displaytext)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the corrected image to the handles structure so it can be used by
%%% subsequent modules.
handles.Pipeline.(CorrectedImageName) = CorrectedImage;