function handles = EnhanceImageForSegmentation(handles)

% Help for the EnhanceImageForSegmentation module:
% Category: Image Processing
% This module give outputs seven distances:
% 
%
%
%

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[~, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the Edge Image?
%infotypeVAR01 = imagegroup
EdgeImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the original outline image?
%infotypeVAR02 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the finel enhanced image?
%defaultVAR03 = EnhancedImage
%infotypeVAR03 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Enter a weigth for the original image. Range is [0,1].
%defaultVAR04 = 0.7
OrigThreshold = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Enter a weigth for the Filtered Image. Range is [0,1].
%defaultVAR05 = 0.285
FiltThreshold = str2double(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Enter a weigth for the Edge Image. Range is [0,1].
%defaultVAR06 = 0.015
EdgeThreshold = str2double(handles.Settings.VariableValues{CurrentModuleNum,6});

%%% get images.
EdgeImage = CPretrieveimage(handles,EdgeImageName,ModuleName,'MustBeGray','DontCheckScale');
OrigImage = CPretrieveimage(handles,OrigImageName,ModuleName,'MustBeGray','DontCheckScale');


maskIm = bwlabel(EdgeImage);
maskIm = imfill(maskIm);
maskIm = maskIm>0; 
%imagesc(maskIm)
maskIm = bwlabel(maskIm);


BackSubstIm = double(OrigImage)-mean(OrigImage(~maskIm));
BackSubstIm(BackSubstIm<0) = 0;
%imagesc(BackSubstIm,[0 quantile(BackSubstIm(:), 0.99)])
%imagesc(OutlineIm,[0 quantile(OutlineIm(:), 0.99)])


[FinalMask,CurrentObjLabels] = bwdist(maskIm);
FinalMask = (FinalMask < 2).*maskIm(CurrentObjLabels);
FinalMask = imfill(FinalMask);

STATS = regionprops(FinalMask,'area');
matArea = [0;struct2array(STATS)'];
numThres = quantile(matArea,0.9);
%matIndex = matArea<numThres;
AreaMask = matArea(FinalMask(:)+1);
AreaMask = reshape(AreaMask,size(FinalMask,1),size(FinalMask,2));
FinalAreaMask = AreaMask>numThres;

BackSubstIm2 = BackSubstIm;
BackSubstIm2(~FinalAreaMask(:)) = 0;

%imagesc(BackSubstIm2,[0 quantile(BackSubstIm(:), 0.99)])

OutputImage = (double(OrigImage).*OrigThreshold + BackSubstIm2.*FiltThreshold + double(EdgeImage).*double(quantile(OrigImage(:), 0.99)).*EdgeThreshold);
% OutputImage = uint16(OutputImage);

% subplot(2,1,1)
% imagesc(OutputImage,[quantile(OutputImage(:), 0.1) quantile(OutputImage(:), 0.99)])
% subplot(2,1,2)
% imagesc(OrigImage,[quantile(OrigImage(:), 0.1) quantile(OrigImage(:), 0.99)])

handles.Pipeline.(OutputName) = OutputImage;










