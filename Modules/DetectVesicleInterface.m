function handles = DetectVesicleInterface(handles)

% Help for the ExtractVesicle module:
% Category: Other
%
% SHORT DESCRIPTION:
% Interface for vesicle detection by Sharif, Chowdhury
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 4904 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = For which object Name?
%defaultVAR01 = Vesicles
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
% disp(ObjectName)

%textVAR02 = Block Size (+/-)(8-64)?
%defaultVAR02 = 32
%infotypeVAR02 = objectgroup indep
strBlockSize = char(handles.Settings.VariableValues{CurrentModuleNum,2});

blkSize = str2num(strBlockSize);


%textVAR03 = Apply Median Filter On Background (+/-)(0/1)?
%defaultVAR03 = 1
%infotypeVAR03 = objectgroup indep
strMedFlag = char(handles.Settings.VariableValues{CurrentModuleNum,3});
% disp(strMedFlag)

medFlag = str2num(strMedFlag);

%textVAR04 = Median Filter Window Length (+/-) (3/5/7)?
%defaultVAR04 = 3
%infotypeVAR04 = objectgroup indep
strMedWindowLen = char(handles.Settings.VariableValues{CurrentModuleNum,4});
% disp(strMedWindowLen)

medWindowLen = 1;
if medFlag > 0.5
    medWindowLen = str2num(strMedWindowLen);
end

%textVAR05 = Replace Median Filter by Min Filter at Threshold level (+) (0-1)?
%defaultVAR05 = 1
%infotypeVAR05 = objectgroup indep
strMFThresh = char(handles.Settings.VariableValues{CurrentModuleNum,5});
% disp(strMFThresh)

MFThresh = 255;

if medFlag > 0.5
    MFThresh  = str2num(strMFThresh)*255;
end





%textVAR06 = Background Level Weight (-) (0-2)?
%defaultVAR06 = 1
%infotypeVAR06 = objectgroup indep
strK1 = char(handles.Settings.VariableValues{CurrentModuleNum,6});
% disp(strK1)

k1 = str2num(strK1);


%textVAR07 = Background Variance Weight (-) (0-3)?
%defaultVAR07 = 1.5
%infotypeVAR07 = objectgroup indep
strK2 = char(handles.Settings.VariableValues{CurrentModuleNum,7});
% disp(strK2)

k2 = str2num(strK2);


%textVAR08 = Maximum Vesicle Size in Pixel (1..)?
%defaultVAR08 = 50
%infotypeVAR08 = objectgroup indep
strVesicleSize = char(handles.Settings.VariableValues{CurrentModuleNum,8});


vesicleSize = str2num(strVesicleSize);


% disp(strVesicleSize)

%textVAR09 = Input Image Name In Matlab?
%defaultVAR09 = OrigGreen
%infotypeVAR09 = objectgroup indep
inputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = Secondary Object Name ( 1- don't use mask)?
%defaultVAR10 = 1
%infotypeVAR10 = objectgroup indep
% disp(inputImageName)
secondaryMaskName = char(handles.Settings.VariableValues{CurrentModuleNum,10});


% disp(secondaryMaskName)


%I = handles.Pipeline.(inputImageName);
I = CPretrieveimage(handles,inputImageName,ModuleName,'MustBeGray','CheckScale');


%secondaryLabelMatrixImage = CPretrieveimage(handles,['SmallRemovedSegmented', secondaryMaskName],ModuleName,'DontCheckColor','DontCheckScale',size(I));

if secondaryMaskName  ~='1'
    fieldname = ['Segmented',secondaryMaskName];
    secondaryMaskImage = CPretrieveimage(handles,fieldname,ModuleName,'DontCheckColor','DontCheckScale',size(I));
    secondaryMaskImage = double(double(secondaryMaskImage)>0.5);
else
    secondaryMaskImage = double( ones( size(I) ) );
end



if max(max(max(I)))>255
    I =  uint8( round(double(I) / 256 ) );
end
maxVal = max(max(I));
if maxVal<2
    I = uint8(double(I)*256);
end

[fI fI2 vcI Ib Is Iflt centr] = detectVesicle(I,1,0,blkSize,medFlag,medWindowLen, k1 ,k2,vesicleSize,MFThresh);
[r c] = size(fI);
l= length(centr(:,1));


fI2 = fI;
for i=1:l
    if centr(i,1)<=r && centr(i,2)<=c && centr(i,1)>0 && centr(i,2)>0
        fI2(centr(i,1), centr(i,2)) =  -fI2(centr(i,1), centr(i,2));
    end
end

fI3 = getconnectedPoints(fI2);
fI3 = double(fI3<0.8).* fI;
fI3  = fI3 .*  secondaryMaskImage;
fieldName = ['Segmented', ObjectName];

usedLabel = unique(fI3);
maxLabel = max(max(fI3))+1;
relabelled = zeros(1, maxLabel);
len = length(usedLabel);

for i=1:len
    currentColor = usedLabel(i)+1;
    relabelled(currentColor) = i-1; 
end
for i= 1:r
    for j=1:c
        fI3(i,j)  = relabelled( fI3(i,j)+1 );
    end
end

mxfI3 = max(max(fI3));
lenfI3 = length(unique(fI3));

handles.Pipeline.(fieldName)= fI3;




ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
 CPfigure(handles,'Image',ThisModuleFigureNumber);
subplot(1,2,1)
imshow(fI3,[])
title('Vesicles');
subplot(1,2,2)
imshow(I,[])
title('Original Image');

linkaxes
colormap('jet')
end
% %end
