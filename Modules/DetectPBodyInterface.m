function handles = DetectPBodyInterface(handles)

% Help for the ExtractVesicle module:
% Category: Other
%
% SHORT DESCRIPTION:
% Interface for P-body detection by Sharif Chowdhury
% 
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

%textVAR01 = For which object?
%defaultVAR01 = PBody
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
disp(ObjectName)

%textVAR02 = P-body radius (-) (2-8)?
%defaultVAR02 = 2
%infotypeVAR02 = objectgroup indep
strRad = char(handles.Settings.VariableValues{CurrentModuleNum,2});
rad = str2num(strRad);


%textVAR03 = Percentange of P-body pixels (+) (0-100)?
%defaultVAR03 = 15
%infotypeVAR03 = objectgroup indep
strPPixel = char(handles.Settings.VariableValues{CurrentModuleNum,3});
disp(strPPixel)
pPixel = str2num(strPPixel);
pPixel = pPixel/100;

%textVAR04 = Score Threshold (-) (0..1)?
%defaultVAR04 = 0.1
%infotypeVAR04 = objectgroup indep
strScoreThresh = char(handles.Settings.VariableValues{CurrentModuleNum,4});
disp(strScoreThresh)
scoreThresh = str2num(strScoreThresh);
%textVAR05 = Input Image Name In Matlab?
%defaultVAR05 = 
%infotypeVAR05 = objectgroup indep




%textVAR06 = Secondary Object Name ( 1- don't use mask)?
%defaultVAR06 = 1
%infotypeVAR06 = objectgroup indep
secondaryMaskName = char(handles.Settings.VariableValues{CurrentModuleNum,6});


%textVAR07 = Max P-Body Size in Pixel ?
%defaultVAR07 = 25
%infotypeVAR07 = objectgroup indep
strMaxSize = char(handles.Settings.VariableValues{CurrentModuleNum,7});

maxSize = str2num(strMaxSize);

%textVAR08 = Min P-Body Size in Pixel?
%defaultVAR08 = 3
%infotypeVAR08 = objectgroup indep
strMinSize = char(handles.Settings.VariableValues{CurrentModuleNum,8});

minSize = str2num(strMinSize);




inputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});
disp(inputImageName)

I = CPretrieveimage(handles,inputImageName,ModuleName,'MustBeGray','CheckScale');
I = double(I);
while max(max(I)) > 1.5 % scale I in 0-1
    I = I/256;
end

if strcmp( secondaryMaskName ,'1') == 0
    fieldname = ['Segmented',secondaryMaskName];
    secondaryMaskImage = CPretrieveimage(handles,fieldname,ModuleName,'DontCheckColor','DontCheckScale',size(I));
    secondaryMaskImage = double(double(secondaryMaskImage)>0.5);
else
    secondaryMaskImage = double( ones( size(I) ) );
end



[pBodies peak] = detect_particles(I,rad,scoreThresh,pPixel,[0 0]);



pBodies = double(pBodies) .* secondaryMaskImage;
I = double(I);
maxVal = max(max(I))

if maxVal<2
    I = I*256;
end


I2 =  uint8( round(I + pBodies*256));

[fI fI2 vcI Ib Is Iflt centr] = detectVesicle(I2,1,0,8,0,3, 0.5 ,1,maxSize ,256);
usedColors = unique(fI.* pBodies);


len = length(usedColors);

maxColor = max(max(fI))+1;

colorFlag = zeros(1, maxColor);


for i= 1:len
    colorFlag(usedColors(i)+1)= 1;
end

[r c] = size(fI);
ar = regionprops(fI, 'area');
for i=1:r
    for j=1:c
        if fI(i,j)>0
            if colorFlag( fI(i,j)+1)<0.5  || ar(fI(i,j)).Area < minSize 
                 fI(i,j)= 0;
            end
        end
    end
end

retImage =getUniqueCompactImageLABEL(pBodies, fI);
fieldName = ['Segmented', ObjectName];
handles.Pipeline.(fieldName)= retImage ;
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    subplot(1,2,1)
    imshow(retImage ,[])
    title('P-Body');
    subplot(1,2,2)
    imshow(I,[])
    title('Original Image');
    linkaxes
    colormap('jet')
end
