function handles = ShiftImage(handles)

% Help for the ShiftImage module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
%
% Authors:
%   Nico Battich
%   
% $Revision: 1808 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images to be shifted?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = x correction in pixels.
%defaultVAR02 = 0
x_correction = str2num(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = y correction in pixels.
%defaultVAR03 = 0
y_correction = str2num(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How do you want to call the corrected images?
%defaultVAR04 = ShiftedGreen
%infotypeVAR04 = imagegroup indep
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = What did you call the images to be be compared with? Only for plotting and control.
%infotypeVAR05 = imagegroup
ComparisonImageName = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

drawnow
%x_correction = 4
%y_correction = 4
%%% "OrigImage".
OrigImage = CPretrieveimage(handles,ImageName,ModuleName);

%%% seconf image
ComparisonImage = CPretrieveimage(handles,ComparisonImageName,ModuleName);



%%% get the output image
ImageOutput = zeros(size(OrigImage));


% define the change
if y_correction > 0
    final_index_one = abs(y_correction)+1:size(OrigImage,1);
    orig_index_one = 1:size(OrigImage,1)-abs(y_correction);
elseif y_correction < 0
    final_index_one = 1:size(OrigImage,1)-abs(y_correction);
    orig_index_one = abs(y_correction-1):size(OrigImage,1);%GG changed +1 to -1
else
    final_index_one = 1:size(OrigImage,1);
    orig_index_one = 1:size(OrigImage,1);
end
    

if x_correction > 0
    final_index_two = abs(x_correction+1):size(OrigImage,2);
    orig_index_two = 1:size(OrigImage,2)-abs(x_correction);
elseif x_correction < 0
    final_index_two = 1:size(OrigImage,2)-abs(x_correction);
    orig_index_two = abs(x_correction-1):size(OrigImage,2);%GG changed +1 to -1
else
    final_index_two = 1:size(OrigImage,2);
    orig_index_two = 1:size(OrigImage,2);
end
    

% change the 
ImageOutput(final_index_one,final_index_two)=OrigImage(orig_index_one,orig_index_two);

%save to handle structure
handles.Pipeline.(OutputName) = ImageOutput;


ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    R = ImageOutput;
    %R = OrigImage;
    q1 = quantile(R(:),0.17); q2 = quantile(R(:),0.92);
    R(R(:)<q1) = q1; R(R(:)>q2) = q2; R = R-q1; R = R / (q2-q1);
    
    G = ComparisonImage;
    q1 = quantile(G(:),0.17); q2 = quantile(G(:),0.92);
    G(G(:)<q1) = q1; G(G(:)>q2) = q2; G = G-q1; G = G/(q2-q1);
    
    B = zeros(size(OrigImage));
    
    PlottingImage = R;
    PlottingImage(:,:,2) = G;
    PlottingImage(:,:,3) = B;
    %PlottingImage=PlottingImage(400:1100,400:1100,:);
    %figure;imshow(PlottingImage)
    
    CPimagesc(PlottingImage,handles); 
end




















