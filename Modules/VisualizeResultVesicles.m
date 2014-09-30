function handles = VisualizeResultVesicles(handles)

% Help for the ExtractVesicle module:
% Category: Other
%
% SHORT DESCRIPTION:
% Visualize the result Image
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

%textVAR01 = Original Image Name In Matlab?
%defaultVAR01 = 
%infotypeVAR01 = objectgroup indep
inputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
disp(inputImageName)

fieldName = ['Segmented', inputImageName];
fieldName 
mask = double(handles.Pipeline.(fieldName));

%count = sum(sum(mask));

org = double(handles.Pipeline.(inputImageName));


subplot(1,2,1)
imshow(org,[]);
title(  'Original Image' );

subplot(1,2,2)
imshow(mask,[]);

count = length( unique(mask) )-1;
title( strcat( 'Result Image-', num2str(count),' particle detected') );
colormap('jet');
linkaxes;

end
