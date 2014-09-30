function HelpChangeColorMap
%   Help for the Change Colormap within the Image Tool Window:
%   Category: Image Tools
%   *****************************************************************
%
%   Change Colormap - Opens a window that allows you to change the 
%   colormap of the selected figure. You can select the default 
%   colormap (which you can set under File > Set Preferences) or 
%   any other predetermined colormap. Note that the colormap selected 
%   will apply to all non-RGBimages in the entire figure, and not 
%   only to the image selected. The Apply To All button will change 
%   the colormap in all module display windows and any other windows 
%   that contain images. If you are running the developer's version 
%   of CellProfiler, you can also open a colormap editor, which 
%   enables you to create personalized colormaps. It will modify the
%   colormap of the last active figure, so be careful if you open it, 
%   click another figure and go back to it, because you might be 
%   changing the colormap of a figure you did not intend to change. 
%   See also Help > General Help > Colormaps.

CPhelpdlg(help('HelpChangeColorMap'))