function [vesmask,CF,colored] = brightspots(C,th_sensitivity,visualize)
% Preliminary version for vesicle detection.

% 13.5.2008 (C) Pekka Ruusuvuori

if nargin < 3
    visualize = 0;
end
if nargin < 2
    th_sensitivity = 1.6;
end
win = 4; %3
fout = ones(2*win+1);
step = 3; %2
fout(win+1-step:win+1+step,win+1-step:win+1+step) = 0;
% these out-->th=1.20
fout(win+1-step,win+1-step) = 1;
fout(win+1-step,win+1+step) = 1;
fout(win+1+step,win+1-step) = 1;
fout(win+1+step,win+1+step) = 1;
fin = ones(2*win+1);
fin = fin - fout;
fout = fout/(sum(fout(:)));
fin = fin/(sum(fin(:)));

fout_res  = conv2(double(C),fout);
fin_res  = conv2(double(C),fin);
CF = fin_res./fout_res;
CF = CF(win+1:end-win,win+1:end-win);
%th = 0.65  % orig
%th = 1.17;  % w/o extra ones in fout
th = median(CF(:)) + th_sensitivity*std(CF(:));
vesmask = (CF>th);
vesmask(1:win,:) = 0; % use only valid area
vesmask(:,1:win) = 0; % use only valid area
vesmask(end-win:end,:) = 0; % use only valid area
vesmask(:,end-win:end) = 0; % use only valid area

% visualization
colored(:,:,1) = C;
colored(vesmask==1) = max(C(:));
colored(:,:,2) = C;
colored(:,:,3) = C;
colored = double(colored);
colored = colored - min(colored(:));
colored = colored/max(colored(:));
if visualize == 1
    figure, imshow(colored,[])
end