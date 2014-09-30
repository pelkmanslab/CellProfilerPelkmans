function [f,fnames] = addfeature(f,fnames,newf,newfname)
% Adds new feature to feature matrix f
% and appends a new feature name.

% 12.3.2008 (C) Pekka Ruusuvuori
if isempty(f)  % first feature to be added
    f = newf;
    runcount = 1;
    for ind = 1:size(newf,2)
        if size(newf,2) > 1
            fnames{ind} = [newfname, num2str(runcount)];
        else
            fnames{ind} = newfname;
        end
        runcount = runcount + 1;
    end
elseif size(newf,1) == size(f,1) % add new one
    f = [f,newf];
    runcount = 1;
    for ind = 1:size(newf,2)
        if size(newf,2) > 1
            fnames{length(fnames)+1} = [newfname, num2str(runcount)];
        else
            fnames{length(fnames)+1} = newfname;
        end
        runcount = runcount + 1;
    end
else
    error(['Feature length for ' newfname ' does not match: this should never happen!'])
end