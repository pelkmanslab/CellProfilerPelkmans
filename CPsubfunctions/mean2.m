function out=mean2(in)
good=find(not(isnan(in)).*not(isinf(in)));
if isempty(good)
	out=NaN;
else
	out=mean(in(good));
end