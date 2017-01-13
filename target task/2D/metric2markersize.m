function [marksize] = metric2markersize(x,allmetrics,minsize_maxsize)

minsize = minsize_maxsize(1);
maxsize = minsize_maxsize(2);
rnge = [min(allmetrics) max(allmetrics)];

marksize = minsize + (maxsize-minsize).*((x-rnge(1))./diff(rnge));