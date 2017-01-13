function [smoothed] = circ_smooth(x,n)

if size(x,2) > size(x,1)
    x = x';
end

xR = repmat(x,3,1);
y = smooth(xR,n);
smoothed = y((length(x)+1):(2*length(x)));