function plot_bounds(lb,ub,xs,varargin)

lb = reshape(lb,1,[]);
ub = reshape(ub,1,[]);

if nargin > 2
    xs = reshape(xs,1,[]);
else
    xs = 1:length(lb);
end

plot([xs fliplr(xs)],[lb fliplr(ub)],'k');