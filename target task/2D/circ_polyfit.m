function[p,resid] = circ_polyfit(x,y)

x2 = circ_dist(x,circ_mean(x));
y2 = y;

[pp] = deal(cell(3,1));
rp = zeros(3,1);
p_init = polyfit(x2,y2,1);
efunc = @(ab) sum(circ_dist(ab(2)+ab(1)*x2,y2).^2);
feval = @(ab) ab(2)+ab(1)*x2;
[pp{1},rp(1)] = fminsearch(efunc,p_init);
[pp{2},rp(2)] = fminsearch(efunc,[0,circ_mean(y2)]);
[pp{3},rp(3)] = fminsearch(efunc,[1,0]);

p = pp{rp==min(rp)};

resid = circ_dist(y2,feval(p));

end