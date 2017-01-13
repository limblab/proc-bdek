function[xplot,yplot] = smooth_circle_hist(angles,nbin)

[n,xout] = hist(angles,nbin);

%normalize n
nnorm = n./sum(n);
nnorm = nnorm./max(nnorm);

pts = [nnorm'.*cos(xout') nnorm'.*sin(xout')];

xplot = pts(:,1);
yplot = pts(:,2);

xplot = [xplot; xplot(1)];
yplot = [yplot; yplot(1)];

