function Nature_circ(centroid_angs,guesses,r,color_plot,newplot)

cp = r.*[cos(centroid_angs) sin(centroid_angs)];

gp = [cos(guesses) sin(guesses)];

if newplot==1
    figure; hold on;
    plot(cos(0:0.1:2*pi),sin(0:.1:2*pi),'k');
    plot(r.*cos(0:0.1:2*pi),r.*sin(0:.1:2*pi),color_plot);
    plot([gp(:,1) cp(:,1)]',[gp(:,2) cp(:,2)]',color_plot);
    plot(gp(:,1), gp(:,2),[color_plot '.']);
else
    plot(r.*cos(0:0.1:2*pi),r.*sin(0:.1:2*pi),color_plot);
    plot([gp(:,1) cp(:,1)]',[gp(:,2) cp(:,2)]',color_plot);
    plot(gp(:,1), gp(:,2),[color_plot '.']);
end