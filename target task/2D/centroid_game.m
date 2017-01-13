h2 = figure; hold on;
ths = 0:0.01:2*pi;

tlocs = 2*pi*rand(200,1);

[x,y] = deal(zeros(200,1));
for i = 1:200

    slice1 = sli{2}(i,:);

    sh = slice1 + (tlocs(i) - circ_mean(slice1'));


    % Set up circle

    p = ang2poly(0,2*pi,7,9);    
    patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','FaceColor','b','FaceAlpha',1);

    plot([7*cos(sh); 9*cos(sh)],...
         [7*sin(sh); 9*sin(sh)],'r','LineWidth',5); 
    title(sprintf('Trial: %d/%d',i,200),'FontSize',16);
    axis equal
     
    [x(i),y(i)] = ginput(1);
    pause(1);
    cla;

end

errs = circ_dist(atan2(y,x),tlocs);

exper_kappa = circ_kappa(errs);

fprintf('Centroid Kappa: %0.1f\n',exper_kappa);
