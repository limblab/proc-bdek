function behavior_scatter(alldays)

llist = flipud(unique(alldays(2).tt(:,3)));
figure; hold on; 
cls = {'b','r','g'};

prioloc = circ_mean(alldays(2).tt(:,2)); if prioloc < 0; prioloc = prioloc + 2*pi; end
centroids = circ_dist(alldays(2).tt(:,9),prioloc)+prioloc;
reaches = circ_dist(alldays(2).tt(:,10),prioloc)+prioloc;
[~,~,kfits,~,fitlines] = behavior_fit_circ([alldays(2).tt(:,2) centroids reaches alldays(2).tt(:,3)]);
for i = 1:length(llist)
    
    is = find(alldays(2).tt(:,3)==llist(i));

    plot(centroids(is),reaches(is),'.','Color',cls{i});
    
    offset = fitlines{i}(2,round(end/2)) - prioloc;
    if  abs(offset) > pi
        plot(fitlines{i}(1,:),fitlines{i}(2,:)-sign(offset)*2*pi,cls{i},'LineWidth',2);
    else
        plot(fitlines{i}(1,:),fitlines{i}(2,:),cls{i},'LineWidth',2);
    end
 
end
    
title(sprintf('%.2f   %.2f',1/(kfits{1}+1),1/(kfits{2}+1)));