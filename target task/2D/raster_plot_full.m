t1 = -200;
tt = alldays(1).tt;
units = alldays(1).PMd_units{6};
figure; hold on;
for tr = 1:length(tt)
    
    tstart = tt(tr,5)-0.2;
    tend = tt(tr,7);
    
    trirast = units(units>=tstart & units<=tend)-tt(tr,5);

    plot(trirast,tt(tr,10)*ones(size(trirast)),'k.','MarkerSize',4);
    %plot(0,tt(tr,10),'r.','MarkerSize',10); 
    plot(tt(tr,6)-tt(tr,5),tt(tr,10),'g.','MarkerSize',10);
    plot(tend-tt(tr,5),tt(tr,10),'r.','MarkerSize',10);
end