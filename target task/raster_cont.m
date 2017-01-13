before = 0;
after = 1.5;
width = 80;

[lowtrain,hightrain,lowind,highind] = raster(bdf,train,trials2,before,after);
lowtrain(lowtrain > 1) = 0;
hightrain(hightrain > 1) = 0;

LOcont_sig = zeros(size(lowtrain));
for i = 1:size(lowtrain,1)
    LOcont_sig(i,:) = train2cont(lowtrain(i,:),width);
end
LO = LOcont_sig;
nonzerolow_ind = find(sum(lowtrain,2)~=0);
LO = LO(nonzerolow_ind,:);
LOav_cont = mean(LOcont_sig,1);

nonzerolow = lowind(nonzerolow_ind);

HIcont_sig = zeros(size(hightrain));
for i = 1:size(hightrain,1)
    HIcont_sig(i,:) = train2cont(hightrain(i,:),width);
end
HI = HIcont_sig;
nonzerohigh_ind = find(sum(hightrain,2)~=0);
HI = HI(nonzerohigh_ind,:);
HIav_cont = mean(HIcont_sig,1);

nonzerohigh = highind(nonzerohigh_ind);

gca;
hold on; 

xlimits = (width:(length(LOcont_sig)+1-width))+before*1000;
ylimits = width:(length(LOcont_sig)+1-width);

scaling = (size(lowtrain,1)+size(hightrain,1)) / max([LOav_cont HIav_cont]);

plot(xlimits,scaling.*LOav_cont(ylimits),'b');

plot(xlimits,scaling.*HIav_cont(ylimits),'r');

figure; hold on;
plot(xlimits,LOav_cont(ylimits),'b');

plot(xlimits,HIav_cont(ylimits),'r');

