PDact = alldays(1).PMd_units;
for i = 1:length(PDS)
    clc; fprintf('%d/%d\n',i,length(PDS));
    if ~isnan(PDnum(i))
        targPD = find(Tnum{1}==PDnum(i));
        targOD = find(Tnum{1}==(mod(PDnum(i)+3,8)+1));
        targOFFPD = find(abs(circ_dist(alldays(1).tt(:,19),PDS(i)))<=pi/4); targOFFPD(ismember(targOFFPD,targPD)) = [];
        targOFFOD = find(abs(circ_dist(alldays(1).tt(:,19)+pi,PDS(i)))<=pi/4); targOFFOD(ismember(targOFFOD,targOD)) = [];
        targORTH = find(Tnum{1}==(mod(PDnum(i)+1,8)+1) | Tnum{1}==(mod(PDnum(i)+5,8)+1));
        
        rP1 = raster_get(alldays(1).PMd_units{i},alldays(1).tt(targPD,:),[0 750]/1000,6);
        rP2 = raster_get(alldays(1).PMd_units{i},alldays(1).tt(targPD,:),[0 250]/1000,7);
                
        rPO1 = raster_get(alldays(1).PMd_units{i},alldays(1).tt(targOFFPD,:),[0 750]/1000,6);
        rPO2 = raster_get(alldays(1).PMd_units{i},alldays(1).tt(targOFFPD,:),[0 250]/1000,7);
               
        
        rORTH1 = raster_get(alldays(1).PMd_units{i},alldays(1).tt(targORTH,:),[0 750]/1000,6);
        rORTH2 = raster_get(alldays(1).PMd_units{i},alldays(1).tt(targORTH,:),[0 250]/1000,7);
        
        % Replace OD
        for j = 1:length(targOD) 
            starttime1 = alldays(1).tt(targOD(j),6);
            stoptime1 = alldays(1).tt(targOD(j),6)+750/1000;
            starttime2 = alldays(1).tt(targOD(j),7);
            stoptime2 = alldays(1).tt(targOD(j),7)+250/1000;
            
            randtimes1 = starttime1+(find(rP1(randperm(size(rP1,1),1),:))./1000)';
            randtimes2 = starttime2+(find(rP2(randperm(size(rP2,1),1),:))./1000)';
            
            PDact{i}(PDact{i}>starttime1 & PDact{i}<stoptime1) = [];
            PDact{i}(PDact{i}>starttime2 & PDact{i}<stoptime2) = [];
            
            PDact{i} = sortrows([PDact{i}; randtimes1; randtimes2]);
        end
        % Replace Off of OD
        for j = 1:length(targOFFOD)
            starttime1 = alldays(1).tt(targOFFOD(j),6);
            stoptime1 = alldays(1).tt(targOFFOD(j),6)+750/1000;
            starttime2 = alldays(1).tt(targOFFOD(j),7);
            stoptime2 = alldays(1).tt(targOFFOD(j),7)+250/1000;
            
            randtimes1 = starttime1+(find(rPO1(randperm(size(rPO1,1),1),:))./1000)';
            randtimes2 = starttime2+(find(rPO2(randperm(size(rPO2,1),1),:))./1000)';
            
            PDact{i}(PDact{i}>starttime1 & PDact{i}<stoptime1) = [];
            PDact{i}(PDact{i}>starttime2 & PDact{i}<stoptime2) = [];
            
            PDact{i} = sortrows([PDact{i}; randtimes1; randtimes2]);
        end
        % Replace PD
        for j = 1:length(targPD)
            starttime1 = alldays(1).tt(targPD(j),6);
            stoptime1 = alldays(1).tt(targPD(j),6)+750/1000;
            starttime2 = alldays(1).tt(targPD(j),7);
            stoptime2 = alldays(1).tt(targPD(j),7)+250/1000;
            
            rP1noj = rP1; rP1noj(j,:) = [];
            rP2noj = rP2; rP2noj(j,:) = [];
            randtimes1 = starttime1+(find(rP1noj(randperm(size(rP1noj,1),1),:))./1000)';
            randtimes2 = starttime2+(find(rP2noj(randperm(size(rP2noj,1),1),:))./1000)';
            
            PDact{i}(PDact{i}>starttime1 & PDact{i}<stoptime1) = [];
            PDact{i}(PDact{i}>starttime2 & PDact{i}<stoptime2) = [];
            
            PDact{i} = sortrows([PDact{i}; randtimes1; randtimes2]);
        end
        % Replace Off PD
        for j = 1:length(targOFFPD)
            starttime1 = alldays(1).tt(targOFFPD(j),6);
            stoptime1 = alldays(1).tt(targOFFPD(j),6)+750/1000;
            starttime2 = alldays(1).tt(targOFFPD(j),7);
            stoptime2 = alldays(1).tt(targOFFPD(j),7)+250/1000;
            
            rPO1noj = rPO1; rPO1noj(j,:) = [];
            rPO2noj = rPO2; rPO2noj(j,:) = [];
            randtimes1 = starttime1+(find(rPO1noj(randperm(size(rPO1noj,1),1),:))./1000)';
            randtimes2 = starttime2+(find(rPO2noj(randperm(size(rPO2noj,1),1),:))./1000)';
            
            PDact{i}(PDact{i}>starttime1 & PDact{i}<stoptime1) = [];
            PDact{i}(PDact{i}>starttime2 & PDact{i}<stoptime2) = [];
            
            PDact{i} = sortrows([PDact{i}; randtimes1; randtimes2]);
        end
        % Replace Orth
        for j = 1:length(targORTH)
            starttime1 = alldays(1).tt(targORTH(j),6);
            stoptime1 = alldays(1).tt(targORTH(j),6)+750/1000;
            starttime2 = alldays(1).tt(targORTH(j),7);
            stoptime2 = alldays(1).tt(targORTH(j),7)+250/1000;
            
            rORTH1noj = rORTH1; rORTH1noj(j,:) = [];
            rORTH2noj = rORTH2; rORTH2noj(j,:) = [];
            randtimes1 = starttime1+(find(rORTH1noj(randperm(size(rORTH1noj,1),1),:))./1000)';
            randtimes2 = starttime2+(find(rORTH2noj(randperm(size(rORTH2noj,1),1),:))./1000)';
            
            PDact{i}(PDact{i}>starttime1 & PDact{i}<stoptime1) = [];
            PDact{i}(PDact{i}>starttime2 & PDact{i}<stoptime2) = [];
            PDact{i} = sortrows([PDact{i}; randtimes1; randtimes2]);
        end
        
    else
       PDact{i} = [];
    end
end

fakePMd_units = cell(length(alldays(1).PMd_units),1);
timeoffsetter = max(alldays(1).tt(:,3))+1;
for i = 1:length(alldays(1).PMd_units)
    if ~isnan(PDnum(i))
        fakePMd_units{i} = [alldays(1).PMd_units{i}; PDact{i}+timeoffsetter];
    else
        fakePMd_units{i} = [];
    end
end

timecols = [3 5 6 7 8 9 10 20];
alldays(2).tt = alldays(1).tt; 
alldays(2).tt(:,timecols) = alldays(2).tt(:,timecols)+timeoffsetter;