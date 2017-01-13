PDact = alldays(1).PMd_units;
windows = {6,[0 750];...
           7,[0 250]};

getPDi = @(targs,pd) find(targs == (mod(pd-1,8)+1));
getODi = @(targs,pd) find(targs == (mod(pd+4-1,8)+1));
getORTHi = @(targs,pd) find(targs == (mod(pd+2-1,8)+1) | targs == (mod(pd+6-1,8)+1));  
getOFFPDi = @(targs,pd) find(targs == (mod(pd,8)+1) | targs == (mod(pd-2,8)+1));
getOFFODi = @(targs,pd) find((mod(targs+4-1,8)+1) == (mod(pd,8)+1) | (mod(targs+4-1,8)+1) == (mod(pd-2,8)+1));
       
tt = alldays(1).tt;
lowlim = min(hist(Tnum{1},8));
throwtrials = cell(8,1);
for i = 1:8
    ttrials = find(Tnum{1}==i);
    throwtrials{i} = ttrials(randperm(length(ttrials),length(ttrials)-lowlim));
end
allthrows = cell2mat(throwtrials);
tt(allthrows,:) = [];
targnums = Tnum{1};
targnums(allthrows) = [];

for i = 1:length(PDS)
    clc; fprintf('%d/%d\n',i,length(PDS));
    if ~isnan(PDnum(i))
        targPD = getPDi(targnums,PDnum(i));
        targOD = getODi(targnums,PDnum(i));
        targOFFPD = getOFFPDi(targnums,PDnum(i));
        targOFFOD = getOFFODi(targnums,PDnum(i));
        targORTH = getORTHi(targnums,PDnum(i));
        
        for k = 1:size(windows,1)
            rP = raster_get(alldays(1).PMd_units{i},tt(targPD,:),windows{k,2}/1000,windows{k,1});
            rPO = raster_get(alldays(1).PMd_units{i},tt(targOFFPD,:),windows{k,2}/1000,windows{k,1});
            rORTH = raster_get(alldays(1).PMd_units{i},tt(targORTH,:),windows{k,2}/1000,windows{k,1});

            % Replace OD
            for j = 1:length(targOD) 
                starttime = tt(targOD(j),windows{k,1})-windows{k,2}(1)/1000;
                stoptime = tt(targOD(j),windows{k,1})+windows{k,2}(2)/1000;

                randtimes = starttime+(find(rP(randperm(size(rP,1),1),:))./1000)';
                randtimes = downsample(randtimes,2);

                PDact{i}(PDact{i}>starttime & PDact{i}<stoptime) = [];
                PDact{i} = sortrows([PDact{i}; randtimes]);
            end
            % Replace Off of OD
            for j = 1:length(targOFFOD)
                starttime = tt(targOFFOD(j),windows{k,1})-windows{k,2}(1)/1000;
                stoptime = tt(targOFFOD(j),windows{k,1})+windows{k,2}(2)/1000;
                
                randtimes = starttime+(find(rPO(randperm(size(rPO,1),1),:))./1000)';
                randtimes = downsample(randtimes,2);

                PDact{i}(PDact{i}>starttime & PDact{i}<stoptime) = [];
                PDact{i} = sortrows([PDact{i}; randtimes]);
            end
            % Replace PD
            for j = 1:length(targPD)
                starttime = tt(targPD(j),windows{k,1})-windows{k,2}(1)/1000;
                stoptime = tt(targPD(j),windows{k,1})+windows{k,2}(2)/1000;

                rPnoj = rP; rPnoj(j,:) = [];
                randtimes = starttime+(find(rPnoj(randperm(size(rPnoj,1),1),:))./1000)';
                randtimes = downsample(randtimes,2);
                     
                PDact{i}(PDact{i}>starttime & PDact{i}<stoptime) = [];
                PDact{i} = sortrows([PDact{i}; randtimes]);
            end
            % Replace Off PD
            for j = 1:length(targOFFPD)
                starttime = tt(targOFFPD(j),windows{k,1})-windows{k,2}(1)/1000;
                stoptime = tt(targOFFPD(j),windows{k,1})+windows{k,2}(2)/1000;

                rPOnoj = rPO; rPOnoj(j,:) = [];
                randtimes = starttime+(find(rPOnoj(randperm(size(rPOnoj,1),1),:))./1000)';
                randtimes = downsample(randtimes,2);
               
                PDact{i}(PDact{i}>starttime & PDact{i}<stoptime) = [];
                PDact{i} = sortrows([PDact{i}; randtimes]);
            end
            % Replace Orth
            for j = 1:length(targORTH)
                starttime = tt(targORTH(j),windows{k,1})-windows{k,2}(1)/1000;
                stoptime = tt(targORTH(j),windows{k,1})+windows{k,2}(2)/1000;

                rORTHnoj = rORTH; rORTHnoj(j,:) = [];
                randtimes = starttime+(find(rORTHnoj(randperm(size(rORTHnoj,1),1),:))./1000)';
                randtimes = downsample(randtimes,2);
                
                PDact{i}(PDact{i}>starttime & PDact{i}<stoptime) = [];
                PDact{i} = sortrows([PDact{i}; randtimes]);
            end
        end
        
    else
       PDact{i} = [];
    end
end

% fakePMd_units = cell(length(alldays(1).PMd_units),1);
% timeoffsetter = max(tt(:,3))+1;
% for i = 1:length(alldays(1).PMd_units)
%     if ~isnan(PDnum(i))
%         fakePMd_units{i} = [alldays(1).PMd_units{i}; PDact{i}+timeoffsetter];
%     else
%         fakePMd_units{i} = [];
%     end
% end
% 
% timecols = [2:10 20];
% alldays(2).tt = tt; 
% alldays(2).tt(:,timecols) = alldays(2).tt(:,timecols)+timeoffsetter;