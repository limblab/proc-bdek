fun1 = @(targs,pd) find(targs == (mod(pd-1,8)+1));
fun2 = @(targs,pd) find(targs == (mod(pd,8)+1) | targs == (mod(pd-2,8)+1));
fun3 = @(targs,pd) find(targs == (mod(pd+2-1,8)+1) | targs == (mod(pd+6-1,8)+1));  
fun4 = @(targs,pd) find((mod(targs+4-1,8)+1) == (mod(pd,8)+1) | (mod(targs+4-1,8)+1) == (mod(pd-2,8)+1));
fun5 = @(targs,pd) find(targs == (mod(pd+4-1,8)+1));

F = {fun1,fun2,fun3,fun4,fun5};

windows = {6, [-200 750];...
           7, [   0 250];...
           8, [   0 250];...
           20,[-100 200]};
           
alldays(1).tunedPMd_units = alldays(1).PMd_units(~isnan(PDS)); 
PDtuned = PDnum(~isnan(PDnum));
tt1 = alldays(1).tt; 
tt1(isnan(tt1(:,20)),20) = tt1(isnan(tt1(:,20)),9); 
TS = cell(length(PDtuned),1);
for i = 1:length(PDtuned) % loop through units
    clc; fprintf('%d/%d\n',i,length(PDtuned));
    for j = 1:length(F) % loop through direction types (PD,OD,etc) 

        trialinds = F{j}(Tnum{1},PDtuned(i));
        
        combrast = cell(1,size(windows,1));
        for k = 1:size(windows,1) % loop through time windows
            rast = raster_get(alldays(1).tunedPMd_units{i},tt1(trialinds,:),windows{k,2}/1000,windows{k,1});
            binrast = bin_array(rast,1,size(rast,2)./50,'mean');
            combrast{k} = reshape(repmat(binrast,50,1),1,[]);
        end
        
        for trep = 1:400
%             TS{i,:}{j}{k}{trep} = find(poissrnd(expandrast)>0)./1000;
            allts = find(poissrnd(cell2mat(combrast))>0)./1000;
            TS{i,:}{j}{1}{trep} = allts(allts<=0.2);
            TS{i,:}{j}{2}{trep} = allts(allts<=1.2 & allts>0.2);
            TS{i,:}{j}{3}{trep} = allts(allts>=1.2);
        end
    end
end

%%
tcols = [6 7 8 9 20 3];

ts = unique(cumsum(reshape(abs(cell2mat(windows(:,2))'),1,[])./1000));
targetlist = Tnum{2};

ACT = cell(length(PDtuned),1);
tt2 = NaN(length(targetlist),20);
tstart = tt1(end,3)+5;


for trial = 1:length(targetlist)
    clc; fprintf('compiling trial %d/%d\n',trial,length(targetlist));
    targnum = targetlist(trial);
    reachdir = mod(targnum*pi/4,2*pi);
    
    tt2(trial,[15 19]) = reachdir;
    tt2(trial,18) = 2;
    
    for j = 1:length(tcols)
        tt2(trial,tcols(j)) = tstart + ts(j);
    end
    
    for i = 1:length(PDtuned)
       
       tunetype = find(cellfun(@(x) ~isempty(x(targnum,PDtuned(i))),F));
       if tunetype == 5
           ACT{i} = [ACT{i}; tstart+TS{i}{5}{1}{1}'; tstart+TS{i}{1}{2}{1}'; tstart + TS{i}{5}{3}{1}'];
           TS{i}{5}{1}(1) = [];
           TS{i}{1}{2}(1) = [];
           TS{i}{5}{3}(1) = [];
       elseif tunetype == 4
           ACT{i} = [ACT{i}; tstart+TS{i}{4}{1}{1}'; tstart+TS{i}{2}{2}{1}'; tstart + TS{i}{4}{3}{1}'];
           TS{i}{4}{1}(1) = [];
           TS{i}{2}{2}(1) = [];
           TS{i}{4}{3}(1) = [];
       else
           ACT{i} = [ACT{i}; tstart+TS{i}{tunetype}{1}{1}'; tstart+TS{i}{tunetype}{2}{1}'; tstart+TS{i}{tunetype}{3}{1}'];
           TS{i}{tunetype}{1}(1) = [];
           TS{i}{tunetype}{2}(1) = [];
           TS{i}{tunetype}{3}(1) = [];
       end
    end
    tstart = tt2(trial,3)+0.5;
end
    
    