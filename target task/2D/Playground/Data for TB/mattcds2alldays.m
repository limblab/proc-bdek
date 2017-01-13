alldays = [];
for blck = 1:length(C);
    cds = C{blck}.cds;
    offsets = [-1 30];
    sorted_i = find(ismember(vertcat(cds.units.ID),1:100));
    arrayInfo = cell(length(sorted_i),5);
    for i = 1:length(sorted_i)
        arrayInfo{i,1} = cds.units(sorted_i(i)).array;
        arrayInfo{i,2} = cds.units(sorted_i(i)).spikes.ts;
        arrayInfo{i,3} = [cds.units(sorted_i(i)).bank, num2str(cds.units(sorted_i(i)).pin)];
        arrayInfo{i,4} = mean(cds.units(sorted_i(i)).spikes.wave);
        arrayInfo{i,5} = double(cds.units(sorted_i(i)).chan)+double(cds.units(sorted_i(i)).ID)*.1;
    end

    alldays(blck).PMd_units = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),2);
    alldays(blck).M1_units = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),2);

    alldays(blck).PMd_locations = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),3);
    alldays(blck).M1_locations = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),3);

    alldays(blck).PMd_waves = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),4);
    alldays(blck).M1_waves = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),4);

    alldays(blck).PMd_ID = cell2mat(arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),5));
    alldays(blck).M1_ID = cell2mat(arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),5));
    
    tt  = cds.trials(cds.trials.result=='R' | cds.trials.result=='F',:);
    res = (cds.trials.result=='R')*32 + (cds.trials.result=='F')*34;
    res(res==0) = [];

    for i = 1:size(cds.trials,2)
        colnm = cds.trials.Properties.VariableNames{i};
        if strcmp(colnm,'result')
            tt.(colnm) = res;
        elseif (~isnumeric(cds.trials.(colnm))) || (size(cds.trials.(colnm),2) > 1)
            tt.(cds.trials.Properties.VariableNames{i}) = [];
        elseif strcmp(colnm,'tgtDir')
            tt.(colnm) = mod(polyval([-2*pi/8 pi/2],tt.tgtID),2*pi);
        end
    end

    labels = tt.Properties.VariableNames;
    alldays(blck).tt = table2array(tt);
    alldays(blck).labels = cell(length(labels),1);
    for i = 1:length(labels)
        alldays(blck).labels{i,:} = [sprintf('%d: ',i) labels{i} sprintf('(%s)',tt.Properties.VariableDescriptions{i})];
    end
    alldays(blck).pos = [cds.kin.t, cds.kin.x+offsets(1), cds.kin.y+offsets(2)];

end
%%
CP = alldays(1).PMd_ID;
CM = alldays(1).M1_ID;
for i = 2:length(alldays)
    CP = intersect(CP,alldays(i).PMd_ID);
    CM = intersect(CM,alldays(i).M1_ID);
end

for i = 1:length(alldays)
    alldays(i).PMd_units(~ismember(alldays(i).PMd_ID,CP)) = [];
    alldays(i).PMd_ID(~ismember(alldays(i).PMd_ID,CP)) = [];
    alldays(i).PMd_locations(~ismember(alldays(i).PMd_ID,CP)) = [];
    alldays(i).PMd_waves(~ismember(alldays(i).PMd_ID,CP)) = [];
    
    alldays(i).M1_units(~ismember(alldays(i).M1_ID,CM)) = [];
    alldays(i).M1_ID(~ismember(alldays(i).M1_ID,CM)) = [];
    alldays(i).M1_locations(~ismember(alldays(i).M1_ID,CM)) = [];
    alldays(i).M1_waves(~ismember(alldays(i).M1_ID,CM)) = [];
end
%%
timeoff = zeros(length(alldays),1); 
for i = 1:length(alldays)
    timeoff(i) = max([cellfun(@(x) max(x),alldays(i).PMd_units);cellfun(@(x) max(x),alldays(i).M1_units);max(alldays(i).pos(:,1))])+1;
end
timeoff = [0 ; cumsum(timeoff(1:end-1))];

timecols = find(cellfun(@(x) ~isempty(strfind(x,'Time')),alldays(1).labels));
for i = 1:length(alldays)
    alldays(i).tt(:,timecols) = alldays(i).tt(:,timecols) + timeoff(i);
    
    alldays(i).PMd_units = cellfun(@(x) x+timeoff(i),alldays(i).PMd_units,'Uni',0);
    alldays(i).M1_units = cellfun(@(x) x+timeoff(i),alldays(i).M1_units,'Uni',0);
    alldays(i).pos(:,1) = alldays(i).pos(:,1)+timeoff(i);
end

for i = 1:length(alldays(1).PMd_units)
    alldays(1).PMd_units{i} = cell2mat(cellfun(@(x) x{i},{alldays.PMd_units},'Uni',0)');
end
for i = 1:length(alldays(1).M1_units)
    alldays(1).M1_units{i} = cell2mat(cellfun(@(x) x{i},{alldays.M1_units},'Uni',0)');
end
alldays(1).kin.pos = cell2mat({alldays.pos}');
alldays(1).reachlength = nanmean(sqrt(C{1}.cds.trials.tgtCtr(:,1).^2 + C{1}.cds.trials.tgtCtr(:,2).^2));
alldays(1).targetsize = nanmean(sqrt((diff(C{1}.cds.trials.tgtCorners(:,[1 3]),[],2).^2 + diff(C{1}.cds.trials.tgtCorners(:,[2 4]),[],2).^2)./2));
alldays = rmfield(alldays,'pos');

    

    