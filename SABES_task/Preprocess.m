lab=3;
ranBy='ranByBrian';
monkey='monkeyMihili';
task='taskUNT2D';
array1='arrayM1';
array2='arrayPMd';
%note the .nev extension is not necessary when providing the file name:
fname  = 'C:/Users/limblab/Desktop/NEV/sorted_Mihili_M1_SAB_10142016_001-01';
fname2 = 'C:/Users/limblab/Desktop/NEV/sorted_Mihili_PMd_SAB_10142016_001-01';
offsets = [-3 33];

% load data into cds:
%make blank cds class:
cds=commonDataStructure();
%load the data:
cds.file2cds(fname,lab,array1,monkey,task,ranBy,'ignoreJumps','mapFileC:/Users/limblab/Desktop/Mihili Left M1 SN 1025-001452.cmp');
cds.file2cds(fname2,lab,array2,monkey,task,ranBy,'ignoreJumps','mapFileC:/Users/limblab/Desktop/Mihili Left PMd SN 6251-001460.cmp');

%%
sorted_i = find(ismember(vertcat(cds.units.ID),1:100));
arrayInfo = cell(length(sorted_i),2);
for i = 1:length(sorted_i)
    arrayInfo{i,1} = cds.units(sorted_i(i)).array;
    arrayInfo{i,2} = cds.units(sorted_i(i)).spikes.ts;
    arrayInfo{i,3} = [cds.units(sorted_i(i)).bank, num2str(cds.units(sorted_i(i)).pin)];
end

alldays(1).PMd_units = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),2);
alldays(1).M1_units = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),2);

alldays(1).PMd_locations = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),3);
alldays(1).M1_locations = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),3);

tt = cds.trials(cds.trials.result=='R',:);
for i = 1:size(cds.trials,2)
    colnm = cds.trials.Properties.VariableNames{i};
    if (~isnumeric(cds.trials.(colnm))) || (size(cds.trials.(colnm),2) > 1)
        tt.(cds.trials.Properties.VariableNames{i}) = [];
    end
end


labels = tt.Properties.VariableNames;
alldays(1).tt = table2array(tt);
alldays(1).labels = cell(length(labels),1);
for i = 1:length(labels)
    alldays(1).labels{i,:} = [sprintf('%d: ',i) labels{i} sprintf('(%s)',tt.Properties.VariableDescriptions{i})];
end
alldays(1).kin.pos = [cds.kin.t, cds.kin.x+offsets(1), cds.kin.y+offsets(2)];

%%
conds = unique(alldays(1).tt(:,end-1:end),'rows');
for i = 1:size(conds,1)
    alldays(i+1).tt = alldays(1).tt(sum(alldays(1).tt(:,end-1:end)==repmat(conds(i,:),size(alldays(1).tt,1),1),2)==2,:);
    alldays(i+1).shifts = conds(i,:);
end

%%
figure; hold on; 
for i = 1:length(alldays)
    subplot(1,length(alldays),i); hold on;
    for j = 1:size(alldays(i).tt,1)
        xy = alldays(1).kin.pos(alldays(1).kin.pos(:,1) > alldays(i).tt(j,5) & alldays(1).kin.pos(:,1) < alldays(i).tt(j,7),2:3);
        plot(xy(:,1),xy(:,2)); 
    end
end

        

