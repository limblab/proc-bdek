switch cds.meta.task
    
    case 'UCK'

        alldays = [];
        offsets = [-1.375 34.65];
        sorted_i = find(ismember(vertcat(cds.units.ID),1:100));
        arrayInfo = cell(length(sorted_i),4);
        for i = 1:length(sorted_i)
            arrayInfo{i,1} = cds.units(sorted_i(i)).array;
            arrayInfo{i,2} = cds.units(sorted_i(i)).spikes.ts;
            arrayInfo{i,3} = [cds.units(sorted_i(i)).bank, num2str(cds.units(sorted_i(i)).pin)];
            arrayInfo{i,4} = mean(cds.units(sorted_i(i)).spikes.wave);
        end

        alldays(1).PMd_units = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),2);
        alldays(1).M1_units = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),2);

        alldays(1).PMd_locations = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),3);
        alldays(1).M1_locations = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),3);

        alldays(1).PMd_waves = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),4);
        alldays(1).M1_waves = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),4);

        tt  = cds.trials(cds.trials.result=='R' | cds.trials.result=='F',:);
        ttA = cds.trials(cds.trials.result=='A' & ~isnan(cds.trials.tgtOnTime),:);

        for i = 1:size(cds.trials,2)
            colnm = cds.trials.Properties.VariableNames{i};
            if (~isnumeric(cds.trials.(colnm))) || (size(cds.trials.(colnm),2) > 1)
                tt.(cds.trials.Properties.VariableNames{i}) = [];
                ttA.(cds.trials.Properties.VariableNames{i}) = [];
            end
        end

        labels = tt.Properties.VariableNames;
        alldays(1).tt = table2array(tt);
        alldays(1).ttA = table2array(ttA);
        alldays(1).labels = cell(length(labels),1);
        for i = 1:length(labels)
            alldays(1).labels{i,:} = [sprintf('%d: ',i) labels{i} sprintf('(%s)',tt.Properties.VariableDescriptions{i})];
        end
        alldays(1).kin.pos = [cds.kin.t, cds.kin.x+offsets(1), cds.kin.y+offsets(2)];

        %% 
        ADTT = alldays(1).tt;
        ADTTA = alldays(1).ttA;
        conds = unique(cds.trials.numTgt);
        conds(conds<0) = [];
        for i = 1:size(conds,1)
            alldays(i).tt = ADTT(ADTT(:,end)==repmat(conds(i,:),size(ADTT,1),1),:);
            alldays(i).numtarg = conds(i);

            alldays(i).ttA = ADTTA(ADTTA(:,end)==repmat(conds(i,:),size(ADTTA,1),1),:);
        end

        %% Add reach directions to completed trials
        for i = 1:length(alldays)
            newcolnum = 19;
            for j = 1:size(alldays(i).tt,1)
                tot = find(alldays(1).kin.pos(:,1)<alldays(i).tt(j,3),1,'last');
                if ~isempty(tot)
                    alldays(i).tt(j,newcolnum) = atan2(alldays(1).kin.pos(tot,3),alldays(1).kin.pos(tot,2));
                else
                    alldays(i).tt(j,newcolnum) = NaN;
                end
            end
            alldays(i).tt(:,15) = alldays(i).tt(:,15)./180*pi;
        end
        alldays(1).labels{19} = '19: Reach angle';

        %% Add reach directions to aborted trials
        for i = 1:length(alldays)
            newcolnum = 19;
            for j = 1:size(alldays(i).ttA,1)
        %         totbeg = find(alldays(1).kin.pos(:,1)>alldays(i).ttA(j,7),1,'first');
                totend = find(alldays(1).kin.pos(:,1)<alldays(i).ttA(j,3),1,'last');
                if ~isempty(tot)
                    alldays(i).ttA(j,newcolnum) = atan2(alldays(1).kin.pos(totend,3),alldays(1).kin.pos(totend,2));
                else
                    alldays(i).ttA(j,newcolnum) = NaN;
                end
            end
            alldays(i).ttA(:,15) = alldays(i).ttA(:,15)./180*pi;
        end
        %% Separate out catch trials
        cueconditions = length(alldays);
        for i = 1:length(alldays)

            nocues = find(isnan(alldays(i).tt(:,8)));
            if ~isempty(nocues)
                alldays(i+cueconditions).tt = alldays(i).tt(nocues,:);
                alldays(i+cueconditions).numtarg = sprintf('%d catch',alldays(i).numtarg);
                alldays(i).tt(nocues,:) = [];
            end
        end

        %% Remove bad trials (missing tgtOnTime)
        for i = 1:length(alldays)

            missingtimes = find(isnan(alldays(i).tt(:,6)));
            if ~isempty(missingtimes)
                alldays(i).tt(missingtimes,:) = [];
            end

            if ~isempty(alldays(i).ttA)
                missingabort = find(isnan(alldays(i).ttA(:,6)));
                if ~isempty(missingabort)
                    alldays(i).tt(missingabort,:) = [];
                end
            end
        end

        %% get Reaction Times
        find_RT;
        for i = 1:length(alldays)
            alldays(i).tt(:,20) = alldays(i).tt(:,9) + RT{i}.*dt;
            alldays(i).ttA(:,20) = NaN;
        end

        alldays(1).labels{20} = '20: Movement onset time';

        % Clean up
        clear arrayInfo tt ADTT RT poscent labels sorted_i TT
        
    case 'UNT'
        
        alldays = [];
        offsets = [-3 33];
        sorted_i = find(ismember(vertcat(cds.units.ID),1:100));
        arrayInfo = cell(length(sorted_i),4);
        for i = 1:length(sorted_i)
            arrayInfo{i,1} = cds.units(sorted_i(i)).array;
            arrayInfo{i,2} = cds.units(sorted_i(i)).spikes.ts;
            arrayInfo{i,3} = [cds.units(sorted_i(i)).bank, num2str(cds.units(sorted_i(i)).pin)];
            arrayInfo{i,4} = mean(cds.units(sorted_i(i)).spikes.wave);
        end

        alldays(1).PMd_units = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),2);
        alldays(1).M1_units = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),2);

        alldays(1).PMd_locations = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),3);
        alldays(1).M1_locations = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),3);

        alldays(1).PMd_waves = arrayInfo(cellfun(@(x) strcmp(x,'PMd'),arrayInfo(:,1)),4);
        alldays(1).M1_waves = arrayInfo(cellfun(@(x) strcmp(x,'M1'),arrayInfo(:,1)),4);

        tt  = cds.trials(cds.trials.result=='R' | cds.trials.result=='F',:);
        ttA = cds.trials(cds.trials.result=='A' & ~isnan(cds.trials.tgtOnTime),:);
        
        SL = tt.cueDir;
        SLA = ttA.cueDir;

        for i = 1:size(cds.trials,2)
            colnm = cds.trials.Properties.VariableNames{i};
            if strcmp(colnm,'result')
                tt.(colnm) = (tt.result=='R')*32 + (tt.result=='F')*34; 
                ttA.(colnm) = (ttA.result=='A')*33;
                
            elseif (~isnumeric(cds.trials.(colnm))) || (size(cds.trials.(colnm),2) > 1)
                tt.(colnm) = [];
                ttA.(colnm) = [];
            elseif ~isempty(strfind(colnm,'Dir'))
                tt.(colnm) = deg2rad(tt.(colnm));
            end
        end
        
        labels = tt.Properties.VariableNames;
%         alldays(1).tt = table2array(tt);
%         alldays(1).ttA = table2array(ttA);
        alldays(1).labels = cell(length(labels),1);
        for i = 1:length(labels)
            alldays(1).labels{i,:} = [sprintf('%d: ',i) labels{i} sprintf('(%s)',tt.Properties.VariableDescriptions{i})];
        end
        alldays(1).kin.pos = [cds.kin.t, cds.kin.x+offsets(1), cds.kin.y+offsets(2)];

        %% Split into different epochs and add slices
        ADTT = table2array(tt);
        ADTTA = table2array(ttA);
        conds = (find(diff(cds.trials.tgtKappa)~=0 | abs(diff(cds.trials.cueKappa))>1000))';
        trialbreaks = [0 cds.trials.number(conds)' cds.trials.number(end)];
        for i = 1:(length(trialbreaks)-1)
            
            blockinds = find(ismember(ADTT(:,1),(trialbreaks(i)+1):trialbreaks(i+1)));
            blockindsA = find(ismember(ADTTA(:,1),(trialbreaks(i)+1):trialbreaks(i+1)));
            
            alldays(i).tt = ADTT(blockinds,:);
            alldays(i).ttA = ADTTA(blockindsA,:);
            
            slcs = SL(blockinds,:);
            slcsA = SLA(blockindsA,:);
            
%             alldays(i).tt = ADTT(ADTT(:,end)==repmat(conds(i,:),size(ADTT,1),1),:);
%             alldays(i).ttA = ADTTA(ADTTA(:,end)==repmat(conds(i,:),size(ADTTA,1),1),:);
%             
%             slcs = SL(ADTT(:,end)==repmat(conds(i,:),size(ADTT,1),1));
%             slcsA = SLA(ADTTA(:,end)==repmat(conds(i,:),size(ADTTA,1),1));
            
            dominant_slnum = mode(cellfun(@length,slcs));
            slcs(cellfun(@(x) length(x)~=dominant_slnum,slcs)) = {nan(1,dominant_slnum)};
            slcsA(cellfun(@(x) length(x)~=dominant_slnum,slcsA)) = {nan(1,dominant_slnum)};
            
            alldays(i).slices = deg2rad(cell2mat(slcs));
            alldays(i).slicesA = deg2rad(cell2mat(slcsA));
        end

        %% Add reach directions to completed trials
        for i = 1:length(alldays)
            newcolnum = 11;
            for j = 1:size(alldays(i).tt,1)
                tot = find(alldays(1).kin.pos(:,1)<alldays(i).tt(j,3),1,'last');
                if ~isempty(tot)
                    alldays(i).tt(j,newcolnum) = atan2(alldays(1).kin.pos(tot,3),alldays(1).kin.pos(tot,2));
                else
                    alldays(i).tt(j,newcolnum) = NaN;
                end
            end
            for j = 1:size(alldays(i).ttA,1)
                totend = find(alldays(1).kin.pos(:,1)<alldays(i).ttA(j,3),1,'last');
                if ~isempty(tot)
                    alldays(i).ttA(j,newcolnum) = atan2(alldays(1).kin.pos(totend,3),alldays(1).kin.pos(totend,2));
                else
                    alldays(i).ttA(j,newcolnum) = NaN;
                end
            end
            
        end
        alldays(1).labels{newcolnum} = sprintf('%d: Reach angle',newcolnum);

        %% Remove bad trials (missing tgtOnTime)
        for i = 1:length(alldays)

            missingtimes = find(isnan(alldays(i).tt(:,6)));
            if ~isempty(missingtimes)
                alldays(i).tt(missingtimes,:) = [];
                alldays(i).slices(missingtimes,:) = [];
            end

            if ~isempty(alldays(i).ttA)
                missingabort = find(isnan(alldays(i).ttA(:,6)));
                if ~isempty(missingabort)
                    alldays(i).tt(missingabort,:) = [];
                    alldays(i).slices(missingabort,:) = [];
                end
            end
        end
        
        %% Rearrange into old order
        alldays(1).labels{12} = '12: Cue Centroid';
        colorder = [2 8 9 5 6 7 3 4 12 11];
        for i = 1:length(alldays)
            centrd = circ_mean(alldays(i).slices,[],2);
            centrdA = circ_mean(alldays(i).slicesA,[],2);
            
            alldays(i).tt(:,12) = centrd;
            alldays(i).ttA(:,12) = centrdA; 
            
            alldays(i).tt = alldays(i).tt(:,colorder);
            alldays(i).ttA = alldays(i).ttA(:,colorder);
        
        end
        labls = alldays(1).labels;
        alldays(1).labels = cell(size(alldays(1).tt,2),1);
        for i = 1:size(alldays(1).tt,2)
            titlestart = strfind(labls{colorder(i)},':')+2;
            alldays(1).labels{i,:} = sprintf('%d: %s',i,labls{colorder(i)}(titlestart:end));
        end
        % Clean up
        clear arrayInfo tt ADTT RT poscent labels sorted_i TT
        
        
end