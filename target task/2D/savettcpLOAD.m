% FLS = {'07', '11'; ...
%        '07', '12'; ...
%        '07', '15'; ...
%        '07', '19'; ...
%        '08', '01'; ...
%        '08', '06'; ...
%        '08', '12'; ...
%        '08', '15'; ...
%        '08', '22'; ...
%        '09', '04'; ...
%        '09', '05'; ...
%        '09', '06'; ...
%        '09', '26'; ...
%        '10', '02'; ...
%        '10', '07'};
%    
% FLS = {'05', '04'; ...
%        '05', '05'; ...
%        '05', '06'; ...
%        '07', '08'};
   
FLS = {'10', '08'; ...
       '10', '27'; ...
       '11', '02'};

for sessnum = 1:size(FLS,1)
   
    MO = FLS{sessnum,1};
    DA = FLS{sessnum,2};
    
    ttcp_LOAD;
    
    flname = sprintf('C:/Users/limblab/Desktop/Mihili_bdf/full alldays/full_alldays_%s_%s',monkey,[MO DA YE]);
    save(flname,'alldays','priors','monkey','MO','DA','YE');

    alldays = alldays_nokin;
    flname2 = sprintf('C:/Users/limblab/Desktop/Mihili_bdf/Processed/%s_%s',monkey,[MO DA YE]);
    save(flname2,'alldays','priors','monkey','MO','DA','YE');
    
    clearvars -except FLS
end

%%
for daynum = 1:size(FLS,1)
    fprintf('%s-%s ',FLS{daynum,1},FLS{daynum,2});
    [alldays, priors, monkey, MO, DA, YE] = load_processed('Mihili',[FLS{daynum,1} FLS{daynum,2} '2015'],0);
    
    badtrials = cell(length(alldays),1);
    for i = 2:length(badtrials)
        badtrials{i} = find(abs(alldays(i).tt(:,3)) < eps | abs(alldays(i).tt(:,3))> 100000 | isnan(alldays(i).tt(:,3)));
    end
    fprintf('(Found %d bad trials)\n',sum(cellfun(@length,badtrials)));
    
    if sum(cellfun(@length,badtrials)) > 0 % If we have some bad trials...
        for i = 2:length(badtrials) % Loop through prior blocks
            if ~isempty(badtrials{i})
                alldays(i).tt(badtrials{i},:) = []; % Eliminate bad trials
            end
        end
        %Save
        flname2 = sprintf('C:/Users/limblab/Desktop/Mihili_bdf/Processed/%s_%s',monkey,[MO DA YE]);
        save(flname2,'alldays','priors','monkey','MO','DA','YE');

        % Load full version
        [alldays, priors, monkey, MO, DA, YE] = load_processed('Mihili',[FLS{daynum,1} FLS{daynum,2} '2015'],1);
        for i = 2:length(badtrials) % Loop through prior blocks
            if ~isempty(badtrials{i}) 
                alldays(i).tt(badtrials{i},:) = []; % Eliminate bad trials
            end
        end
        %Save
        flname = sprintf('C:/Users/limblab/Desktop/Mihili_bdf/full alldays/full_alldays_%s_%s',monkey,[MO DA YE]);
        save(flname,'alldays','priors','monkey','MO','DA','YE');
        
    end
end
    