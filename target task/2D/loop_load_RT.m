%% Files to Load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
FileName = {'Mihili','07112013',  2,  [1, 2]    ;...
    
            'Mihili','07152013',  2,  [1, 2]    ;...
            'Mihili','07152013',  3,  [1, 2]    ;...
            
            'Mihili','07192013',  2,  [1, 2]    ;...
            'Mihili','07192013',  3,  [1, 2]    ;...
            'Mihili','07192013',  5,  [1, 2]    ;...
            
            'Mihili','08062013',  2,  [1, 2]    ;...
            'Mihili','08062013',  3,  [1, 2]    ;...
            
            'Mihili','08122013',  2,  [1, 2]    ;...
            'Mihili','08122013',  3,  [1, 2]    ;...
            
            'Mihili','08152013',  2,  [1, 2]    ;...
            
            'Mihili','07122013',  2,  [1, 2]    ;...
    
            'Mihili','08012013',  2,  [1, 2]    ;...
            'Mihili','08012013',  3,  [1, 2]    ;...
            'Mihili','08012013',  4,  [1, 2]    ;...
            
            'Mihili','08222013',  2,  [1, 2]    ;...
            'Mihili','08222013',  3,  [1, 2]    ;...
            
            'Mihili','09042013',  2,  [1, 2]    ;...
            
            'Mihili','09052013',  2,  [1, 2]    ;...
            
            'Mihili','09062013',  2,  [1, 2]    ;...
            
            'Mihili','09262013',  2,  [1, 2]    ;...

            'Mihili','10022013',  2,  [1, 2]    ;...


            'Mihili','10072013',  2,  [1, 2]    ;...
            
            'MrT'   ,'05042013',  2,  [1, 2]    ;...
            'MrT'   ,'05052013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 3]    ;...
            'MrT'   ,'07082013',  2,  [1, 2]    ...

};    

session_limit = 10e10;

%% Do activity and Behavior

[REACTTIME,PEAKSPEED,TIME2TARG,REWARDRATE] = ...
    deal(cell(size(FileName,1),1)); %initialize

counter = 1;
for daynum = 1:size(FileName,1)
   
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2},1); 
    alldays(1).tt(isnan(alldays(1).tt(:,3)),3) = alldays(1).tt(find(isnan(alldays(1).tt(:,3)))-1,3);
    if isfield(alldays,'bdfM')
        BDF = alldays(1).bdfM;
    elseif isfield(alldays,'kin')
        BDF = alldays(1).kin;
    end
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    
    alldays(2).tt(alldays(2).tt(:,3)>100,:) = [];

%     alldays(2).tt = alldays(2).tt(ISNew{daynum},:);
    
    llist = flipud(unique(alldays(2).tt(:,3))); 
    alldays(2).tt(~ismember(alldays(2).tt(:,3),llist(FileName{daynum,4})),:) = [];
% %     
    pdi = 1; speed_script; alldays(1).tt(:,12) = alldays(1).tt(:,6)+react_time./1000;
    pdi = 2; speed_script; alldays(2).tt(:,12) = alldays(2).tt(:,6)+react_time./1000;
   
    REACTTIME{daynum,:} = [alldays(2).tt(:,12)-alldays(2).tt(:,6) alldays(2).tt(:,3)];
    PEAKSPEED{daynum,:} = [topspeed, alldays(2).tt(:,3)];
    TIME2TARG{daynum,:} = [alldays(2).tt(:,7)-alldays(2).tt(:,6) alldays(2).tt(:,3)];
    REWARDRATE{daynum,:} = [double(alldays(2).tt(:,8)==32) alldays(2).tt(:,3)];
    
end
G = 1:size(FileName,1);
[~,G2] = sortrows(cellfun(@(x) str2double(x),FileName(:,2)));

%%

figure; hold on; 
for i = 1:length(GM)
    is1 = find(REACTTIME{i}(:,2)==max(REACTTIME{i}(:,2)));
    is2 = find(REACTTIME{i}(:,2)==min(REACTTIME{i}(:,2)));
    
    [y1,x1] = ecdf(REACTTIME{i}(is1,1));
    [y2,x2] = ecdf(REACTTIME{i}(is2,1));
    
    plot(x1,y1,'b','LineWidth',1); 
    plot(x2,y2,'r','LineWidth',1);
    
    title('React Time');    
end

%%
figure; hold on; 
for i = 1:length(GM)
    is1 = find(PEAKSPEED{i}(:,2)==max(PEAKSPEED{i}(:,2)));
    is2 = find(PEAKSPEED{i}(:,2)==min(PEAKSPEED{i}(:,2)));
    
    [y1,x1] = ecdf(PEAKSPEED{i}(is1,1));
    [y2,x2] = ecdf(PEAKSPEED{i}(is2,1));
    
    plot(x1,y1,'b','LineWidth',1); 
    plot(x2,y2,'r','LineWidth',1);
    
    title('Peak Speed');    
end

%%
figure; hold on; 
for i = 1:length(GM)
    is1 = find(TIME2TARG{i}(:,2)==max(TIME2TARG{i}(:,2)));
    is2 = find(TIME2TARG{i}(:,2)==min(TIME2TARG{i}(:,2)));
    
    [y1,x1] = ecdf(TIME2TARG{i}(is1,1)-REACTTIME{i}(is1,1));
    [y2,x2] = ecdf(TIME2TARG{i}(is2,1)-REACTTIME{i}(is2,1));
    
    plot(x1,y1,'b','LineWidth',1); 
    plot(x2,y2,'r','LineWidth',1);
    
    title('TTT');    
end


%%
figure; hold on; 
for i = 1:length(GM)
    is1 = find(PEAKSPEED{i}(:,2)==max(PEAKSPEED{i}(:,2)));
    is2 = find(PEAKSPEED{i}(:,2)==min(PEAKSPEED{i}(:,2)));
    
    peaksttest(i,:) = ttest2(PEAKSPEED{i}(is1,1),PEAKSPEED{i}(is2,1));
    
end

%%
rp = zeros(length(GM),2);
for i = 1:length(GM)
    is1 = find(REWARDRATE{i}(:,2)==max(REWARDRATE{i}(:,2)));
    is2 = find(REWARDRATE{i}(:,2)==min(REWARDRATE{i}(:,2)));
    
    rp(i,1) = round(100*mean(REWARDRATE{i}(is1,1)));
    rp(i,2) = round(100*mean(REWARDRATE{i}(is2,1)));
    
end
