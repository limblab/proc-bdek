%% Files to Load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
BRAIN_AREA = 'PMd';
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
            
%             'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             
%             'MrT'   ,'05062013',  2,  [1, 3]    ;...
%             
%             'MrT'   ,'07082013',  2,  [1, 2]    ...

};

% FileName = {'Mihili','08062013',  2,  [1, 2]    };...
% FileName = {'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 3]    ;...
%             'MrT'   ,'07082013',  2,  [1, 2]    ...
% };        

session_limit = 10e10;

%% Do activity and Behavior

[split_indices] = ...
    deal(cell(size(FileName,1),1)); %initialize
[LI,PEAKSPD,REACTTIME,TIMEPEAK,TT]=deal(cell(100,1));
[av_ts,av_rt,dRT_bnd,dTS_bnd,dTT_bnd] = deal(nan(10000,2));
[dRT,dTS,dTT] = deal(nan(10000,1));
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

    %
    pdi =2;
   
    split_indices{daynum} = 1:session_limit:size(alldays(2).tt,1);
    split_indices{daynum}(end) = size(alldays(2).tt,1);
    if(length(split_indices{daynum})==1); split_indices{daynum} = [1 size(alldays(2).tt,1)]; end
    
    ALLDAYS = alldays;
    if strcmp(FileName{daynum,1},'MrT'); numslcs = 10; else numslcs = 5; end
    for section = 1:(length(split_indices{daynum})-1)

        indis = split_indices{daynum}(section):split_indices{daynum}(section+1);
        alldays(2).tt = ALLDAYS(2).tt(indis,:);
                
        liklist = flipud(unique(alldays(2).tt(:,3)));
        for lks = 1:length(liklist)
            LI{counter}{lks} = find(alldays(2).tt(:,3)==liklist(lks));
        end
        
        speed_script;
        
%         react_time(abs(react_time)>500) = nan;
 
        tts = [time_topspeed alldays(2).tt(:,3)];
        rts = [react_time alldays(2).tt(:,3)];
        tss = [topspeed alldays(2).tt(:,3)];
        
        PEAKSPD{counter} = tss;
        REACTTIME{counter} = rts;
        TIMEPEAK{counter} = tts;
   
        dmeanfunc = @(x) nanmean(x(x(:,2)==min(x(:,2)),1))-nanmean(x(x(:,2)==max(x(:,2)),1));
        
        av_rt(counter) = nanmean(react_time); 
        av_ts(counter) = nanmean(topspeed);
        
        dRT(counter,:) = dmeanfunc(rts);
        dTS(counter,:) = dmeanfunc(tss);
        dTT(counter,:) = dmeanfunc(tts);
        
        [dRT_bnd(counter,1), dRT_bnd(counter,2)] = boot_bounds(1000,dmeanfunc,rts,2.5,97.5);
        [dTS_bnd(counter,1), dTS_bnd(counter,2)] = boot_bounds(1000,dmeanfunc,tss,2.5,97.5);
        [dTT_bnd(counter,1), dTT_bnd(counter,2)] = boot_bounds(1000,dmeanfunc,tts,2.5,97.5);
        
%         dRTmed(counter,:) = dmedfunc(rts);
%         dTSmed(counter,:) = dmedfunc(tss);
%         
%         [dRTmed_bnd(counter,1), dRTmed_bnd(counter,2)] = boot_bounds(1000,dmedfunc,rts,2.5,97.5);
%         [dTSmed_bnd(counter,1), dTSmed_bnd(counter,2)] = boot_bounds(1000,dmedfunc,tss,2.5,97.5);

        TTT{counter} = [alldays(2).tt(:,7)-alldays(2).tt(:,6) alldays(2).tt(:,3)];

        counter = counter + 1;
    end
    clc; fprintf('%d/%d\n',daynum,size(FileName,1)); 

    clearvars -except FileName BRAIN_AREA G counter ...
        split_indices LI session_limit loop_ranges av_rt av_ts dRT dTS dRT_bnd ...
        dTS_bnd dTT dTT_bnd PEAKSPD REACTTIME TIMEPEAK ISNew TTT
    
end
av_rt(counter:end) = [];
av_ts(counter:end) = [];

dRT(counter:end,:) = [];
dRT_bnd(counter:end,:) = [];

dTS(counter:end,:) = [];
dTS_bnd(counter:end,:) = [];

dTT(counter:end,:) = [];
dTT_bnd(counter:end,:) = [];

PEAKSPD(counter:end) = [];
REACTTIME(counter:end) = [];
TIMEPEAK(counter:end) = [];

TTT(counter:end) = [];

dRES = dRT;
dRES_bnd = dRT_bnd;


% G = 1:length(PEAKSPD);
% [~,G2] = sortrows(cellfun(@(x) str2double(x),FileName(:,2)));


splitfunc = @(x) {x(x(:,2)==max(x(:,2)),1), x(x(:,2)==min(x(:,2)),1)};

cellPS = cellfun(splitfunc,PEAKSPD,'UniformOutput',0);
cellRT = cellfun(splitfunc,REACTTIME,'UniformOutput',0);
cellTP = cellfun(splitfunc,TIMEPEAK,'UniformOutput',0);
cellttt = cellfun(splitfunc,TTT,'UniformOutput',0)';


[RT,PS,TP] = deal(cell(length(cellPS),1));
 for i = 1:length(cellPS)

     for j = 1:2
         [RT{i}{j}(1), RT{i}{j}(2)] = boot_bounds(1000,@nanmean,cellRT{i}{j},2.5,97.5); 
         [PS{i}{j}(1), PS{i}{j}(2)] = boot_bounds(1000,@nanmean,cellPS{i}{j},2.5,97.5); 
         [TP{i}{j}(1), TP{i}{j}(2)] = boot_bounds(1000,@nanmean,cellTP{i}{j},2.5,97.5); 
         
         RT{i}{j}(3) = nanmean(cellRT{i}{j});
         PS{i}{j}(3) = nanmean(cellPS{i}{j});
         TP{i}{j}(3) = nanmean(cellTP{i}{j});
     end 
 end
%%
% figure; subplot(2,3,1); hold on; 
% for i = 1:length(G2)
%     plot(i*[1 1],RT{G2(i)}{1},'b'); 
%     plot(i*[1 1]+.2,RT{G2(i)}{2},'r'); 
% end
% 
% subplot(2,3,2); hold on; 
% for i = 1:length(G2)
%     plot(i*[1 1],PS{G2(i)}{1},'b'); 
%     plot(i*[1 1]+.2,PS{G2(i)}{2},'r'); 
% end

% %%
% Trialnums = cellfun(@(x) size(x,1),REACTTIME);
% figure; hold on; 
% for i = 1:length(RT)
%     
%     s = metric2markersize(Trialnums(i),Trialnums,[10 40]);
%     plot(RT{i}{1}(3),RT{i}{2}(3),'.','MarkerSize',s);
%     
%     plot(RT{i}{1}(1:2),[1 1]*RT{i}{2}(3));
%     plot([1 1]*RT{i}{1}(3),RT{i}{2}(1:2));
% end
% xl = xlim; yl = ylim;
% xlim([min([xl yl]) max([xl yl])]);
% ylim([min([xl yl]) max([xl yl])]);
% text(mean(xl),mean(yl),sprintf('N = %d -- %d',min(Trialnums),max(Trialnums)),'FontSize',14);
% plot(xl,xl,'k--');
% title('Reaction Time','FontSize',18);

% %%
% Trialnums = cellfun(@(x) size(x,1),REACTTIME);
% figure; hold on; 
% for i = 1:length(dRT)
%    
%     s = metric2markersize(Trialnums(i),Trialnums,[10 40]);
%     plot(dRES_posts(i),dRT(i),'.','MarkerSize',s)
%     
%     plot([1 1]*dRES_posts(i),dRT_bnd(i,:));
%     plot(dRES_bnd_posts(i,:),[1 1]*dRT(i));
% 
% end
% % xl = xlim; yl = ylim;
% % xlim([min([xl yl]) max([xl yl])]);
% % ylim([min([xl yl]) max([xl yl])]);
% % text(mean(xl),mean(yl),
% title(sprintf('N = %d -- %d',min(Trialnums),max(Trialnums)),'FontSize',14);
% % plot(xl,xl,'k--');
% xlabel('\Delta behavioral uncertainty','FontSize',14);
% ylabel('\Delta Reaction Time','FontSize',14);
% title('Reaction Time','FontSize',18);

%%
for i = 1:length(REACTTIME); 
    REACTTIME{i}(REACTTIME{i}(:,1)<0,:) = nan;
end

figure; hold on; 
for i = 1:length(REACTTIME)
    is1 = find(REACTTIME{i}(:,2)==max(REACTTIME{i}(:,2)));
    is2 = find(REACTTIME{i}(:,2)==min(REACTTIME{i}(:,2)));
    
    [y1,x1] = ecdf(REACTTIME{i}(is1,1));
    [y2,x2] = ecdf(REACTTIME{i}(is2,1));
    
    plot(x1,y1,'b','LineWidth',1); 
    plot(x2,y2,'r','LineWidth',1);
    
%     patch([x1' fliplr(x2')],[y1' fliplr(y2')],'k','FaceAlpha',0.2,'EdgeColor','none');
    
end





