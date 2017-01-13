%% Files to Load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
BRAIN_AREA = 'PMd';
% FileName = {'Mihili','07112013',  2,  [1, 2]    ;...
%     
%             'Mihili','07152013',  2,  [1, 2]    ;...
%             'Mihili','07152013',  3,  [1, 2]    ;...
%             
%             'Mihili','07192013',  2,  [1, 2]    ;...
%             'Mihili','07192013',  3,  [1, 2]    ;...
%             'Mihili','07192013',  5,  [1, 2]    ;...
%             
%             'Mihili','08062013',  2,  [1, 2]    ;...
%             'Mihili','08062013',  3,  [1, 2]    ;...
%             
%             'Mihili','08122013',  2,  [1, 2]    ;...
%             'Mihili','08122013',  3,  [1, 2]    ;...
%             
%             'Mihili','08152013',  2,  [1, 2]    ;...
%             
%             'Mihili','07122013',  2,  [1, 2]    ;...
%     
%             'Mihili','08012013',  2,  [1, 2]    ;...
%             'Mihili','08012013',  3,  [1, 2]    ;...
%             'Mihili','08012013',  4,  [1, 2]    ;...
%             
%             'Mihili','08222013',  2,  [1, 2]    ;...
%             'Mihili','08222013',  3,  [1, 2]    ;...
%             
%             'Mihili','09042013',  2,  [1, 2]    ;...
%             
%             'Mihili','09052013',  2,  [1, 2]    ;...
%             
%             'Mihili','09062013',  2,  [1, 2]    ;...
%             
%             'Mihili','09262013',  2,  [1, 2]    ;...
% 
%             'Mihili','10022013',  2,  [1, 2]    ;...
% 
% 
%             'Mihili','10072013',  2,  [1, 2]    ;...
% 
% };

% FileName = {'Mihili'   ,'10082015',  2, [1, 2];...           
%             'Mihili'   ,'10272015',  2, [1, 2];...         
%             'Mihili'   ,'11022015',  2, [1, 2] ...
% };


% 
% % FileName = {'Mihili','08062013',  2,  [1, 2]    };...
FileName = {'MrT'   ,'05042013',  2,  [1, 2]    ;...
            'MrT'   ,'05052013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 3]    ;...
%             'MrT'   ,'05062013',  2,  [2, 3]    ;...
            'MrT'   ,'07082013',  2,  [1, 2]    ...
};        

session_limit = 100;

%% Do activity and Behavior
Act_aligns = {'target','target','target',12};
% Act_wins = {-100:50:250, 300:50:700, -100:50:200};
Act_wins = {[50 250], [300 500], [500 700], [-100 100]};

[FL, FU, FM, RB, BPDS, split_indices,KRATS, PRIS, DDIR, DSTD, DPRC, DVAF] = ...
    deal(cell(size(FileName,1),1)); %initialize
[PLR,VA,LI,PLRb]=deal(cell(100,1));
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

    llist = flipud(unique(alldays(2).tt(:,3))); 
    alldays(2).tt(~ismember(alldays(2).tt(:,3),llist(FileName{daynum,4})),:) = [];
    alldays(1).tt = alldays(1).tt(1,:); alldays(1).tt(:,12) = 1;
%     pdi = 1; speed_script; alldays(1).tt(:,12) = alldays(1).tt(:,6)+react_time./1000;
    pdi = 2; speed_script; alldays(2).tt(:,12) = alldays(2).tt(:,6)+react_time./1000;
%     load('C:\Users\limblab\Desktop\Figure list\Revis_figs\kinematics_new.mat');
%     alldays(2).tt(:,12) = alldays(2).tt(:,6)+RTS{daynum}./1000;

    load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_V-D-RT\\%s_PD90_%s_%s_V-D-RT.mat',...
         BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}));

    PDS = get_PDS_fulldat(fulldat,pi);
    PDmat = [PDS(:,1), PDS(:,2), PDS(:,2), PDS(:,3)];
    
    split_indices{daynum} = 1:session_limit:size(alldays(2).tt,1);
    split_indices{daynum}(end) = size(alldays(2).tt,1);
    if(length(split_indices{daynum})==1); split_indices{daynum} = [1 size(alldays(2).tt,1)]; end
    
    ALLDAYS = alldays;
    for section = 1:(length(split_indices{daynum})-1)
        indis = split_indices{daynum}(section):split_indices{daynum}(section+1);
        alldays(2).tt = ALLDAYS(2).tt(indis,:);

        do_baselines = 1;
        for TIM = 1:length(Act_aligns)
            %%%% DO Cisek analysis and decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%
            fire_align = Act_aligns(TIM);
            fire_ranges = Act_wins(TIM);
            best_PDS = PDmat(:,TIM);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum(~isnan(best_PDS))>10
                Cisek_plot_dotuning_loop;
                prccor = cell2mat(cellfun(@(x) mean(abs(x)<=pi/8),decode_dir,'UniformOutput',0));
                do_baselines = 0;
            else
                decode_dir = {nan; nan};
                decode_std = {nan; nan};
                decodeVAF = {nan; nan};
                prccor = {nan; nan};
                do_baselines = 1;
            end
                
            DDIR{counter}(:,TIM) = decode_dir; 
            DSTD{counter}(:,TIM) = decode_std;
            DVAF{counter}(:,TIM) = decodeVAF;
            DPRC{counter}{1,TIM} = prccor;

            
        end
       counter = counter+1; 
    end

    clearvars -except FL FU FM FileName RB VA BPDS BRAIN_AREA G counter PLR ...
        split_indices std_priors std_likes std_posts LI PLRb session_limit ...
        Indall loop_ranges KRATS PRIS dRES_PMd dRES_M1 ISNew REACTTIME DDIR ...
        DSTD DPRC Act_aligns Act_wins DVAF RTS;
    
end

%%
vaf = cellfun(@(x) cell2mat(x),DVAF,'UniformOutput',0);
figure; hold on; 
[VE,VE1,VE2] = deal(zeros(length(vaf),size(vaf{1},2)));
for i = 1:size(vaf{1},2)    
    subplot(1,size(vaf{1},2),i);
    hold on;
    for j = 1:length(vaf)
        
          ve1 = circ_var(DDIR{j}{1,i});
          ve2 = circ_var(DDIR{j}{2,i});
%           vaf1 = (circ_var(RDS{j}{1})-ve1)./circ_var(RDS{j}{1});
%           vaf2 = (circ_var(RDS{j}{2})-ve2)./circ_var(RDS{j}{2});
%         vaf1 = (1-vaf{j}(1,i)).*TVAR{j}(1);
%         vaf2 = (1-vaf{j}(2,i)).*TVAR{j}(2);
%         
%         maxv = max(TVAR{j});
%         
%         vafnew{j}(1,i) = (maxv - vaf1)./maxv;
%         vafnew{j}(2,i) = (maxv - vaf2)./maxv;
%         
%         plot(vaf{j}(1,i),vaf{j}(2,i),'b.','MarkerSize',20); 
%         plot(vafnew{j}(1,i),vafnew{j}(2,i),'b.','MarkerSize',20); 
          plot(1-ve1,1-ve2,'b.','MarkerSize',20);
%           plot(vaf1,vaf2,'b.','MarkerSize',20);
            
          VE1(j,i) = 1-ve1; 
          VE2(j,i) = 1-ve2; 
          VE(j,i) = (1-ve1)-(1-ve2);

    end
    VE(VE1==1) = NaN;
    VE(VE2==1) = NaN;
    ylim([0 1]); xlim([0 1]);
    plot([0 1],[0 1],'k--');
end

%%
% for daynum = 1:size(FileName,1)
%     fprintf('%d/%d\n',daynum,size(FileName,1));
%     load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_V-D-RT\\%s_PD90_%s_%s_V-D-RT.mat',...
%              'PMd',FileName{daynum,1},FileName{daynum,2}));
%     PDS = get_PDS_fulldat(fulldat,pi/2);
%     numneurs_PMD(daynum,:) = sum(~isnan(PDS(:,end)));
%     
%     load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_V-D-RT\\%s_PD90_%s_%s_V-D-RT.mat',...
%              'M1',FileName{daynum,1},FileName{daynum,2}));
%     PDS = get_PDS_fulldat(fulldat,pi/2);
%     numneurs_M1(daynum,:) = sum(~isnan(PDS(:,end)));
% end

%%
% figure; hold on; 
% [DSTDcat,DPRCcat] = deal(cell(length(DSTD),1));
% for i = goodhighs'%1:length(DSTD)%highnuminds'
%     
%     DSTDcat{i}{1,:} = cell2mat(DSTD{i}(1,:));
%     DSTDcat{i}{2,:} = cell2mat(DSTD{i}(2,:));
%     
%     DPRCcat{i} = cell2mat(DPRC{i});
%     
%     plot(DSTDcat{i}{1}(1,:),'b'); 
%     plot(DSTDcat{i}{2}(1,:),'r');
%     
%     patch([1:21 21:-1:1],[DSTDcat{i}{1}(1,:) fliplr(DSTDcat{i}{2}(1,:))],'k','FaceAlpha',1,'EdgeAlpha',1);
% %         plot(DSTDcat{i}{2}(1,:)-DSTDcat{i}{1}(1,:),'b');
% 
% %     plot(DPRCcat{i}(2,:)./DPRCcat{i}(1,:));
% %     
% %     plot(DPRCcat{i}(1,:),'b');
% %     plot(DPRCcat{i}(2,:),'r');
% % pause
% end

%%
% figure; hold on; 
% for i = 1:size(prcpmd,1)
%     
%     plot(prcpmd(i,1),prcpmd(i,2),'.','Color','b','MarkerSize',...
%         metric2markersize(numneurs_PMD(i),[numneurs_PMD; numneurs_M1],[10 40]));
%     
%     plot(prcm1(i,1),prcm1(i,2),'.','Color','r','MarkerSize',...
%         metric2markersize(numneurs_M1(i),[numneurs_M1; numneurs_PMD],[10 40]));
%     
%     
% end
% plot([0 1],[0 1],'k--');

%%
% figure; hold on; 
% for i = 1:size(prcpmd,1)
%     
%     plot(prcpmd(i,1),prcpmd(i,2),'.','Color','b','MarkerSize',...
%         metric2markersize(dRES(i),dRES,[10 40]));
%     
%     plot(prcm1(i,1),prcm1(i,2),'.','Color','r','MarkerSize',...
%         metric2markersize(dRES(i),dRES,[10 40]));
%     
%     
% end
% plot([0 1],[0 1],'k--');

%%
% figure; hold on; 
% for i = 1:size(prcpmd,1)
%     
%     plot(stdpmd(i,1),stdpmd(i,2),'.','Color','b','MarkerSize',...
%         metric2markersize(numneurs_PMD(i),[numneurs_PMD; numneurs_M1],[10 40]));
%     
%     plot(stdm1(i,1),stdm1(i,2),'.','Color','r','MarkerSize',...
%         metric2markersize(numneurs_M1(i),[numneurs_M1; numneurs_PMD],[10 40]));
%     
%     
% end
% plot([0.1 1.5],[0.1 1.5],'k--');

%%
% figure; hold on; 
% for i = 1:size(prcpmd,1)
%     
% %     if ismember(i,over100PMd)
%     if ismember(i,over50ppmd)
%     plot(stdpmd(i,1),stdpmd(i,2),'.','Color','b','MarkerSize',...
%         20);%metric2markersize(numneurs_PMD(i),[numneurs_PMD; numneurs_M1],[10 40]));
%     end
% %     if ismember(i,over100M1)
%     if ismember(i,over50pm1)
%     plot(stdm1(i,1),stdm1(i,2),'.','Color','r','MarkerSize',...
%          20);%metric2markersize(numneurs_M1(i),[numneurs_M1; numneurs_PMD],[10 40]));
%     end
% end
% plot([0.1 .9],[0.1 0.9],'k--');