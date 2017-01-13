%%
% Set up list of all days for inclusion
% FileName = {'Mihili','07112013',  2        ;...
%             'Mihili','07152013', [2,3]     ;...
%             'Mihili','07192013', [2,3,4]   ;...
%             'Mihili','08062013', [2,3]     ;...
%             'Mihili','08122013', [2,3]     ;...
%             'Mihili','08152013',  2        ;...
%             'MrT'   ,'05042013',  2        ;...
%             'MrT'   ,'05052013',  2        ;...
%             'MrT'   ,'05062013',  2        };
        
% FileName = {'Mihili','09282015', 2};
%   
% FileName = {'Mihili','07112013',  2,  [1, 2]    ;...
%             'Mihili','07152013',  2,  [1, 2]    ;...
%             'Mihili','07192013',  2,  [1, 2]    ;...
%             'Mihili','08062013',  2,  [1, 2]    ;...
%             'Mihili','08122013',  2,  [1, 2]    ;...
%             'Mihili','08152013',  2,  [1, 2]    ;...
%             'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 3]    };
        
        
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
            
            'MrT'   ,'05042013',  2,  [1, 2]    ;...
            'MrT'   ,'05052013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 3]    };
%   
% FileName = {'Mihili','07192013',  2};
% 
% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10122015',  2        ;...
%             'Mihili','10272015',  2        };
        
% CONTROLS
% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10272015',  2        ;...
%             'Mihili','11022015',  2        };

% FileName = {'Mihili','05062014', 2};
        
BRAIN_AREA = 'PMd';
[FL, FU, FM, RB, VA, BPDS] = deal(cell(size(FileName,1),1)); %initialize
for daynum = 1:size(FileName,1)
    
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2},0); 
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

    if exist(sprintf('%s_PD90_%s_%s_full.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}),'file')
        load(sprintf('%s_PD90_%s_%s_full.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}));
    else
        brain_area = BRAIN_AREA;
        CO_index = 1;
        reachdir_col = 10;
        loop_alignmentsPD = {'target','go'};
        loop_rangesPD = {[600 800],[50 250]};
        [~,pref,alltuned] = tuning_types_PD_func(alldays,brain_area,CO_index,reachdir_col,loop_alignmentsPD,loop_rangesPD);
        save(sprintf('%s_PD90_%s_%s_full',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}),'pref','alltuned');
    end
    
    prs = nan(size(alltuned)); 
    bpds = nan(size(alltuned,1),1);
    for i = 1:size(prs,1);
        if ~isnan(pref(i))
            prs(i,pref(i)) = alltuned(i,pref(i)); 
            bpds(i) = alltuned(i,pref(i));
        end
    end
    best_PDS = bpds;%prs(:,2);
    pretarget_plot_PDOD; close;
    
    % Fill variables
    FL{daynum} = lows;
    FU{daynum} = highs;
    FM{daynum} = mids;
    
    RB{daynum} = rbs;
    VA{daynum} = voif_all;
    BPDS{daynum} = best_PDS;
    
    clearvars -except FL FU FM FileName RB VA BPDS BRAIN_AREA;
    
end
loop_load_behavior; close;

%%
if 0
figure; hold on; 
col = {'r','b'};

for t = 1:3
subplot(1,3,t); hold on; 
    for PO = 1:2

        for i = 1:length(FM)

            plot_cross(FM{i}{PO}(1,t+1),[FL{i}{PO}(1,t+1) FU{i}{PO}(1,t+1)],...
                       FM{i}{PO}(2,t+1),[FL{i}{PO}(2,t+1) FU{i}{PO}(2,t+1)],col{PO});
                   
        end
    end
    
    for i = 1:length(FM) 
        plot([FM{i}{1}(1,t+1) FM{i}{2}(1,t+1)],[FM{i}{1}(2,t+1) FM{i}{2}(2,t+1)],'k');
    end
        
        
    xl = xlim;
    plot(xl,xl,'k--');
    axis square; 
    xlabel('deltaFR (Low Unc)','FontSize',16);
    ylabel('deltaFR (High Unc)','FontSize',16);
    title(sprintf('Time bin: %d',t),'FontSize',18);
end
%%
figure; hold on; 
col = {'r','b'};

for t = 1:3
subplot(1,3,t); hold on; 
    for like = 1:2

        for i = 1:length(FM)

            plot_cross(FM{i}{1}(like,t+1),[FL{i}{1}(like,t+1) FU{i}{1}(like,t+1)],...
                       FM{i}{2}(like,t+1),[FL{i}{2}(like,t+1) FU{i}{2}(like,t+1)],col{like});
                   
        end
    end
    xl = xlim;
    plot(xl,xl,'k--');
    axis square; 
end
%%
figure; hold on; 
col = {'r','b'};

[lb,ub,ms] = deal(zeros(1,2));
for t = 1:3
subplot(1,3,t); hold on; 

    for i = 1:length(FM)
        for PO = 1:2

            % RB{session}{PD/OD}{time window}{uncertainty condition}
            unc_diff = RB{i}{PO}{t+1}{2}-RB{i}{PO}{t+1}{1};
            unc_diff = unc_diff./(mean(RB{i}{PO}{t+1}{1}));
            
            
%             unc_diff = RB{i}{PO}{t+1}{1}./RB{i}{PO}{t+1}{2};

            
            
            
            lb(PO) = prctile(unc_diff,2.5);
            ub(PO) = prctile(unc_diff,97.5);
            ms(PO) = mean(unc_diff);
                   
        end
        plot([i i],[lb(1) ub(1)],'b','LineWidth',2); plot([i i]+0.1,[lb(2) ub(2)],'r','LineWidth',2);
        plot(i,ms(1),'b.','MarkerSize',10); plot(i+.1,ms(2),'r.','MarkerSize',10);
%         plot_cross(ms(1),[lb(1) ub(1)],ms(2),[lb(2) ub(2)],'b');
    end 
        
    xl = xlim;
    plot(xl,xl,'k--');
    axis square; 
    xlabel('Uncertainty effect (PD)','FontSize',16);
    ylabel('Uncertainty effect (OD)','FontSize',16);
    title(sprintf('Time bin: %d',t),'FontSize',18);
end
%%
figure; hold on; 
for j = 1:9
sess = j;
subplot(9,2,1+(j-1)*2); hold on; 
linehist(RB{sess}{1}{3}{1},20,'b'); linehist(RB{sess}{1}{3}{2},20,'r'); 
subplot(9,2,2*j); hold on; 
linehist(RB{sess}{2}{3}{1},20,'b'); linehist(RB{sess}{2}{3}{2},20,'r'); 
end
%%
% figure; hold on;
[SL_bnds,SH_bnds] = deal(cell(2,1));
for i = 1:length(VA)
    
    for j = 1:length(VA{i})
        
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            clc; fprintf('day: %d/%d\nunit type:%d/%d\ntime: %d/%d\n',...
                    i,length(VA),j,length(VA{i}),j2,length(VA{i}{j})-1);
            
                [SL_bnds{j}{i}{j2}(:,1), SL_bnds{j}{i}{j2}(:,2)] = ...
                      boot_bounds(1000,@nanmean,VA{i}{j}{k}(LI{i}{1},:),2.5,97.5);
                  
                [SH_bnds{j}{i}{j2}(:,1), SH_bnds{j}{i}{j2}(:,2)] = ...
                      boot_bounds(1000,@nanmean,VA{i}{j}{k}(LI{i}{2},:),2.5,97.5);
                   
        end
    end
end
%% T-test between uncertainty conditions
[CIs] = deal(cell(3,1));
figure; hold on;
pdodmaxes{1} = cellfun(@(x) max(abs(x{1}(1,:))),FM);
pdodmaxes{2} = cellfun(@(x) max(abs(x{2}(1,:))),FM);
for i = 1:length(VA)
    
    for j = 1:length(VA{i})
        
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            clc; fprintf('day: %d/%d\nunit type:%d/%d\ntime: %d/%d\n',...
                    i,length(VA),j,length(VA{i}),j2,length(VA{i}{j})-1);
            
            d1 = reshape(VA{i}{j}{k}(LI{i}{1},:),[],1); d1(isnan(d1))=[];
            d2 = reshape(VA{i}{j}{k}(LI{i}{2},:),[],1); d2(isnan(d2))=[];  
            
            [h,p,CIs{j2}{j}(i,:)] = ttest2(d2,d1,0.05,'both','unequal');
            CIs{j2}{j}(i,:) = CIs{j2}{j}(i,:)./max([pdodmaxes{1}(i) pdodmaxes{2}(i)]);
        end        
    end
end
%

for i = 1:length(CIs)
    
    subplot(1,length(CIs),i);
    for j = 1:size(CIs{i}{1},1)
        
        if j<=6; cfc = 'b'; else cfc = 'r'; end
        plot_cross(mean(CIs{i}{1}(j,:)),CIs{i}{1}(j,:),mean(CIs{i}{2}(j,:)),CIs{i}{2}(j,:),cfc)
        
    end
    axis equal;
    axis square;
    cura = axis;
    plot([-max(abs(cura)) max(abs(cura))],[0 0],'k--')
    plot([0 0],[-max(abs(cura)) max(abs(cura))],'k--')
    xlabel('uncertainty change in FR (PD)','FontSize',16);
    ylabel('uncertainty change in FR (OD)','FontSize',16);
    box off;
end
end
if 0
%% Percentage of neurons
[PON,DS,PONBOOT] = deal(cell(3,1));
bnum = 100;
cutoff_count = 1;
for i = 1:length(VA)
    for j = 1:length(VA{i}) 
        for j2 = 1:(length(VA{i}{j})-1)
            clc; fprintf('%d/%d\n',i,length(VA));
            k = j2 +1;
     
            cs1 = VA{i}{j}{k}(LI{i}{1},:);
            cs2 = VA{i}{j}{k}(LI{i}{2},:);
            
            d1 = nanmean(cs1);
            d2 = nanmean(cs2);
            
            [~,b1ind] = bootstrp(bnum,[],LI{i}{1});
            [~,b2ind] = bootstrp(bnum,[],LI{i}{2});
            
            pover_boot = zeros(1,bnum);
            for b = 1:bnum
                
                good1 = sum(~isnan(cs1(b1ind(:,b),:)))>cutoff_count;
                good2 = sum(~isnan(cs2(b2ind(:,b),:)))>cutoff_count;
                
                good_all = good1&good2;
                
                d1boot = nanmean(cs1(b1ind(:,b),:));
                d2boot = nanmean(cs2(b2ind(:,b),:));

%                 pover_boot(b) = sum(d2boot>d1boot)./(sum(~isnan(d1boot+d2boot)));
                pover_boot(b) = sum(d2boot(good_all)>d1boot(good_all))./sum(good_all);
%                 pover_boot(b) = sum(d2boot(good_all)>d1boot(good_all))./sum(good_all);
%                 
            end
            PONBOOT{j2}{j}(i,:) = pover_boot;

            DS{j2}{j}{i} = [d1;d2];
            
            good_tot = sum(~isnan(cs1))>10 & sum(~isnan(cs2))>10;
            pover = sum(d2(good_tot)>d1(good_tot))./sum(good_tot);
            punder = sum(d1(good_tot)>d2(good_tot))./sum(good_tot);
%             pover = sum(d2>d1)./(sum(~isnan(d1+d2)));
%             punder = sum(d2<d1)./(sum(~isnan(d1+d2)));
            
            PON{j2}{j}(i,:) = [pover punder];
                        
        end
    end
end
%
figure; hold on;
for i = 1:length(PON)
    
    subplot(1,length(PON),i); hold on; 
    
    totalsPD = horzcat(DS{i}{1}{:});
    totalsOD = horzcat(DS{i}{2}{:});
    
    totalpercPD = sum(totalsPD(2,:)>totalsPD(1,:))./(sum(~isnan(sum(totalsPD))));
    totalpercOD = sum(totalsOD(2,:)>totalsOD(1,:))./(sum(~isnan(sum(totalsOD))));
    
    for j = 1:size(PON{i}{1},1)
        
        if j<=6 
            fprintf('%d - M1\n',j);
            cfc = 'b.'; 
        else
            fprintf('%d - M2\n',j);
            cfc = 'bo';
        end
        
%         plot(PON{i}{1}(j,1),PON{i}{2}(j,1),cfc,'MarkerSize',8);
        
        plot(mean(PONBOOT{i}{1}(j,:)),mean(PONBOOT{i}{2}(j,:)),cfc,'MarkerSize',10);
        plot_cross(mean(PONBOOT{i}{1}(j,:)),prctile(PONBOOT{i}{1}(j,:),[2.5,97.5]),...
                   mean(PONBOOT{i}{2}(j,:)),prctile(PONBOOT{i}{2}(j,:),[2.5,97.5]),'k');
%         plot_cross(PON{i}{1}(j,1),prctile(PONBOOT{i}{1}(j,:),[2.5,97.5]),...
%                    PON{i}{2}(j,1),prctile(PONBOOT{i}{2}(j,:),[2.5,97.5]),'k');
%         pause; 
        
    end
    
    plot(totalpercPD,totalpercOD,'r.','MarkerSize',5);
    
    xlim([0 1]); ylim([0 1]);
    plot([0 1],[0 1],'k--');
    plot([0 1],[.5 .5],'k-');
    plot([.5 .5],[0 1],'k-');
    axis square;
    xlabel('% \DeltaFR_h_i_g_h > \DeltaFR_l_o_w (PD)','FontSize',16);
    ylabel('% \DeltaFR_h_i_g_h > \DeltaFR_l_o_w (OD)','FontSize',16);
    box off;
end
end
%%
if 0
%% T-test of neurons
[PON,H,CI] = deal(cell(3,1));
for i = 1:length(VA)
    for j = 1:length(VA{i}) 
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            
            d1 = VA{i}{j}{k}(LI{i}{1},:);
            d2 = VA{i}{j}{k}(LI{i}{2},:);
            
            nd1 = sum(~isnan(d1));
            nd2 = sum(~isnan(d2));

            for u = 1:size(d1,2);
                
                if nd1(u)>5 && nd2(u)>5
                    [H{j2}{j}{i}(u,:),~,CI{j2}{j}{i}(u,:)] = ttest2(d2(:,u),d1(:,u),0.05,'both','unequal');
                else
                    H{j2}{j}{i}(u,:) = NaN;
                end
%                 if nd1(u)>10 && nd2(u)>10
%                     [~,hreturn] = ranksum(d2(:,u),d1(:,u));
%                     H{j2}{j}{i}(u,:) = double(hreturn);
%                 else
%                     H{j2}{j}{i}(u,:) = NaN;
%                 end
            end
            validunits = ~isnan(H{j2}{j}{i});
            nvalidunits = sum(validunits);
            
            incs = nanmean(d2)>nanmean(d1);
            decs = nanmean(d2)<nanmean(d1);
                       
            pover = nansum(H{j2}{j}{i}==1 & incs'==1);%./nvalidunits;
            punder = nansum(H{j2}{j}{i}==1 & decs'==1);%./nvalidunits;
            pnon = sum(nd1>5 & nd2>5);%nansum(H{j2}{j}{i}==0)./nvalidunits;
 
            PON{j2}{j}(i,:) = [pover punder pnon];
                        
        end
    end
end
%%
figure; hold on;
    
ns = PON{1}{1}(:,3)+PON{1}{2}(:,3);
num2size = @(x) 5+floor(25*x/max(ns));
msizes = num2size(ns);
[SDPLOT, ODPLOT] = deal(zeros(size(PON{1}{1},1),length(PON)));
for i = 1:length(PON)
    subplot(1,length(PON),i); hold on; 
    for j = 1:size(PON{i}{1},1)
        if j<=6 
            cfc = '.'; 
            colr = 'b';
        else
            cfc = '.';
            colr = 'r';
        end
        
        SDPLOT(j,i) = (PON{i}{1}(j,1)-PON{i}{1}(j,2))./PON{i}{1}(j,3);
        ODPLOT(j,i) = (PON{i}{2}(j,1)-PON{i}{2}(j,2))./PON{i}{2}(j,3);

        plot(SDPLOT(j,i),ODPLOT(j,i),cfc,'Color',colr,'MarkerSize',msizes(j));
%         plot(PON{i}{2}(j,1),PON{i}{2}(j,2),cfc,'Color','r','MarkerSize',5);
    end

%     xlim([0 .5]); ylim([0 .5]);
    plot([-.5 .5],[-.5 .5],'k--');
    plot([-.5 .5],[0 0],'k-');
    plot([0 0],[-.5 .5],'k-');
    axis square;
    xlabel('%increase - %decrease (SD)','FontSize',16);
    ylabel('%increase - %decrease (OD)','FontSize',16);
    box off;
end
%%
figure; hold on;
    
colrs = {'b','r'};
[SDPLOT, ODPLOT] = deal(zeros(size(PON{1}{1},1),length(PON)));


for i = 1:length(PON)
    subplot(1,length(PON),i); hold on; 
    
    MMSD(1,:) = sum(PON{i}{1}(1:6,:));
    MMOD(1,:) = sum(PON{i}{2}(1:6,:));
    
    MMSD(2,:) = sum(PON{i}{1}(7:9,:));
    MMOD(2,:) = sum(PON{i}{2}(7:9,:));
    
    for j = 1:size(MMSD,1)

        SDPLOT(j,i) = (MMSD(j,1)-MMSD(j,2))./MMSD(j,3);
        ODPLOT(j,i) = (MMOD(j,1)-MMOD(j,2))./MMOD(j,3);

        plot(SDPLOT(j,i),ODPLOT(j,i),'.','Color',colrs{j});
    end

%     xlim([0 .5]); ylim([0 .5]);
    plot([-.5 .5],[-.5 .5],'k--');
    plot([-.5 .5],[0 0],'k-');
    plot([0 0],[-.5 .5],'k-');
    axis square;
    xlabel('%increase - %decrease (SD)','FontSize',16);
    ylabel('%increase - %decrease (OD)','FontSize',16);
    box off;
end
%% Wilcoxon Test
[PON,DS,H,CI] = deal(cell(3,1));
figure; hold on;
for i = 1:length(VA)
    for j = 1:length(VA{i}) 
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            
            d1 = VA{i}{j}{k}(LI{i}{1},:);
            d2 = VA{i}{j}{k}(LI{i}{2},:);

            for u = 1:size(d1,2);
                [H{j2}{j}{i}(u,:),~,CI{j2}{j}{i}(u,:)] = kstest2(d2(:,u),d1(:,u),0.05,'both','unequal');
            end
            validunits = ~isnan(H{j2}{j}{i});
            
%             pover = nansum(H{j2}{j}{i}==1 & CI{j2}{j}{i}(:,1)>0)./sum(validunits);
            pover = nansum(CI{j2}{j}{i}(:,1) >= 0)./sum(validunits);
%             punder = nansum(H{j2}{j}{i}==1 & CI{j2}{j}{i}(:,2)<0)./sum(validunits);
            punder = nansum(CI{j2}{j}{i}(:,2) <= 0)./sum(validunits);
%             pnon = nansum(H{j2}{j}{i}(validunits)==0)./sum(validunits);
            pnon = nansum(prod(CI{j2}{j}{i},2)<0)./sum(validunits);
            
            PON{j2}{j}(i,:) = [pover+pnon punder pnon];
                        
        end
    end
end
%
for i = 1:length(PON)
    
    subplot(1,length(PON),i); hold on; 
    for j = 1:size(PON{i}{1},1)
        
        if j<=6 
%             fprintf('%d - M1\n',j);
            cfc = 'k.'; 
        else
%             fprintf('%d - M2\n',j);
            cfc = 'ko';
        end
        
        plot(PON{i}{1}(j,1),PON{i}{2}(j,1),cfc,'MarkerSize',5);
%         plot(PON{i}{1}(j,1),PON{i}{1}(j,2),'b.','MarkerSize',5);
%         plot(PON{i}{2}(j,1),PON{i}{2}(j,2),'r.','MarkerSize',5);
%         plot(PON{i}{1}(j,1),PON{i}{2}(j,1),cfc,'MarkerSize',5);
        
    end
    
%     plot(totalpercPD,totalpercOD,'r.','MarkerSize',5);
    
    xlim([0 1]); ylim([0 1]);
    plot([0 1],[0 1],'k--');
    plot([0 1],[.5 .5],'k-');
    plot([.5 .5],[0 1],'k-');
    axis square;
    xlabel('% \DeltaFR_h_i_g_h > \DeltaFR_l_o_w (PD)','FontSize',16);
    ylabel('% \DeltaFR_h_i_g_h > \DeltaFR_l_o_w (OD)','FontSize',16);
    box off;
end
%% Percentage of trials
[CON,PON_t] = deal(cell(3,1));
figure; hold on;
for i = 1:length(VA)
    for j = 1:length(VA{i}) 
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            clc; fprintf('day: %d/%d\nunit type:%d/%d\ntime: %d/%d\n',...
                    i,length(VA),j,length(VA{i}),j2,length(VA{i}{j})-1);
            
            d1 = nanmean(VA{i}{j}{k}(LI{i}{1},:),2)';
            d2 = nanmean(VA{i}{j}{k}(LI{i}{2},:),2)';

            CON{j2}{i}{1}(j,:) = d1;
            CON{j2}{i}{2}(j,:) = d2;
                        
        end
    end
end
%
for i = 1:length(CON)
    
    subplot(1,length(CON),i); hold on; 
    
    for j = 1:size(CON{i}{1},2)
        
        if j<=6 
            fprintf('%d - M1\n',j);
            cfc = 'k.'; 
        else
            fprintf('%d - M2\n',j);
            cfc = 'ko';
        end
        
        plot(CON{i}{j}{1}(1,:),CON{i}{j}{1}(2,:),'b.','MarkerSize',5);
        plot(CON{i}{j}{2}(1,:),CON{i}{j}{2}(2,:),'r.','MarkerSize',5);
        
    end
    
%     plot(totalpercPD,totalpercOD,'r.','MarkerSize',5);
    
%     xlim([0 1]); ylim([0 1]);
%     plot([0 1],[0 1],'k--');
%     plot([0 1],[.5 .5],'k-');
%     plot([.5 .5],[0 1],'k-');
    axis square;
    xlabel('% \DeltaFR_h_i_g_h > \DeltaFR_l_o_w (PD)','FontSize',16);
    ylabel('% \DeltaFR_h_i_g_h > \DeltaFR_l_o_w (OD)','FontSize',16);
    box off;
end
%% Uncertainty change for all neurons
[FON,quad,quads,perces] = deal(cell(3,1));
figure; hold on;
for i = 1:length(VA)
    for j = 1:length(VA{i}) 
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            clc; fprintf('day: %d/%d\nunit type:%d/%d\ntime: %d/%d\n',...
                    i,length(VA),j,length(VA{i}),j2,length(VA{i}{j})-1);
            
            d1 = nanmean(VA{i}{j}{k}(LI{i}{1},:));
            d2 = nanmean(VA{i}{j}{k}(LI{i}{2},:));
            
            fchange = d2-d1;
            
            FON{j2}{j}{i} = fchange;
            
                        
        end
    end
end
%
for i = 1:length(FON)
    
    subplot(1,length(FON),i); hold on; 
    for j = 1:size(FON{i}{1},2)
        
        if j<=6 
            fprintf('%d - M1\n',j);
            cfc = 'k.'; 
        else
            fprintf('%d - M2\n',j);
            cfc = 'kx';
        end
        
%         pdodangs = atan2(FON{i}{2}{j},FON{i}{1}{j});
%         quad{i}{j} = (ceil((mod(pdodangs+2*pi,2*pi))./pi*2));
   
%         plot(FON{i}{1}{j},FON{i}{2}{j},cfc,'MarkerSize',4);  
    end
    
    PDmonkey1 = horzcat(FON{i}{1}{1:6});
    PDmonkey2 = horzcat(FON{i}{1}{7:9});

    ODmonkey1 = horzcat(FON{i}{2}{1:6});
    ODmonkey2 = horzcat(FON{i}{2}{7:9});

    [x(1,:), y(1,:)] = linehist(PDmonkey1,200,'b');
    [x(2,:), y(2,:)] = linehist(PDmonkey2,200,'c');
    [x(3,:), y(3,:)] = linehist(ODmonkey1,200,'r'); 
    [x(4,:), y(4,:)] = linehist(ODmonkey2,200,'m');
    
    cla; 
    
    plot(x(1,:),y(1,:)./sum(y(1,:)),'b');
    plot(x(2,:),y(2,:)./sum(y(2,:)),'c');
    plot(x(3,:),y(3,:)./sum(y(3,:)),'r');
    plot(x(4,:),y(4,:)./sum(y(4,:)),'m');
%     quads{i} = horzcat(quad{i}{:});
    
%     for k = 1:4
%         perces{i}(k) = sum(quads{i}==k)./sum(~isnan(quads{i}));
%     end
    
    cura = axis;
%     plot([-max(abs(cura)) max(abs(cura))],[0 0],'k--')
%     plot([0 0],[-max(abs(cura)) max(abs(cura))],'k--')
    xlabel('uncertainty change in FR (PD)','FontSize',16);
    ylabel('uncertainty change in FR (OD)','FontSize',16);
    box off;
    axis square;
    xlabel('PD','FontSize',16);
    ylabel('OD','FontSize',16);
%     title('%d | %d | %d | %d',perces{i}(1),perces{i}(2),perces{i}(3),perces{i}(4));
    box off;
end
%%
Monkey1_PD = {vertcat(SL_bnds{1}{1:6}), vertcat(SH_bnds{1}{1:6})};
Monkey1_OD = {vertcat(SL_bnds{2}{1:6}), vertcat(SH_bnds{2}{1:6})};
        
Monkey2_PD = {vertcat(SL_bnds{1}{7:9}), vertcat(SH_bnds{1}{7:9})};
Monkey2_OD = {vertcat(SL_bnds{2}{7:9}), vertcat(SH_bnds{2}{7:9})};
%%
Monkey1_PD = {cell2mat(vertcat(SL_bnds{1}{1:6})), cell2mat(vertcat(SH_bnds{1}{1:6}))};
Monkey1_OD = {cell2mat(vertcat(SL_bnds{2}{1:6})), cell2mat(vertcat(SH_bnds{2}{1:6}))};
        
Monkey2_PD = {cell2mat(vertcat(SL_bnds{1}{7:9})), cell2mat(vertcat(SH_bnds{1}{7:9}))};
Monkey2_OD = {cell2mat(vertcat(SL_bnds{2}{7:9})), cell2mat(vertcat(SH_bnds{2}{7:9}))};

figure; hold on;
for i = 1:length(Monkey1_PD{1})

    plot_cross(mean(Monkey1_PD{1}(i,3:4),2),[Monkey1_PD{1}(i,3:4)],...
               mean(Monkey1_PD{2}(i,3:4),2),[Monkey1_PD{2}(i,3:4)],'b');

    plot_cross(mean(Monkey1_OD{1}(i,3:4),2),[Monkey1_OD{1}(i,3:4)],...
               mean(Monkey1_OD{2}(i,3:4),2),[Monkey1_OD{2}(i,3:4)],'r');
end
%% Low Vs High discriminability
trace = cell(length(VA),1);
for d = 1:length(VA)
    for i = 1:length(VA{d}{1})
%         trace{d}(i,1) = nanmean(nanmean(VA{d}{1}{i}(LI{d}{1},:),2)-nanmean(VA{d}{2}{i}(LI{d}{1},:),2));
%         trace{d}(i,2) = nanmean(nanmean(VA{d}{1}{i}(LI{d}{2},:),2)-nanmean(VA{d}{2}{i}(LI{d}{2},:),2));
        
        trace{d}(i,1) = nanmean(VA{d}{1}{i}(LI{d}{1}));
        trace{d}(i,2) = nanmean(VA{d}{1}{i}(LI{d}{2}));
       
    end
end
figure; hold on; 
% plot([-4 10],[-4 10],'k--');
plot([0 length(VA{1}{1})+1],[0 0],'k--');
for d = 1:length(VA)
    plot(trace{d}(:,1),trace{d}(:,2),'k');
%     plot(trace{d}(:,1)-trace{d}(:,2),'k');
%     plot(trace{d}(:,1),'b');
%     plot(trace{d}(:,2),'r');
    
    for t = 1:length(VA{d}{1})
        plot(trace{d}(t,1),trace{d}(t,2),'.','Color',1-[1 1 1]*t*(1/length(VA{d}{1})),'MarkerSize',20);
%         plot(t,trace{d}(t,1)-trace{d}(t,2),'.','MarkerSize',20);
    end
%     pause; 
    
end
plot([-5 15],[-5 15],'k--')
%% 
trace = cell(length(VA),1);
for d = 1:length(VA)
    for i = 1:4
        trace{d}(i,1) = nanmean(nanmean(VA{d}{1}{i}(LI{d}{1},:),2)-nanmean(VA{d}{2}{i}(LI{d}{1},:),2));
        trace{d}(i,2) = nanmean(nanmean(VA{d}{1}{i}(LI{d}{2},:),2)-nanmean(VA{d}{2}{i}(LI{d}{2},:),2));
    end
end
figure; hold on; 
% plot([-4 10],[-4 10],'k--');
for d = 1:length(VA)
%     plot(trace{d}(2:end,1),'b');
    plot(trace{d}(1:end,2)-trace{d}(1:end,1));
    for t = 1:3
        plot(trace{d}(t+1,1),trace{d}(t+1,2),'.','Color',1-[1 1 1]*t*(1/3),'MarkerSize',20);
    end
%     pause; cla;
end
end
%% Effect mag
[EF,H,CI,NN,AV] = deal(cell(1,1));
av_estim = cell(length(VA),1);
for i = 1:length(VA) % Days
    fprintf('%d/%d\n',i,length(VA));
    for j = 1:length(VA{i}) % PD/OD
        for j2 = 1:(length(VA{i}{j})) % Time Bin
           
            k = j2;
            
            d1 = VA{i}{j}{k}(LI{i}{1},:);
            d2 = VA{i}{j}{k}(LI{i}{2},:);
%             e1 = VA{i}{j+2}{k}(LI{i}{1},:);
%             e2 = VA{i}{j+2}{k}(LI{i}{2},:);

            nd1 = sum(~isnan(d1));
            nd2 = sum(~isnan(d2));

            drop1 = nan(1,size(d1,2)); drop1(nd1>5) = 1; drop1 = repmat(drop1,size(d1,1),1);
            drop2 = nan(1,size(d2,2)); drop2(nd2>5) = 1; drop2 = repmat(drop2,size(d2,1),1);
            
            d1_cor = d1.*drop1;
            d2_cor = d2.*drop2;
%             e1_cor = e1.*drop1;
%             e2_cor = e2.*drop2;
            
            HS = nan(1,size(d1_cor,2));
            for u = 1:size(d1_cor,2)
                hret = ttest2(d2_cor(:,u),d1_cor(:,u),0.05,'both','unequal');
                if hret == 1; HS(u) = 1; 
                elseif hret == 0; HS(u) = 0;
                end
            end
            
            incs = HS == 1 & (nanmean(d2_cor)>nanmean(d1_cor)); 
            decs = HS == 1 & (nanmean(d2_cor)<nanmean(d1_cor)); 
            nneurs = sum(~isnan(HS));
            
%             d1_cor = d1_cor.*repmat(HS,size(d1_cor,1),1);
%             d2_cor = d2_cor.*repmat(HS,size(d2_cor,1),1);
            
            d1d2 = [d1_cor ones(size(d1_cor,1),1); d2_cor 2*ones(size(d2_cor,1),1)];
%             e1e2 = [e1_cor 3*ones(size(e1_cor,1),1); e2_cor 4*ones(size(e2_cor,1),1)];
            
%             ed = [d1d2; e1e2];
            
            e = (nanmean(d2_cor)-nanmean(d1_cor));
            e(isinf(e)) = NaN;
            
%             [l,h] = boot_bounds(1000,@nanmean,e,2.5,97.5);

            difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1))-nanmean(x(x(:,end)==1,1:end-1)));
            
%             difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1))-nanmean(x(x(:,end)==1,1:end-1)))./...
%                             nanmean(nanmean(x(x(:,end)==3,1:end-1))-nanmean(x(x(:,end)==4,1:end-1)));
                            
                        
%             av_estim{i}(j,j2) = abs(nanmean(nanmean(ed(ismember(ed(:,end),[3 4]),1:end-1))));
            
            
            [l,h,rb] = boot_bounds(1000,difffunc,d1d2,2.5,97.5);
            
%             m = nanmean(e);
            m = nanmean(rb);
            
            EF{j2}{j}(i,:) = [l,h,m];
            %EF{TIMEBIN}{PDOD}(DAY,:)
            NN{j2}{j}(i,:) = sum(nd1>5 | nd2>5);
            AV{j2}{j}(i,:) = [nanmean(d1_cor(:)) nanmean(d2_cor(:))];
                        
        end
    end
end
%%
figure; hold on;

ns = NN{1}{1}+NN{1}{2};
num2size = @(x) 5+floor(25*x/max(ns));
msizes = num2size(ns);

[SDPLOT, ODPLOT] = deal(zeros(size(EF{1}{1},1),length(EF)));
for i = 1:length(EF)
    subplot(1,length(EF),i); hold on; 
    for j = 1:size(EF{i}{1},1)
        if j<=6 
            cfc = '.'; 
            colr = 'b';
        else
            cfc = '.';
            colr = 'r';
        end
        
%         SDPLOT(j,i) = EF{i}{1}(j,3);
%         ODPLOT(j,i) = EF{i}{2}(j,3);

        X_ax = EF{i}{1}(j,:);%./max(av_estim{j}(1,:));
        Y_ax = EF{i}{2}(j,:);%./max(av_estim{j}(2,:));

        plot_cross(X_ax(3),X_ax(1:2),Y_ax(3),Y_ax(1:2),colr);
%         
%         plot(X_ax(3),Y_ax(3),'.','MarkerSize',msizes(j));
        
%         plot([j j],EF{i}{1}(j,1:2),'b');
%         plot([j j]+0.1,EF{i}{2}(j,1:2),'r');
%         plot(SDPLOT(j,i),ODPLOT(j,i),cfc,'Color',colr,'MarkerSize',msizes(j));
%         plot(PON{i}{2}(j,1),PON{i}{2}(j,2),cfc,'Color','r','MarkerSize',5);

%         ax_mult = 5;
%         plot([-1 1].*ax_mult,[0 0],'k-');
%         plot([0 0],[-1 1].*ax_mult,'k-');
%         pause; 
        
    end

%     xlim([-5 5]); ylim([-5 5]);
    ax_mult = 5;
    plot([-1 1].*ax_mult,[-1 1].*ax_mult,'k--');
    plot([-1 1].*ax_mult,[0 0],'k-');
    plot([0 0],[-1 1].*ax_mult,'k-');
    axis square;
    xlabel('FR_H - FR_L (SD)','FontSize',16);
    ylabel('FR_H - FR_L (OD)','FontSize',16);
    box off;
end

%%
dpths = zeros(length(VA),2);
for i = 1:length(VA)
    
    f = @(x) nanmean(reshape(x,[],1));
    
    dH = f(VA{i}{1}{2}(LI{i}{2},:)) - f(VA{i}{2}{2}(LI{i}{2},:));
    dL = f(VA{i}{1}{2}(LI{i}{1},:)) - f(VA{i}{2}{2}(LI{i}{1},:));
        
    dpths(i,:) = [nanmean(dL),nanmean(dH)];
end

    