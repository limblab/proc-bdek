%% Files to Load
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
%             'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 3]    };

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
            'Mihili','09262013',  2,  [1, 3]    ;...

            'Mihili','10022013',  2,  [1, 2]    ;...
            'Mihili','10022013',  2,  [1, 3]    ;...

            'Mihili','10072013',  2,  [1, 2]    ;...
            'Mihili','10072013',  2,  [1, 3]    ;...
            
            'MrT'   ,'05042013',  2,  [1, 2]    ;...
            'MrT'   ,'05052013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 3]    ;...
};
G = 1:29; G([22,24,26]) = []; G(end-2:end) = [];     


% FileName = {'Mihili','05062014', 2, [1 2]};
% G =1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%   
% FileName = {'Mihili','07192013',  2};
% 
% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10122015',  2        ;...
%             'Mihili','10272015',  2        };
        
% CONTROLS
% FileName = {'Mihili','05062014',  2, [1 2]    ;...
%             'Mihili','10082015',  2, [1 2]    ;...
%             'Mihili','10272015',  2, [1 2]    ;...
%             'Mihili','11022015',  2, [1 2]    };
% G = 1:4;
% FileName = {'Mihili','05062014', 2};
        
BRAIN_AREA = 'M1';
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
    
    clearvars -except FL FU FM FileName RB VA BPDS BRAIN_AREA G;
    
end
days = str2num(cell2mat(FileName(G,2)));
[sortdays,sortinds] = sortrows(days);
G = G(sortinds);

loop_load_behavior; close;
%% Effect mag
[EF,H,CI,NN,AV,AV_bnd,AV_bnd_low,AV_bnd_high,diff_significance,D1D2] = deal(cell(1,1));
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
            
            D1D2{j2}{j}{i} = d1d2;
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
            [~,~,booted_1] = boot_bounds(1000,@(x) nanmean(x(:)),d1_cor,2.5,97.5);
            [~,~,booted_2] = boot_bounds(1000,@(x) nanmean(x(:)),d2_cor,2.5,97.5);
            
            AV_bnd{j2}{j}(i,:) = prctile(booted_2-booted_1,[2.5,97.5]);
            AV_bnd_low{j2}{j}(i,:) = [prctile(booted_1,[2.5,97.5]) mean(booted_1)];
            AV_bnd_high{j2}{j}(i,:) = [prctile(booted_2,[2.5,97.5]),mean(booted_2)];
            
            
            diff_significance{j2}(i,j) = ttest2(booted_1,booted_2,.05,'both');
            
        end
    end
end
%%
for tp = 1:length(AV)
    

av_PD = AV{tp}{1};
av_OD = AV{tp}{2};
av_ORTH = AV{tp}{3};

d_PD = EF{tp}{1}(:,3);%av_PD(:,2)-av_PD(:,1);
d_OD = EF{tp}{2}(:,3);%av_OD(:,2)-av_OD(:,1);
d_ORTH = EF{tp}{3}(:,3);%av_ORTH(:,2)-av_ORTH(:,1);

bnd_PD = EF{tp}{1}(:,1:2);%AV_bnd{1}{1};
bnd_OD = EF{tp}{2}(:,1:2);%AV_bnd{1}{2};
bnd_ORTH = EF{tp}{3}(:,1:2);%AV_bnd{1}{3};


%%
xlm = [floor(10*min(min(dRES_bnd(G,:))))/10 ceil(10*max(max(dRES_bnd(G,:))))/10];
% xlm = [0,0.3];
% xlm = [-3 1.5];
xtic = xlm(1):0.01:xlm(2);
pclr = 'b';

[pPD,S_PD] = polyfit(dRES(G),d_PD(G),1); 
[pOD,S_OD] = polyfit(dRES(G),d_OD(G),1);
[pORTH,S_ORTH] = polyfit(dRES(G),d_ORTH(G),1);

[~,pdes,pdrs] = regress(d_PD(G),[ones(length(G),1) dRES(G)]);
[~,odes,odrs] = regress(d_OD(G),[ones(length(G),1) dRES(G)]);
[~,orthes,orthrs]=  regress(d_ORTH(G),[ones(length(G),1) dRES(G)]);

Q = [zeros(length(G),1); ones(length(G),1)];
combx = [d_PD(G); d_OD(G)];
combdRES = [dRES(G); dRES(G)];

[~,combint] = regress(combx,[combdRES, combdRES.*Q, ones(length(Q),1), Q]);

[Y_PD,DEL_PD] = polyconf(pPD,xtic,S_PD,'predopt','curve');
[Y_ORTH,DEL_ORTH] = polyconf(pORTH,xtic,S_ORTH,'predopt','curve');
[Y_OD,DEL_OD] = polyconf(pOD,xtic,S_OD,'predopt','curve');

R2_PD = 1 - sum((polyval(pPD,dRES(G)) - d_PD(G)).^2)./sum((mean(d_PD(G)) - d_PD(G)).^2);
R2_OD = 1 - sum((polyval(pOD,dRES(G)) - d_OD(G)).^2)./sum((mean(d_OD(G)) - d_OD(G)).^2);
R2_ORTH = 1 - sum((polyval(pORTH,dRES(G)) - d_ORTH(G)).^2)./sum((mean(d_ORTH(G)) - d_ORTH(G)).^2);

upperlim = ceil(max(reshape([bnd_PD(G,:);bnd_OD(G,:);bnd_ORTH(G,:)],[],1)));
lowerlim = floor(min(reshape([bnd_PD(G,:);bnd_OD(G,:);bnd_ORTH(G,:)],[],1)));

clrs = distinguishable_colors(6); clrs(4,:) = [0 0 0]; clrs(1,:) = [0 1 1];
figure;
subplot(1,3,1); hold on;
plot(dRES(G),d_PD(G),'.','Color',pclr,'MarkerSize',18); 
for i = 1:length(G)
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_PD(G(i)),bnd_PD(G(i),:),pclr,.5);
end
% for i = 1:length(G); plot(dRES(G(i)),d_PD(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pPD,xlm),'k','LineWidth',2);
% plot(xtic,Y_PD+DEL_PD,'b');
% plot(xtic,Y_PD-DEL_PD,'b');
patch([xtic fliplr(xtic)],[Y_PD+DEL_PD fliplr(Y_PD-DEL_PD)],'k','FaceAlpha',0.1,'EdgeAlpha',0)

plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pPD(1),R2_PD),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pPD(1),R2_PD),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('SD Neurons','FontSize',18);
xlabel('\Delta SE_r_e_s','FontSize',14)
ylabel('\Delta FR','FontSize',18);

subplot(1,3,2); hold on;
plot(dRES(G),d_ORTH(G),'.','Color',pclr,'MarkerSize',18); 
for i = 1:length(G)
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_ORTH(G(i)),bnd_ORTH(G(i),:),pclr,.5);
end
% for i = 1:length(G); plot(dRES(G(i)),d_ORTH(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pORTH,xlm),'k','LineWidth',2);
patch([xtic fliplr(xtic)],[Y_ORTH+DEL_ORTH fliplr(Y_ORTH-DEL_ORTH)],'k','FaceAlpha',0.1,'EdgeAlpha',0)
plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pORTH(1),R2_ORTH),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pORTH(1),R2_ORTH),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('ORTH Neurons','FontSize',18);

subplot(1,3,3); hold on;
plot(dRES(G),d_OD(G),'.','Color',pclr,'MarkerSize',20); 
for i = 1:length(G)
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_OD(G(i)),bnd_OD(G(i),:),pclr,.5);
%     plot(dRES(G(i)),d_OD(G(i)),'.','Color',c(round(2.78*i),:),'MarkerSize',18);
end

% for i = 1:length(G); plot(dRES(G(i)),d_OD(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pOD,xlm),'k','LineWidth',2);
patch([xtic fliplr(xtic)],[Y_OD+DEL_OD fliplr(Y_OD-DEL_OD)],'k','FaceAlpha',0.1,'EdgeAlpha',0)
plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pOD(1),R2_OD),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pOD(1),R2_OD),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('OD Neurons','FontSize',18);

fitE_PD = sqrt(diag((S_PD.R)\inv(S_PD.R'))./S_PD.normr.^2./S_PD.df);
fitE_ORTH = sqrt(diag((S_ORTH.R)\inv(S_ORTH.R'))./S_ORTH.normr.^2./S_ORTH.df);
fitE_OD = sqrt(diag((S_OD.R)\inv(S_OD.R'))./S_OD.normr.^2./S_OD.df);

slopes_TYPES(tp,:) = [pPD(1) pORTH(1) pOD(1)];
Rsquareds(tp,:) = [R2_PD R2_ORTH R2_OD];
slope_LOWS(tp,:) = [pdes(2,1) orthes(2,1) odes(2,1)];
slope_HIGHS(tp,:) = [pdes(2,2) orthes(2,2) odes(2,2)];
slopes_different(tp,1) = combint(2,1) > 0 | combint(2,2) < 0;
%[fitE_PD(1) fitE_ORTH(1) fitE_OD(1)];
end
%%
figure; hold on; 
pd_order = [1 3 2];
for j = 1:length(pd_order)
    subplot(1,3,j); hold on; 
    for i = 1:length(G)
        plot([i i],AV_bnd_low{1}{pd_order(j)}(G(i),1:2),'b','LineWidth',2);
        plot([i i]+.1,AV_bnd_high{1}{pd_order(j)}(G(i),1:2),'r','LineWidth',2);
        plot(i,AV_bnd_low{1}{pd_order(j)}(G(i),3),'b.','MarkerSize',7);
        plot(i+0.1,AV_bnd_high{1}{pd_order(j)}(G(i),3),'r.','MarkerSize',7);
    end
end

%%
figure; hold on; 

subplot(1,3,1); hold on; 
plot(1:length(G),d_PD(G),'b.','MarkerSize',18);
plot(repmat(1:length(G),2,1),bnd_PD(G,:)','b','LineWidth',2);
% plot(G(notrend),d_PD(G(notrend)),'r.','MarkerSize',18);
title('SD','FontSize',18);

subplot(1,3,2); hold on; 
plot(1:length(G),d_ORTH(G),'b.','MarkerSize',18);
plot(repmat(1:length(G),2,1),bnd_ORTH(G,:)','b','LineWidth',2);
% plot(G(notrend),d_ORTH(G(notrend)),'r.','MarkerSize',18);
title('ORTH','FontSize',18);

subplot(1,3,3); hold on; 
plot(1:length(G),d_OD(G),'b.','MarkerSize',18);
plot(repmat(1:length(G),2,1),bnd_OD(G,:)','b','LineWidth',2);
% plot(G(notrend),d_OD(G(notrend)),'r.','MarkerSize',18);
title('OD','FontSize',18);    


%%
allses = 1:length(G);
notrend = [9:17];%[6:14];
figure; hold on; 

subplot(1,3,1); hold on; 
plot(allses,d_PD(G),'b.','MarkerSize',18);
plot(repmat(allses,2,1),bnd_PD(G,:)','b','LineWidth',2);
plot(allses(notrend),d_PD(G(notrend)),'r.','MarkerSize',18);
title('SD','FontSize',18);

subplot(1,3,2); hold on; 
plot(allses,d_ORTH(G),'b.','MarkerSize',18);
plot(repmat(allses,2,1),bnd_ORTH(G,:)','b','LineWidth',2);
plot(allses(notrend),d_ORTH(G(notrend)),'r.','MarkerSize',18);
title('ORTH','FontSize',18);

subplot(1,3,3); hold on; 
plot(allses,d_OD(G),'b.','MarkerSize',18);
plot(repmat(allses,2,1),bnd_OD(G,:)','b','LineWidth',2);
plot(allses(notrend),d_OD(G(notrend)),'r.','MarkerSize',18);
title('OD','FontSize',18);    
        
%%
allses = 1:length(G);
notrend = nan;%[6:14];
figure; hold on; 
plot(allses,dRES(G),'b.','MarkerSize',18);
plot(repmat(allses,2,1),dRES_bnd(G,:)','b','LineWidth',2);
plot(allses(notrend),dRES(G(notrend)),'r.','MarkerSize',18);

%%
PDODORTH = cell(3,1);
for i = 1:length(D1D2) %Time Bin
    
    for j = 1:length(D1D2{i}) % PD OD ORTH
        
        totpad = cell(length(G),1);
        starti = 1;
        for k = 1:length(G)
            padded = nan(size(D1D2{i}{j}{G(k)},1),sum(cellfun(@(x) size(x,2),D1D2{i}{j}(G))));

            padded(:,starti:(starti+size(D1D2{i}{j}{G(k)},2)-2)) = D1D2{i}{j}{G(k)}(:,1:(end-1));
            padded(:,end) = D1D2{i}{j}{G(k)}(:,end);
            
            totpad{k} = padded; 
            starti = starti+size(D1D2{i}{j}{G(k)},2)-1;
        end
            
        catpad = vertcat(totpad{:});
        
        [l,h,rb] = boot_bounds(1000,difffunc,catpad,2.5,97.5);
        m = nanmean(rb);
        
        PDODORTH{j}(:,i) = [m;l;h];
    end
end
      %%  
figure; hold on;
c2p = {'b','r','g'};
xlis = [0:100:800 1000:100:1400];
x2p = (xlis(1:end-1)+xlis(2:end))./2; x2p(9) = [];
for i = 1:3
    plot(x2p,PDODORTH{i}(1,:),'.-','Color',c2p{i});
    for j = 1:length(x2p)
%     patch([x2p fliplr(x2p)],[PDODORTH{i}(2,:) fliplr(PDODORTH{i}(3,:))],c2p{i},'EdgeColor','none');
    plot(x2p(j)*[1 1]+.1*[i i],[PDODORTH{i}(2,j); PDODORTH{i}(3,j)],c2p{i});
    end
end
        