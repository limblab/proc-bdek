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


% FileName = {'Mihili','08062013',  2,  [1, 2]    };...
FileName = {'MrT'   ,'05042013',  2,  [1, 2]    ;...
            'MrT'   ,'05052013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 3]    ;...
%             'MrT'   ,'05062013',  2,  [2, 3]    ;...
            'MrT'   ,'07082013',  2,  [1, 2]    ...
};        

session_limit = 10e10;

%% Do activity and Behavior

[FL, FU, FM, RB, BPDS, split_indices,PLRb,KRATS, PRIS] = ...
    deal(cell(size(FileName,1),1)); %initialize
[PLR,VA,LI,PLRb]=deal(cell(100,1));
[std_priors,std_likes,std_posts] = deal(nan(10000,2));
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
%     pdi = 1; speed_script; alldays(1).tt(:,12) = alldays(1).tt(:,6)+react_time./1000;
    pdi = 2; speed_script; alldays(2).tt(:,12) = alldays(2).tt(:,6)+react_time./1000;
   
%     addreact = REACTTIME{daynum}(:,1)./1000; addreact(isnan(addreact)) = 0;
%     alldays(2).tt(:,6) = alldays(2).tt(:,6) + addreact;
%     alldays(2).tt(:,3) = Indall{daynum};
    
    if 1 %exist(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\%s_PD90_%s_%s_full.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}),'file')
        load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_V-D-RT\\%s_PD90_%s_%s_V-D-RT.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}));
%         load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_LD\\%s_PD90_%s_%s_LD.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}));

        PDS = get_PDS_fulldat(fulldat,pi);
        alltuned = PDS;
        pref = 2*ones(size(pref));
    else
        brain_area = BRAIN_AREA;
        CO_index = 1;
        reachdir_col = 10;
        loop_alignmentsPD = {12};
%         loop_alignmentsPD = {12};
        loop_rangesPD = {[0 200]};%,[-100 100]};
%         loop_rangesPD = {[-200 0]};
        [~,pref,alltuned] = tuning_types_PD_func(alldays,brain_area,CO_index,reachdir_col,loop_alignmentsPD,loop_rangesPD);
        save(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\%s_PD90_%s_%s_full',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}),'pref','alltuned');
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
    pretarget_plot_PDOD_limited; close;
    
    % Fill variables
    FL{daynum} = lows;
    FU{daynum} = highs;
    FM{daynum} = mids;
    
    RB{daynum} = rbs;
%     VA{daynum} = voif_all;
    BPDS{daynum} = best_PDS;
    voif_all(cellfun(@isempty,voif_all)) = [];
    
    split_indices{daynum} = 1:session_limit:size(voif_all{1}{1},1);
    split_indices{daynum}(end) = size(voif_all{1}{1},1);
    if(length(split_indices{daynum})==1); split_indices{daynum} = [1 size(voif_all{1}{1},1)]; end
    
    [~,~,fitkrats,fitpris] = behavior_fit_circ(alldays(2).tt(:,[2 9 10 3]));
    
    ALLDAYS = alldays;
    if strcmp(FileName{daynum,1},'MrT'); numslcs = 10; else numslcs = 5; end
    for section = 1:(length(split_indices{daynum})-1)

        indis = split_indices{daynum}(section):split_indices{daynum}(section+1);
        alldays(2).tt = ALLDAYS(2).tt(indis,:);
                
        liklist = flipud(unique(alldays(2).tt(:,3)));
        for lks = 1:length(liklist)
            LI{counter}{lks} = find(alldays(2).tt(:,3)==liklist(lks));
        end
        [PLR{counter},PLRb{counter}] = behavior_fit_circ(alldays(2).tt(:,[2 9 10 3]));%,fitkrats,fitpris);

        std_priors(counter,:) = 1./sqrt([PLR{counter}{1}(1) PLR{counter}{2}(1)]);
        std_likes(counter,:)  = 1./sqrt([PLR{counter}{1}(2) PLR{counter}{2}(2)]);
        std_posts(counter,:)  = 1./sqrt([PLR{counter}{1}(3) PLR{counter}{2}(3)]);
        
        for j = 1:length(voif_all)
            for k = 1:length(voif_all{j})  
                VA{counter}{j}{k} = voif_all{j}{k}(indis,:);
            end
        end
        PRIS{counter} = fitpris;
        KRATS{counter} = fitkrats;
        counter = counter + 1;
    end
    clc; fprintf('%d/%d\n',daynum,size(FileName,1)); 

    clearvars -except FL FU FM FileName RB VA BPDS BRAIN_AREA G counter PLR ...
        split_indices std_priors std_likes std_posts LI PLRb session_limit ...
        Indall loop_ranges KRATS PRIS dRES_PMd dRES_M1 ISNew REACTTIME;
    
end
VA(cellfun(@isempty,VA)) = [];
PLR(cellfun(@isempty,PLR)) = [];
PLRb(cellfun(@isempty,PLRb)) = [];
LI(cellfun(@isempty,LI)) = [];
std_priors(isnan(sum(std_priors,2)),:) = [];
std_likes(isnan(sum(std_likes,2)),:) = [];
std_posts(isnan(sum(std_posts,2)),:) = [];

G = 1:length(VA);
[~,G2] = sortrows(cellfun(@(x) str2double(x),FileName(:,2)));

dRES_prior = diff(std_priors,[],2);
dRES_likes = diff(std_likes,[],2);
dRES_posts = diff(std_posts,[],2);

dRES_bnd_prior = cell2mat(cellfun(@(x) prctile(1./sqrt(x{2}(:,1))-1./sqrt(x{1}(:,1)),[2.5 97.5]),PLRb,'UniformOutput',0));
dRES_bnd_likes = cell2mat(cellfun(@(x) prctile(1./sqrt(x{2}(:,2))-1./sqrt(x{1}(:,2)),[2.5 97.5]),PLRb,'UniformOutput',0));
dRES_bnd_posts = cell2mat(cellfun(@(x) prctile(1./sqrt(x{2}(:,3))-1./sqrt(x{1}(:,3)),[2.5 97.5]),PLRb,'UniformOutput',0));

tanslopes = cell2mat(cellfun(@(x) [x{1}(:,2)./sum(x{1}(:,[1 2]),2) x{2}(:,2)./sum(x{2}(:,[1 2]),2)],PLR,'UniformOutput',0));
tanslopes_bnd = cellfun(@(x) [x{1}(:,2)./sum(x{1}(:,[1 2]),2) x{2}(:,2)./sum(x{2}(:,[1 2]),2)],PLRb,'UniformOutput',0);
dslopes = diff(tanslopes,[],2);
dslopes_bnd = cell2mat(cellfun(@(x) prctile(x(:,2)-x(:,1),[2.5 97.5]),tanslopes_bnd,'UniformOutput',0));

dRES = dRES_posts;
dRES_bnd = dRES_bnd_posts;

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
            
            D1D2{j2}{j}{i} = d1d2; %D1D2{TIME}{PD_OD_ORTH}{DAYS}
%             ed = [d1d2; e1e2];
            
            e = (nanmean(d2_cor)-nanmean(d1_cor));
            e(isinf(e)) = NaN;
            
%             [l,h] = boot_bounds(1000,@nanmean,e,2.5,97.5);

            difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1))-nanmean(x(x(:,end)==1,1:end-1)));
%             difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1)))-nanmean(nanmean(x(x(:,end)==1,1:end-1)));
            
%             difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1))-nanmean(x(x(:,end)==1,1:end-1)))./...
%                             nanmean(nanmean(x(x(:,end)==3,1:end-1))-nanmean(x(x(:,end)==4,1:end-1)));
                            
                        
%             av_estim{i}(j,j2) = abs(nanmean(nanmean(ed(ismember(ed(:,end),[3 4]),1:end-1))));
            
            
            [l,h,~,rb] = boot_bounds(1000,difffunc,d1d2,2.5,97.5);
            
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

%% Alternate Method
if 0
dfunc = @(x) nanmean(x(x(:,end)==2,1))-nanmean(x(x(:,end)==1,1));
difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1))-nanmean(x(x(:,end)==1,1:end-1)));
splitfunc = @(x) {find(x(:,end)==1),find(x(:,end)==2)};
[Ld,Md,Hd] = deal(zeros(3,length(D1D2)));
adiffs = cell(3,1);
[DDmeanall] = deal(zeros(3,length(D1D2)));
E = cell(1,length(D1D2));
for tim = 1:length(D1D2)
    fprintf('%d/%d\n',tim,length(D1D2));
    for pdodorth = 1:3
        split_ses = D1D2{tim}{pdodorth}';

        S = cell(length(split_ses),1);
        for i = 1:length(split_ses)
            inds = splitfunc(D1D2{tim}{pdodorth}{i});
            
            [A,B] = meshgrid(inds{1},inds{2});
            c = cat(2,A',B');
            d = reshape(c,[],2);

            DD = D1D2{tim}{pdodorth}{i}(d(:,2),:) - D1D2{tim}{pdodorth}{i}(d(:,1),:);
%             S{i} = typecast(reshape(DD(~isnan(DD)),[],1),'int8');
            S{i} = sparse(reshape(DD(~isnan(DD)),[],1));
            adiffs{pdodorth}(i,tim) = full(nanmean(S{i}));
            E{tim}{pdodorth}(i,3) = full(nanmean(S{i}));
            se = full(nanmean(S{i})./sqrt(length(S{i})));
            E{tim}{pdodorth}(i,1:2) = [E{tim}{pdodorth}(i,3)-se E{tim}{pdodorth}(i,3)+se];
%             DDcat{i}{pdodorth,tim} = reshape(DD(~isnan(DD)),[],1);
%             dcheck{i}(pdodorth,tim) = difffunc2(split_ses{i});
        end
        DDmeanall(pdodorth,tim) = mean(vertcat(S{:}));
        
    end
end
clear S 
end
%% Plot Correlations
avdiffs = cell(3,1);
[statsPD,statsOD,statsORTH,statsODORTH] = deal(zeros(length(AV),4));
for tp = 1:length(AV)
    
avdiffs{1}(:,tp) = EF{tp}{1}(:,3);
avdiffs{2}(:,tp) = EF{tp}{3}(:,3);
avdiffs{3}(:,tp) = EF{tp}{2}(:,3);

% avdiffs{1}(:,tp) = diff(AV{tp}{1},[],2);
% avdiffs{2}(:,tp) = diff(AV{tp}{3},[],2);
% avdiffs{3}(:,tp) = diff(AV{tp}{2},[],2);

av_PD = AV{tp}{1};
av_OD = AV{tp}{2};
av_ORTH = AV{tp}{3};

d_PD = EF{tp}{1}(:,3);
d_OD = EF{tp}{2}(:,3);
d_ORTH = EF{tp}{3}(:,3);

d_ODORTH = [EF{tp}{2}(:,3); EF{tp}{3}(:,3)];

% missings = find(isnan(sum([d_PD d_OD d_ORTH],2)));
% dRES(missings) = [];
% d_PD(missings) = [];
% d_OD(missings) = [];
% d_ORTH(missings) = [];


bnd_PD = EF{tp}{1}(:,1:2);
bnd_OD = EF{tp}{2}(:,1:2);
bnd_ORTH = EF{tp}{3}(:,1:2);
bnd_ODORTH = [EF{tp}{2}(:,1:2); EF{tp}{3}(:,1:2)];

% d_PD = av_PD(:,2)-av_PD(:,1);
% d_OD = av_OD(:,2)-av_OD(:,1);
% d_ORTH = av_ORTH(:,2)-av_ORTH(:,1);
% 
% bnd_PD = AV_bnd{1}{1};
% bnd_OD = AV_bnd{1}{2};
% bnd_ORTH = AV_bnd{1}{3};

N_PD = NN{tp}{1};
N_OD = NN{tp}{2};
N_ORTH = NN{tp}{3};
N_ODORTH = [NN{tp}{2}; NN{tp}{3}];

totN = vertcat(NN{1}{:});
size_range = [10 40];%[10 40];

% s_PD = metric2markersize(N_PD,totN,size_range);
% s_OD = metric2markersize(N_OD,totN,size_range);
% s_ORTH = metric2markersize(N_ORTH,totN,size_range);
% s_ODORTH = metric2markersize(N_ODORTH,totN,size_range);

Gcat = [G (G+length(dRES))];
dREScat = [dRES; dRES];
dRES_bndcat = [dRES_bnd; dRES_bnd];

%%
xlm = [floor(10*min(min(dRES_bnd(G,:))))/10 ceil(10*max(max(dRES_bnd(G,:))))/10];
% xlm = [0,0.3];
% xlm = [-3 1.5];
xtic = xlm(1):0.01:xlm(2);
pclr = 'b';

[pPD,S_PD] = polyfit(dRES(G),d_PD(G),1); 
[pOD,S_OD] = polyfit(dRES(G),d_OD(G),1);
[pORTH,S_ORTH] = polyfit(dRES(G),d_ORTH(G),1);
[pODORTH,S_ODORTH] = polyfit(dREScat(Gcat),d_ODORTH(Gcat),1);

[~,pdes,pdrs,~,statsPD(tp,:)] = regress(d_PD(G),[ones(length(G),1) dRES(G)]);
[~,odes,odrs,~,statsOD(tp,:)] = regress(d_OD(G),[ones(length(G),1) dRES(G)]);
[~,orthes,orthrs,~,statsORTH(tp,:)]=  regress(d_ORTH(G),[ones(length(G),1) dRES(G)]);
[~,odorthes,odorthrs,~,statsODORTH(tp,:)]=  regress(d_ODORTH(Gcat),[ones(length(Gcat),1) dREScat(Gcat)]);

Q = [zeros(length(G),1); ones(length(G),1)];
combx = [d_PD(G); d_OD(G)];
combdRES = [dRES(G); dRES(G)];

[~,combint] = regress(combx,[combdRES, combdRES.*Q, ones(length(Q),1), Q]);

[Y_PD,DEL_PD] = polyconf(pPD,xtic,S_PD,'predopt','curve');
[Y_ORTH,DEL_ORTH] = polyconf(pORTH,xtic,S_ORTH,'predopt','curve');
[Y_OD,DEL_OD] = polyconf(pOD,xtic,S_OD,'predopt','curve');
[Y_ODORTH,DEL_ODORTH] = polyconf(pODORTH,xtic,S_ODORTH,'predopt','curve');

R2_PD = 1 - sum((polyval(pPD,dRES(G)) - d_PD(G)).^2)./sum((mean(d_PD(G)) - d_PD(G)).^2);
R2_OD = 1 - sum((polyval(pOD,dRES(G)) - d_OD(G)).^2)./sum((mean(d_OD(G)) - d_OD(G)).^2);
R2_ORTH = 1 - sum((polyval(pORTH,dRES(G)) - d_ORTH(G)).^2)./sum((mean(d_ORTH(G)) - d_ORTH(G)).^2);
R2_ODORTH = 1 - sum((polyval(pODORTH,dREScat(Gcat)) - d_ODORTH(Gcat)).^2)./...
                sum((mean(d_ODORTH(Gcat)) - d_ODORTH(Gcat)).^2);

upperlim = ceil(max(reshape([bnd_PD(G,:);bnd_OD(G,:);bnd_ORTH(G,:)],[],1)));
lowerlim = floor(min(reshape([bnd_PD(G,:);bnd_OD(G,:);bnd_ORTH(G,:)],[],1)));
% OD/ORTH
figure; subplot(1,3,3); hold on; 
for i = 1:length(Gcat)
    if i > length(G); dotcolor = 'b'; else dotcolor = 'g'; end
    plot(dREScat(Gcat(i)),d_ODORTH(Gcat(i)),'.','Color',dotcolor,'MarkerSize',...
            metric2markersize(N_ODORTH(Gcat(i)),N_ODORTH,size_range)); 
    plot_cross(dREScat(Gcat(i)),dRES_bndcat(Gcat(i),:),d_ODORTH(Gcat(i)),bnd_ODORTH(Gcat(i),:),pclr,.5);
end
plot(xlm,polyval(pODORTH,xlm),'k','LineWidth',2);
patch([xtic fliplr(xtic)],[Y_ODORTH+DEL_ODORTH fliplr(Y_ODORTH-DEL_ODORTH)],'k','FaceAlpha',0.1,'EdgeAlpha',0)
plot(xlm,[0 0],'k');
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pODORTH(1),R2_ODORTH),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('ORTH+OD Neurons','FontSize',18);
xlabel(sprintf('N = [%d %d]',min(N_ODORTH(G)),max(N_ODORTH(G))),'FontSize',18);
ylabel('\Delta FR','FontSize',18);


%

clrs = distinguishable_colors(6); clrs(4,:) = [0 0 0]; clrs(1,:) = [0 1 1];
figure;
subplot(1,3,1); hold on;
for i = 1:length(G)
    plot(dRES(G(i)),d_PD(G(i)),'.','Color',pclr,'MarkerSize',metric2markersize(N_PD(G(i)),N_PD,size_range)); 
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
xlabel(sprintf('N = [%d %d]',min(N_PD(G)),max(N_PD(G))),'FontSize',18);
%xlabel('\Delta SE_r_e_s','FontSize',14)
ylabel('\Delta FR','FontSize',18);

subplot(1,3,2); hold on;
for i = 1:length(G)
    plot(dRES(G(i)),d_ORTH(G(i)),'.','Color',pclr,'MarkerSize',metric2markersize(N_ORTH(G(i)),N_ORTH,size_range)); 
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
xlabel(sprintf('N = [%d %d]',min(N_ORTH(G)),max(N_ORTH(G))),'FontSize',18);

subplot(1,3,3); hold on;
for i = 1:length(G)
    plot(dRES(G(i)),d_OD(G(i)),'.','Color',pclr,'MarkerSize',metric2markersize(N_OD(G(i)),N_OD,size_range)); 
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_OD(G(i)),bnd_OD(G(i),:),pclr,.5);
%     plot(dRES(G(i)),d_OD(G(i)),'.','Color',c(round(2.78*i),:),'MarkerSize',18);
end
plot(mean(xlm),0,'r.','MarkerSize',size_range(1));
plot(mean(xlm),0,'r.','MarkerSize',size_range(2));

% for i = 1:length(G); plot(dRES(G(i)),d_OD(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pOD,xlm),'k','LineWidth',2);
patch([xtic fliplr(xtic)],[Y_OD+DEL_OD fliplr(Y_OD-DEL_OD)],'k','FaceAlpha',0.1,'EdgeAlpha',0)
plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pOD(1),R2_OD),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pOD(1),R2_OD),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('OD Neurons','FontSize',18);
xlabel(sprintf('N = [%d %d]',min(N_OD(G)),max(N_OD(G))),'FontSize',18);

fitE_PD = sqrt(diag((S_PD.R)\inv(S_PD.R'))./S_PD.normr.^2./S_PD.df);
fitE_ORTH = sqrt(diag((S_ORTH.R)\inv(S_ORTH.R'))./S_ORTH.normr.^2./S_ORTH.df);
fitE_OD = sqrt(diag((S_OD.R)\inv(S_OD.R'))./S_OD.normr.^2./S_OD.df);

slopes_TYPES(tp,:) = [pPD(1) pORTH(1) pOD(1)];
Rsquareds(tp,:) = [R2_PD R2_ORTH R2_OD];
slope_LOWS(tp,:) = [pdes(2,1) orthes(2,1) odes(2,1)];
slope_HIGHS(tp,:) = [pdes(2,2) orthes(2,2) odes(2,2)];
slopes_different(tp,1) = combint(2,1) > 0 | combint(2,2) < 0;
%[fitE_PD(1) fitE_ORTH(1) fitE_OD(1)];
if length(AV)>1; close; end
% 
%%
if length(AV)> 1 && tp==length(AV)
    figure; hold on;
    xranges = {-100:100:200};%, 1000+(-100:100:400)};
    xcents = cell2mat(cellfun(@(x) .5*(x(1:end-1) + x(2:end)),xranges,'UniformOutput',0));
    pmarks = {'.-','s-','o-'};
    pols = {'b','g','r'};
    for j = 1:3
        plot(xcents+2*j,slopes_TYPES(:,j),pmarks{j},'Color',pols{j});
        for i = 1:length(xcents)
            plot([1 1]*xcents(i)+2*j,[slope_LOWS(i,j) slope_HIGHS(i,j)],'Color',pols{j});
        end
    end
    
    figure; hold on; 
    for j = 1:3
        plot(xcents+2*j,Rsquareds(:,j),pmarks{j},'Color',pols{j});
    end
end

end
%    
