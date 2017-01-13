%% Files to Load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
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
%             'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 3]    ;...
%             'MrT'   ,'07082013',  2,  [1, 2]    ...
% };

FileName = {'Mihili'   ,'10082015',  2, [1, 2];...           
            'Mihili'   ,'10272015',  2, [1, 2];...         
            'Mihili'   ,'11022015',  2, [1, 2] ...
};
G = 1:size(FileName,1);
%%
[T_CO, T_LOW, T_HIGH, UN_k, UN_ls, UN_mu, RB, VA, BPDS, N_PMD, N_M1] = ...
    deal(cell(size(FileName,1),1)); %initialize
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

    kept_llist = flipud(unique(alldays(2).tt(:,3)));
    LIlowhigh = zeros(1,2);
    for lks = 1:length(kept_llist)
        LIlowhigh(lks) = sum(alldays(2).tt(:,3)==kept_llist(lks));
    end
    
    load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_V-D-RT\\%s_PD90_%s_%s_V-D-RT.mat','PMd',FileName{daynum,1},FileName{daynum,2}));
    PDS_pmd = get_PDS_fulldat(fulldat,pi-(pi/2)*(strcmp(FileName{daynum,1},'Mihili')));
%     load(sprintf('C:\\Users\\limblab\\Desktop\\dropbox_overfill\\PD_V-D-RT\\%s_PD90_%s_%s_V-D-RT.mat','M1',FileName{daynum,1},FileName{daynum,2}));
%     PDS_m1 = get_PDS_fulldat(fulldat,pi-(pi/2)*(strcmp(FileName{daynum,1},'Mihili')));
    PDS_m1 = [nan nan nan];
    
    % Fill variables
    T_CO{daynum} = size(alldays(1).tt,1); % Number of CO trials
    T_LOW{daynum} = LIlowhigh(1); % Number of UN trials
    T_HIGH{daynum} = LIlowhigh(2); % Number of UN trials
    UN_ls{daynum} = kept_llist; % likelihoods used
    UN_mu{daynum} = round(circ_mean(alldays(2).tt(:,2))./pi*180); % Prior mean
    UN_k{daynum} = round(circ_kappa(alldays(2).tt(:,2))); % Prior kappa
    N_PMD{daynum} = [sum(~isnan(PDS_pmd(:,1))), sum(~isnan(PDS_pmd(:,2))), sum(~isnan(PDS_pmd(:,3)))];
    N_M1{daynum} = [sum(~isnan(PDS_m1(:,1))), sum(~isnan(PDS_m1(:,2))), sum(~isnan(PDS_m1(:,3)))];
    
    clearvars -except T_CO T_LOW T_HIGH UN_ls UN_mu UN_k N_PMD N_M1 FileName RB VA BPDS BRAIN_AREA G;
    
end
%%
[~,G2] = sortrows(cellfun(@(x) str2double(x),FileName(:,2)));

task_table = zeros(length(G2),14);
for i = 1:length(G2)
    daynum = G2(i);
    task_table(i,:) = round([double(strcmp(FileName{daynum,1},'Mihili')), T_CO{daynum}, T_LOW{daynum}, T_HIGH{daynum}...
                       UN_ls{daynum}(1), UN_ls{daynum}(2), UN_mu{daynum}, UN_k{daynum}...
                       N_PMD{daynum}, N_M1{daynum}]);
end

