% % Set up list of all days for inclusion
% FileName = {'Mihili','07112013',  2        ;...
%             'Mihili','07152013', [2,3]     ;...
%             'Mihili','07192013', [2,3,4]   ;...
%             'Mihili','08062013', [2,3]     ;...
%             'Mihili','08122013', [2,3]     ;...
%             'Mihili','08152013',  2        ;...
%             'MrT'   ,'05042013',  2        ;...
%             'MrT'   ,'05052013',  2        ;...
%             'MrT'   ,'05062013',  2        };
        
% FileName = {'Mihili','07112013',  2,  [1, 2]    ;...
%             'Mihili','07152013',  2,  [1, 2]    ;...
%             'Mihili','07192013',  2,  [1, 2]    ;...
%             'Mihili','08062013',  2,  [1, 2]    ;...
%             'Mihili','08122013',  2,  [1, 2]    ;...
%             'Mihili','08152013',  2,  [1, 2]    ;...
%             'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 3]    };
    
if ~exist('FileName','var')
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
end
%         
% FileName = {'Mihili','09282015', 2};%

% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10122015',  2        ;... 
%             'Mihili','10272015',  2        };

        
% FileName = {'Mihili','05062014', 2};
% FileName = {'Mihili','05062014',  2        ;...
%             'Mihili','11212013',  2        };

% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10272015',  2        ;...
%             'Mihili','11022015',  2         };
%%
% FileName = {'Mihili','08062013',  2 };
% FileName = {'Mihili','07192013', [2 3 4]};
    
[slopes,slope_bounds,LI,Tcount,resid_var,resid, slope_boot,centroids,reaches,targets] = ...
    deal(cell(size(FileName,1),1)); %initialize
figure; hold on; 
cols2p = {'b','r','g'};
for daynum = 1:size(FileName,1)
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2},0); 
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    alldays(2).tt(alldays(2).tt(:,3)>100,:) = [];
    alldays(2).tt(isnan(alldays(2).tt(:,9)),:) = [];
    alldays(2).tt(isnan(alldays(2).tt(:,10)),:) = [];
    
    llist = flipud(unique(alldays(2).tt(:,3)));
    alldays(2).tt(~ismember(alldays(2).tt(:,3),llist(FileName{daynum,4})),:) = [];
        
    like_conds = flipud(unique(alldays(2).tt(:,3)));
    Tcount{daynum}(1) = size(alldays(1).tt,1);
    for i = 1:length(like_conds)
        is = find(alldays(2).tt(:,3)==like_conds(i));
        LI{daynum}{i} = is;
%         [b,bint] = regress(alldays(2).tt(is,10),[ones(length(is),1) alldays(2).tt(is,9)]);
        
        Tcount{daynum}(i+1) = length(is);
        [lll, hhh ,bbb] = boot_bounds(1000,@(x) circ_polyfit(x(:,1),x(:,2)),alldays(2).tt(is,[9 10]),2.5,97.5);
        [b(1:2),resid{daynum}{i}] = circ_polyfit(alldays(2).tt(is,9),alldays(2).tt(is,10));
        bint(1,:) = [lll(1) hhh(1)];
        slope_boot{daynum}(:,i) = bbb(:,1);
        resid_var{daynum}(i) = circ_var(resid{daynum}{i});
        slope_bounds{daynum}(i,:) = bint(1,:); 
        slopes{daynum}(i) = b(1);
        centroids{daynum}{i} = alldays(2).tt(is,9);
        reaches{daynum}{i} = alldays(2).tt(is,10);
        targets{daynum}{i} = alldays(2).tt(is,2);
%         plot((daynum+0.05*(i-1)),slopes{daynum}(i),'.','Color',cols2p{i});
%         plot((daynum+0.05*(i-1))*[1 1],slope_bounds{daynum}(i,:),'Color',cols2p{i});
            
    end
    plot_cross(slopes{daynum}(1),slope_bounds{daynum}(1,:),slopes{daynum}(2),slope_bounds{daynum}(2,:),cols2p{(daynum>6)+1},1);
    
    clc; fprintf('%d/%d\n',daynum,size(FileName,1)); 
end
slope_std = [cellfun(@(x) std(x(:,1)),slope_boot), cellfun(@(x) std(x(:,2)),slope_boot)];

RES = [cellfun(@(x) circ_std(x{1}),resid), cellfun(@(x) circ_std(x{2}),resid)];
dRES = RES(:,2)-RES(:,1);
%%
dRES_bnd = zeros(size(RES,1),2);

for i = 1:size(RES,1)

res_label = [resid{i}{1} ones(size(resid{i}{1},1),1); resid{i}{2} 2*ones(size(resid{i}{2},1),1)];
difffunction = @(x) circ_std(x(x(:,2)==2))-circ_std(x(x(:,2)==1)); 

[l,h] = boot_bounds(1000,difffunction,res_label,2.5,97.5);

dRES_bnd(i,:) = [l,h];
end
%%
MT = 27:29;
figure; hold on; 
[res_1,res_2] = deal(zeros(size(RES,1),2));
for i = 1:size(RES,1) 
    [res_1(i,1),res_1(i,2)] = boot_bounds(1000,@circ_std,resid{i}{1},2.5,97.5);
    [res_2(i,1),res_2(i,2)] = boot_bounds(1000,@circ_std,resid{i}{2},2.5,97.5);
    if ismember(i,G)
        plot_cross(RES(i,1),res_1(i,:),RES(i,2),res_2(i,:),'b',1);
        plot(RES(i,1),RES(i,2),'b.'); 
    elseif ismember(i,MT)
        plot_cross(RES(i,1),res_1(i,:),RES(i,2),res_2(i,:),'b',1);
        plot(RES(i,1),RES(i,2),'bo'); 
    end
end
plot([0 .15],[0 0.15],'k--');
xlim([0.05 0.14]);
ylim([0.05 0.8]); 

%%
figure; hold on; 
for daynum = 1:size(FileName,1)
    
    if ismember(daynum,G)
        
        plot_cross(slopes{daynum}(1),slope_bounds{daynum}(1,:),...
                   slopes{daynum}(2),slope_bounds{daynum}(2,:),...
                   'b',1);
               
        plot(slopes{daynum}(1),slopes{daynum}(2),'b.');
    elseif ismember(daynum,MT)
        plot_cross(slopes{daynum}(1),slope_bounds{daynum}(1,:),...
                   slopes{daynum}(2),slope_bounds{daynum}(2,:),...
                   'r',1);
               
        plot(slopes{daynum}(1),slopes{daynum}(2),'ro');
    end
end
        

