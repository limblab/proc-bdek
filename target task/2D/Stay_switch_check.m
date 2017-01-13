dircol = 19;
align_col = 7;
window_bounds = [0 250];

[~,~,~,fulldatPDS] = tuning_types_PD_func(alldays,'PMd',1,dircol,{align_col},{[0 750]});
PDS = get_PDS_fulldat(fulldatPDS,pi/2);

PDnum = mod(round((PDS+2*pi)./(pi/4))-1,8)+1;
Tnum{1} = mod(round((alldays(1).tt(:,dircol)+2*pi)./(pi/4))-1,8)+1;
Tnum{2} = mod(round((alldays(2).tt(:,dircol)+2*pi)./(pi/4))-1,8)+1;

rast = cell(length(PDS),2);
for i = 1:length(PDS)
    clc; fprintf('getting raster %d/%d\n',i,length(PDS));
    if ~isnan(PDS(i))
        r1 = raster_get(alldays(1).PMd_units{i},alldays(1).tt,window_bounds/1000,align_col);
        rast{i,1} = bin_array(r1,size(r1,1),size(r1,2)./10,'sum');
        r2 = raster_get(alldays(1).PMd_units{i},alldays(2).tt,window_bounds/1000,align_col);
        rast{i,2} = bin_array(r2,size(r2,1),size(r2,2)./10,'sum');
        
        clear r1 r2
    end
end
%%
lenwin = 25;
avwind = ones(1,lenwin)./(lenwin/100);
sizewin = max(cellfun(@(x) size(x,2),rast(:,1)))-lenwin+1;
%%
aved = cell(size(rast,1),size(rast,2));
for i = 1:size(rast,1)
    clc; fprintf('%d/%d\n',i,size(rast,1));
    if ~isempty(rast{i,1})
        for j = 1:size(rast{i,1},1)
            aved{i,1}(j,:) = conv(rast{i,1}(j,:),avwind,'valid');
        end

        for j = 1:size(rast{i,2},1)
            aved{i,2}(j,:) = conv(rast{i,2}(j,:),avwind,'valid');
        end
    end
end

%%
getPD = @(pd,targs) pd == (mod(targs-1,8)+1);
getOD = @(pd,targs) pd == (mod(targs+4-1,8)+1);
getORTH = @(pd,targs) pd == (mod(targs+2-1,8)+1) | pd == (mod(targs+6-1,8)+1);

[T1,T2] = deal(cell(size(aved,1),3));
for i = 1:size(aved,1)
    
    if ~isnan(PDnum(i))
        T1{i,1} = aved{i,1}(getPD(PDnum(i),Tnum{1}),:);
        T1{i,2} = aved{i,1}(getOD(PDnum(i),Tnum{1}),:);
        T1{i,3} = aved{i,1}(getORTH(PDnum(i),Tnum{1}),:);
        
        T2{i,1} = aved{i,2}(getPD(PDnum(i),Tnum{2}),:);
        T2{i,2} = aved{i,2}(getOD(PDnum(i),Tnum{2}),:);
        T2{i,3} = aved{i,2}(getORTH(PDnum(i),Tnum{2}),:);
    end
end
T1full = T1(cellfun(@(x) ~isempty(x),T1(:,1)),:);
T2full = T2(cellfun(@(x) ~isempty(x),T2(:,1)),:);

% %%
medfrs1 = cellfun(@(x) nanmedian(x(:)),T1(:,3));
medfrs2 = cellfun(@(x) nanmedian(x(:)),T2(:,3)); 

percover = [];
[percover.pd{1},percover.od{1},percover.pd{2},percover.od{2}] = deal(nan(size(T1,1),2));
for i = 1:size(T1,1)

    if ~isempty(T1{i,1})
        
        if medfrs1(i) > 0
            percover.pd{1}(i,1) = nanmean(reshape(T1{i,1} > medfrs1(i),[],1))*100;
            percover.pd{1}(i,2) = nanmean(reshape(T1{i,1} == medfrs1(i),[],1))*100;

            percover.od{1}(i,1) = nanmean(reshape(T1{i,2} > medfrs1(i),[],1))*100;
            percover.od{1}(i,2) = nanmean(reshape(T1{i,2} == medfrs1(i),[],1))*100;       

            percover.pd{2}(i,1) = nanmean(reshape(T2{i,1} > medfrs2(i),[],1))*100;
            percover.pd{2}(i,2) = nanmean(reshape(T2{i,1} == medfrs2(i),[],1))*100;

            percover.od{2}(i,1) = nanmean(reshape(T2{i,2} > medfrs2(i),[],1))*100;
            percover.od{2}(i,2) = nanmean(reshape(T2{i,2} == medfrs2(i),[],1))*100;
        end
        
    end
        
end
%  
% %%
% medbytime{1} = cellfun(@(x) nanmedian(x,1),T1(:,3),'Uni',0);
% medbytime{2} = cellfun(@(x) nanmedian(x,1),T2(:,3),'Uni',0);
% 
% for i = 1:size(T1,1)
%     
%     if ~isempty(PDnum(i))
%         
%         nonzeroinds = find(medbytime{1}{i} > 0);
%         povermed.pd{1}(i,:) = 100*nanmean(reshape(T1{i,1}(:,nonzeroinds)>repmat(medbytime{1}{i}(:,nonzeroinds),size(T1{i,1},1),1),[],1));
%         povermed.od{1}(i,:) = 100*nanmean(reshape(T1{i,2}(:,nonzeroinds)>repmat(medbytime{1}{i}(:,nonzeroinds),size(T1{i,2},1),1),[],1));
%         
%         nonzeroinds = find(medbytime{2}{i} > 0);
%         povermed.pd{2}(i,:) = 100*nanmean(reshape(T2{i,1}(:,nonzeroinds)>repmat(medbytime{2}{i}(:,nonzeroinds),size(T2{i,1},1),1),[],1));
%         povermed.od{2}(i,:) = 100*nanmean(reshape(T2{i,2}(:,nonzeroinds)>repmat(medbytime{2}{i}(:,nonzeroinds),size(T2{i,2},1),1),[],1));
%     end
% end
% 
% %%
% numtrials = max(cellfun(@(x) size(x,1),aved(:,2)));
% medians = {cellfun(@(x) nanmedian(x(:)),T1(:,3)), cellfun(@(x) nanmedian(x(:)),T2(:,3))};
% medians{1}(medians{1}<1) = NaN; medians{2}(medians{2}<1) = NaN;
% pover_pdod = zeros(numtrials,3);
% for i = 1:numtrials % 2-T trials
%     clc; fprintf('%d/%d\n',i,numtrials);
%     pdinds   = getPD(Tnum{2}(i),PDnum);
%     odinds   = getOD(Tnum{2}(i),PDnum);
%     orthinds = getORTH(Tnum{2}(i),PDnum);
% %     
% %     pover_pdod(i,1) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,2),'Uni',0)) >= repmat(medians{2}(pdinds),1,sizewin),2)-...
% %                nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,2),'Uni',0)) <= repmat(medians{2}(pdinds),1,sizewin),2));
% %     
% %     pover_pdod(i,2) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,2),'Uni',0)) >= repmat(medians{2}(odinds),1,sizewin),2)-...
% %                nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,2),'Uni',0)) <= repmat(medians{2}(odinds),1,sizewin),2));
% %            
% %     pover_pdod(i,3) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) >= repmat(medians{2}(orthinds),1,sizewin),2)-...
% %                nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) <= repmat(medians{2}(orthinds),1,sizewin),2));
%            
%    pover_pdod(i,1) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,2),'Uni',0)) > repmat(medians{2}(pdinds),1,sizewin),2));
%     
%     pover_pdod(i,2) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,2),'Uni',0)) > repmat(medians{2}(odinds),1,sizewin),2));
%            
%     pover_pdod(i,3) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) > repmat(medians{2}(orthinds),1,sizewin),2));
% %            
% %     pover_pdod(i,1) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,2),'Uni',0)) > repmat(medians{2}(pdinds),1,sizewin),2)+...
% %                .5*nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,2),'Uni',0)) == repmat(medians{2}(pdinds),1,sizewin),2));
% %     
% %     pover_pdod(i,2) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,2),'Uni',0)) > repmat(medians{2}(odinds),1,sizewin),2)+...
% %                .5*nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,2),'Uni',0)) == repmat(medians{2}(odinds),1,sizewin),2));
% %            
% %     pover_pdod(i,3) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) > repmat(medians{2}(orthinds),1,sizewin),2)+...
% %                .5*nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) == repmat(medians{2}(orthinds),1,sizewin),2));
% end
%    
% %%
% numtrials1 = max(cellfun(@(x) size(x,1),aved(:,1)));
% pover_pdod1 = zeros(numtrials1,3);
% for i = 1:numtrials1 % 1-T trials
%     clc; fprintf('%d/%d\n',i,numtrials1);
%     pdinds = getPD(Tnum{1}(i),PDnum);
%     odinds = getOD(Tnum{1}(i),PDnum);
%     orthinds = getORTH(Tnum{1}(i),PDnum);
% 
%     pover_pdod1(i,1) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,1),'Uni',0)) >= repmat(medians{1}(pdinds),1,sizewin),2)-...
%                nanmean(cell2mat(cellfun(@(x) x(i,:), aved(pdinds,1),'Uni',0)) <= repmat(medians{1}(pdinds),1,sizewin),2));
%     
%     pover_pdod1(i,2) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,1),'Uni',0)) >= repmat(medians{1}(odinds),1,sizewin),2)-...
%                nanmean(cell2mat(cellfun(@(x) x(i,:), aved(odinds,1),'Uni',0)) <= repmat(medians{1}(odinds),1,sizewin),2));
%            
%     pover_pdod1(i,3) = nanmean(nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) >= repmat(medians{1}(orthinds),1,sizewin),2)-...
%                nanmean(cell2mat(cellfun(@(x) x(i,:), aved(orthinds,2),'Uni',0)) <= repmat(medians{1}(orthinds),1,sizewin),2));
% end
% 
% %%
% if 0
% numtrials = max(cellfun(@(x) size(x,1),aved(:,2)));
% pover_pdod = zeros(numtrials,3);
% 
% totav = cellfun(@(x) nanmean(x,2),aved,'Uni',0);
% 
% [t1_pdod,t2_pdod] = deal(zeros(size(totav,1),1));
% for i = 1:size(totav,1)
%     
%     orthinds1 = getORTH(PDnum(i),Tnum{1});
%     pdinds1   = getPD(PDnum(i),Tnum{1});
%     odinds1   = getOD(PDnum(i),Tnum{1});
%     
%     orthinds2 = getORTH(PDnum(i),Tnum{2});
%     pdinds2   = getPD(PDnum(i),Tnum{2});
%     odinds2   = getOD(PDnum(i),Tnum{2});
%     
% %     t1_pdod(i,1) = nanmean(totav{i,1}(pdinds1)>nanmedian(totav{i,1}(orthinds1)));
% %     t1_pdod(i,2) = nanmean(totav{i,1}(odinds1)>nanmedian(totav{i,1}(orthinds1)));
% %     t1_pdod(i,3) = nanmean(totav{i,1}(odinds1|pdinds1)>nanmedian(totav{i,1}(orthinds1)));
%     
%     t1_pdod(i,1) = nanmean(totav{i,1}(pdinds1)>nanmedian(totav{i,1}(orthinds1)))+.5*nanmean(totav{i,1}(pdinds1)==nanmedian(totav{i,1}(orthinds1)));
%     t1_pdod(i,2) = nanmean(totav{i,1}(odinds1)>nanmedian(totav{i,1}(orthinds1)))+.5*nanmean(totav{i,1}(odinds1)==nanmedian(totav{i,1}(orthinds1)));
%     t1_pdod(i,3) = nanmean(totav{i,1}(odinds1|pdinds1)>nanmedian(totav{i,1}(orthinds1)))+.5*nanmean(totav{i,1}(odinds1|pdinds1)==nanmedian(totav{i,1}(orthinds1)));
%     
% %     t2_pdod(i,1) = nanmean(totav{i,2}(pdinds2)>nanmedian(totav{i,2}(orthinds2)));
% %     t2_pdod(i,2) = nanmean(totav{i,2}(odinds2)>nanmedian(totav{i,2}(orthinds2)));
% %     t2_pdod(i,3) = nanmean(totav{i,2}(odinds2|pdinds2)>nanmedian(totav{i,2}(orthinds2)));
%     
%     t2_pdod(i,1) = nanmean(totav{i,2}(pdinds2)>nanmedian(totav{i,2}(orthinds2)))+.5*nanmean(totav{i,2}(pdinds2)==nanmedian(totav{i,2}(orthinds2)));
%     t2_pdod(i,2) = nanmean(totav{i,2}(odinds2)>nanmedian(totav{i,2}(orthinds2)))+.5*nanmean(totav{i,2}(odinds2)==nanmedian(totav{i,2}(orthinds2)));
%     t2_pdod(i,3) = nanmean(totav{i,2}(odinds2|pdinds2)>nanmedian(totav{i,2}(orthinds2)))+.5*nanmean(totav{i,2}(odinds2|pdinds2)==nanmedian(totav{i,2}(orthinds2)));
%     
%     
% end
% %%
% mediansbytime = {cellfun(@(x) nanmedian(x,1),T1(:,3),'Uni',0), cellfun(@(x) nanmedian(x,1),T2(:,3),'Uni',0)};
% [PD_over,OD_over,ORTH_over] = deal(cell(2,1));
% 
% oneTavs = cell(8,3);
% for i = 1:8
%     targtrials = find(Tnum{1}==i);
%     pdinds = getPD(i,PDnum);
%     odinds = getOD(i,PDnum);
%     orthinds = getORTH(i,PDnum);
%     oneTavs{i,1} = cell2mat(cellfun(@(x) nanmean(x(targtrials,:),1),aved(pdinds,1),'Uni',0));
%     oneTavs{i,2} = cell2mat(cellfun(@(x) nanmean(x(targtrials,:),1),aved(odinds,1),'Uni',0));
%     oneTavs{i,3} = cell2mat(cellfun(@(x) nanmean(x(targtrials,:),1),aved(orthinds,1),'Uni',0));
%     
% end
% 
% 
% for numtrg = 1:2
%     numtrials = max(cellfun(@(x) size(x,1),aved(:,numtrg))); 
%     for i = 1:numtrials % 2-T trials
%         clc; fprintf('%d/%d\n',i,numtrials);
%         
%         t1num = Tnum{numtrg}(i);
%         t2num = mod(t1num-1+4,8)+1;
%         
%         pdinds   = getPD(t1num,PDnum);
%         odinds   = getOD(t1num,PDnum);
%         orthinds = getORTH(t1num,PDnum);
% 
%         pdact = cell2mat(cellfun(@(x) x(i,:), aved(pdinds,numtrg),'Uni',0));
%         odact = cell2mat(cellfun(@(x) x(i,:), aved(odinds,numtrg),'Uni',0));
%         orthact = cell2mat(cellfun(@(x) x(i,:),aved(orthinds,numtrg),'Uni',0));
%         
%         pdmeds = cell2mat(mediansbytime{numtrg}(pdinds));
%         odmeds = cell2mat(mediansbytime{numtrg}(odinds));
%         orthmeds=cell2mat(mediansbytime{numtrg}(orthinds));
%         
% %         pdmeds(pdmeds <1) = nan;
% %         odmeds(odmeds <1) = nan;
% %         orthmeds(orthmeds<1) = nan;
% % 
% %         PD_over{numtrg}(i,:) = nanmean(reshape(pdact>=pdmeds,[],1));
% %         OD_over{numtrg}(i,:) = nanmean(reshape(odact>=odmeds,[],1));
% %         ORTH_over{numtrg}(i,:) = nanmean(reshape(orthact>=orthmeds,[],1));
% 
% 
%         PD_over{numtrg}(i,:) = nanmean(reshape(pdact>pdmeds,[],1)) + ...
%                                .5*nanmean(reshape(pdact==pdmeds,[],1));
%         OD_over{numtrg}(i,:) = nanmean(reshape(odact>odmeds,[],1)) + ...
%                                .5*nanmean(reshape(odact==odmeds,[],1));
%         ORTH_over{numtrg}(i,:) = nanmean(reshape(orthact>orthmeds,[],1)) + ...
%                                .5*nanmean(reshape(orthact==orthmeds,[],1));
% %               
% %         PD_over{numtrg}(i,:) = nanmean(pdact(:)-oneTavs{t1num,1}(:));
% %         OD_over{numtrg}(i,:) = nanmean(odact(:)-oneTavs{t2num,1}(:));
% % %         ORTH_over{numtrg}(i,:) = nanmean(orthact(:)-oneTavs{Tnum{numtrg}(i),3}(:));
% 
% %         PD_over{numtrg}{i,:} = nanmean(pdact-pdmeds,1);
% %         OD_over{numtrg}{i,:} = nanmean(odact-odmeds,1);
% %         ORTH_over{numtrg}{i,:} = nanmean(orthact-orthmeds,1);
% %         
%         
% %         PD_over{numtrg}{i,:} = nanmean(sqrt(pdact)-sqrt(pdmeds),1);
% %         OD_over{numtrg}{i,:} = nanmean(sqrt(odact)-sqrt(odmeds),1);
% %         ORTH_over{numtrg}{i,:} = nanmean(sqrt(orthact)-sqrt(orthmeds),1);
%         
%         
%         
%     end
% end
% end
