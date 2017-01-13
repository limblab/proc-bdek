%D1D2{TIME}{PD_OD_ORTH}{DAYS} 
dfunc = @(x) nanmean(x(x(:,end)==2,1))-nanmean(x(x(:,end)==1,1));
difffunc = @(x) nanmean(nanmean(x(x(:,end)==2,1:end-1))-nanmean(x(x(:,end)==1,1:end-1)));
splitfunc = @(x) {find(x(:,end)==1),find(x(:,end)==2)};
[Ld,Md,Hd] = deal(zeros(3,length(D1D2)));
[adiffs,neuravs,navs] = deal(cell(3,1));
[DDmeanall,DDSE] = deal(zeros(3,length(D1D2)));
for tim = 1:length(D1D2)
    fprintf('%d/%d\n',tim,length(D1D2));
    for pdodorth = 1:3
        split_ses = D1D2{tim}{pdodorth}';

        S = cell(length(split_ses),1);
        ML = zeros(length(split_ses),2);
        for i = 1:length(split_ses)
            inds = splitfunc(D1D2{tim}{pdodorth}{i});
            
            [A,B] = meshgrid(inds{1},inds{2});
            c = cat(2,A',B');
            d = reshape(c,[],2);

            DD = D1D2{tim}{pdodorth}{i}(d(:,2),:) - D1D2{tim}{pdodorth}{i}(d(:,1),:);
%             S{i} = typecast(reshape(DD(~isnan(DD)),[],1),'int8');
            ML(i,:) = [nansum(DD(:)) sum(~isnan(DD(:)))];

            neuravs{pdodorth}{i}(tim,:) = nanmean(DD); 
%             S{i} = sparse(reshape(DD(~isnan(DD)),[],1));
            adiffs{pdodorth}(i,tim) = nanmean(DD(:));
%             DDcat{i}{pdodorth,tim} = reshape(DD(~isnan(DD)),[],1);
%             dcheck{i}(pdodorth,tim) = difffunc2(split_ses{i});
        end
%         DDmeanall(pdodorth,tim) = mean(vertcat(S{:}));
        DDmeanall(pdodorth,tim) = sum(ML(:,1))./sum(ML(:,2));
%         DDSE(pdodorth,tim) = std(vertcat(S{:}))./sqrt(length(vertcat(S{:})));
%         
    end
end
for pdodorth = 1:3
    navs{pdodorth} = horzcat(neuravs{pdodorth}{:});
end
% DDmeanall([2 3],:) = DDmeanall([3 2],:);

% %% Plot adiffs
% ord = [1 3 2];
% figure; hold on; 
% for ss = 1:length(ord)
%     i = ord(ss);
%     subplot(1,3,ss); hold on; 
%     plot(adiffs{i}','Color',[.75 .75 .75]); 
%     plot(DDmeanall(i,:),'k','LineWidth',5); 
%     ylim([-6 6]); 
%     plot([1 size(adiffs{1},2)],[0 0],'k--');
% end

%% Plot avdiffs
ord = [1 3 2];
figure; hold on; 

for ss = 1:length(ord)
    i = ord(ss);
    subplot(1,3,ss); hold on; 

    if size(avdiffs{1},2)>1 
        xp = 1:size(avdiffs{1},2); 
        yp1 = avdiffs{ss}';
        yp2 = DDmeanall(i,:);
    else
        xp = [0 1];
        yp1 = [avdiffs{ss} avdiffs{ss}];
        yp2 = [DDmeanall(i,:) DDmeanall(i,:)];
    end
    plot(xp,yp1,'Color',[.75 .75 .75]); 
    plot(xp,yp2,'Color','k','LineWidth',5); 
    ylim([-3 3]); 
    plot([1 size(avdiffs{ss},2)],[0 0],'k--');
end