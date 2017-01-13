%% concatenated trajs for each direction
us = cell(length(seqlikes),1);
for i = 1:length(seqlikes)  
    us{i} = [];
    for j = 1:length(seqlikes{i})  
        us{i} = [us{i}; seqlikes{i}(j).xorth(1:end,:)'];
    end
end
uscat = vertcat(us{:});
p0 = mean(uscat,1);
%%

cost_helper = @(x) plane_fit_costfunc_kmeans(us,x,p0);

[sol, fval] = fminsearch(cost_helper,[0.5774 0.5774 0.5774]);
% [sol, fval] = fminsearch(cost_helper,[1 0 0]);
 
proj_func = @(P,n,p0) [P(1) - n(1).*((sum(n.*P)-sum(n.*p0))./(sum(n.^2))),...
                          P(2) - n(2).*((sum(n.*P)-sum(n.*p0))./(sum(n.^2))),...
                          P(3) - n(3).*((sum(n.*P)-sum(n.*p0))./(sum(n.^2)))];
                     
solrep = repmat(sol,size(uscat,1),1);
p0rep = repmat(p0,size(uscat,1),1);
proj_func_all = @(P,nrep,p0rep) [P(:,1) - nrep(:,1).*((sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2))),...
                   P(:,2) - nrep(:,2).*((sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2))),...
                   P(:,3) - nrep(:,3).*((sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2)))];

uscat_proj = proj_func_all(uscat,solrep,p0rep);

us_proj = cell(length(us),1);
us_projXX = uscat_proj;
for i = 1:length(us)
    us_proj{i} = uscat_proj(1:length(us{i}),:);
    us_projXX(1:length(us{i}),:) = [];
end

%%
N = size(uscat_proj,1);
origin = uscat_proj(1,:);

localz = cross(uscat_proj(2,:)-origin, uscat_proj(3,:)-origin);
unitz = localz/norm(localz,2);

localx = uscat_proj(2,:)-origin;
unitx = localx/norm(localx,2);

localy = cross(localz, localx);
unity = localy/norm(localy,2);

T = [unitx(:), unity(:), unitz(:), origin(:); 0 0 0 1];
C = [uscat_proj, ones(N,1)];
uscat_proj2 = T' \ C';
uscat_proj2 = uscat_proj2(1:2,:)';

us_cell = cell(length(us),1);
uscat_projXX = uscat_proj2;
for i = 1:length(us)
    us_cell{i} = uscat_projXX(1:length(us{i}),:);
    uscat_projXX(1:length(us{i}),:) = [];
end
%us_proj2 = reshape(uscat_proj2,size(us,1),size(us,2),2);

%%
usproj_sep = cell(length(us_cell),1);
uscell_XX = us_cell;
for i = 1:length(uscell_XX)
    numTrials = length(seqlikes{i});
    for j= 1:numTrials
        lengthTrial = seqlikes{i}(j).T;
        
        usproj_sep{i}{j} = uscell_XX{i}(1:lengthTrial,:);
        
        uscell_XX{i}(1:lengthTrial,:) = [];
    end
end

%%
figure; hold on; 
% subplot(1,2,1); hold on;
% for i = 1:size(us,1)
%     plot3(us(i,:,1),us(i,:,2),us(i,:,3),'.-','color',cols2plot(i,:));
% end
% subplot(1,2,2); hold on; 
for i = 1:size(usproj_sep,1)
    for j = 1:length(usproj_sep{i})
        %plot(usproj_sep{i}{j}(:,1),usproj_sep{i}{j}(:,2),'.-','color',cols2plot{i});
        plot(usproj_sep{i}{j}((end-20):end,1),usproj_sep{i}{j}((end-20):end,2),'.','color',cols2plot{i});

    end
end

