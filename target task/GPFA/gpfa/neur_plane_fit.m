%% concatenated trajs for each direction
ds = cell(length(seqdirs),1);
for i = 1:length(seqdirs)  
    ds{i} = [];
    for j = 1:length(seqdirs{i})  
        ds{i} = [ds{i}; seqdirs{i}(j).xorth(1:3,:)'];
    end
end

%% averaged trajs (cut to shortest)
dscut = cell(length(seqdirs),1);

% find shortest trajectory
seqcell = struct2cell(seqTrain);
minL = min(vertcat(seqcell{2,:,:}));

for i = 1:length(seqdirs)
    
    dircell = struct2cell(seqdirs{i});
    xorthcell = vertcat(dircell(end,:,:));
    xorths = vertcat(xorthcell(:));
    
    for dim = 1:3
        
        cutdown = cellfun(@(x) x(dim,1:minL),xorths,'UniformOutput',0);
        cutarray = vertcat(cutdown{:});
        
        dscut{i}(dim,:) = mean(cutarray,1);
    end
end

% create trajectory array (direction x time x dimension)
dar = zeros(length(dscut),size(dscut{1},2),size(dscut{1},1));
for dim = 1:3
    dimavtrace = cellfun(@(x) x(dim,:),dscut,'UniformOutput',0);
    dar(:,:,dim) = vertcat(dimavtrace{:});
end
darcat = reshape(dar,size(dar,1)*size(dar,2),3);
p0 = mean(darcat,1);
%%

cost_helper = @(x) plane_fit_costfunc(dar,x,p0);

[sol,fval] = fminsearch(cost_helper,[1 1 1]);

proj_func = @(P,n,p0) [P(1) - n(1).*((sum(n.*P)-sum(n.*p0))./(sum(n.^2))),...
                          P(2) - n(2).*((sum(n.*P)-sum(n.*p0))./(sum(n.^2))),...
                          P(3) - n(3).*((sum(n.*P)-sum(n.*p0))./(sum(n.^2)))];
                     
solrep = repmat(sol,size(darcat,1),1);
p0rep = repmat(p0,size(darcat,1),1);
proj_func_all = @(P,nrep,p0rep) [P(:,1) - nrep(:,1).*((sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2))),...
                   P(:,2) - nrep(:,2).*((sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2))),...
                   P(:,3) - nrep(:,3).*((sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2)))];

darcat_proj = proj_func_all(darcat,solrep,p0rep);
dar_proj = reshape(darcat_proj,size(dar));

%%
N = size(darcat_proj,1);
origin = darcat_proj(1,:);
localz = cross(darcat_proj(2,:)-origin, darcat_proj(3,:)-origin);
unitz = localz/norm(localz,2);
localx = darcat_proj(2,:)-origin;
unitx = localx/norm(localx,2);
localy = cross(localz, localx);
unity = localy/norm(localy,2);
T = [localx(:), localy(:), localz(:), origin(:); 0 0 0 1];
C = [darcat_proj, ones(N,1)];
darcat_proj2 = T \ C';
darcat_proj2 = darcat_proj2(1:2,:)';

dar_proj2 = reshape(darcat_proj2,size(dar,1),size(dar,2),2);

%%
figure; hold on; 
subplot(1,2,1); hold on;
for i = 1:size(dar,1)
    plot3(dar(i,:,1),dar(i,:,2),dar(i,:,3),'.-','color',cols2plot(i,:));
end
subplot(1,2,2); hold on; 
for i = 1:size(darcat_proj,1)
    plot(dar_proj2(i,:,1),dar_proj2(i,:,2),'.-','color',cols2plot(i,:));
end

