Mdl = cell(4,1);
for i = 1:4
    
    G1 = cell2mat(DIRREG(i,:)');
    G2 = cell2mat(DIRREG(i+4,:)');
    
    G = [G1;G2];
    
    Lab1 = ones(size(G1,1),1);
    Lab2 = 2*ones(size(G2,1),1);
    Lab = [Lab1; Lab2];
    
    T = [array2table(G), table(Lab)];
    
    Mdl{i} = fitcsvm(T,'Lab');
    clc; fprintf('%d/%d\n',i,4);
end
% 
cls = cell(4,1);
for bl = 1:4
    stateG = permute(cell2mat(cellfun(@(x) permute(x',[1 3 2]),EPOCH{bl}','Uni',0)),[3 2 1]);
    for i = 1:size(alldays(bl).tt,1)
        tgt = targids{bl}(i);
        if tgt > 4
            tgt = tgt-4;
            signv = -1;
        else
            signv = 1;
        end

        [~,score] = predict(Mdl{tgt},squeeze(stateG(i,:,:))');

        cls{bl}(i,:) = score(:,1)'.*signv;

    end
end

%% Built classifier for EVERY target
Mdl_targ = cell(8,1);
targnums = 1:8;
orths = @(x) mod(x+[2 6]-1,8)+1;
for i = 1:8
    clc; fprintf('%d/%d\n',i,4);
    
%     targother = targnums(targnums~=i); 
    targother = orths(i);
    
    G1 = cell2mat(DIRREG(i,:)');
    G2 = cell2mat(reshape(DIRREG(targother,:),[],1));
    G = [G1;G2];
    
    Lab1 = ones(size(G1,1),1);
    Lab2 = 2*ones(size(G2,1),1);
    Lab = [Lab1; Lab2];
    
    T = [array2table(G), table(Lab)];
    
    Mdl_targ{i} = fitcsvm(T,'Lab');
    
end
% 
cls_targ = cell(4,1);
for bl = 1:4
    stateG = permute(cell2mat(cellfun(@(x) permute(x',[1 3 2]),EPOCH{bl}','Uni',0)),[3 2 1]);
    for i = 1:size(alldays(bl).tt,1)
        tgt = targids{bl}(i);

        allscores = zeros(8,size(stateG,3));
        for k = 1:8
            [~,score] = predict(Mdl_targ{k},squeeze(stateG(i,:,:))');
            allscores(k,:) = score(:,1)';
        end
        cls_targ{bl}(:,:,i) = allscores;
    end
end

%%
ABrep = M.pmd{2}.d([1 5],:,:)./repmat(nanmean(M.pmd{2}.d([3 7],:,:),1),2,1,1);
AB_mah{1} = reshape(ABrep(:,delay1,:),2,[])';
AB_mah{2} = reshape(ABrep(:,delay2,:),2,[])';
AB_mah{3} = reshape(ABrep(:,mem,:),2,[])';
figure; hold on; 
for i = 1:3
    subplot(1,3,i); hold on; 
    plot(AB_mah{i}(:,1),AB_mah{i}(:,2),'.'); 
end
