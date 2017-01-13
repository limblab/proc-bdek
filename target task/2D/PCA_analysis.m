FC = trial_counts{1};
clear FCC co_inds PC;
for i = 1:length(FC); FCC(:,:,i) = FC{i}'; end

for i = 1:8; co_inds{i} = find(circ_dist(alldays(1).tt(:,10),pi/4*(i-1)) < pi/8); end

Ftot = reshape(FCC,size(FCC,1),numel(FCC(1,:,:)));

colors = {'r','g','b','m','k','c','m','b'};
figure; hold on; 
for targ = 1:length(co_inds)

    FCC_cur = FCC(:,co_inds{targ},:);
    F_all = reshape(FCC_cur,size(FCC_cur,1),numel(FCC_cur(1,:)));
    F_norm = (F_all - repmat(min(Ftot,[],2),1,size(F_all,2)))./repmat((max(Ftot,[],2) - min(Ftot,[],2)),1,size(F_all,2));

    coeff = princomp(F_norm);

    pc1_m = mean(reshape(coeff(:,1),size(FCC_cur,2),size(FCC_cur,3)));
    pc2_m = mean(reshape(coeff(:,2),size(FCC_cur,2),size(FCC_cur,3)));
    pc3_m = mean(reshape(coeff(:,3),size(FCC_cur,2),size(FCC_cur,3)));
    pc4_m = mean(reshape(coeff(:,4),size(FCC_cur,2),size(FCC_cur,3)));
    
    PC{1}(targ,:) = pc1_m; 
    PC{2}(targ,:) = pc2_m;
    PC{3}(targ,:) = pc3_m;
    PC{4}(targ,:) = pc4_m;
    
    %plot3(smooth(pc1_m),smooth(pc2_m),smooth(pc3_m),'Color',colors{targ});
    plot(smooth(pc1_m),smooth(pc2_m),colors{targ});
end

%%
figure; hold on; 
[coeffs,scores] = deal(cell(length(co_inds),1));
%for targ = 1:length(co_inds)

    for timebin = 1:size(FCC,3)
        FCC_cur = FCC(:,co_inds{targ},timebin);
        F_norm = (FCC_cur - repmat(min(Ftot,[],2),1,size(FCC_cur,2)))./repmat((max(Ftot,[],2) - min(Ftot,[],2)),1,size(FCC_cur,2));
        F_norm(isnan(F_norm))= 0;

        [coeffs{targ}{timebin},scores{targ}{timebin}] = princomp(F_norm');

        plot(scores{targ}{timebin}(:,1),scores{targ}{timebin}(:,2),'.','Color',colors{targ});

%         PC{1}(targ,:) = pc1_m; 
%         PC{2}(targ,:) = pc2_m;
%         PC{3}(targ,:) = pc3_m;
%         PC{4}(targ,:) = pc4_m;
% 
%         plot3(smooth(pc1_m),smooth(pc2_m),smooth(pc3_m),'Color',colors{targ});
    end
%end
%%
[coeff,score] = princomp(zscore(Ftot'));
sc = score';
scr = reshape(sc,size(FCC));

[scrdir,avscrdir] = deal(cell(length(co_inds),1));
for i = 1:length(co_inds)
    scrdir{i} = scr(:,co_inds{i},:);
    avscrdir{i} = reshape(mean(scrdir{i},2),size(FCC,1),size(scrdir{i},3));
end

figure; hold on; 
for i = 2:24 
    for j = 1:8 
        plot3(avscrdir{j}(1,1:i),avscrdir{j}(2,1:i),avscrdir{j}(3,1:i),'-','Color',colors{j}); 
        xlim([-3 7]); ylim([-3 7]);  
    end 
    pause(0.25); 
end

%%
% FT_m = Ftot - repmat(mean(Ftot,2),1,size(Ftot,2));
% FT_s = FT_m./repmat(max(abs(FT_m),[],2),1,size(Ftot,2));
% FT_s(isnan(FT_s)) = 0;

% FT_trial = reshape(FT_s,size(FCC));

% FC_fact = FCC_factor';
% FCF_m = FC_fact - repmat(mean(FC_fact,2),1,size(FC_fact,2));
% FCF = FCF_m./repmat(max(abs(FCF_m),[],2),1,size(FCF_m,2));
% FCF(isnan(FCF)) = 0;
% 
% nonzerneurs = find(max(FCF,[],2)~=0);
% lambdanz = factoran(FCF(:,nonzerneurs),3);
[FCC,FCU] = deal(zeros(size(trial_counts{1}{1},1),size(trial_counts{1}{1},2),length(trial_counts{1})));
for i = 1:length(trial_counts{1})
    FCC(:,:,i) = trial_counts{1}{i};
end

for i = 1:length(trial_counts{1})
    FCU(:,:,i) = trial_counts{1}{i};
end
%%

Factor_base = FCC_factor;
nonzerneurs = find(max(Factor_base,[],1)~=0);
lambdanz = factoran(Factor_base(:,nonzerneurs),3);

lambda = zeros(size(Factor_base,2),3);
lambda(nonzerneurs,:) = lambdanz;
%% Center Out
[f1,f2,f3] = deal(zeros(size(Factor_base,1),3));
for timeb = 1:size(FCC,3)
    proj = squeeze(FCC(:,:,timeb)) * lambda;
    f1(:,timeb) = proj(:,1);
    f2(:,timeb) = proj(:,2); 
    f3(:,timeb) = proj(:,3);
end

[cof1,cof2,cof3,avf1,avf2,avf3] = deal(cell(length(co_inds),1));
for i = 1:length(co_inds)

    cof1{i} = f1(co_inds{i},:);
    cof2{i} = f2(co_inds{i},:);
    cof3{i} = f3(co_inds{i},:);
    
    avf1{i} = mean(cof1{i});
    avf2{i} = mean(cof2{i});
    avf3{i} = mean(cof3{i});
end
%% Uncertainty
[fu1,fu2,fu3] = deal(zeros(size(FCU,1),3));
for timeb = 1:size(FCU,3)
    proj = squeeze(FCU(:,:,timeb)) * lambda;
    fu1(:,timeb) = proj(:,1);
    fu2(:,timeb) = proj(:,2); 
    fu3(:,timeb) = proj(:,3);
end

[unf1,unf2,unf3,avuf1,avuf2,avuf3] = deal(cell(length(like_ind{1}),1));
for i = 1:length(like_ind{1})

    unf1{i} = fu1(like_ind{1}{i},:);
    unf2{i} = fu2(like_ind{1}{i},:);
    unf3{i} = fu3(like_ind{1}{i},:);
    
    avuf1{i} = mean(unf1{i});
    avuf2{i} = mean(unf2{i});
    avuf3{i} = mean(unf3{i});
end

figure; hold on; 
for i = 2:size(FCC,3) 
    for j = 1:length(co_inds) 
        plot3(avf1{j}(1:i),avf2{j}(1:i),avf3{j}(1:i),'-','Color',colors{j}); 
        %xlim([-3 7]); ylim([-3 7]);  
    end 
    pause(0.25); 
end

figure; hold on;
uncols = {'b','g','r'};
for i = 2:size(FCC,3)
    for j = 1:length(like_ind{1})
        plot3(avuf1{j}(1:i),avuf2{j}(1:i),avuf3{j}(1:i),'-','Color',colors{j},'LineWidth',3); 
    end
    pause(0.25);
end
        
%figure; hold on;
%uncols = {'b','g','r'};
for i = 2:size(FCC,3)
    for j = 1:length(like_ind{1})
        for k = 1:size(unf1{j},1)
            plot3(unf1{j}(k,1:i),unf2{j}(k,1:i),unf3{j}(k,1:i),'-','Color',colors{j},'LineWidth',0.5); 
        end
    end
    pause(0.25);
end
        
        
        