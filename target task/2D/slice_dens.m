

around_cent = slices - repmat(comp_tt(:, 9),1,numslices);
around_est  = slices - repmat(comp_tt(:,10),1,numslices);

cent_L = reshape(around_cent(Linds,:),1,numslices*length(Linds));
cent_H = reshape(around_cent(Hinds,:),1,numslices*length(Hinds));

est_L = reshape(around_est(Linds,:),1,numslices*length(Linds));
est_H = reshape(around_est(Hinds,:),1,numslices*length(Hinds));

L_inds = (1:length(Linds))';
H_inds = (1:length(Hinds))';

[junk, lorder] = sortrows(comp_tt(Linds,:));
[junk, horder] = sortrows(comp_tt(Hinds,:));

sort_centroids = 1;

if sort_centroids == 1
    LORDER = L_inds(lorder);
    HORDER = H_inds(horder);
else 
    LORDER = L_inds;
    HORDER = H_inds;
end

mup = pi/2;
kappa = 10;
mu_hats = comp_tt(:,10);
sig_hats = std(slices,0,2);
sigp = 1 - besseli(1,kappa)./besseli(0,kappa);

posBayes = (mu_hats.*sigp.^2 + mup.*sig_hats.^2)./(sigp.^2+sig_hats.^2);

aroundBayes = slices - repmat(posBayes,1,numslices);

Bayes_L = reshape(aroundBayes(Linds,:),1,numslices*length(Linds));
Bayes_H = reshape(aroundBayes(Hinds,:),1,numslices*length(Hinds));
    
figure; subplot(1,2,1);
plot(around_est(Linds,:),repmat(LORDER,1,numslices),'b.');
title('Centered on estimate');
subplot(1,2,2);
plot(around_est(Hinds,:),repmat(HORDER,1,numslices),'r.');


figure; subplot(1,2,1); 
plot(around_cent(Linds,:),repmat(LORDER,1,numslices),'b.');
title('Centered on centroid');
subplot(1,2,2); 
plot(around_cent(Hinds,:),repmat(HORDER,1,numslices),'r.');


figure; subplot(1,2,1); 
plot(aroundBayes(Linds,:),repmat(LORDER,1,numslices),'b.');
title('Centered on Bayes');
subplot(1,2,2);
plot(aroundBayes(Hinds,:),repmat(HORDER,1,numslices),'r.');


[nL_est, xout1] = hist(est_L,50);
[nH_est, xout2] = hist(est_H,50);

[nL_cent,xout3] = hist(cent_L,50);
[nH_cent,xout4] = hist(cent_H,50);

[nL_Bayes,xout5]= hist(Bayes_L,50);
[nH_Bayes,xout6]= hist(Bayes_H,50);

figure; hold on;
plot(xout5,nL_Bayes,'g');
plot(xout1,nL_est,'b.-');

figure; hold on;
plot(xout6,nH_Bayes,'g');
plot(xout2,nH_est,'r.-');


%%
stepplot = 0;
if stepplot == 1
    figure; hold on;
    inout = [7.5 8.5];

    for i = 1:length(comp_tt)
        plot((inout(1)-0.01)*cos(0:0.01:2*pi),(inout(1)-0.01)*sin(0:0.01:2*pi),'k');
        plot((inout(2)+0.01)*cos(0:0.01:2*pi),(inout(2)+0.01)*sin(0:0.01:2*pi),'k');
        for j = 1:numslices
            plot(inout.*cos(slices(i,j)),inout.*sin(slices(i,j)),'Linewidth',2,'Color','b');
        end

        plot(inout.*cos(comp_tt(i,10)),inout.*sin(comp_tt(i,10)),'Linewidth',2,'Color','r');

        pause; cla;
    end
end     
        
%%
for kp = 1:1:100
    vp = 1 - besseli(1,kp)./besseli(0,kp); 
    for kf = 1:1:100
        vf = 1 - besseli(1,kf)./besseli(0,kf); 
        
        bayesest(kp,kf) = vp/(vp+vf/10);
        
    end
end
figure; surf(1:100,1:100,bayesest);
xlabel('Prior'); ylabel('Feedback');
           
           
           
           
           
           
           
           
           
        