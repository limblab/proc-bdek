function[valatend,percBayes,distfromopt,spread_pos] = departure_from_Bayes(tt)
% ii = i;

% From "Cue combination on the circle and the sphere" (Richard F. Murray
% and Yaniv Morgenstern)
prodvm = @(k1,k2,k3,mu3,theta) besseli(0,k3)/(2*pi*besseli(0,k1)*besseli(0,k2))*circ_vmpdf(theta,mu3,k3);
getmu3 = @(mu1,mu2,k1,k2) mu1 + atan2(sin(mu2-mu1),(k1/k2)+cos(mu2-mu1));
getk3 = @(k1,k2,mu1,mu2) sqrt(k1^2 + k2^2 + 2*k1*k2*cos(mu2-mu1));

%sample kappa
%khat = @(angs,ii) 2*abs(1/length(angs)*sum(exp(ii*angs)))^2/(1 - abs(1/length(angs)*sum((exp(ii*angs)).^2)));

% % variance from std
% cvar = @(sd) 1 - exp(-sd^2/2);

% numerically compute k3
kappas = 0:0.001:700;
model_vars = 1 - besseli(1,kappas)./besseli(0,kappas);

mup_true = circ_mean(tt(:,2));
mup_cent = circ_mean(tt(:,9));
kappap_true = kappa_calc(tt(:,2),0.05);
kappap_cent = kappa_calc(tt(:,9),0.05);

valatend = zeros(size(tt,1),1);
percBayes = zeros(size(tt,1),1);
distfromopt = zeros(size(tt,1),1);
spread_pos = zeros(size(tt,1),1);
for i = 1:size(tt,1)
    
    mu_l = tt(i,9);
    var_l = tt(i,11);
    kappa_l = interp1(model_vars,kappas,var_l/5,'nearest','extrap');
    
    mu_p = mup_cent;
    kappa_p = kappap_true;
    
    k3 = getk3(kappa_l,kappa_p,mu_l,mu_p);
    mu3 = getmu3(mu_l,mu_p,kappa_l,kappa_p);
    
    endpoint = tt(i,10);
    
    xs = 0:0.01:2*pi;
    totdist = prodvm(kappa_l,kappa_p,k3,mu3,xs);
    totdist = totdist/(sum(totdist)*0.01);
    
    if sum(isnan(totdist)) > 0 
        totdist = circ_vmpdf(xs,mu3,700); 
        valatend(i) = circ_vmpdf(endpoint,mu3,700);
        spread_pos(i) = 700;
    else
        valatend(i) = prodvm(kappa_l,kappa_p,k3,mu3,endpoint);     
        spread_pos(i) = k3;
    end
    
    percBayes(i) = valatend(i)/max(totdist);
    distfromopt(i) = circ_dist(endpoint,mu3);
    
end


 %%
% figure; hold on; plot(cos(0:0.1:2*pi),sin(0:0.1:2*pi),'k');
% plot(cos(alldays(2).slices(1,:)),sin(alldays(2).slices(1,:)),'ro')
% 
% plot(cos(ll),sin(ll),'bo')
% plot(cos(ul),sin(ul),'bo')
% 
% %%
% mu = mean(alldays(2).slices(1,:));
% k = circ_kappa(alldays(2).slices(1,:));
% 
% figure; hold on;
% plot(alldays(2).slices(1,:),ones(1,5),'ro')
% plot(0:0.01:2*pi,circ_vmpdf(0:0.01:2*pi,mu,k),'b');

    
% figure; hold on;
% plot(bin_array(spinv(linds),5,1),'b'); 
% plot(bin_array(spinv(minds),5,1),'g'); 
% plot(bin_array(spinv(hinds),5,1),'r'); 
%     
    
% figure; hold on;
% plot(bin_array(abs(lfromp(linds)),5,1),'b'); 
% plot(bin_array(abs(lfromp(minds)),5,1),'g'); 
% plot(bin_array(abs(lfromp(hinds)),5,1),'r'); 
% 
% figure; hold on;
% plot(bin_array(abs(efromp(linds)),5,1),'b'); 
% plot(bin_array(abs(efromp(minds)),5,1),'g'); 
% plot(bin_array(abs(efromp(hinds)),5,1),'r'); 
% 
% figure; hold on;
% plot(bin_array(alldays(2).tt(linds,11),5,1),'b'); 
% plot(bin_array(alldays(2).tt(minds,11),5,1),'g'); 
% plot(bin_array(alldays(2).tt(hinds,11),5,1),'r'); 
%     
% figure; hold on;
% plot(bin_array(abs(dfo(linds)),5,1),'b'); 
% plot(bin_array(abs(dfo(minds)),5,1),'g'); 
% plot(bin_array(abs(dfo(hinds)),5,1),'r'); 
% 
% tblocks = round(linspace(1,700,6));
% for i = 1:5
%     
%     dfob = dfo(tblocks(i):tblocks(i+1));
%     ctt = alldays(2).tt(tblocks(i):tblocks(i+1),:);
%     
%     [dfol1(i), dfol2(i)] = boot_bounds(1000,@nanmean,abs(dfob(ctt(:,3) == 50)),2.5,97.5);
%     [dfom1(i), dfom2(i)] = boot_bounds(1000,@nanmean,abs(dfob(ctt(:,3) == 5)),2.5,97.5);
%     [dfoh1(i), dfoh2(i)] = boot_bounds(1000,@nanmean,abs(dfob(ctt(:,3) == 1)),2.5,97.5);
%     
%     dfolm(i) = nanmean(abs(dfob(ctt(:,3) == 50)));
%     dfomm(i) = nanmean(abs(dfob(ctt(:,3) == 5)));
%     dfohm(i) = nanmean(abs(dfob(ctt(:,3) == 1)));
%     
%     slope_l(i,:) = polyfit(ctt(ctt(:,3) == 50,9),ctt(ctt(:,3) == 50,10),1);
%     slope_m(i,:) = polyfit(ctt(ctt(:,3) == 5,9),ctt(ctt(:,3) == 5,10),1);
%     slope_h(i,:) = polyfit(ctt(ctt(:,3) == 1,9),ctt(ctt(:,3) == 1,10),1);
% end
% 
% figure; hold on;
% plot(1:5,dfolm,'b');
% patch([1:5 fliplr(1:5)],[dfol1' fliplr(dfol2')],'b','FaceAlpha',.5,'EdgeAlpha',.5);
% 
% plot(1:5,dfomm,'g');
% patch([1:5 fliplr(1:5)],[dfom1 fliplr(dfom2)],'g','FaceAlpha',.5,'EdgeAlpha',.5);
% 
% plot(1:5,dfohm,'r');
% patch([1:5 fliplr(1:5)],[dfoh1 fliplr(dfoh2)],'r','FaceAlpha',.5,'EdgeAlpha',.5);
% 
% figure; hold on; 
% plot(1-slope_l(:,1),'b');
% plot(1-slope_m(:,1),'g');
% plot(1-slope_h(:,1),'r');




