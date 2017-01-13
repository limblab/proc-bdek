%% Plot Correlations
avdiffs = cell(3,1);
for tp = 1:length(AV)
    
avdiffs{1}(:,tp) = EF{tp}{1}(:,3);
avdiffs{2}(:,tp) = EF{tp}{3}(:,3);
avdiffs{3}(:,tp) = EF{tp}{2}(:,3);

% avdiffs{1}(:,tp) = diff(AV{tp}{1},[],2);
% avdiffs{2}(:,tp) = diff(AV{tp}{3},[],2);
% avdiffs{3}(:,tp) = diff(AV{tp}{2},[],2);

av_PD = AV{tp}{1};
av_OD = AV{tp}{2};
av_ORTH = AV{tp}{3};

d_PD = EF{tp}{1}(:,3);
d_OD = EF{tp}{2}(:,3);
d_ORTH = EF{tp}{3}(:,3);

bnd_PD = EF{tp}{1}(:,1:2);
bnd_OD = EF{tp}{2}(:,1:2);
bnd_ORTH = EF{tp}{3}(:,1:2);

% d_PD = av_PD(:,2)-av_PD(:,1);
% d_OD = av_OD(:,2)-av_OD(:,1);
% d_ORTH = av_ORTH(:,2)-av_ORTH(:,1);
% 
% bnd_PD = AV_bnd{1}{1};
% bnd_OD = AV_bnd{1}{2};
% bnd_ORTH = AV_bnd{1}{3};

N_PD = NN{tp}{1};
N_OD = NN{tp}{2};
N_ORTH = NN{tp}{3};

totN = vertcat(NN{1}{:});
size_range = [10 40];%[10 40];

s_PD = metric2markersize(N_PD,totN,size_range);
s_OD = metric2markersize(N_OD,totN,size_range);
s_ORTH = metric2markersize(N_ORTH,totN,size_range);

%%
xlm = [floor(10*min(min(dRES_bnd(G,:))))/10 ceil(10*max(max(dRES_bnd(G,:))))/10];
% xlm = [0,0.3];
% xlm = [-3 1.5];
xtic = xlm(1):0.01:xlm(2);
pclr = 'b';

[pPD,S_PD] = polyfit(dRES(G),d_PD(G),1); 
[pOD,S_OD] = polyfit(dRES(G),d_OD(G),1);
[pORTH,S_ORTH] = polyfit(dRES(G),d_ORTH(G),1);

[~,pdes,pdrs] = regress(d_PD(G),[ones(length(G),1) dRES(G)]);
[~,odes,odrs] = regress(d_OD(G),[ones(length(G),1) dRES(G)]);
[~,orthes,orthrs]=  regress(d_ORTH(G),[ones(length(G),1) dRES(G)]);

Q = [zeros(length(G),1); ones(length(G),1)];
combx = [d_PD(G); d_OD(G)];
combdRES = [dRES(G); dRES(G)];

[~,combint] = regress(combx,[combdRES, combdRES.*Q, ones(length(Q),1), Q]);

[Y_PD,DEL_PD] = polyconf(pPD,xtic,S_PD,'predopt','curve');
[Y_ORTH,DEL_ORTH] = polyconf(pORTH,xtic,S_ORTH,'predopt','curve');
[Y_OD,DEL_OD] = polyconf(pOD,xtic,S_OD,'predopt','curve');

R2_PD = 1 - sum((polyval(pPD,dRES(G)) - d_PD(G)).^2)./sum((mean(d_PD(G)) - d_PD(G)).^2);
R2_OD = 1 - sum((polyval(pOD,dRES(G)) - d_OD(G)).^2)./sum((mean(d_OD(G)) - d_OD(G)).^2);
R2_ORTH = 1 - sum((polyval(pORTH,dRES(G)) - d_ORTH(G)).^2)./sum((mean(d_ORTH(G)) - d_ORTH(G)).^2);

upperlim = ceil(max(reshape([bnd_PD(G,:);bnd_OD(G,:);bnd_ORTH(G,:)],[],1)));
lowerlim = floor(min(reshape([bnd_PD(G,:);bnd_OD(G,:);bnd_ORTH(G,:)],[],1)));

clrs = distinguishable_colors(6); clrs(4,:) = [0 0 0]; clrs(1,:) = [0 1 1];
figure;
subplot(1,3,1); hold on;
for i = 1:length(G)
    plot(dRES(G(i)),d_PD(G(i)),'.','Color',pclr,'MarkerSize',metric2markersize(N_PD(G(i)),N_PD,size_range)); 
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_PD(G(i)),bnd_PD(G(i),:),pclr,.5);
end
% for i = 1:length(G); plot(dRES(G(i)),d_PD(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pPD,xlm),'k','LineWidth',2);
% plot(xtic,Y_PD+DEL_PD,'b');
% plot(xtic,Y_PD-DEL_PD,'b');
patch([xtic fliplr(xtic)],[Y_PD+DEL_PD fliplr(Y_PD-DEL_PD)],'k','FaceAlpha',0.1,'EdgeAlpha',0)

plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pPD(1),R2_PD),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pPD(1),R2_PD),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('SD Neurons','FontSize',18);
xlabel(sprintf('N = [%d %d]',min(N_PD(G)),max(N_PD(G))),'FontSize',18);
%xlabel('\Delta SE_r_e_s','FontSize',14)
ylabel('\Delta FR','FontSize',18);

subplot(1,3,2); hold on;
for i = 1:length(G)
    plot(dRES(G(i)),d_ORTH(G(i)),'.','Color',pclr,'MarkerSize',metric2markersize(N_ORTH(G(i)),N_ORTH,size_range)); 
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_ORTH(G(i)),bnd_ORTH(G(i),:),pclr,.5);
end
% for i = 1:length(G); plot(dRES(G(i)),d_ORTH(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pORTH,xlm),'k','LineWidth',2);
patch([xtic fliplr(xtic)],[Y_ORTH+DEL_ORTH fliplr(Y_ORTH-DEL_ORTH)],'k','FaceAlpha',0.1,'EdgeAlpha',0)
plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pORTH(1),R2_ORTH),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pORTH(1),R2_ORTH),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('ORTH Neurons','FontSize',18);
xlabel(sprintf('N = [%d %d]',min(N_ORTH(G)),max(N_ORTH(G))),'FontSize',18);

subplot(1,3,3); hold on;
for i = 1:length(G)
    plot(dRES(G(i)),d_OD(G(i)),'.','Color',pclr,'MarkerSize',metric2markersize(N_OD(G(i)),N_OD,size_range)); 
    plot_cross(dRES(G(i)),dRES_bnd(G(i),:),d_OD(G(i)),bnd_OD(G(i),:),pclr,.5);
%     plot(dRES(G(i)),d_OD(G(i)),'.','Color',c(round(2.78*i),:),'MarkerSize',18);
end
plot(mean(xlm),0,'r.','MarkerSize',size_range(1));
plot(mean(xlm),0,'r.','MarkerSize',size_range(2));

% for i = 1:length(G); plot(dRES(G(i)),d_OD(G(i)),'.','Color',clrs(Kgrps_G(i),:),'MarkerSize',18); end
plot(xlm,polyval(pOD,xlm),'k','LineWidth',2);
patch([xtic fliplr(xtic)],[Y_OD+DEL_OD fliplr(Y_OD-DEL_OD)],'k','FaceAlpha',0.1,'EdgeAlpha',0)
plot(xlm,[0 0],'k');
% text(0.01,2.5,sprintf('slope: %.2f (R2=%.2f)',pOD(1),R2_OD),'FontSize',16);
text(xlm(1),lowerlim+.2,sprintf('slope: %.1f (R^2=%.2f)',pOD(1),R2_OD),'FontSize',16);
xlim(xlm);
ylim([lowerlim upperlim]);
title('OD Neurons','FontSize',18);
xlabel(sprintf('N = [%d %d]',min(N_OD(G)),max(N_OD(G))),'FontSize',18);

fitE_PD = sqrt(diag((S_PD.R)\inv(S_PD.R'))./S_PD.normr.^2./S_PD.df);
fitE_ORTH = sqrt(diag((S_ORTH.R)\inv(S_ORTH.R'))./S_ORTH.normr.^2./S_ORTH.df);
fitE_OD = sqrt(diag((S_OD.R)\inv(S_OD.R'))./S_OD.normr.^2./S_OD.df);

slopes_TYPES(tp,:) = [pPD(1) pORTH(1) pOD(1)];
Rsquareds(tp,:) = [R2_PD R2_ORTH R2_OD];
slope_LOWS(tp,:) = [pdes(2,1) orthes(2,1) odes(2,1)];
slope_HIGHS(tp,:) = [pdes(2,2) orthes(2,2) odes(2,2)];
slopes_different(tp,1) = combint(2,1) > 0 | combint(2,2) < 0;
%[fitE_PD(1) fitE_ORTH(1) fitE_OD(1)];
if length(AV)>1; close; end
% 
%%
if length(AV)> 1 && tp==length(AV)
    figure; hold on;
    xranges = {-100:100:200};%, 1000+(-100:100:400)};
    xcents = cell2mat(cellfun(@(x) .5*(x(1:end-1) + x(2:end)),xranges,'UniformOutput',0));
    pmarks = {'.-','s-','o-'};
    pols = {'b','g','r'};
    for j = 1:3
        plot(xcents+2*j,slopes_TYPES(:,j),pmarks{j},'Color',pols{j});
        for i = 1:length(xcents)
            plot([1 1]*xcents(i)+2*j,[slope_LOWS(i,j) slope_HIGHS(i,j)],'Color',pols{j});
        end
    end
    
    figure; hold on; 
    for j = 1:3
        plot(xcents+2*j,Rsquareds(:,j),pmarks{j},'Color',pols{j});
    end
end

end
%    