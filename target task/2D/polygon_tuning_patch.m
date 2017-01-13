function[extent] = polygon_tuning_patch(thetas,rates,colr,style,falpha,unitname,varargin)


hold on; 
plotcol = colr;

extent = 10*ceil(max(rates(:))/10);

plot([-extent extent],[0 0],'k','LineWidth',.5);
plot([0 0],[-extent extent],'k','LineWidth',.5);

datax = [thetas, thetas(1)];
datay = [rates(:,1)', rates(1,1)];
datar = [datay.*cos(datax); datay.*sin(datax)];

plot(datar(1,:), datar(2,:),'LineStyle',style,'Color',plotcol,'LineWidth',2);

thetas2 = thetas;

rates2 = rates(:,2)';
rates3 = rates(:,3)';

ratesL = rates2;
ratesH = rates3;

ratesL(isnan(rates2+rates3)') = [];
ratesH(isnan(rates2+rates3)') = [];
thetas2(isnan(rates2+rates3)') = [];

if length(thetas2)==length(thetas)
    datayL = [ratesL ratesL(1)];
    datayH = [ratesH ratesH(1)];
    datax2 = [thetas2 thetas2(1)];    
else
    datayL = ratesL;% ratesL(1)];
    datayH = ratesH;% ratesH(1)];
    datax2 = thetas2;% thetas2(1)];
end


datarL = [datayL.*cos(datax2); datayL.*sin(datax2)]; 
datarH = [datayH.*cos(datax2); datayH.*sin(datax2)];

patch([datarL(1,:) fliplr(datarH(1,:))],[datarL(2,:) fliplr(datarH(2,:))],colr,'FaceAlpha',falpha);

title(sprintf('Unit %.1f',unitname),'FontSize',14);
end