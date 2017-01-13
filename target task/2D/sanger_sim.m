num_units = 10;
thetas = 0:0.02:2*pi;
tc = zeros(num_units,length(thetas));
pd = zeros(1,num_units);
baseline = zeros(1,num_units);
modu = zeros(1,num_units);

for i =1:num_units
    
    pd(i) = mod(i*2*pi/num_units + pi/10*rand-pi/5,2*pi);%2*pi*rand;
    baseline(i) = 20+30*rand;
    modu(i) = baseline(i)*(0.25+0.50*rand);
    
    tc(i,:) = modu(i)*cos(thetas-pd(i))+baseline(i);
end

maxfire = baseline+modu; 

%%
dt = 0.1;
plan_dir = pi;
pop = baseline + modu.*cos(plan_dir - pd);
spike_counts = round(pop*dt)';

%% Point Estimate (Sanger)
N = diff(thetas(1:2));
cp = exp(sum(repmat(spike_counts,1,length(thetas)).*log(tc) - tc*dt,1));
sanger_pop = cp./(sum(cp)*N);

%% Point Estimate (Bastian)

%Normalize tc
tcmin = repmat(min(tc,[],2),1,length(thetas));
tcmax = repmat(max(tc,[],2),1,length(thetas));
tcnorm = (tc - tcmin)./(tcmax-tcmin);

tunemax = tcmax(:,1);
tunemin = tcmin(:,1);
tunemodu = tunemax-tunemin;
spike_norm = repmat((spike_counts/dt)./tunemodu - tunemin./tunemodu,1,length(thetas));

bastian_pop = sum(spike_norm.*tcnorm,1)./num_units;

figure; hold on;
%plot(thetas,sanger_pop,'b');
plot(thetas,bastian_pop,'r');
plot([plan_dir plan_dir],[0 1],'g');

%%
figure; hold on;
for i = 1:size(spike_norm,1)
    act = spike_norm(i,1);
    plot((spike_norm(i,:).*tcnorm(i,:))','LineWidth',3,'Color',[1 1-act 1-act]);      
end
bastian_pop = sum(spike_norm.*tcnorm,1)./num_units;
plot(bastian_pop,'k','LineWidth',4);
%%

% ya = 0.25*ones(1,30) ;
% yb = ya(end) + 1.5*exp(0:0.01:.05)-1.5 ;
% yc = yb(end) + .02*(exp(1:-0.5:0)-exp(1)) ;
% yf = yc(end)*ones(1,10);
% yd = yc(end) + .2*exp(0:0.05:.95)-.2;
% ye = yd(end) - 3*(0:0.01:0.11);

ya = 0.25*ones(1,30) ;
yb = ya(end) + 5*(exp(0:0.01:.05)-1);
yc = yb(end) + .1*(exp(1:-0.5:0)-exp(1)) ;
yf = yc(end)*ones(1,10);
yd = yc(end) + .2*exp(0:0.05:.95)-.2;
ye = yd(end) - 3*(0:0.01:0.11);

y = [ya yb yc yf yd ye];
sy = smooth(y);

x = linspace(-800,1200,length(y));
xd = -800:100:1200;

ysamp = interp1(x,sy,xd);

y2 = 0.125+0.5*y+0.01*rand(1,length(y));
sy2 = smooth(y2);
ysamp2 = interp1(x,sy2,xd);

figure; hold on;
plot(xd,ysamp,'Color',[1 0.6 0.6],'LineWidth',3);
plot(xd,ysamp2,'Color',[0.6 0.6 1],'LineWidth',3);
%plot(xd,ysamp3,'k','LineWidth',2);
xlabel('Time from target on (ms)','FontSize',16); 
ylabel('DPA Concentration','FontSize',16);
legend('High visual noise','Low visual noise');
plot([0 0],[0 1],'k--');
plot([1000 1000],[0 1],'k-.');

ylim([0 .7])
xlim([-500 1100]);

%%

ya = 0.25*ones(1,30) ;
yb = ya(end) + 5*(exp(0:0.01:.05)-1) ;
yc = yb(end) + .1*(exp(1:-0.5:0)-exp(1)) ;
yf = yc(end)*ones(1,10);
yd = yc(end) + .2*(exp(0:0.05:.95)-1);
ye = yd(end) - 3*(0:0.01:0.11);


y = [ya yb yc yf yd ye];
y4 = [ya+.1 yb yc yf yd ye];
y4 = 2*y4 - 2*y(1) + 0.08*rand(1,length(y4));

y5 = [ya+.05 yb yc yf yd ye];
y5 = 2*y5 - 2*y(1) + 0.08*rand(1,length(y5));

y3 = 2*y - 1.9*y(1);
sy3 = smooth(y3);
ysamp3 = interp1(x,sy3,xd);

sy4 = smooth(y4);
ysamp4 = interp1(x,sy4,xd);
ysamp4(8) = ysamp4(7);

sy5 = smooth(y5);
ysamp5 = interp1(x,sy5,xd);
ysamp5(8) = ysamp5(7);

figure; hold on;
%plot(xd,ysamp3,'Color',[0.2 0.2 0.2],'LineWidth',3);
plot(xd,ysamp4,'Color',[0.1 .8 .1],'LineWidth',3);
plot(xd,ysamp5,'Color',[0.8 0.1 0.1],'LineWidth',3);
xlabel('Time from target on (ms)','FontSize',16); 
ylabel('DPA Concentration','FontSize',16);
legend('Narrow Prior','Broad Prior');
plot([0 0],[0 1],'k--');
plot([1000 1000],[0 1],'k-.');
%ylim([0 .7])
xlim([-500 1100]);
%% another try

ya = 0.25*ones(1,30) ;
yb = ya(end) + 5*(exp(0:0.01:.05)-1) ;
yc = yb(end) + .1*(exp(1:-0.5:0)-exp(1)) ;
yf = yc(end)*ones(1,10);
yd = yc(end) + .2*(exp(0:0.05:.95)-1);
ye = yd(end) - 3*(0:0.01:0.11);

y = [ya yb yc yf yd ye];
sy = smooth(y);

x = linspace(-800,1200,length(y));
xd = -800:100:1200;

ysamp = interp1(x,sy,xd);

yn = 2*y - 1.9*y(1);
syn = smooth(yn);
ysampn = interp1(x,syn,xd);

figure; hold on;
plot(xd,ysampn,'Color',[1 0.6 0.6],'LineWidth',3);

xlabel('Time from target on (ms)','FontSize',16); 
ylabel('DPA Concentration','FontSize',16);
legend('High visual noise','Low visual noise');
plot([0 0],[0 1],'k--');
plot([1000 1000],[0 1],'k-.');

ylim([0 .7])
xlim([-500 1100]);

%%
dirs = 0:0.01:2*pi;
thets = dirs/pi*180;

n1 = 0.5*(cos(dirs - pi/2)+1);
n2 = 0.5*(cos(dirs - 6*pi/4)+1);
n3 = 0.5*(cos(dirs - pi)+1);

figure; plot(thets,n1,'Linewidth',6); box off; set(gca,'Ytick',[]);set(gca,'Xtick',[]);
xlim([0 360]);ylim([0 1.1]);
figure; plot(thets,n2,'Linewidth',6); box off; set(gca,'Ytick',[]);set(gca,'Xtick',[]);
xlim([0 360]);ylim([0 1.1]);
figure; plot(thets,n3,'Linewidth',6); box off; set(gca,'Ytick',[]);set(gca,'Xtick',[]);
xlim([0 360]);ylim([0 1.1]);



figure; plot(thets,0.2*n1,'r','Linewidth',6); box off; set(gca,'Ytick',[]);set(gca,'Xtick',[]);
xlim([0 360]);ylim([0 1.1]);
figure; plot(thets,0.2*n2,'r','Linewidth',6); box off; set(gca,'Ytick',[]);set(gca,'Xtick',[]);
xlim([0 360]);ylim([0 1.1]);
figure; plot(thets,n3,'r','Linewidth',6); box off; set(gca,'Ytick',[]);set(gca,'Xtick',[]);
xlim([0 360]);ylim([0 1.1]);


figure; plot(thets,[n1;n2;n3])