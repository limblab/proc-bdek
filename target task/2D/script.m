%%

num_units = 10;
thetas = 0:0.02:2*pi;
tc = zeros(num_units,length(thetas));
pd = zeros(1,num_units);
pd2 = pd;
baseline = zeros(1,num_units);
modu = zeros(1,num_units);

for i =1:num_units
    
    pd(i) = mod(i*2*pi/num_units + pi/5*rand-pi/2.5,2*pi);%2*pi*rand;
    pd2(i) = mod(i*2*pi/num_units + pi/5*rand-pi/2.5,2*pi);
       
    tc(i,:) = 0.5*cos(thetas-pd(i))+0.5;
end
plan_dir = pi;

pop = 0.5*cos(plan_dir - pd)+0.5;
pop2 = 1.5*(cos(plan_dir - pd2)+1);

figure; hold on;
for i = 1:num_units
    
    plot(thetas/pi*180,pop(i)*tc(i,:),'LineWidth',3,'Color',[0.9 1-pop(i)/3 1-pop(i)/3]);
end

minv = min(sum(repmat(pop',1,size(tc,2)).*tc,1)./2);
minv2 = min(sum(repmat(pop2',1,size(tc,2)).*tc,1)./2);

plot(thetas/pi*180,sum(repmat(pop',1,size(tc,2)).*tc,1)./2 + 4-minv,'k','LineWidth',4);
ylim([0 8.5]);
xlim([0 360]);
xlabel('Direction (degrees)','FontSize',16);
set(gca,'ytick',[])

figure; hold on;
for i = 1:num_units
    
    plot(thetas/pi*180,pop2(i)*tc(i,:),'LineWidth',3,'Color',[1-pop2(i)/3 1-pop2(i)/3 0.9]);
end
plot(thetas/pi*180,sum(repmat(pop2',1,size(tc,2)).*tc,1)./2 + 4-minv2,'k','LineWidth',4);
ylim([0 8.5]);
xlim([0 360]);
xlabel('Direction (degrees)','FontSize',16);
set(gca,'ytick',[]);
