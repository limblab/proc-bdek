muf = [-10 10];
slices = 10;
targetsize = 2;

mup = 2;
sigp = 1;

sigf1 = 1/sqrt(10);
sigf2 = 4/sqrt(10);

s1 = sigf1^2./(sigf1^2 + sigp^2);

pos1 = (muf.*sigp^2 + mup*sigf1^2)./(sigp.^2+sigf1.^2);
pos2 = (muf.*sigp^2 + mup*sigf2^2)./(sigp.^2+sigf2.^2);

sig1 = sqrt((sigp.^2*sigf1.^2)./(sigp.^2+sigf1.^2));
sig2 = sqrt((sigp.^2*sigf2.^2)./(sigp.^2+sigf2.^2));

figure; hold on; 
plot(muf,pos1,'b','LineWidth',3); 
plot(muf,pos2,'r','LineWidth',3);
plot(muf,pos1+sig1,'b--'); plot(muf,pos1-sig1,'b--','LineWidth',1);
plot(muf,pos2+sig2,'r--'); plot(muf,pos2-sig2,'r--','LineWidth',1);

legend(sprintf('slope: %.3f',diff(pos1)/diff(muf)),...
       sprintf('slope: %.3f',diff(pos2)/diff(muf)));
   
dt = 0.05;
priordist = normpdf(-10:dt:10,mup,sigp);
   
plot(-10:dt:10,priordist-5,'k');   
plot(priordist-10,-10:dt:10,'k');  

axis([-10 10 -5 5]);

target = ones(1,targetsize./dt);

broadtarg = conv(priordist,target);
broadtarg = broadtarg(targetsize/dt/2:end-targetsize/dt/2)./(targetsize./dt);
broadstd = 1./(max(broadtarg)*sqrt(2*pi));

finalstd = (sigp^2*broadstd^2)./(sigp^2+broadstd^2);

