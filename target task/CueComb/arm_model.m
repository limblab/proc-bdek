l1 = 35;%20; %
l2 = 30;%25; 

theta1 = 0:0.01:pi/2;
theta2 = 0:0.01:pi;

%zeroangs = [1.06 0.77]; %
zeroangs = [1.28 1.57];

[THETA1,THETA2] = meshgrid(theta1,theta2);

X = l1*cos(THETA1) + l2*cos(THETA1+THETA2);
Y = l1*sin(THETA1) + l2*sin(THETA1+THETA2);

data1 = [X(:) Y(:) THETA1(:)];
data2 = [X(:) Y(:) THETA2(:)];

zeropoint = data1(data1(:,3)==zeroangs(1) & data2(:,3)==zeroangs(2),1:2);

%%
angs = .01:0.01:2*pi;
[TS,posits] = deal(zeros(length(angs),2));
for i = 1:length(angs)
    
    XY = zeropoint + [cos(angs(i)) sin(angs(i))];
    
    posits(i,:) = XY;
    
    X = XY(1);
    Y = XY(2);

    c2 = (X.^2 + Y.^2 - l1^2 - l2^2)/(2*l1*l2);
    s2 = sqrt(1-c2.^2);
    THETA2D = atan2(s2,c2);

    k1 = l1 + l2.*c2;
    k2 = l2*s2;
    THETA1D = atan2(Y,X) - atan2(k2,k1);
    
    TS(i,:) = [THETA1D, THETA2D];
    
end
figure; hold on; subplot(1,2,1);
plot(angs,TS(:,1),'b'); title('Shoulder','FontSize',18); xlim([0 2*pi]);
subplot(1,2,2);
plot(angs,TS(:,2),'r'); title('Elbow','FontSize',18); xlim([0 2*pi]);
%%
zerosrep = repmat(zeroangs,length(TS),1);
deltaT = TS-zerosrep;

relativeElbow = abs(deltaT(:,2))-abs(deltaT(:,1));

ang2relElbow = @(x) relativeElbow(ceil(100*x));
%% Check
THETA1 = TS(:,1);
THETA2 = TS(:,2);

X = l1*cos(THETA1) + l2*cos(THETA1+THETA2);
Y = l1*sin(THETA1) + l2*sin(THETA1+THETA2);

ps = [X(:), Y(:)];

if 1
    %%
    inds = bumps{1};
    [Delta_reach, Delta_bump] = deal(zeros(length(inds),2));
    for i = 1:length(inds)
        
        % bump
        bumpdir = bump_ang(inds(i));
        bumpmag = bump_mag(inds(i));
        chosendir = tt(inds(i),12);

        % calculate actual joint changes (Inverse Kinematics)
        X = zeropoint(1) + bumpmag*cos(bumpdir);
        Y = zeropoint(2) + bumpmag*sin(bumpdir);

        c2 = (X.^2 + Y.^2 - l1^2 - l2^2)/(2*l1*l2);
        s2 = sqrt(1-c2.^2);
        THETA2D = atan2(s2,c2);

        k1 = l1 + l2.*c2;
        k2 = l2*s2;
        THETA1D = atan2(Y,X) - atan2(k2,k1);
    
    	ANGS = [THETA1D, THETA2D]; 
        Delta_bump(i,:) = ANGS - zeroangs;
    
        % calculate perceived joint changes (Inverse Kinematics)
        X = zeropoint(1) + bumpmag*cos(chosendir);
        Y = zeropoint(2) + bumpmag*sin(chosendir);

        c2 = (X.^2 + Y.^2 - l1^2 - l2^2)/(2*l1*l2);
        s2 = sqrt(1-c2.^2);
        THETA2D = atan2(s2,c2);

        k1 = l1 + l2.*c2;
        k2 = l2*s2;
        THETA1D = atan2(Y,X) - atan2(k2,k1);
    
    	ANGS = [THETA1D, THETA2D]; 
        Delta_reach(i,:) = ANGS - zeroangs;
        Rang(i,:) = ANGS;
    
    end
    %%
    
end
Rdir = atan2(Delta_reach(:,2),Delta_reach(:,1));
Bdir = atan2(Delta_bump(:,2),Delta_bump(:,1));

%% Angle changes
figure; hold on;

badeg = bump_ang(inds)/pi*180;
DBdeg = Delta_bump/pi*180;
DRdeg = Delta_reach/pi*180;
Edeg = DRdeg-DBdeg;

subplot(1,2,1); hold on; 
plot(badeg,DBdeg(:,1),'k.'); 
plot(badeg,DRdeg(:,1),'b.'); 
xlabel('Bump Direction (deg)','FontSize',16);
ylabel('Change in Joint (deg)','FontSize',16);
set(gca,'FontSize',14);
xlim([-180 180]);
ylim([-2.5 2.5]);
title('Shoulder','FontSize',18);

subplot(1,2,2); hold on; 
plot(badeg,DBdeg(:,2),'k.'); 
plot(badeg,DRdeg(:,2),'r.'); 
xlabel('Bump Direction (deg)','FontSize',16);
ylabel('Change in Joint (deg)','FontSize',16);
set(gca,'FontSize',14);
xlim([-180 180]);
ylim([-2.5 2.5]);
title('Elbow','FontSize',18);

%% -% Error
figure; hold on; 
subplot(1,2,1); hold on; 
plot(badeg,Edeg(:,1),'b.'); 
xlabel('Bump Direction (deg)','FontSize',16);
ylabel('Joint Error (deg)','FontSize',16);
set(gca,'FontSize',14);
xlim([-180 180]);
ylim([-1.5 1.5]);
title('Shoulder','FontSize',18);

subplot(1,2,2); hold on; 
plot(badeg,Edeg(:,2),'r.'); 
xlabel('Bump Direction (deg)','FontSize',16);
ylabel('Joint Error (deg)','FontSize',16);
set(gca,'FontSize',14);
xlim([-180 180]);
ylim([-1.5 1.5]);
title('Elbow','FontSize',18);

%% -% Error v Joint Change Magnitude
figure; hold on; 
subplot(1,2,1); hold on; 
plot(abs(DBdeg(:,1)),Edeg(:,1),'b.'); 
xlabel('|Change in Joint| (deg)','FontSize',16);
ylabel('Joint Error (deg)','FontSize',16);
set(gca,'FontSize',14);
title('Shoulder','FontSize',18);

subplot(1,2,2); hold on; 
plot(abs(DBdeg(:,2)),Edeg(:,2),'r.'); 
xlabel('|Change in Joint| (deg)','FontSize',16);
ylabel('Joint Error (deg)','FontSize',16);
set(gca,'FontSize',14);
title('Elbow','FontSize',18);

%% -% Error v Other Joint
figure; hold on; 
subplot(1,2,1); hold on; 
plot(DBdeg(:,2),Edeg(:,1),'b.'); 
xlabel('Change in Elbow Joint (deg)','FontSize',16);
ylabel('Joint Error (deg)','FontSize',16);
set(gca,'FontSize',14);
title('Shoulder','FontSize',18);

subplot(1,2,2); hold on; 
plot(DBdeg(:,1),Edeg(:,2),'r.'); 
xlabel('Change in Shoulder Joint (deg)','FontSize',16);
ylabel('Joint Error (deg)','FontSize',16);
set(gca,'FontSize',14);
title('Elbow','FontSize',18);

%%

%err_color = @(x) floor((x-min(x))./(max(x)-min(x))*63)+1;
err_color = @(x) floor(x.*(32/max(abs(x)))+32);

get_sens = @(x) interp1(angs,dRdD,x);
