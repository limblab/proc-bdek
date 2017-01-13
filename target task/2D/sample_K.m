numslices = 5;

% Draw samples with specific SAMPLE concentrations
%ks = exp(0:0.1:4)' - 0.9;
ks = [25];

samps = zeros(1,numslices);
[SLICE, SLICE_N] = deal(cell(length(ks),1));
for i = 1:length(ks)
    for j = 1:1000
        fprintf('k: %d/%d (%d/1000)\n',i,length(ks),j);
        
        keep = false;
        while ~keep
            
            for nums = 1:numslices
                samps(nums) = vonmisrand(ks(i));
            end
            
            if abs(circ_kappa(samps) - ks(i)) < 0.05*ks(i)
                keep = true;
            end
        end
        clc;
        SLICE{i}(j,:) = samps;
        
        for nums2 = 1:numslices
            samps(nums2) = vonmisrand(ks(i));
        end
        SLICE_N{i}(j,:) = samps;
    end
end

%%
[c, cn, kc, kn] = deal(zeros(length(ks),length(SLICE{1})));
[ck, cnk] = deal(zeros(length(ks),1));
for i = 1:length(ks)   
    for j = 1:length(SLICE{i})     
        c(i,j) = circ_mean(SLICE{i}(j,:)');
        cn(i,j) = circ_mean(SLICE_N{i}(j,:)');
        
        kc(i,j) = circ_kappa(SLICE{i}(j,:));
        kn(i,j) = circ_kappa(SLICE_N{i}(j,:));
    end
    ck(i) = circ_kappa(c(i,:));   
    cnk(i) = circ_kappa(cn(i,:));
end

%% Plot low examples
h1 = figure; hold on;
ths = 0:0.01:2*pi;

for i = 1:length(SLICE{1})
  
    plot(7*cos(ths),7*sin(ths),'k'); 
    plot(9*cos(ths),9*sin(ths),'k');
    axis('equal');
    
    sl = SLICE{2}(i,:) + pi/2;
    
    plot([7*cos(sl); 9*cos(sl)],...
         [7*sin(sl); 9*sin(sl)],'b','LineWidth',3); 
     
    plot([7*cos(pi/2); 9*cos(pi/2)],...
         [7*sin(pi/2); 9*sin(pi/2)],'k','LineWidth',2); 
    
    pause; cla;
end

%% plot high examples
h2 = figure; hold on;
ths = 0:0.01:2*pi;

plot(7*cos(ths),7*sin(ths),'k'); axis('equal');
plot(9*cos(ths),9*sin(ths),'k');

for i = 1:length(SLICE{2})
     
    plot(7*cos(ths),7*sin(ths),'k'); axis('equal');
    plot(9*cos(ths),9*sin(ths),'k');
    
    sh = SLICE{1}(i,:) + pi/2;

    plot([7*cos(sh); 9*cos(sh)],...
         [7*sin(sh); 9*sin(sh)],'r','LineWidth',3); 

    plot([7*cos(pi/2); 9*cos(pi/2)],...
    [7*sin(pi/2); 9*sin(pi/2)],'k','LineWidth',2); 
    
    pause; cla;
end

%%
slice1 = sli{2}(6,:);

sh = slice1 + (pi/2 - circ_mean(slice1'));


%[m,mh,ml] = circ_mean(sh');
[ml,mh] = boot_bounds(10000,@circ_mean,sh',2.5,97.5);
m = circ_mean(sh');

h2 = figure; hold on;
ths = 0:0.01:2*pi;

plot(7*cos(ths),7*sin(ths),'k'); axis('equal');
plot(9*cos(ths),9*sin(ths),'k');
     
plot(7*cos(ths),7*sin(ths),'k');
plot(9*cos(ths),9*sin(ths),'k'); axis off;

 p = ang2poly(circ_mean(sh'),2*pi,7,9);    
 patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','FaceColor','b','FaceAlpha',1);
%plot(8*cos(ths),8*sin(ths),'b','LineWidth',45); 


% pm = ang2poly(circ_mean(sh'),abs(circ_dist(ml,mh)),7,9);
% patch(pm(:,1),pm(:,2),pm(:,3),'EdgeColor','none','FaceColor','g','FaceAlpha',0.9);


plot([7*cos(sh); 9*cos(sh)],...
     [7*sin(sh); 9*sin(sh)],'r','LineWidth',5); 

plot([7*cos(pi/2); 9*cos(pi/2)],...
[7*sin(pi/2); 9*sin(pi/2)],'k','LineWidth',2); 

%
uncK = circ_kappa(sh)*5;

unc_s = zeros(10000,1);
for i = 1:10000
    unc_s(i) = pi/2 + vonmisrand(uncK);
end
unc_ss = sortrows(unc_s);

lowidx = round(2.5/100*10000);
highidx = round(97.5/100*10000);

b95_l = unc_ss(lowidx);
b95_h = unc_ss(highidx);

plot(6.9*cos(b95_l:0.01:b95_h),6.9*sin(b95_l:0.01:b95_h),'r.-','LineWidth',6);

plot(6.9*cos(ml:0.01:mh),6.9*sin(ml:0.01:mh),'g.-','LineWidth',6);
plot(6.9*cos(m),6.9*sin(m),'k.','MarkerSize',20);

%%
slice1 = sli{11}(1,:);

sh = slice1 + (pi/2 - circ_mean(slice1'));


%[m,mh,ml] = circ_mean(sh');
[ml,mh] = boot_bounds(10000,@circ_mean,sh',2.5,97.5);
m = circ_mean(sh');

h2 = figure; hold on;
ths = 0:0.01:2*pi;

plot(7*cos(ths),7*sin(ths),'k'); axis('equal');
plot(9*cos(ths),9*sin(ths),'k');
     
plot(7*cos(ths),7*sin(ths),'k');
plot(9*cos(ths),9*sin(ths),'k'); axis off;

 p = ang2poly(circ_mean(sh'),2*pi,7,9);    
 patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','FaceColor','b','FaceAlpha',1);
%plot(8*cos(ths),8*sin(ths),'b','LineWidth',45); 


% pm = ang2poly(circ_mean(sh'),abs(circ_dist(ml,mh)),7,9);
% patch(pm(:,1),pm(:,2),pm(:,3),'EdgeColor','none','FaceColor','g','FaceAlpha',0.9);


plot(8*cos(sh), 8*sin(sh),'r.','MarkerSize',50); 

plot([7*cos(pi/2); 9*cos(pi/2)],...
[7*sin(pi/2); 9*sin(pi/2)],'k','LineWidth',2); 





