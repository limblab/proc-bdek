
clear Fr Fr2 Fr3 Fr4

for i = 1:size(F,1); Fr(i,:) = interp1(1:size(F,2),F(i,:),1:.5:size(F,2)); end
for i = 1:size(Fr,2); Fr2(:,i) = interp1(1:size(F,1),Fr(:,i)',1:.5:size(F,1)); end
for i = 1:size(Fr2,1); Fr3(i,:) = smooth(Fr2(i,:),20); end
for i = 1:size(Fr3,2); Fr4(:,i) = circ_smooth(Fr3(:,i),10); end


met2col = @(X) floor(63*(X - min(X(:)))./(max(X(:))-min(X(:))))+1;
met2col2 = @(X,S) floor(63*(X - min(S(:)))./(max(S(:))-min(S(:))))+1;
% FA = F{3}; FA2 = reshape(smooth(FA),size(FA,1),size(FA,2));
% FA3 = reshape(smooth(FA2'),size(FA2,2),size(FA2,1))';

% G = met2col(FA3); G(end,:) = []; c = colormap;%('jet');

% G = met2col(Fr4); G(end,:) = []; c = colormap('jet');
G = met2col2(Fr4,SC); G(end,:) = []; c = colormap('jet');
rstep = 1;
figure; hold on;

ds = 2*pi/size(G,1);
scent = linspace(-pi,pi,size(G,1)+1); scent(end) = [];
for i = 1:size(G,1)
    
%     ths = [spatial_cents(i)-ds/2 spatial_cents(i)+ds/2];
    ths = [scent(i)-ds/2 scent(i)+ds/2];
    
    for ts = 1:size(G,2)
        
        rs = rstep*[ts ts+1];
        
        patch([rs(1).*cos(ths(1):pi/100:ths(2)) fliplr(rs(2).*cos(ths(1):pi/100:ths(2)))],...
              [rs(1).*sin(ths(1):pi/100:ths(2)) fliplr(rs(2).*sin(ths(1):pi/100:ths(2)))],...
              c(G(i,ts),:),'EdgeColor','none');
          
        patch([rs(1).*cos([ths(1),ths(2)]) fliplr(rs(2).*cos([ths(1),ths(2)]))],...
            [rs(1).*sin([ths(1),ths(2)]) fliplr(rs(2).*sin([ths(1),ths(2)]))],...
            c(G(i,ts),:),'EdgeColor','none');
    end
end