met2col = @(X) floor(63*(X - min(X(:)))./(max(X(:))-min(X(:))))+1;

% FA = F{3}; FA2 = reshape(smooth(FA),size(FA,1),size(FA,2));
% FA3 = reshape(smooth(FA2'),size(FA2,2),size(FA2,1))';

% G = met2col(FA3); G(end,:) = []; c = colormap;%('jet');

G = met2col(FA3); G(end,:) = []; c = colormap;%('jet');

rstep = 1;
figure; hold on;

ds = 2*pi/size(G,1);
for i = 1:size(G,1)
    
    ths = [spatial_cents(i)-ds/2 spatial_cents(i)+ds/2];
    
    for ts = 1:size(G,2)
        
        rs = rstep*[ts ts+1];
        
        patch([rs(1).*cos(ths(1):pi/80:ths(2)) fliplr(rs(2).*cos(ths(1):pi/80:ths(2)))],...
              [rs(1).*sin(ths(1):pi/80:ths(2)) fliplr(rs(2).*sin(ths(1):pi/80:ths(2)))],...
              c(G(i,ts),:),'EdgeAlpha',0);
    end
end