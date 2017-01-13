traj = cell(length(seqlikes{1}),1);
for i = 1:length(seqlikes{1})

    traj{i} = seqlikes{1}(i).xorth(1:3,:);

end

alltrajs = horzcat(traj{:})';

[n,V,p] = affine_fit(alltrajs);

nrep = repmat(n',size(alltrajs,1),1);
p0rep = repmat(p,size(alltrajs,1),1);
proj_func_all = @(P,nrep,p0rep) [P(:,1) - nrep(:,1).*(sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2)),...
                   P(:,2) - nrep(:,2).*(sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2)),...
                   P(:,3) - nrep(:,3).*(sum(nrep.*P,2)-sum(nrep.*p0rep,2))./(sum(nrep.^2,2))];

                      
proj_alltraj = proj_func_all(alltrajs,nrep,p0rep);

u_traj = cell(length(seqlikes{2}),1);
dh2l = nan(length(seqlikes{2}),100);

d = -sum(n'.*p);
for j = 1:length(seqlikes{2})
    
    u_traj{j} = seqlikes{2}(j).xorth(1:3,:);
    
    for k = 1:size(u_traj{j},2);
        
        dh2l(j,k) = abs(n(1)*u_traj{j}(1,k) + n(2)*u_traj{j}(2,k) + n(3)*u_traj{j}(3,k) + d)/sqrt(sum(n.^2));
    end
    
end

avd = nanmean(dh2l,2);
sumd = nansum(dh2l,2);
maxd = max(dh2l,[],2);
%%
figure; hold on;
for j = 1:length(seqlikes{2})
    
    if ismember(j,closeH)
        plot3(mean(seqlikes{2}(j).xorth(1,:)),mean(seqlikes{2}(j).xorth(2,:)),mean(seqlikes{2}(j).xorth(3,:)),'.-','color','r'); 
    else
        plot3(mean(seqlikes{2}(j).xorth(1,:)),mean(seqlikes{2}(j).xorth(2,:)),mean(seqlikes{2}(j).xorth(3,:)),'.-','color','b'); 
    end
end

%%
CFD = @(x) -0.5*x(:,1:(end-2)) + 0.5*x(:,3:end);
CFD2 = @(x) x(:,1:(end-2)) -2*x(:,2:(end-1)) + x(:,3:end);

[d,dprime,ddubprime,dfin,f,y,pks,pkloc] = deal(cell(length(seqlikes),1));
figure; hold on; 
for unclev = 1:length(seqlikes)
    for tri = 1:length(seqlikes{unclev})

        x = seqlikes{unclev}(tri).xorth;
        
%         X1 = x(1:3,1); 
%         X2 = x(1:3,end);
        X1 = x(:,1);
        X2 = x(:,end);
        
        for pt = 1:(size(x,2)-1)
            
%             Xc1 = x(
            
%             X0 = x(1:3,pt);
            X0 = x(:,pt);
            %d{unclev}{tri}(:,pt) = sqrt(sum((cross((X0 - X1),(X0 - X2))).^2))/sqrt(sum((X2-X1).^2));
            dfin{unclev}{tri}(:,pt) = sqrt(sum((X0-X2).^2));
        end
%         dprime{unclev}{tri} = CFD(d{unclev}{tri});
%         ddubprime{unclev}{tri} = CFD2(d{unclev}{tri});
        
        [pks{unclev}{tri},pkloc{unclev}{tri}] = findpeaks(dfin{unclev}{tri});
%         
%         L = length(ddubprime{unclev}{tri});
%         NFFT = 2^nextpow2(L);
%         Y = fft(ddubprime{unclev}{tri},NFFT)/L;
%         f{unclev}{tri} = (1/20)/2*linspace(0,1,NFFT/2+1);
%         y{unclev}{tri} = 2*abs(Y(1:NFFT/2+1));

        plot(dfin{unclev}{tri},'Color',cols2plot{unclev});
        %plot(dprime{unclev}{tri},dfin{unclev}{tri}(2:(end-1)),'.','Color',cols2plot{unclev});
    end
end

%%
[projang] = deal(cell(length(seqlikes),1));
figure; hold on; 
r2 = 100;
for unclev = 1:length(seqlikes)
    for tri = 1:length(seqlikes{unclev})

        x = seqlikes{unclev}(tri).xorth;
        
        X1 = x(1:3,1); 
        X2 = x(1:3,end);
        
        X2cent = X2 - X1;
        
        Co = [sum((X2cent).^2), 0 , -r2];
        
        to = roots(Co);
        top = to(to>=0);
        
        Tsphere = top*X2cent;
        
        for pt = 1:(size(x,2)-1)
            
            X01 = x(1:3,pt);
            X01cent = X01 - X1;
            X02 = x(1:3,pt+1);
            X02cent = X02 - X1;
            
            C = [sum((X02cent-X01cent).^2), 2*sum(X01cent.*(X02cent-X01cent)), sum(X01cent.^2)-r2];
            
            t = roots(C);
            tp = t(t>=0);
            
            Psphere = X01cent + tp*(X02cent-X01cent);
            
            projang{unclev}{tri}(:,pt) = atan2(norm(cross(Psphere,Tsphere)),dot(Psphere,Tsphere));
        end
        
        plot(projang{unclev}{tri},'Color',cols2plot{unclev});
    end
end


%%
CFD = @(x) -0.5*x(:,1:(end-2)) + 0.5*x(:,3:end);
CFD2 = @(x) x(:,1:(end-2)) -2*x(:,2:(end-1)) + x(:,3:end);

[d,dprime,ddubprime,dfin,f,y,pks,pkloc] = deal(cell(length(seqlikes),1));
figure; hold on; 
for unclev = 1:length(seqlikes)
    for tri = 1:length(seqlikes{unclev})

        x = seqlikes{unclev}(tri).xorth;

        for pt = 1:(size(x,2)-1)
            
            X1 = x(:,pt);
            XE = x(:,end);
            X2 = x(:,pt+1);
            
            Xtan = X2(2:3) - X1(2:3); Xtan = [Xtan; 0];
            Xend = XE(2:3) - X1(2:3); Xend = [Xend; 0];
            
            d{unclev}{tri}(:,pt) = atan2(norm(cross(Xtan,Xend)), dot(Xtan,Xend));
        end

        plot(d{unclev}{tri},'Color',cols2plot{unclev});
        %plot(dprime{unclev}{tri},dfin{unclev}{tri}(2:(end-1)),'.','Color',cols2plot{unclev});
    end
end

%%
figure; hold on; 
for unclev = 1:length(seqlikes)
    for tri = 1:length(seqlikes{unclev})

        x = seqlikes{unclev}(tri).xorth;

        plot(x(2,end)-x(2,1),x(3,end)-x(3,1),'.','Color',cols2plot{unclev});
        %plot(dprime{unclev}{tri},dfin{unclev}{tri}(2:(end-1)),'.','Color',cols2plot{unclev});
    end
end
