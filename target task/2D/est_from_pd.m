lowind = find(comp_tt(:,3)==max(comp_tt(:,3)));
highind = find(comp_tt(:,3)==min(comp_tt(:,3)));

figure; hold on; plot([0 60],[0 60],'k--');
for i = 1:length(lowind)
    
    target = comp_tt(lowind(i),2);
    
    estimate = mean(pd(:,2:3),2) + 0.5.*(pd(:,3)-pd(:,2)).*cos(target-pd(:,1));
    
    estimate = estimate';
    
    actual = rastbin(lowind(i),:);
    
    error = (estimate - actual)./pd(:,3)';
    
    eL(i) = sum(error.^2);
    esL(i,:) = error;
    
    plot(actual,estimate,'b.')
    drawnow; %pause;
    
end

for i = 1:length(highind)
    
    target = comp_tt(highind(i),2);
    
    estimate = mean(pd(:,2:3),2) + 0.5.*(pd(:,3)-pd(:,2)).*cos(target-pd(:,1));
    
    estimate = estimate';
    
    actual = rastbin(highind(i),:);
    
    error = (estimate - actual)./pd(:,3)';
    
    eH(i) = sum(error.^2);
    esH(i,:) = error;
end
    
figure; errorbar(1,mean(eL),std(eL),'b.'); hold on;
errorbar(2,mean(eH),std(eH),'r.');

    