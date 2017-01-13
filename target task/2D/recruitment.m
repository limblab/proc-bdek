pad_base = [base_recruit zeros(1,length(base_recruit))];


Cbase = fft(pad_base,length(pad_base));

[Crecov,recov] = deal(cell(length(av_gain),1));
for i = 1:length(av_gain)
    
    pad_rec = [av_gain{i}, zeros(1,length(base_recruit))];
    
    Crecruit = fft(pad_rec,length(pad_rec));
    
    Crecov{i} = Crecruit./Cbase;
    
    rec = ifft(Crecov{i},length(base_recruit));
    recov{i} = rec(1:length(base_recruit));
end
%%

figure; hold on; 
for i = 1:size(X,1)
    
    plot(repmat(alldays(2).slices(i,:),2,1),[0 max(FUNC(i,:))],'r'); 
    plot([1 1].*alldays(2).tt(i,10),[0 max(FUNC(i,:))],'g');
    
    plot(wrapped_cents,FUNC(i,:),'k');
    plot(wrapped_cents,EFUNC(i,:),'k--');
    
    pause;
    cla;
end