for i = 1:(length(alldays)-1)

    linds = find(alldays(i+1).tt(:,3)==max(alldays(i+1).tt(:,3)));
    hinds = find(alldays(i+1).tt(:,3)==min(alldays(i+1).tt(:,3)));

    pl = polyfit(alldays(i+1).tt(linds,9),alldays(i+1).tt(linds,10),1);
    ph = polyfit(alldays(i+1).tt(hinds,9),alldays(i+1).tt(hinds,10),1);

    fprintf('block %d: slopes = [%.2f %.2f] (diff %.2f)\n',i+1,pl(1),ph(1),pl(1)-ph(1));
end