times = {[-100 800],[0 500]};
aligns = [5 6];

pdi = 2;
TT = alldays(pdi).tt;

cop = {'b','g','r'};
likes = flipud(unique(TT(:,3)));
figure; hold on;
sB = cell(length(times),1);
for ti = 1:length(times)
    for li = 1:length(likes)
        
        l_inds = find(TT(:,3)==likes(li)); 
        speeds= kin_exam(alldays(1).bdfM,TT(l_inds,:),times{ti}(1),times{ti}(2),aligns(ti));
        
        xplots = (times{ti}(1):times{ti}(2)) + (times{1}(end)+200)*(ti-1);
        
        [sB{ti}(1,:),sB{ti}(2,:)] = boot_bounds(1000,@nanmean,speeds,2.5,97.5);
        
        plot(xplots,nanmean(speeds,1),cop{li});
        patch([xplots fliplr(xplots)],[sB{ti}(1,:) sB{ti}(2,end:-1:1)],cop{li},'FaceAlpha',0.25,'EdgeAlpha',0);
        
        clc; fprintf('%d/%d - %d/%d\n',ti,length(times),li,length(likes));
    end
end
%%   
figure; hold on;
trace_sampy = cell(length(likes),1);
for li = 1:length(likes)

    l_inds = find(TT(:,3)==likes(li)); 
    
    for tri = 1:length(l_inds)
        
        i1 = find(alldays(1).bdfM.pos(:,1) > TT(l_inds(tri),6),1,'first');
        i2 = find(alldays(1).bdfM.pos(:,1) < TT(l_inds(tri),7),1,'last');
        
        xtrace = alldays(1).bdfM.pos(i1:i2,2)';
        ytrace = alldays(1).bdfM.pos(i1:i2,3)';
        
        xtr = xtrace - xtrace(1);
        ytr = ytrace - ytrace(1);
        
        movang = atan2(ytr(end),xtr(end));
        
        rotmax = [cos(-movang), -sin(-movang); sin(-movang), cos(-movang)];
        
        trace_rot = rotmax*[xtr; ytr];
        trace_rotsx = smooth(trace_rot(1,:))';
        trace_rotsy = smooth(trace_rot(2,:))';

        [n,bin] = histc(trace_rotsx,unique(trace_rotsx));
        dups = find(n>1);
        nonduploc = find(~ismember(bin,dups));
        trace_rotsu = [trace_rotsx(nonduploc);trace_rotsy(nonduploc)];

        trace_sampy{li}(tri,:) = interp1(trace_rotsu(1,:),trace_rotsu(2,:),linspace(trace_rotsu(1,1),trace_rotsu(1,end),1000));
        clc; fprintf('%d/%d - %d/%d\n',ti,length(times),li,length(likes));

    end
    
    badtri = sum(abs(diff(trace_sampy{li}'))) > 10;
    trace_sampy{li}(badtri,:) = [];
    
    plot(1:1000,nanmean(trace_sampy{li},1),cop{li});
    [ltrace,utrace] = boot_bounds(1000,@nanmean,trace_sampy{li},2.5,97.5);
    
    patch([1:1000 1000:-1:1],[ltrace' fliplr(utrace')],cop{li},'FaceAlpha',0.25,'EdgeAlpha',0);
    title('Deviation');
end

%%
figure; hold on;
[av_dev,av_absdev] = deal(cell(length(likes),1));

for li = 1:length(likes)
    
    plot(1:1000,nanmean(abs(trace_sampy{li}),1),cop{li});
    [ltrace,utrace] = boot_bounds(1000,@nanmean,abs(trace_sampy{li}),2.5,97.5);

    patch([1:1000 1000:-1:1],[ltrace' fliplr(utrace')],cop{li},'FaceAlpha',0.25,'EdgeAlpha',0);
    
    title('Absolute deviation');
    
    av_absdev{li} = nanmean(abs(trace_sampy{li}),2);
    av_dev{li} = nanmean(trace_sampy{li},2);
    clc; fprintf('%d/%d\n',li,length(likes));

end


        
                