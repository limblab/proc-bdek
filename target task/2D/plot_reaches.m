sess = alldays(2);
sessbdf = alldays(1);
Bayes_estimates = alldays(2).tt(:,2);

figure; hold on;
depart_fromB = circ_dist(sess.tt(:,10),Bayes_estimates);

trials2use = 1:length(sess.tt);
%trials2use = misfits;
for i = 1:length(trials2use)
    
    subplot(1,2,1); hold on;
    
    movei = trials2use(i);
    
    ths = 0:0.01:2*pi;
    sl = sess.slices(movei,:);
    
    plot(7*cos(ths),7*sin(ths),'k'); axis('equal');
    plot(9*cos(ths),9*sin(ths),'k');
    
    if sess.tt(movei,3) == 5
        cop = 'r';
    elseif sess.tt(movei,3) > 50
        cop = 'g';
    else
        cop = 'b';
    end
    
    plot([7*cos(sl); 9*cos(sl)],[7*sin(sl); 9*sin(sl)],cop,'LineWidth',4); 
    
    if abs(circ_dist(sess.tt(movei,10),sess.tt(movei,2))) < (7.5/180*pi)
        plot(7*cos(sess.tt(movei,10)),7*sin(sess.tt(movei,10)),'y.','MarkerSize',25);
    else
        plot(7*cos(sess.tt(movei,10)),7*sin(sess.tt(movei,10)),'yo');
    end
    
    
    bdfstart = find(BDF.pos(:,1) > sess.tt(movei,6),1,'first');
    %bdfend = find(BDF.pos(:,1) > sess.tt(i,7),1,'first');
    bdfend = find(BDF.pos(:,1) > sess.tt(movei,6) & ...
        sqrt((BDF.pos(:,2)-3).^2 + (BDF.pos(:,3)+33).^2) > 6.9,1,'first');
    %bdfend = bdfmove(end);
    
    plot(BDF.pos(bdfstart:bdfend,2)-3,BDF.pos(bdfstart:bdfend,3)+33,'b');
    
    if ~isempty('Bayes_estimates')
        plot([7*cos(Bayes_estimates(movei)); 9*cos(Bayes_estimates(movei))],...
             [7*sin(Bayes_estimates(movei)); 9*sin(Bayes_estimates(movei))],'k','LineWidth',3); 
    end
    if i > 1
        plot([7*cos(sess.tt(movei-1,2)); 9*cos(sess.tt(movei-1,2))],...
             [7*sin(sess.tt(movei-1,2)); 9*sin(sess.tt(movei-1,2))],'Color',[.5 .5 .5],'LineWidth',2); 
    end
    axis off;
%     pause;
%     cla;
    
    subplot(1,2,2); hold on;
    
    plot(1:movei,depart_fromB(1:movei),'k.');
	plot(movei,depart_fromB(movei),'b.','MarkerSize',20);
    plot([0 1].*length(sess.tt),[0 0],'k--');
%     plot([-pi pi],[0 0],'k');
%     plot([0 0],[-pi pi],'k');
    
    pause; 
    subplot(1,2,1); cla;
    subplot(1,2,2); cla;
    
end