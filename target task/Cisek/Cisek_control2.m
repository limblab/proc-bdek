task = 'two';

tt_2T = alldays(2).tt;
tt_1T = alldays(1).tt;

targ_dirs = tt_2T(:,13);
targ_dirs1 = tt_1T(:,13);

columnalign = 6;
%    figure; hold on;
[perc_overPDOD,perc_overPD,perc_overOD,...
    perc_overPDOD1,perc_overPD1,perc_overOD1,po_PD,po_OD,po_PDOD,numzers] = deal(zeros(length(alldays(1).PMd_units),1));
[ysPD,ysOD,ysPDOD,ysO] = deal(cell(length(alldays(1).PMd_units),1));
for n = 1:length(alldays(1).PMd_units)
    clc; fprintf('Unit: %d/%d\n',n,length(alldays(1).PMd_units));
    unit = alldays(1).PMd_units{n};
    pd = PD1{1}(n);

    if strcmp(task,'both') || strcmp(task,'two')
        [raster,allinds] = raster_plotCK(unit(2:end),tt_2T,[0 1],columnalign,0,'none',[],20);
        %rast_bin = bin_array(raster{1},size(raster{1},1),size(raster{1},2)/10,'sum');
        rast_bin = bin_array(raster{1},size(raster{1},1),2,'sum');
        wts = ones(1,25);
        zertrials = sum(rast_bin,2)==0;
        numzers(n) = sum(zertrials);
        
        T_frompd = circ_dist(targ_dirs,pd); 
        O_frompd = circ_dist(targ_dirs+pi,pd);
        %% 2 target 
        %-% PD , OD or ORTH
        pdonly = find(abs(T_frompd)<(pi/4) );
        odonly = find(abs(O_frompd)<(pi/4) );
        orth = find(abs(T_frompd)>(pi/4) & abs(O_frompd)>(pi/4) );
        
        %-% Loop through PD and OD trials
        pd_trials = rast_bin(pdonly,:);
%         yS_pdonly = zeros(size(pd_trials,1),76);
%         for i = 1:size(pd_trials,1)
%             yS_pdonly(i,:) = conv(pd_trials(i,:),wts,'valid');
%         end
        yS_pdonly = pd_trials;

        od_trials = rast_bin(odonly,:);
%         yS_odonly = zeros(size(od_trials,1),76);
%         for i = 1:size(od_trials,1)
%             yS_odonly(i,:) = conv(od_trials(i,:),wts,'valid');
%         end
        yS_odonly = od_trials;
        
        yS_pdod = [yS_pdonly; yS_odonly];

        %-% Orthogonal

        %-% Loop through orthogonal trials
        o_trials = rast_bin(orth,:);
%         yS_o = zeros(size(o_trials,1),76);
%         for i = 1:size(o_trials,1)
%             yS_o(i,:) = conv(o_trials(i,:),wts,'valid');
%         end
        
        yS_o = o_trials;
        
        %%
        
        ysPD{n} = yS_pdonly;
        ysOD{n} = yS_odonly;
        ysPDOD{n} = yS_pdod; 
        ysO{n} = yS_o;

        Counts_O = sqrt(reshape(yS_o,[],1));
        Counts_PDOD = sqrt(reshape(yS_pdod,[],1));
        Counts_PD = sqrt(reshape(yS_pdonly,[],1));
        Counts_OD = sqrt(reshape(yS_odonly,[],1));

        median_O = median(Counts_O);
        perc_overPDOD(n) = sum(Counts_PDOD > median_O)/length(Counts_PDOD);
        perc_overPD(n) = sum(Counts_PD > median_O)/length(Counts_PD);
        perc_overOD(n) = sum(Counts_OD > median_O)/length(Counts_OD);
        
        
%         po_PDnum = sum(sum((yS_pdonly) > repmat(mean((yS_o)),size((yS_pdonly),1),1)))
%         po_ODnum = sum(sum((yS_odonly) > repmat(mean((yS_o)),size((yS_odonly),1),1)))
%         po_PDODnum = sum(sum((yS_pdod) > repmat(mean((yS_o)),size((yS_pdod),1),1)))
   
        po_PD(n) = sum(sum((yS_pdonly) > repmat(mean((yS_o)),size((yS_pdonly),1),1)))/numel(yS_pdonly);
        po_OD(n) = sum(sum((yS_odonly) > repmat(mean((yS_o)),size((yS_odonly),1),1)))/numel(yS_odonly);
        po_PDOD(n) = sum(sum((yS_pdod) > repmat(mean((yS_o)),size((yS_pdod),1),1)))/numel(yS_pdod);
        
        
%         plot(mean(yS_pdonly)); plot(mean(yS_odonly),'r plot(mean(yS_o),'g');
%         title(sprintf('PD: %.0f  OD: %.0f',100*perc_overPD(n),100*perc_overOD(n)));
%         pause; cla;
    end
    
    %%
    if strcmp(task,'both') || strcmp(task,'one')
        [raster1,allinds1] = raster_plotCK(unit(2:end),tt_1T,[0 1],columnalign,0,'none',[],20);
        rast_bin1 = bin_array(raster1{1},size(raster1{1},1),size(raster1{1},2)/10,'sum');

        T_frompd1 = circ_dist(targ_dirs1,pd); 
        O_frompd1 = circ_dist(targ_dirs1+pi,pd);
        % 1 target 
        %-% PD or OD
        pdonly1 = find(abs(T_frompd1)<(pi/4));
        odonly1 = find(abs(O_frompd1)<(pi/4));
        %-% Loop through PD and OD trials
        pd_trials1 = rast_bin1(pdonly1,:);
        yS_pdonly1 = zeros(size(pd_trials1,1),51);
        for i = 1:size(pd_trials1,1)
            yS_pdonly1(i,:) = conv(pd_trials1(i,:),wts,'valid');
        end

        od_trials1 = rast_bin(odonly1,:);
        yS_odonly1 = zeros(size(od_trials1,1),51);
        for i = 1:size(od_trials1,1)
            yS_odonly1(i,:) = conv(od_trials1(i,:),wts,'valid');
        end

        yS_pdod1 = [yS_pdonly1; yS_odonly1];

        %-% Orthogonal
        orth1 = find(abs(T_frompd1)>(pi/4) & abs(O_frompd1)>(pi/4));
        %-% Loop through orthogonal trials
        o_trials1 = rast_bin(orth1,:);
        yS_o1= zeros(size(o_trials1,1),51);
        for i = 1:size(o_trials1,1)
            yS_o1(i,:) = conv(o_trials1(i,:),wts,'valid');
        end

        %%
        Counts_O1 = sqrt(reshape(yS_o1,[],1));
        Counts_PDOD1 = sqrt(reshape(yS_pdod1,[],1));
        Counts_PD1 = sqrt(reshape(yS_pdonly1,[],1));
        Counts_OD1 = sqrt(reshape(yS_odonly1,[],1));

        median_O1 = median(Counts_O1);
        perc_overPDOD1(n) = sum(Counts_PDOD1 >= median_O1)/length(Counts_PDOD1);
        perc_overPD1(n) = sum(Counts_PD1 >= median_O1)/length(Counts_PD1);
        perc_overOD1(n) = sum(Counts_OD1 >= median_O1)/length(Counts_OD1);


    end
end
