function[Recruit] = PPC_sanger(alldays,PDI,TA,good_neurs,tune_ranges,spike_count_raw,like_ind) 

sanger_func_notnorm = @(spikes,prior,tuning,deltaT) ...
    prior.*exp(nansum(bsxfun(@times,spikes,log(tuning)) - (tuning)*deltaT));

% sanger_func_notnorm = @(spikes,prior,tuning,deltaT) ...
%     prior.*prod(exp(-tuning*deltaT),1).*prod(tuning.^repmat(spikes,1,size(tuning,2)),1);

dts = diff(tune_ranges)./1000;
for i = 1:(length(TA)-1)
    clc; fprintf('time bin: %d/%d\n',i,length(TA));
    
    %% Get UN recruitments
    [AF] = deal(zeros(size(alldays(PDI).tt,1),629));
    for trial = 1:size(alldays(PDI).tt,1)

        spikes = spike_count_raw{1}{i}(trial,:)';
        prior = ones(1,629);
        tuning = TA{i}(good_neurs{i},:);
        deltaT = dts(i);

        PDF = sanger_func_notnorm(spikes,prior,tuning,deltaT);
 
        N = 0.01*sum(PDF(1:end-1)) + PDF(end);
        
        PDFn = PDF./N;
        
        PDFmove = PDFn(1:(end-1));
        PDFstate = PDFn(end);

        moveloc = alldays(PDI).tt(trial,10);
        move_ind = round(100*moveloc);


        if move_ind == 0 
            atzero = PDFmove; 
        else
            atzero = [PDFmove(move_ind:end) PDFmove(1:(move_ind-1))];
        end
        
        AF(trial,:) = [atzero(end/2:end) atzero(1:(end/2)-1) PDFstate];
    end

    % Clean up

    [bL,bU,av_rec,REC] = deal(cell(length(like_ind{1}),1));
    for lik = 1:length(like_ind{1})

        [bL{lik}, bU{lik}] = boot_bounds(1000,@nanmean,AF(like_ind{1}{lik},:),2.5,97.5);

        REC{lik} = AF(like_ind{1}{lik},:);
        av_rec{lik} = nanmean(REC{lik})';
    end
        
    Recruit.low(:,i) = bL;
    Recruit.high(:,i) = bU; 
    Recruit.av(:,i) = av_rec;
    Recruit.all{i} = REC;
end

 