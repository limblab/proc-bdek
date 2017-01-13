pks = [25];
lks = [1 5 50];
numslices = 50;

sli = cell(length(lks),1);
samps = zeros(1,numslices);
realk = zeros(length(lks));
for j = 1:length(lks)   
    for s = 1:1000        
        keep = false;
        while ~keep       
            for nums = 1:numslices
                samps(nums) = vonmisrand(lks(j));
            end         
            if abs(circ_kappa(samps) - lks(j)) < 0.05*lks(j)
                keep = true;
            end
        end
        sli{j}(s,:) = samps;
        realk(j) = circ_kappa(circ_mean(sli{j},[],2));
        clc; fprintf('generating slices %d/%d (%d/%d)\n',j,length(lks),s,1000);
    end

end

%%
slopes = zeros(length(pks),length(lks));
slofunc = @(p,l,n) 1/(p/(l*n)+1);
for i = 1:length(pks) 
    for j = 1:length(ls); 
        slopes(i,j) = slofunc(pks(i),ls(j),numslices); 
    end
end

conf_kap = zeros(size(sli{j},1),length(lks));
for j = 1:length(lks)
    for s = 1:size(sli{j},1)
        mean_ests = bootstrp(1000,@circ_mean,sli{j}(s,:)');
        
        conf_kap(s,j) = circ_kappa(mean_ests);
        clc; fprintf('%d/%d (%d/%d)\n',j,length(lks),s,size(sli{j},1));
    end
end

        
        
        
        