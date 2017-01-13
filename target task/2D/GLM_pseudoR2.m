[d_null,d_part,d_full,p_val,pass] =deal(nan(size(VA{1},2),length(VA)));
[b_part,b_full] = deal(cell(size(VA{1},2),length(VA)));
for i =1:length(VA)
    
    for j = 1:size(VA{i},2);
        
        Y = VA{i}(:,j)./10;
        
        X_part = [sin(alldays(2).tt(:,10)) cos(alldays(2).tt(:,10))];
        X_full = [X_part alldays(2).tt(:,3)];
        X_null = ones(size(alldays(2).tt(:,2)));
        
        if sum(~isnan(Y))>50
        [b_part{j,i},d_part(j,i)] = glmfit(X_part,Y); % ALL
        [b_full{j,i},d_full(j,i)] = glmfit(X_full,Y);
        [~,d_null(j,i)] = glmfit(X_null,Y,'normal','constant','off');
        
%         pseudoR2.part(j,i) = 100*(1 - d_part(j,i)./d_null(j,i));
%         pseudoR2.full(j,i) = 100*(1 - d_full(j,i)./d_null(j,i));
%         
        p_val(j,i) = 1-chi2cdf(d_null(j,i)-d_full(j,i),2);
%         if p_val(j,i) < 0.05
%             pass(j,i) = 1;
%         end
%         
%         else
%             
%         b_part{j,i} = [nan; nan; nan];
%         b_full{j,i} = [nan; nan; nan; nan];
%             
%         pseudoR2.part(j,i) = NaN;
%         pseudoR2.full(j,i) = NaN;
%         
%         p_val(j,i) = NaN;
%         pass(j,i) = NaN;
        
        end
        pseudoR2.part = 100*(1 - d_part./d_null);
        pseudoR2.full = 100*(1 - d_full./d_null);

    end
    clc; fprintf('%d/%d\n',i,length(VA));

end
%

rR2 = 100*(pseudoR2.full-pseudoR2.part)./pseudoR2.full;