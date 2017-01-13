[LOW,HIGH,LOWc, HIGHc] = deal(cell(3,1));
for t = 1:length(VA{1}{1})
    
    for ses = 1:length(VA)

        linds = LI{ses}{1};
        hinds = LI{ses}{2};

        for PDODORTH = 1:3

            LOW{PDODORTH}{ses,t} = reshape(VA{ses}{PDODORTH}{t}(linds,:),[],1);
            HIGH{PDODORTH}{ses,t} = reshape(VA{ses}{PDODORTH}{t}(hinds,:),[],1);
            
%             CHANGE{PDODORTH}{ses,t} = nanmean(HIGH{PDODORTH}{ses,t}
        end

    end
end

CHANGE = zeros(3,size(LOW{1},2));
for PDODORTH = 1:3
    
    for t = 1:size(LOW{1},2)
        
        cfbL = vertcat(LOW{PDODORTH}{:,t}); 
        cfbH = vertcat(HIGH{PDODORTH}{:,t});
        
        CHANGE(PDODORTH,t) = nanmean(cfbH)-nanmean(cfbL);
    end
end
