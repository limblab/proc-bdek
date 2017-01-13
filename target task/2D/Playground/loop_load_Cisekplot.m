files2run = {'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10262016_alldays.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10282016_alldays.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10312016_alldays.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Mihili_07082014_alldays_transformed.mat'};
         
FD = cell(length(files2run),1);
 
for fi = 1:length(files2run)
    
    
    load(files2run{fi});

    for i = 1:length(alldays);  
        alldays(i).tt(isnan(alldays(i).tt(:,20)),20) = alldays(i).tt(isnan(alldays(i).tt(:,20)),9);
    end
    Cisek_plot_dotuning; 
    save([files2run{fi}(1:end-11) 'PDS'],'best_PDS');
    clear best_PDS;
    
    FD{fi} = Fire_dir{1};
end