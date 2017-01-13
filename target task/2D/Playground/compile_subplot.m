function[h] = compile_subplot(SUBPLOT_SIZE,NAMES,cleanup,varargin)

fp = cell(SUBPLOT_SIZE(1),SUBPLOT_SIZE(2));
[ax,s] = deal(cell(length(NAMES),1));
for i = 1:length(NAMES)
    
    fp{ceil(i/SUBPLOT_SIZE(2)),mod(i-1,SUBPLOT_SIZE(2))+1} = get(NAMES(i),'Position');

end

totwidth = sum(mean(cellfun(@(x) x(3),fp)))*1.15;
totheight = sum(mean(cellfun(@(x) x(4),fp)))*1.25;

h = figure('Position',[-827-0.5*totwidth,415-0.5*totheight,totwidth,totheight]);

for i = 1:length(NAMES)

    ax_i = gca(NAMES(i));

    s{i} = subplot(SUBPLOT_SIZE(1),SUBPLOT_SIZE(2),i); 
    hold on; 
    
    ax{i} = get(ax_i,'children');
    copyobj(ax{i},s{i});
    set(s{i},'XLim',ax_i.XLim,'YLim',ax_i.YLim);
    if nargin > 2 && cleanup == 1
        close(NAMES(i));
    end
end
