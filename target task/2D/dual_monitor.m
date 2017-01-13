function dual_monitor(type)
% dual_monitor(type) switches default figure location to second monitor. 
%   type:
%        'small' - small window on second monitor
%        'large' - large window
%        'full'  - figure occupies entire second monitor
%        'off'   - returns settings to default
switch type
    
    case 'small'
        set(0,'defaultfigureposition',[-700 135 660 520]);    
    case 'large'
        set(0,'defaultfigureposition',[-1021 -49 1010 684]);  
    case 'full'
        set(0,'defaultfigureposition',[-1599 1 1595 818]);   
    case 'micro'
        set(0,'defaultfigureposition',[-528 253 518 205]);
    case 'off'
        set(0,'defaultfigureposition','default');
    case 'on'
        set(0,'defaultfigureposition',[-1018 135 1010 684]);     
    case 'wide'
        set(0,'defaultfigureposition',[-1184 135 1118 350]);
    case 'super_wide'
        set(0,'defaultfigureposition',[-1601 354 1595 382]);
end
