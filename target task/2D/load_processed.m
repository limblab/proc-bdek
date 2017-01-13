function[alldaysout, priorsout, monkeyout, MO, DA, YE] = load_processed(monkey,mmddyy,with_kin,varargin)

if nargin < 2
    with_kin = 0;
end

orig_place = cd;

thiscd = which('load_processed'); 
thiscd(end-15:end) = [];

if strcmp(orig_place(1),'/')
    cd('/Users/brian/Desktop/Northwestern/Data/Mihili/Processed');
else
   
    if with_kin
        cd('C:\Users\limblab\Desktop\Mihili_bdf\full alldays');
%         cd('C:\Users\limblab\Desktop\Original Data\full alldays');
    else
        cd('C:\Users\limblab\Desktop\Mihili_bdf\Processed');
%         cd('C:\Users\limblab\Desktop\Original Data\Processed');
    end
end

if with_kin
    filename = ['full_alldays_' monkey '_' mmddyy '.mat'];
else
    filename = [monkey '_' mmddyy '.mat'];
end
load(filename);

if exist('alldays_nokin','var')
    alldaysout = alldays_nokin; 
elseif exist('alldays','var')
    alldaysout = alldays;    
else alldaysout = []; fprintf('NO alldays found');
end

if exist('priors','var'); priorsout = priors; else priorsout = []; fprintf('NO priors found\n'); end
monkeyout = monkey;

MO = mmddyy(1:2);
DA = mmddyy(3:4);
YE = mmddyy(5:8);

cd(orig_place);