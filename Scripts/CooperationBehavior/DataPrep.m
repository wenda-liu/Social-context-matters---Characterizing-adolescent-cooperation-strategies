clear all;
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load data_se.mat
% make response consistant 
datafile_se((datafile_se(:,9) == 3),9) = 0;
datafile_se(:,11) = datafile_se(:,11) - 1;
datafile_se((datafile_se(:,11) == 2),11) = 0;

load data_cc.mat
load data_sc.mat
load data_ci.mat
load data_si.mat
load data_30.mat
load data_70.mat


data_all = [datafile_se;datafile_cc;datafile_sc;datafile_ci;datafile_si;datafile_30;datafile_70]; %conbine all conditions
data_all(isnan(data_all(:,9)),:) = [];
data_all((data_all(:,9 )== 0),:) = []; %exclude no response trial
% calculate data exclusion percent
cond_code = unique(data_all(:,2),"stable");
for cond_i = 1:length(cond_code)
    sub_id = unique(data_all(data_all(:,2) == cond_code(cond_i),1),"stable");
    for sub_i = 1:length(sub_id)
        cond_data = data_all(data_all(:,2) == cond_code(cond_i),:);
        data_size(sub_i,cond_i) = sum(cond_data(:,1) == sub_id(sub_i));
    end
    clearvars cond_data
end

[exclude_index(:,1),exclude_index(:,2)] = find(data_size<15 & data_size>0);
for exclude_i = 1:size(exclude_index,1)
    sub_id = unique(data_all(data_all(:,2) == cond_code(exclude_index(exclude_i,2)),1),"stable");
    data_all(data_all(:,2) == cond_code(exclude_index(exclude_i,2)) & data_all(:,1) == sub_id(exclude_index(exclude_i,2)),:) = [];
end
data_all(isnan(data_all(:,10)),10) = 0;
save data_cleaned_all.mat data_all cond_code
%%
%% after keep|after share
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load data_cleaned_all.mat
cond_code = unique(data_all(:,2),"stable");
a = data_all;
% a(isnan(a(:,11)),11) = 0;
% a(a(:,10) == 0,11) = nan;
% a(a(:,10) == 0,10) = nan;
lme(:,1) = a(:,1); %subj
lme(:,2) = a(:,9)-1; %pshare
lme(:,3) = a(:,10); %ashare
lme(:,4) = a(:,11); %ptrustee
for cond_i = 1:length(cond_code) 
    subid = unique(a(a(:,2) == cond_code(cond_i),1),"stable");
    for i = 1:length(subid)
        if cond_i == 1 && i ==1
            t4t(1,1) = 0.5;
            tempo(:,1) = a(a(:,2) == cond_code(cond_i) & a(:,1) == subid(i),11);
            t4t(end+1:end+length(tempo)-1,1) = tempo(1:end-1,1);
            clear tempo
        else
            t4t(end+1,1) = 0.5;
            tempo(:,1) = a(a(:,2) == cond_code(cond_i) & a(:,1) == subid(i),11);
            t4t(end+1:end+length(tempo)-1,1) = tempo(1:end-1,1);
            clear tempo
        end
    end
    clear subid
end

lme(:,5) = t4t;
gt4t = t4t;
for i = 0:length(t4t)-2
    j = length(t4t) - i;
    if gt4t(j) == 0 && t4t(j-1) ~= 0
        gt4t(j) = 1;
    end
end
lme(:,6) = gt4t;
lme(:,7) = a(:,2); %cond
%% for plot
for cond_i = 1:length(cond_code) 
    subid = unique(a(a(:,2) == cond_code(cond_i),1),"stable");
    for i = 1:length(subid)
        sub_data = lme(lme(:,1) == subid(i) & lme(:,7) == cond_code(cond_i),:);
        afterkeep = sub_data(sub_data(:,5) == 0,3); %after keep
        afternan = [sub_data(isnan(sub_data(:,5)),3)]; %after no choice
        aftershare = sub_data(sub_data(:,5) == 1,3); %after share
        a_k_s{cond_i}(i,3) = mean(afterkeep);
        a_k_s{cond_i}(i,2) = mean(aftershare);
        a_k_s{cond_i}(i,1) = mean(afternan);
        clear sub_data afterkeep aftershare

%         a_k_s{cond_i}(i,1) = mean(lme(lme(:,1) == subid(i) & lme(:,7) == cond_code(cond_i) & lme(:,5) == 0,2)); %2p 3ashare after keep
%         a_k_s{cond_i}(i,2) = mean(lme(lme(:,1) == subid(i) & lme(:,7) == cond_code(cond_i) & lme(:,5) == 1,2)); %2p a3share after share
%         a_k_s{cond_i}(i,3) = mean(lme(lme(:,1) == subid(i) & lme(:,7) == cond_code(cond_i) & lme(:,5) == .5,2)); %2p ashare after share
    end
    clear subid
end
%%
a_k_s_30 = a_k_s{6};
a_k_s_70 = a_k_s{7};
a_k_s_ns = [a_k_s{2};a_k_s{3}];
a_k_s_s = [a_k_s{4};a_k_s{5}];
a_k_s_se = a_k_s{1};

to_plot_violin = [a_k_s_30(:,1);a_k_s_30(:,2);a_k_s_70(:,1);a_k_s_70(:,2);a_k_s_ns(:,1);a_k_s_ns(:,2);...
    a_k_s_s(:,1);a_k_s_s(:,2);a_k_s_se(:,1);a_k_s_se(:,2);];
%% violin plot
% to_plot_violin = [keep(1:70);share(1:70);keep(71:132);share(71:132);keep(133:end);share(133:end)];
name{1} = 'fix30 keep';
name{2} = 'fix30 share';
name{3} = 'fix70 keep';
name{4} = 'fix70 share'; 
name{5} = 'nonsocial keep';
name{6} = 'nonsocial share';
name{7} = 'social keep';
name{8} = 'social share'; 
name{9} = 'real keep';
name{10} = 'real share';

% size_condition = [size(invest_cc,2)+size(invest_sc,2) size(invest_ci,2)+size(invest_si,2)];%1nonsocial 2social
size_condition = [sum(~isnan(a_k_s_30(:,1))) sum(~isnan(a_k_s_30(:,2))) sum(~isnan(a_k_s_70(:,1)))...
    sum(~isnan(a_k_s_70(:,2))) sum(~isnan(a_k_s_ns(:,1))) sum(~isnan(a_k_s_ns(:,2)))...
    sum(~isnan(a_k_s_s(:,1))) sum(~isnan(a_k_s_s(:,2)))...
    sum(~isnan(a_k_s_se(:,1))) sum(~isnan(a_k_s_se(:,2)))];%1nonsocial 2social
to_plot_violin(isnan(to_plot_violin)) = [];
cat_keepshare(1:size_condition(1),1) = repmat(name(1),size_condition(1),1);
cat_keepshare(end+1:end+size_condition(2),1) = repmat(name(2),size_condition(2),1);
cat_keepshare(end+1:end+size_condition(3),1) = repmat(name(3),size_condition(3),1);
cat_keepshare(end+1:end+size_condition(4),1) = repmat(name(4),size_condition(4),1);
cat_keepshare(end+1:end+size_condition(5),1) = repmat(name(5),size_condition(5),1);
cat_keepshare(end+1:end+size_condition(6),1) = repmat(name(6),size_condition(6),1);
cat_keepshare(end+1:end+size_condition(7),1) = repmat(name(7),size_condition(7),1);
cat_keepshare(end+1:end+size_condition(8),1) = repmat(name(8),size_condition(8),1);
cat_keepshare(end+1:end+size_condition(9),1) = repmat(name(9),size_condition(9),1);
cat_keepshare(end+1:end+size_condition(10),1) = repmat(name(10),size_condition(10),1);

% to_plot_violin = to_plot_violin +1;
figure
C = colororder;
C(4,:) = [152;245;255]./255;
C(5,:) = [122;197;205]./255;

n=1;
grouporder={'fix30 keep','fix30 share','fix70 keep', 'fix70 share',...
    'nonsocial keep','nonsocial share','social keep','social share','real keep','real share'};
% left violin
vsl1 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.05,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(4,:)});  

% Same for the right violin
vsr1 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.05,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(4,:)});  
% left violin
n=2;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.05,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(5,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.05,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(5,:)});  
% left violin
n=3;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.05,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(1,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.05,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(1,:)});  
% left violin
n=4;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.05,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(2,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.05,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(2,:)});  % left violin
n=5;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.05,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(3,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.05,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(3,:)});  

%  ylim([0.9 2.1]);
%   ylim([-0.1 1.1]);
%    yticks([0 .2 .4 .6 .8 1])

%  yticks([1 1.2 1.4 1.6 1.8 2])
%  yticklabels({'0','0.2','0.4','0.6','0.8','1'})
ylim([-0.2 5.2]);
yticks([0 1 2 3 4 5])
%  yticklabels({'0','0.2','0.4','0.6','0.8','1'})
ylabel('share probability')
xticks([1 2 3 4 5])
%  xticklabels({sprintf('fix30\\newlinekeep/share'),sprintf('fix70\\newlinekeep/share'),...
%      sprintf('nonsocial\\newlinekeep/share'),sprintf('social\\newlinekeep/share'),sprintf('real\\newlinekeep/share')})
 xticklabels({sprintf('Fix30'),sprintf('Fix70'),...
     sprintf('Nonsocial'),sprintf('Social'),sprintf('Real')})
 hold off
%%
%% investor|trustee pshare
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load data_cleaned_all.mat
cond_code = unique(data_all(:,2),"stable");
a = data_all;
a(a(:,10) == 0,11) = nan;
% a(a(:,10) == 0,10) = nan;

lme(:,1) = a(:,1); %subj
lme(:,2) = a(:,9)-1; %pshare
lme(:,3) = a(:,10); %ashare
lme(:,4) = a(:,11); %ptrustee
lme(:,7) = a(:,2); %cond

%% for plot
for cond_i = 1:length(cond_code) 
    subid = unique(a(a(:,2) == cond_code(cond_i),1),"stable");
    for i = 1:length(subid)
        p_i_t{cond_i}(i,1) = mean(lme(lme(:,1) == subid(i) & lme(:,7) == cond_code(cond_i),2));
        p_i_t{cond_i}(i,2) = nanmean(lme(lme(:,1) == subid(i) & lme(:,7) == cond_code(cond_i),4));
    end
    clear subid
end
%%
fix30 = p_i_t{6};
fix70 = p_i_t{7};
nonsocial = [p_i_t{2};p_i_t{3}];
social = [p_i_t{4};p_i_t{5}];
real = p_i_t{1};

to_plot_violin = [fix30(:,1); fix30(:,2); fix70(:,1); fix70(:,2); nonsocial(:,1);...
    nonsocial(:,2); social(:,1); social(:,2); real(:,1); real(:,2)];
%% violin plot
name{1} = 'fix30 keep';
name{2} = 'fix30 share';
name{3} = 'fix70 keep';
name{4} = 'fix70 share'; 
name{5} = 'nonsocial keep';
name{6} = 'nonsocial share';
name{7} = 'social keep';
name{8} = 'social share'; 
name{9} = 'real keep';
name{10} = 'real share';

% size_condition = [size(invest_cc,2)+size(invest_sc,2) size(invest_ci,2)+size(invest_si,2)];%1nonsocial 2social
size_condition = [sum(~isnan(fix30(:,1))) sum(~isnan(fix30(:,2))) sum(~isnan(fix70(:,1)))...
    sum(~isnan(fix70(:,2))) sum(~isnan(nonsocial(:,1))) sum(~isnan(nonsocial(:,2)))...
    sum(~isnan(social(:,1))) sum(~isnan(social(:,2)))...
    sum(~isnan(real(:,1))) sum(~isnan(real(:,2)))];%1nonsocial 2social
to_plot_violin(isnan(to_plot_violin)) = [];
cat_keepshare(1:size_condition(1),1) = repmat(name(1),size_condition(1),1);
cat_keepshare(end+1:end+size_condition(2),1) = repmat(name(2),size_condition(2),1);
cat_keepshare(end+1:end+size_condition(3),1) = repmat(name(3),size_condition(3),1);
cat_keepshare(end+1:end+size_condition(4),1) = repmat(name(4),size_condition(4),1);
cat_keepshare(end+1:end+size_condition(5),1) = repmat(name(5),size_condition(5),1);
cat_keepshare(end+1:end+size_condition(6),1) = repmat(name(6),size_condition(6),1);
cat_keepshare(end+1:end+size_condition(7),1) = repmat(name(7),size_condition(7),1);
cat_keepshare(end+1:end+size_condition(8),1) = repmat(name(8),size_condition(8),1);
cat_keepshare(end+1:end+size_condition(9),1) = repmat(name(9),size_condition(9),1);
cat_keepshare(end+1:end+size_condition(10),1) = repmat(name(10),size_condition(10),1);

% to_plot_violin = to_plot_violin +1;
figure(10)
C = colororder;
C(4,:) = [152;245;255]./255;
C(5,:) = [122;197;205]./255;
n=1;
grouporder={'fix30 keep','fix30 share','fix70 keep', 'fix70 share',...
    'nonsocial keep','nonsocial share','social keep','social share','real keep','real share'};
% left violin
vsl1 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.105,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(4,:)});  

% Same for the right violin
vsr1 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.105,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(4,:)});  
% left violin
n=2;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.105,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(5,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.105,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(5,:)});  
% left violin
n=3;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.105,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(1,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.105,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(1,:)});  
% left violin
n=4;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.105,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(2,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.105,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(2,:)});  % left violin
n=5;
vsl2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n-1}))},...         %Your data
    n-0.105,...                                                                          % 0.05 is the offset
    'HalfViolin','left',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(3,:)});                                               

% Same for the right violin
vsr2 = Violin({to_plot_violin(strcmp(cat_keepshare, grouporder{2*n}))},...
    n+0.105,... % Note here that I add the offset
    'HalfViolin','right',...                                                           % just left side violin plot
    'QuartileStyle','shadow',...
    'DataStyle', 'scatter',...
    'ShowNotches', false,...
    'ShowMean', true,...
    'ShowMedian', false,...
    'ViolinColor',{C(3,:)});  

%  ylim([0.9 2.1]);
  ylim([-0.1 1.1]);
   yticks([0 .2 .4 .6 .8 1])

%  yticks([1 1.2 1.4 1.6 1.8 2])
 ylabel('share probability')
 yticklabels({'0','0.2','0.4','0.6','0.8','1'})
 xticks([1 2 3 4 5])
 xticklabels({sprintf('Fix30'),sprintf('Fix70'),...
     sprintf('Nonsocial'),sprintf('Social'),sprintf('Real')})
hold off
%% connect dots
fix30(:,3) = .895;
fix30(:,4) = 1.105;
line([fix30(:,3) fix30(:,4)]',[fix30(:,1) fix30(:,2)]','Color',C(4,:))

fix70(:,3) = 1.895;
fix70(:,4) = 2.105;
line([fix70(:,3) fix70(:,4)]',[fix70(:,1) fix70(:,2)]','Color',C(5,:))

nonsocial(:,3) = 2.895;
nonsocial(:,4) = 3.105;
line([nonsocial(:,3) nonsocial(:,4)]',[nonsocial(:,1) nonsocial(:,2)]','Color',C(1,:))

social(:,3) = 3.895;
social(:,4) = 4.105;
line([social(:,3) social(:,4)]',[social(:,1) social(:,2)]','Color',C(2,:))

real(:,3) = 4.895;
real(:,4) = 5.105;
line([real(:,3) real(:,4)]',[real(:,1) real(:,2)]','Color',C(3,:))
%%

%% lme data glm
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load data_cleaned_all.mat
cond_code = unique(data_all(:,2),"stable");
a = data_all;
% a(isnan(a(:,11)),11) = 0;
a(a(:,10) == 0,11) = nan;
% a(a(:,10) == 0,10) = nan;
lme(:,1) = a(:,1); %subj
lme(:,2) = a(:,9)-1; %pshare
lme(:,3) = a(:,10); %ashare
lme(isnan(lme(:,3)),3)= 0;
lme(:,4) = a(:,11); %ptrustee
for cond_i = 1:length(cond_code) 
    subid = unique(a(a(:,2) == cond_code(cond_i),1),"stable");
    for i = 1:length(subid)
        if cond_i == 1 && i ==1
            t4t(1,1) = 0.5;
            tempo(:,1) = a(a(:,2) == cond_code(cond_i) & a(:,1) == subid(i),11);
            t4t(end+1:end+length(tempo)-1,1) = tempo(1:end-1,1);
            clear tempo
        else
            t4t(end+1,1) = 0.5;
            tempo(:,1) = a(a(:,2) == cond_code(cond_i) & a(:,1) == subid(i),11);
            t4t(end+1:end+length(tempo)-1,1) = tempo(1:end-1,1);
            clear tempo
        end
    end
    clear subid
end
lme(:,5) = t4t;
gt4t = t4t;
for i = 0:length(t4t)-2
    j = length(t4t) - i;
    if gt4t(j) == 0 && t4t(j-1) ~= 0
        gt4t(j) = 1;
    end
end
lme(:,6) = gt4t;
lme(lme(:,5)==0.5,7) = 0.5; %keep model
lme(lme(:,5)~=0.5,7) = 0;
lme(:,8) = lme(:,7); %share model
lme(lme(:,7)==0,8) = 1;
% nan = 0.5 t4t gt4t
% lme(isnan(lme(:,5)),5) = 0.5;
% lme(isnan(lme(:,6)),6) = 0.5;
lme(:,10) = a(:,2); %cond
lme(:,11) = a(:,3); %trial
% game
lme(lme(:,10)==11 | lme(:,10)==12,12) = 3; %cc sc
lme(lme(:,10)==30,12) = 1; %30
lme(lme(:,10)==70,12) = 2; %70
lme(lme(:,10)==13 | lme(:,10)==14,12) = 4; %ci si
lme(lme(:,10)==21,12) = 5; %real
% social
lme(lme(:,10)==30 | lme(:,10)==70,13) = NaN; %fix
lme(lme(:,10)==11 | lme(:,10)==12,13) = 0; %cc sc
lme(lme(:,10)==13 | lme(:,10)==14,13) = 1; %ci si
lme(lme(:,10)==21,13) = 1; %real
% AI
lme(lme(:,10)==30 | lme(:,10)==70,14) = 0; %fix
lme(lme(:,10)==11 | lme(:,10)==12,14) = 1; %cc sc
lme(lme(:,10)==13 | lme(:,10)==14,14) = 1; %ci si
lme(lme(:,10)==21,14) = NaN; %real
load demographic.mat
for j = 1:size(demographic,1)
    lme(lme(:,1) == demographic(j,1),15) = demographic(j,2);
    lme(lme(:,1) == demographic(j,1),16) = demographic(j,3);
    lme(lme(:,1) == demographic(j,1),17) = demographic(j,4);
    lme(lme(:,1) == demographic(j,1),18) = demographic(j,5);
end
header = {"ID","pshare","ashare","ptrustee","t4t","gt4t","keepm","sharem","nan","cond",...
    "trial","cond1-5","sociao/nonsocial","AI"};

%% mean lme
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load lme_g.mat
%
% lme(isnan(lme(:,3)),3)= 0;
lme(isnan(lme(:,5)),5)= 0.5;

%
cond_i = unique(lme(:,10),'stable');
n_count = 0;
for i = 1:length(cond_i)
    sub_j = unique(lme(lme(:,10)==cond_i(i),1),'stable');
    for j = 1:length(sub_j)
        n_count = n_count+1;
        m_lme(n_count,1) = sub_j(j);
        m_lme(n_count,2) = cond_i(i);
        m_lme(n_count,3) = nanmean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i),2)); %pshare
        m_lme(n_count,4) = nanmean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i),3)); %ashare
        m_lme(n_count,5) = nanmean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i),4)); %ptrustee
        m_lme(n_count,14) = mean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i) & lme(:,5) == 0,3)); %2p 3ashare after keep
        m_lme(n_count,15) = mean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i) & lme(:,5) == 1,3)); %2p a3share after share
        m_lme(n_count,16) = mean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i) & lme(:,5) == .5,3)); %2p ashare after share
    end
end
m_lme(:,17) = m_lme(:,15) - m_lme(:,14);

m_lme(:,6) = m_lme(:,3)./m_lme(:,5);
% game
m_lme(m_lme(:,2)==11 | m_lme(:,2)==12,11) = 3; %cc sc
m_lme(m_lme(:,2)==30,11) = 1; %30
m_lme(m_lme(:,2)==70,11) = 2; %70
m_lme(m_lme(:,2)==13 | m_lme(:,2)==14,11) = 4; %ci si
m_lme(m_lme(:,2)==21,11) = 5; %real
% social
m_lme(m_lme(:,2)==30 | m_lme(:,2)==70,12) = NaN; %fix
m_lme(m_lme(:,2)==11 | m_lme(:,2)==12,12) = 0; %cc sc
m_lme(m_lme(:,2)==13 | m_lme(:,2)==14,12) = 1; %ci si
m_lme(m_lme(:,2)==21,12) = 1; %real
% AI
m_lme(m_lme(:,2)==30 | m_lme(:,2)==70,13) = 0; %fix
m_lme(m_lme(:,2)==11 | m_lme(:,2)==12,13) = 1; %cc sc
m_lme(m_lme(:,2)==13 | m_lme(:,2)==14,13) = 1; %ci si
m_lme(m_lme(:,2)==21,13) = NaN; %real

load demographic.mat
for j = 1:size(demographic,1)
    m_lme(m_lme(:,1) == demographic(j,1),7) = demographic(j,2);
    m_lme(m_lme(:,1) == demographic(j,1),8) = demographic(j,3);
    m_lme(m_lme(:,1) == demographic(j,1),9) = demographic(j,4);
    m_lme(m_lme(:,1) == demographic(j,1),10) = demographic(j,5);
end

% within subject
m_lme(~ismember(m_lme(:,1),subj4cond),:) = [];
%% tiral  plot
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load lme_g.mat
AI = lme(lme(:,14)==1,:);
Fix = lme(lme(:,14)==0,:);
trial(1,:) = [1:20];
% to_plot_AI = NaN(200,20);
% to_plot_Fix = NaN(200,20);
for i = 1:20
    l1 = length(AI(AI(:,11)==i,2));
%     to_plot_AI(1:l1,i) = AI(AI(:,11)==i,2);
    to_plot_AI(:,i) = nanmean(AI(AI(:,11)==i,2));
    errorbar_AI(:,i) = nanstd(AI(AI(:,11)==i,2))/sqrt(l1);

    l2 = length(Fix(Fix(:,11)==i,2));
%     to_plot_Fix(1:l2,i) = Fix(Fix(:,11)==i,2);
    to_plot_Fix(:,i) = nanmean(Fix(Fix(:,11)==i,2));
    errorbar_Fix(:,i) = nanstd(Fix(Fix(:,11)==i,2))/sqrt(l2);
end
clf
shadedErrorBar(trial,to_plot_AI,errorbar_AI,'lineprops','-r','patchSaturation',0.33)
hold on
shadedErrorBar(trial,to_plot_Fix,errorbar_Fix,'lineprops','-b','patchSaturation',0.33)

  ylim([0 1]);
   yticks([0 .2 .4 .6 .8 1])
  xlim([0 20]);
   xticks([1 5 10 15 20])

%  yticks([1 1.2 1.4 1.6 1.8 2])
 ylabel('share probability')
 yticklabels({'0','0.2','0.4','0.6','0.8','1'})
 xlabel('Trial')

hold off
%% tiral  plot 4 cond
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load lme_g.mat
social = lme(lme(:,12)==4,:);
fix30 = lme(lme(:,12)==1,:);
nonsocial = lme(lme(:,12)==3,:);
fix70 = lme(lme(:,12)==2,:);
trial(1,:) = [1:20];
% to_plot_AI = NaN(200,20);
% to_plot_Fix = NaN(200,20);
for i = 1:20
    l1 = length(social(social(:,11)==i,2));
%     to_plot_AI(1:l1,i) = AI(AI(:,11)==i,2);
    to_plot_social(:,i) = nanmean(social(social(:,11)==i,2));
    errorbar_social(:,i) = nanstd(social(social(:,11)==i,2))/sqrt(l1);

    l2 = length(fix30(fix30(:,11)==i,2));
%     to_plot_Fix(1:l2,i) = Fix(Fix(:,11)==i,2);
    to_plot_30(:,i) = nanmean(fix30(fix30(:,11)==i,2));
    errorbar_30(:,i) = nanstd(fix30(fix30(:,11)==i,2))/sqrt(l2);

    l3 = length(fix70(fix70(:,11)==i,2));
%     to_plot_Fix(1:l2,i) = Fix(Fix(:,11)==i,2);
    to_plot_70(:,i) = nanmean(fix70(fix70(:,11)==i,2));
    errorbar_70(:,i) = nanstd(fix70(fix70(:,11)==i,2))/sqrt(l2);

    l4 = length(nonsocial(nonsocial(:,11)==i,2));
%     to_plot_Fix(1:l2,i) = Fix(Fix(:,11)==i,2);
    to_plot_non(:,i) = nanmean(nonsocial(nonsocial(:,11)==i,2));
    errorbar_non(:,i) = nanstd(nonsocial(nonsocial(:,11)==i,2))/sqrt(l2);
end
clf
shadedErrorBar(trial,to_plot_social,errorbar_social,'lineprops','-r','patchSaturation',0.33)
hold on
shadedErrorBar(trial,to_plot_30,errorbar_30,'lineprops','-b','patchSaturation',0.33)
hold on
shadedErrorBar(trial,to_plot_70,errorbar_70,'lineprops','-g','patchSaturation',0.33)
hold on
shadedErrorBar(trial,to_plot_non,errorbar_non,'lineprops','-y','patchSaturation',0.33)
  ylim([0 1]);
   yticks([0 .2 .4 .6 .8 1])
  xlim([0 20]);
   xticks([1 5 10 15 20])

%  yticks([1 1.2 1.4 1.6 1.8 2])
 ylabel('share probability')
 yticklabels({'0','0.2','0.4','0.6','0.8','1'})
 xlabel('Trial')

hold off
%% find within subject data
clear all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load lme_m.mat

% m_lme(m_lme(:,11)==5,:) = [];
subid = unique(m_lme(:,1),'stable');
for i_sub = 1:length(subid)
    subdata = m_lme(m_lme(:,1) == subid(i_sub),:);
    sub_cond(i_sub,1) = subid(i_sub);
    cond = unique(subdata(:,11),'sorted');
    for i = 1:length(cond)
        switch cond(i)
            case 1
            sub_cond(i_sub,2) = 1;
            case 2
            sub_cond(i_sub,3) = 1;
            case 3
            sub_cond(i_sub,4) = 1;
            case 4
            sub_cond(i_sub,5) = 1;   
            case 5
            sub_cond(i_sub,6) = 1;              
        end
    end
    clear subdata
end
%
% subj4cond = sub_cond((sum(sub_cond(:,2:5),2) == 4),1);
% sum(sub_cond(:,4)==1 & sub_cond(:,5)==1)
%
subj4cond = sub_cond((sum(sub_cond(:,2:3),2) > 0 & sum(sub_cond(:,4:5),2) > 0),1);
sum(sub_cond(:,4)==1 & sub_cond(:,5)==1)
%% yale gw
clear all
close all
load lme_m.mat

suball = unique(m_lme(:,1));
GW = suball(round(suball(:)) == (suball(:)));
Yale = suball(round(suball(:)) ~= (suball(:)));

for i = 1:length(GW)
    GWdemo(i,:) = mean(m_lme(m_lme(:,1)==GW(i),7:10),1);
end
GWdemo(:,1) = GWdemo(:,1)/12;
mean(GWdemo)
std(GWdemo)
sum(GWdemo(:,2)==1)
for i = 1:length(Yale)
    Yaledemo(i,:) = mean(m_lme(m_lme(:,1)==Yale(i),7:10),1);
end
Yaledemo(:,1) = Yaledemo(:,1)/12;
nanmean(Yaledemo)
nanstd(Yaledemo)
sum(Yaledemo(:,2)==1)
sum(isnan(Yaledemo(:,4)))
%% adaptivity/non
clear all
close all
load lme_m.mat

GW = unique(m_lme(m_lme(:,11)==1 | m_lme(:,11)==2,1));
Yale = unique(m_lme(m_lme(:,11)==3 | m_lme(:,11)==4,1));

for i = 1:length(GW)
    GWdemo(i,:) = mean(m_lme(m_lme(:,1)==GW(i),7:10),1);
end
GWdemo(:,1) = GWdemo(:,1)/12;
nanmean(GWdemo)
nanstd(GWdemo)
sum(GWdemo(:,2)==1)
sum(isnan(GWdemo(:,3)))
for i = 1:length(Yale)
    Yaledemo(i,:) = mean(m_lme(m_lme(:,1)==Yale(i),7:10),1);
end
Yaledemo(:,1) = Yaledemo(:,1)/12;
nanmean(Yaledemo)
nanstd(Yaledemo)
sum(Yaledemo(:,2)==1)
sum(isnan(Yaledemo(:,3)))
%% social/non
clear all
close all
load lme_m.mat

GW = unique(m_lme(m_lme(:,11)==3,1));
Yale = unique(m_lme(m_lme(:,11)==4,1));

for i = 1:length(GW)
    GWdemo(i,:) = mean(m_lme(m_lme(:,1)==GW(i),7:10),1);
end
GWdemo(:,1) = GWdemo(:,1)/12;
nanmean(GWdemo)
nanstd(GWdemo)
sum(GWdemo(:,2)==1)
sum(isnan(GWdemo(:,4)))
for i = 1:length(Yale)
    Yaledemo(i,:) = mean(m_lme(m_lme(:,1)==Yale(i),7:10),1);
end
Yaledemo(:,1) = Yaledemo(:,1)/12;
nanmean(Yaledemo)
nanstd(Yaledemo)
sum(Yaledemo(:,2)==1)
sum(isnan(Yaledemo(:,4)))
%%
clear all
close all
load lme_m.mat
GWdata = m_lme((round(m_lme(:,1)) ~= (m_lme(:,1))),:);

GW = unique(GWdata(GWdata(:,11)==2,1));
Yale = unique(GWdata(GWdata(:,11)==4,1));
sum(ismember(GW,Yale))

for i = 1:length(GW)
    GWdemo(i,:) = mean(m_lme(m_lme(:,1)==GW(i),7:10),1);
end
GWdemo(:,1) = GWdemo(:,1)/12;
nanmean(GWdemo)
nanstd(GWdemo)
sum(GWdemo(:,2)==1)
sum(isnan(GWdemo(:,4)))
for i = 1:length(Yale)
    Yaledemo(i,:) = mean(m_lme(m_lme(:,1)==Yale(i),7:10),1);
end
Yaledemo(:,1) = Yaledemo(:,1)/12;
nanmean(Yaledemo)
nanstd(Yaledemo)
sum(Yaledemo(:,2)==1)
sum(isnan(Yaledemo(:,4)))
%%
%%
clear all
load matlab.mat
%% gw vs yale
% GW = demographic(demographic_cond(:,9)==1,1);
% Yale = demographic(demographic_cond(:,10)==1,1);

% fix vs adap
GW = demographic(demographic_cond(:,2)==1|demographic_cond(:,3)==1,1);
Yale = demographic(demographic_cond(:,4)==1|demographic_cond(:,5)==1,1);
%% social vs non-social
GW = demographic(demographic_cond(:,4)==1,1);
Yale = demographic(demographic_cond(:,5)==1,1);
for i = 1:size(GW)
    GW(i,2:5) = demographic(demographic(:,1)==GW(i,1),2:5);
end

for i = 1:size(Yale)
    Yale(i,2:5) = demographic(demographic(:,1)==Yale(i,1),2:5);
end

% test1 = ismember(GW(:,1),Yale(:,1));
% test2 = ismember(Yale(:,1),GW(:,1));
% GW(test1,:)=[];
% Yale(test2,:)=[];

[H,P,CI,STATS] = ttest2(GW(:,2),Yale(:,2))
[H,P,CI,STATS] = ttest2(GW(:,4),Yale(:,4))
[H,P,CI,STATS] = ttest2(GW(:,5),Yale(:,5))

       n1 = sum(GW(:,3)-1); N1 = size(GW(:,3),1);
       n2 = sum(Yale(:,3)-1); N2 = size(Yale(:,3),1);
       x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       [tbl,chi2stat,pval] = crosstab(x1,x2)