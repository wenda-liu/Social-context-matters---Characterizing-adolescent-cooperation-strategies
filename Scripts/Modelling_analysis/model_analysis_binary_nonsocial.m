%% read data
clear all
close all
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame
load lme_g.mat

%2 for psahre 3 for ashare
model_ori = lme(:,[1 2 4 10 12]); 

%2 for psahre 3 for ashare
% model_ori = lme(:,[1 3 4 10 12 4]); 
% model_ori(isnan(model_ori(:,2)),2) = 0;
% % a outcome
% model_ori(~isnan(model_ori(:,3)),3) = model_ori(~isnan(model_ori(:,3)),3).*model_ori(~isnan(model_ori(:,3)),2).*2 + 5-model_ori(~isnan(model_ori(:,3)),2);
% model_ori(~isnan(model_ori(:,6)),6) = model_ori(~isnan(model_ori(:,6)),6).*model_ori(~isnan(model_ori(:,6)),2).*2;
% model_ori(isnan(model_ori(:,3)),3) = 5;
% model_ori(~isnan(model_ori(:,3)),3) = model_ori(~isnan(model_ori(:,3)),3).*model_ori(~isnan(model_ori(:,3)),2).*2;

%% 90% 10%
model_all = model_ori;
cond_i = unique(lme(:,10),'stable');
n_count = 0;
for i = 1:length(cond_i)
    sub_j = unique(lme(lme(:,10)==cond_i(i),1),'stable');
    for j = 1:length(sub_j)
        n_count = n_count+1;
        pshare(n_count,1) = nanmean(lme(lme(:,1) == sub_j(j) & lme(:,10) == cond_i(i),2)); %pshare
        pshare(n_count,2) = cond_i(i);
        pshare(n_count,3) = sub_j(j);
    end
end

pshare((pshare(:,1)>0.9 | pshare(:,1)<0.1),1) = nan;
to_exclude = pshare(isnan(pshare(:,1)),:);

for i = 1:size(to_exclude,1)
    model_all(model_all(:,1) == to_exclude(i,3) & model_all(:,4) == to_exclude(i,2),:) = [];
end
clearvars -except model_ori model_all to_exclude
%%
sub_id = unique(model_all(:,1),"stable");
for i = 1:length(sub_id)
    model_sub{i} = model_all(model_all(:,1)==sub_id(i),:);
%     model_sub{i}(:,9) = model_all(model_all(:,1)==sub_id(i),6);
    model_sub{i}(1,6) = 1;
    for j = 2:size(model_sub{i},1)
        if model_sub{i}(j,5) == model_sub{i}(j-1,5)
           model_sub{i}(j,6) = 0;
        else
           model_sub{i}(j,6) = 1;
        end
    end
end

for i = 1:length(sub_id)
    model_sub{i}(1,7) = 1;
    for j = 2:size(model_sub{i},1)
        if model_sub{i}(j,4) == model_sub{i}(j-1,4)
           model_sub{i}(j,7) = 0;
        else
           model_sub{i}(j,7) = 1;
        end
    end
end
header_model = {"ID","pshare","ptrustee","cond","cond1-5","reset1","reset2"};
%% get social only
for i = 1:length(sub_id)
    model_sub_social{i,1} = model_sub{i}(model_sub{i}(:,5)==4,:);
%     model_sub_social{i,1} = model_sub{i}(model_sub{i}(:,5)==4 | model_sub{i}(:,5)==3,:);

end
data_con = model_sub_social(~cellfun('isempty',model_sub_social));
for i = 1:size(data_con,1)
    index = isnan(data_con{i}(:,3));
    data_con{i,1}(index,3) = 0;
    clear index;
end
%% find start point from real
load sub_pair.mat
for i = 1:size(data_con,1)
    id = data_con{i}(1,1);
    sub_data = model_ori(model_ori(:,1)==id,2:5);%2:5 6trustee share
%     social_id(i,1) = id;
    %ci si order
    switch social_id(i,2)
        case 13
        data_con{i} = sortrows(data_con{i},4,'ascend');
        case 14
        data_con{i} = sortrows(data_con{i},4,'descend');
    end
    %find partner
    switch isnan(nanmean(sub_data(sub_data(:,4)==5,1)))
%        data_con{i}(1,8) = nanmean(sub_data(sub_data(:,4)==5,1));
        case 1
       id2 = to_find_subid(i,2);
       sub_data2 = model_ori(model_ori(:,1)==id2,2:5);%2:5 6trustee share
           if isnan(nanmean(sub_data2(sub_data2(:,4)==5,2)))
           data_con{i}(1,8) = .5;%0.5 or 2.5
           nodata(i,1) = id;
           nodata(i,2) = id2;
           else 
           data_con{i}(1,8) = nanmean(sub_data2(sub_data2(:,4)==5,2))/2;%,2))outcome  ,5))trustee share
           end
        case 0
    data_con{i}(1,8) = nanmean(sub_data(sub_data(:,4)==5,2))/2;%,2))outcome  ,5))trustee share
%     clear sub_data sub_data2
    end    
    data_con{i}(:,6) = 0;
    data_con{i}(1,6) = 1;
end
%% nonsocial
for i = 1:length(sub_id)
    model_sub_non{i,1} = model_sub{i}(model_sub{i}(:,5)==3,:);
end
data_con_non = model_sub_non(~cellfun('isempty',model_sub_non));
for i = 1:size(data_con_non,1)
    index = isnan(data_con_non{i}(:,3));
    data_con_non{i,1}(index,3) = 0;
    clear index;
    data_con_non{i}(1,8) = .5;
    nonsocial_id(i,1) = data_con_non{i}(1,1);
end
%% everyone
% for i = 1:size(data_con_non,1)
%     [LIA,LOCB] = ismember(nonsocial_id(i,1),social_id);
%     if LIA == 1
%        data_con{LOCB} = [data_con{LOCB};data_con_non{i}];
%     else
%         data_con{size(data_con,1)+1,1} = data_con_non{i}; 
%     end
% end
%% have both
% for i = 1:size(data_con_non,1)
%     [LIA,LOCB] = ismember(nonsocial_id(i,1),social_id);
%     if LIA == 1
%        data_con2{LOCB,1} = [data_con{LOCB};data_con_non{i}];
%     end
% end
% clear data_con
% data_con = data_con2(~cellfun('isempty',data_con2));

%% model
cd /Users/wl/Documents/DSN_Lab/Inperson_trustgame/modeling_scripts/selectedmodels/binary/
n1 = 8; % for plots
n2 = 8;
which_kids = [1:size(data_con_non,1)];

for i_mod = [1]
% for i_mod = [1 2 3 4 5 6 7 8 9] 
 figure
for i_sub = 1:length(which_kids)
%  for i_sub=1
   % for i_sub = which_kids
        
        A.investor = data_con_non{ which_kids(i_sub), 1 }( :, 2); % my decision
        A.trustee = data_con_non{ which_kids(i_sub), 1 }( :, 3); % trustee 3outcome/9trustee share
        A.reset_sub = data_con_non{ which_kids(i_sub), 1 }( :, 6); % Reset social/non etc 6;cc/sc etc 7
        A.condition = data_con_non{ which_kids(i_sub), 1 }( :, 5); % condition
        A.prior = data_con_non{ which_kids(i_sub), 1 }( 1, 8 ); % prior
%         A.prior = 0.5; % reset

         if i_mod == 1
            % v23
            %% simple RL model with binarized ratings and answers
            % adaptation of i_mod == 1
            n = size(A.investor,1);
            % note: reset is always at 0.5, which is the natural midpoint
            [estimates_m{ i_mod }(i_sub,:), model, fval, exitflag(i_sub,i_mod), A.v_out, A.delta_out] = fitmodel_mle_rl_BIN_bound_onepar( A.trustee, A.reset_sub, A.investor, A.condition, A.prior  ); % xdata1 is outcome
            % Note: minimization not for SSE (sum of squared errors) but
            % for "LLE" negative log-likelhood (basically difference bt
            % "normal" regression & logistic regression
            % This changes formulas for AIC/BIC
            [LLE(i_sub,i_mod), FittedCurve2] = model(estimates_m{ i_mod }(i_sub,:)); 
            AIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 3*1; % Akaike Information Criterion  % k1 = 1; % number of free parameters
            BIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 3*log(n)/2; % Bayesian Information Criterion
            
          
%             to_plot = fitmodel_mle_rl_BIN_bound_to_plot(A.trustee, A.reset_sub, A.investor, A.condition, A.prior, estimates_m{ i_mod }(i_sub,:));
            to_plot = A.v_out; 
              elseif i_mod == 2
            % v23
            %% simple RL model with binarized ratings and answers
            % adaptation of i_mod == 1
            n = size(A.investor,1);
            % note: reset is always at 0.5, which is the natural midpoint
            [estimates_m{ i_mod }(i_sub,:), model, fval, exitflag(i_sub,i_mod), A.v_out, A.delta_out] = fitmodel_mle_rl_BIN_bound_onepar_nodecay( A.trustee, A.reset_sub, A.investor, A.condition, A.prior  ); % xdata1 is outcome
            % Note: minimization not for SSE (sum of squared errors) but
            % for "LLE" negative log-likelhood (basically difference bt
            % "normal" regression & logistic regression
            % This changes formulas for AIC/BIC
            [LLE(i_sub,i_mod), FittedCurve2] = model(estimates_m{ i_mod }(i_sub,:)); 
            AIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 3*1; % Akaike Information Criterion  % k1 = 1; % number of free parameters
            BIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 3*log(n)/2; % Bayesian Information Criterion
            
          
%             to_plot = fitmodel_mle_rl_BIN_bound_noreset_to_plot(A.trustee, A.reset_sub, A.investor, A.condition, A.prior, estimates_m{ i_mod }(i_sub,:));
            to_plot = A.v_out;
             elseif i_mod == 3
            % v23
            %% simple RL model with binarized ratings and answers
            % adaptation of i_mod == 1
            n = size(A.investor,1);
            % note: reset is always at 0.5, which is the natural midpoint
            [estimates_m{ i_mod }(i_sub,:), model, fval, exitflag(i_sub,i_mod), A.v_out] =  fitmodel_mle_tit_f_tat_onepar( A.trustee, A.reset_sub, A.investor, A.condition, A.prior); % xdata1 is outcome
            % Note: minimization not for SSE (sum of squared errors) but
            % for "LLE" negative log-likelhood (basically difference bt
            % "normal" regression & logistic regression
            % This changes formulas for AIC/BIC
            [LLE(i_sub,i_mod), FittedCurve2] = model(estimates_m{ i_mod }(i_sub,:)); 
            AIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 2*1; % Akaike Information Criterion  % k1 = 1; % number of free parameters
            BIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 2*log(n)/2; % Bayesian Information Criterion
            
%             to_plot = fitmodel_mle_tit_f_tat_fixed_to_plot(A.trustee, A.reset_sub, A.investor, estimates_m{ i_mod }(i_sub,:));
         
                     to_plot = A.v_out;


          elseif i_mod == 4
            n = size(A.investor,1);
            [estimates_m{ i_mod }(i_sub,:), model, fval, exitflag(i_sub,i_mod), A.v_out] =  fitmodel_mle_tit_f_tat_onepar_nop_nor_nodecay( A.trustee, A.reset_sub, A.investor); % xdata1 is outcome
            [LLE(i_sub,i_mod), FittedCurve2] = model(estimates_m{ i_mod }(i_sub,:));
            AIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 2*1; % Akaike Information Criterion  % k1 = 1; % number of free parameters
            BIC(i_sub,i_mod) = LLE(i_sub,i_mod) + 2*log(n)/2; % Bayesian Information Criterion
%             to_plot = fitmodel_mle_tit_f_tat_fixed_to_plot(A.trustee, A.reset_sub, A.investor, estimates_m{ i_mod }(i_sub,:));
          to_plot = A.v_out;
         end
            
            %to_plot = [];
            
            subplot(n1,n2,i_sub)

%             A.investor(A.investor>0) = 1;
%             to_plot(to_plot>0) = 1;

            plot( A.investor )
            hold on
            plot( to_plot, 'r' )
%             plot( A.trustee, 'k' )
            %stem( to_plot_stem, 'g' )
            title(['m ' num2str(i_mod) ' p ' num2str(i_sub)])
            if i_sub == 1
            legend({'ratings','model','outcome','reset'})
            end
        
            check_size( i_sub, 1 ) = size( A.trustee, 1);
            
       % clear A to_plot
            
            
 end
end

% BIC(sum(exitflag,2)<size(exitflag,2),:)=[];
% AIC(sum(exitflag,2)<size(exitflag,2),:)=[];

% re_arra2 = [2,4,5]; % [1,7]; % 
BIC_sum = sum( BIC,1);
BIC_diff = BIC_sum-max(BIC_sum);
figure
bar( BIC_diff, 'k' )
%  ylim([1000 1800]);
 ylabel('Delta BIC')
 xticklabels({sprintf('rl'),sprintf('rl\\newlinenodecay'),...
     sprintf('tit4tat'),sprintf('tit4tat\\newlinenodecay')})
xtickangle(45)
%%
% clear BMS
% re_arra3=[1 4];
% i_sit = 1;
% % [BMS.alpha( i_sit, :), BMS.exp_r( i_sit, :), BMS.xp( i_sit, :), BMS.pxp( i_sit, :), BMS.bor( i_sit, :)] = spm_BMS( -BIC(:,re_arra2) ); % spm_BMS( -BIC );
% 
% %[BMS.alpha( i_sit, :), BMS.exp_r( i_sit, :), BMS.xp( i_sit, :), BMS.pxp( i_sit, :), BMS.bor( i_sit, :)] = spm_BMS( -AIC(which_kids,re_arra2) ); % spm_BMS( -BIC );
% % [BMS.alpha( i_sit, :), BMS.exp_r( i_sit, :), BMS.xp( i_sit, :), BMS.pxp( i_sit, :), BMS.bor( i_sit, :)] = spm_BMS( -BIC(:,:) ); % spm_BMS( -BIC );
% [BMS.alpha( i_sit, :), BMS.exp_r( i_sit, :), BMS.xp( i_sit, :), BMS.pxp( i_sit, :), BMS.bor( i_sit, :)] = spm_BMS( -BIC(:,re_arra3) ); % spm_BMS( -BIC );
% 
% figure
% bar( BMS.pxp, 'k' )
% ylim([0,1])
%  xticklabels({sprintf('rl'),sprintf('rl\\newline1par'),sprintf('rl\\newline1parnoprior'),sprintf('rl\\newlinenoresetprior'),...
%      sprintf('tit4tat'),sprintf('tit4tat\\newline1par'),sprintf('tit4tat\\newline1parnoprior'),sprintf('tit4tat\\newlinenoresetprior'),...
%      sprintf('nodecay'),sprintf('nodecayreset')})
% xtickangle(45)
%% BIC bootstrapped 95% CI of the summed difference
% clear BIC_bs
% BIC_bs = BIC(:,4)-BIC(:,1); 
% nboot = 1000000;
% sum_boot = NaN(nboot,1);
% for i = 1:nboot
%    % Draw bootstrap sample
%    BIC_boot = randsample(BIC_bs, length(BIC_bs), 1);
%    sum_boot(i) = sum(BIC_boot);
% end
% sumBIC = sum(BIC_bs)
% CI_lowerBIC = quantile(sum_boot, 0.025);
% CI_upperBIC = quantile(sum_boot, 0.975);
% disp('CI_Low_BIC')
% disp(num2str(CI_lowerBIC,'%.2f'))
% disp('CI_Up_BIC')
% disp(num2str(CI_upperBIC,'%.2f'))