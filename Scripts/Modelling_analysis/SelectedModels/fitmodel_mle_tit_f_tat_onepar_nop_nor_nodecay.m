function [estimates, model, fval, exitflag, v_vec] = fitmodel_mle_tit_f_tat_onepar_nop_nor_nodecay(trustee, reset_sub, ydata) % xdata1 is outcome, xdat2 is for resetting
% Call fminsearch with a random starting point.
% start_point = rand(1);
start_point = [0.1, -12]; % i.e., HERE p(1) = alpha % p(2) = beta (inverse temperature)
model = @expfun;
options = optimset('Display', 'final', 'MaxIter', 10000000);
[estimates,fval,exitflag] = fminsearch(model, start_point, options);
% expfun accepts model parameters as inputs, and outputs sse
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [lle, le] = expfun(p) % BIN changes to lle, le from sse, se (lle = log-likelihood)
        
        %% RL
        v_vec = zeros( length(trustee), 1 );
        %         v_vec(1,1) = 5; % set initial value to 0.5
        for i = 1:length(trustee)
            if i == 1
               v_vec(i,1) = 0.5; % keep rest fixed at 0.5 for BIN version            
            else
            
                if ydata(i-1,1) > 0 %if person shares just do your rl
                    %t4t 0/1
                    if trustee(i-1,1) == 0
                        v_vec(i,1) = 0;
                    else
                        v_vec(i,1) = 1;
                    end
%                     %t4t prob
%                     v_vec(i,1) = trustee(i-1,1)*ydata(i-1,1);
                    
                elseif ydata(i-1,1) == 0
                  
                    %% DECAY version that mimics RL
%                     % Basically we just assume that the investor acts as if
%                     % the trustee has shared (although this has of course 
%                     % not happend). We do the very same RL as above with a
%                     % new learning rate. This makes sure that v_vec cannot
%                     % exceed 1.
%                     % scheme.
%                     delta_decay = 1-v_vec(i-1,1); % here "1" instead of "xdata1(i,1)" 
%                     % because we assume that trustee has shared and that
%                     % the investor calculates the delta_decay (i.e., PE)
%                     % for such a scenario
% 
%                     v_vec(i,1) = v_vec(i-1,1) + p(1)*delta_decay;


                      delta_decay = 0-v_vec(i-1,1);
                      v_vec(i,1) = v_vec(i-1,1) + p(1)*delta_decay;

%                     v_vec(i,1) = p(1);
                    %% here add some probability
                   % if y(xdata3(i,1))==0
                  % v_vec(i+1,1) = p(1)*(xdata3(i,1)+1);
                  %  else
                   % end
                    
                     %% insert decay
                %else
                  
                % here additionally decay the calculated v_vec(i+1,1)
                %decay_help = v_vec(i,1) + p(1)*delta(i,1);                
                %for i_d = 1:decay_tria
                %    decay_help = p(3)*(decay_help-reset) + reset;
               % end
               % v_vec(i+1,1) = decay_help;
                                
                end
            end
            
        end
        % BIN: do the sigmoid
        
%         v_vec(end) = [];
       % v_vec(v_vec>1)=1;
        % careful delta has all trials
        % BIN: do the sigmoid
        
        %% replace prob
%         vec = [0:1];
% %         v_vec(v_vec>5) = 5;
%         for j = 1:length(v_vec)
%             diff_vec = abs(v_vec(j) - vec);
%             exp_vec  = exp( -diff_vec*p(2));
%             sigmoid_vec = exp_vec / sum( exp_vec );
%             le(j,1) = sigmoid_vec(ydata(j)+1);
%         end
% %         le = ydata .* (1 - prob)  + (1 - ydata) .* (prob);
%         lle = -sum(log(le));
        %% 
        prob = sigmoid( v_vec, p(2) );

%         BIN: do le, lle instead of sse, se
%         distanz proband siugmoid 1 oder distanz zu 0
        
        le = ydata .* (1 - prob)  + (1 - ydata) .* (prob);
        lle = -sum(log(le));
         if p(1) <= 0 || p(1) >= 1
            lle = 10^10;
         end

    end
end