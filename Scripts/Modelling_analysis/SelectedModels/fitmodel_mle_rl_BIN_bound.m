function [estimates, model, fval, exitflag, v_vec, delta] = fitmodel_mle_rl_BIN_bound(trustee, reset_sub, ydata, condition, prior) % trustee is outcome, xdat2 is for resetting
% Call fminsearch with a random starting point.
% start_point = rand(1);
start_point = [0.5, 0.5, 0.5, 0.1]; % i.e., HERE p(1) = alpha % p(2) = beta (inverse temperature) 
model = @expfun;
options = optimset('Display', 'final', 'MaxIter', 10000000);
[estimates,fval,exitflag] = fminsearch(model, start_point, options);
% expfun accepts model parameters as inputs, and outputs sse
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [lle, le] = expfun(p) % BIN changes to lle, le from sse, se (lle = log-likelihood)
        
        %% RL 
        v_vec = zeros( length(trustee), 1 );
        for i = 1:length(trustee)
            if reset_sub(i,1) == 1
                switch  condition(i,1) %3nonsocial 4 social
                    case 4
                        v_vec(i,1) = prior; % prior
                    otherwise
                        v_vec(i,1) = .5; % keep rest fixed at 0.5 for BIN version
                end
            else
                if ydata(i-1,1) > 0 %if person shares just do your rl
                    delta(i,1) = trustee(i-1,1)-v_vec(i-1,1);
%                     delta(i,1) = trustee(i-1,1)-v_vec(i-1,1);
                    
                    switch condition(i,1) %3nonsocial 4 social
                        case 3 
                             v_vec(i,1) = v_vec(i-1,1) + p(1)*delta(i,1);
                        case 4
                             v_vec(i,1) = v_vec(i-1,1) + p(2)*delta(i,1);
                    end            
                elseif ydata(i-1,1) == 0
                %% here add some probability
                
                
                    %% DECAY version that mimics RL
                    % Basically we just assume that the investor acts as if
                    % the trustee has shared (although this has of course 
                    % not happend). We do the very same RL as above with a
                    % new learning rate. This makes sure that v_vec cannot
                    % exceed 1.
                    % scheme.
                    delta_decay = 1-v_vec(i-1,1); % here "1" instead of "trustee(i,1)" 
                    % because we assume that trustee has shared and that
                    % the investor calculates the delta_decay (i.e., PE)
                    % for such a scenario
                    v_vec(i,1) = v_vec(i-1,1) + p(3)*delta_decay; 
                    %%
%                     v_vec(i,1) = p(3);

                end
            end
              
%             end
        end
%%
%         % BIN: do the sigmoid    
%         vec = [0:1];
% %         v_vec(v_vec>5) = 5;
%         for j = 1:length(v_vec)
%             diff_vec = abs(v_vec(j) - vec);
%             exp_vec  = exp( -diff_vec*p(4));
%             sigmoid_vec = exp_vec / sum( exp_vec );
%             le(j,1) = sigmoid_vec(ydata(j)+1);
%         end
%         lle = -sum(log(le));
%% 
        prob = sigmoid( v_vec, p(4) );

%         BIN: do le, lle instead of sse, se
%         distanz proband siugmoid 1 oder distanz zu 0
        
        le = ydata .* (1 - prob)  + (1 - ydata) .* (prob);
        lle = -sum(log(le));
        %%
         if p(1) <= 0 || p(1) >= 1
            lle = 10^10;
         end
         
         if p(2) <= 0 || p(2) >= 1
            lle = 10^10;
         end
         
         if p(3) <= 0 || p(3) >= 1
            lle = 10^10;
         end
        end
end