%% tit for tat
clear all
%model random share 0~5 on first trial
for ite = 1:1000
    trialnum=20; %number of trials
    alpha=0.5; %free parameter
    model(1,1) = randi([1 5]);
    %model keep on first trial 
    if model(1,1) == 0
        game(1,1) = 0;
        %model share on first trial
        else game(1,1) = 1;
    end
    % from second trial
    for trial = 2:trialnum
        %model share on previous
        if model(trial-1,1) > 0 
            % tit for tat
            %game keep on previous
            if game(trial-1,1) == 0
            model(trial,1) = 0;
            game(trial,1) = 0;
            %game share on previous
            else
                vec = [1:5];
                diff_vec = abs(model(trial-1,1)-vec);
    %             exp_vec  = exp( -diff_vec*p(3));
                exp_vec  = exp( -diff_vec*rand*10);
                sigmoid_vec = exp_vec / sum( exp_vec );
                prob = rand;
                for j = 1:5
                    if prob<sum(sigmoid_vec(1:j))
                        model(trial,1) = vec(j);
                        break
                    else
                        continue
                    end
                end
                         % game
            % CooperativeAgent_Logistic_PIML(s_t,s_t1,alpha);
            delta_t = (model(trial,1) - model(trial-1,1));
            Pshare = 1/(1+exp(-1*alpha*(model(trial,1) + delta_t)));
            game(trial,1) = +(rand < Pshare);  
            end
        else 
           %model keep on previous
           model(trial,1) = 0;
           game(trial,1) = 0;
        end  
    end
total_model(:,ite) = model;
total_game(:,ite) = game;
clearvars -except total_game total_model
end
toplot(:,1) = mean(total_model,2);
toplot(:,2) = mean(total_game,2);
ttd(:,1) = std(total_model')%/sqrt(1000);
ttd(:,2) = std(total_game')%/sqrt(1000);

figure
histogram(total_model)

trial(1,:) = [1:20];
figure
plot(toplot)
hold on
  ylim([0 5]);
  shadedErrorBar(trial,toplot(:,1),ttd(:,1),'lineprops','-b','patchSaturation',0.33)
  hold on
  shadedErrorBar(trial,toplot(:,2),ttd(:,2),'lineprops','-y','patchSaturation',0.33)

%% generous tit for tat
clear all
%model random share 0~5 on first trial
for ite = 1:1000
    trialnum=20; %number of trials
    alpha=0.5; %free parameter
    model(1,1) = randi([1 5]);
    %model keep on first trial 
    if model(1,1) == 0
        game(1,1) = 0;
        %model share on first trial
        else game(1,1) = 1;
    end
    % from second trial
    for trial = 2:trialnum
        %model share on previous
        if model(trial-1,1) > 0 
            % tit for tat
            %game keep on previous
            if game(trial-1,1) == 0
            model(trial,1) = 0;
            game(trial,1) = 0;
            %game share on previous
            else
                vec = [0:5];
                diff_vec = abs(model(trial-1,1)-vec);
    %             exp_vec  = exp( -diff_vec*p(3));
                exp_vec  = exp( -diff_vec*rand);
                sigmoid_vec = exp_vec / sum( exp_vec );
                prob = rand;
                for j = 1:6
                    if prob<sum(sigmoid_vec(1:j))
                        model(trial,1) = vec(j);
                        break
                    end
                end
            % game
            % CooperativeAgent_Logistic_PIML(s_t,s_t1,alpha);
            delta_t = (model(trial,1) - model(trial-1,1));
            Pshare = 1/(1+exp(-1*alpha*(model(trial,1) + delta_t)));
            game(trial,1) = +(rand < Pshare);  
            end
        else 
           %model keep on previous
            delta_decay = 5-model(trial-1); 
            tempo = round(model(trial-1,1) + rand*delta_decay);
            vec = [0:5];
            diff_vec = abs(tempo-vec);
    %       exp_vec  = exp( -diff_vec*p(3));
            exp_vec  = exp( -diff_vec*rand);
            sigmoid_vec = exp_vec / sum( exp_vec );
            prob = rand;
                for j = 1:6
                    if prob<sum(sigmoid_vec(1:j))
                        model(trial,1) = vec(j);
                        break
                    end
                end
                if model(trial,1) == 0
                    game(trial,1) = 0;
                else
                % game
                % CooperativeAgent_Logistic_PIML(s_t,s_t1,alpha);
                delta_t = (model(trial,1) - model(trial-1,1));
                Pshare = 1/(1+exp(-1*alpha*(model(trial,1) + delta_t)));
                game(trial,1) = +(rand < Pshare);  
                end
        end  
    end
total_model(:,ite) = model;
total_game(:,ite) = game;
clearvars -except total_game total_model
end
toplot(:,1) = mean(total_model,2);
toplot(:,2) = mean(total_game,2);
ttd(:,1) = std(total_model')/sqrt(1000);
ttd(:,2) = std(total_game')/sqrt(1000);

figure
histogram(total_model)

trial(1,:) = [1:20];
figure
plot(toplot)
hold on
  ylim([0 5]);
  shadedErrorBar(trial,toplot(:,1),ttd(:,1),'lineprops','-b','patchSaturation',0.33)
  hold on
  shadedErrorBar(trial,toplot(:,2),ttd(:,2),'lineprops','-y','patchSaturation',0.33)

%% RL
clear all
%model random share 0~5 on first trial
for ite = 1:1000
    trialnum=20; %number of trials
    alpha=0.5; %free parameter
    model(1,1) = randi([0 5]);
    %model keep on first trial 
    if model(1,1) == 0
        game(1,1) = 0;
        %model share on first trial
        else game(1,1) = 1;
    end
    % from second trial
    for trial = 2:trialnum
        if model(trial-1,1) > 0 
            delta = game(trial-1,1) - model(trial-1); 
            tempo = round(model(trial-1,1) + rand*delta);
            vec = [0:5];
            diff_vec = abs(tempo-vec);
    %       exp_vec  = exp( -diff_vec*p(3));
            exp_vec  = exp( -diff_vec*rand);
            sigmoid_vec = exp_vec / sum( exp_vec );
            prob = rand;
                for j = 1:6
                    if prob<sum(sigmoid_vec(1:j))
                        model(trial,1) = vec(j);
                        break
                    end
                end
                if model(trial,1) == 0
                    game(trial,1) = 0;
                else
                % game
                % CooperativeAgent_Logistic_PIML(s_t,s_t1,alpha);
                delta_t = (model(trial,1) - model(trial-1,1));
                Pshare = 1/(1+exp(-1*alpha*(model(trial,1) + delta_t)));
                game(trial,1) = +(rand < Pshare);  
                end
        else
           %model keep on previous
            delta_decay = 5-model(trial-1); 
            tempo = round(model(trial-1,1) + rand*delta_decay);
            vec = [0:5];
            diff_vec = abs(tempo-vec);
    %       exp_vec  = exp( -diff_vec*p(3));
            exp_vec  = exp( -diff_vec*rand);
            sigmoid_vec = exp_vec / sum( exp_vec );
            prob = rand;
                for j = 1:6
                    if prob<sum(sigmoid_vec(1:j))
                        model(trial,1) = vec(j);
                        break
                    end
                end
                if model(trial,1) == 0
                    game(trial,1) = 0;
                else
                % game
                % CooperativeAgent_Logistic_PIML(s_t,s_t1,alpha);
                delta_t = (model(trial,1) - model(trial-1,1));
                Pshare = 1/(1+exp(-1*alpha*(model(trial,1) + delta_t)));
                game(trial,1) = +(rand < Pshare);  
                end 
        end
    end 
total_model(:,ite) = model;
total_game(:,ite) = game;
clearvars -except total_game total_model
end
toplot(:,1) = mean(total_model,2);
toplot(:,2) = mean(total_game,2);
figure
histogram(total_model)
%%
histogram(total_model)