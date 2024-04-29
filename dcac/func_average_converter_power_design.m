function OptAvg = func_average_converter_power_design(CL,ENV)
%Input the diffferetial power delivery connection
trial_num = ENV.Avg_Conv.trial_num;
converter_num = ENV.Avg_Conv.Num;
Node_num = CL.Stat.Bat_num;

comb_explosion_lim = 1e6;
%[tempX,tempY] = ndgrid([CL.Bat,CL.Cap],[CL.Bat,CL.Cap]);
%all_diff_conn = [reshape(tempX,[size(tempX,1)*size(tempX,2),1]),reshape(tempY,[size(tempY,1)*size(tempY,2),1])];
all_diff_conn = nchoosek([CL.Bat],2);
linked_conn = {};
exhaust_num = nchoosek(size(all_diff_conn,1),converter_num);

% for i = 1:size(CL.Bat,2)  % use the averaged value of q for primary converter design
%     CL.Bat{i}.qlim = CL.Bat{i}.qlim_mu ;
% end
CL = Flatten_battery_stat(CL); % use the flatened distributions of q for primary converter design

Bat_power_rating = zeros(CL.Stat.Bat_num,1);
for i = 1:CL.Stat.Bat_num 
    Bat_power_rating(i) = CL.Bat{i}.volt*CL.Bat{i}.curlim;
end

if (exhaust_num<=trial_num)&&(exhaust_num<=comb_explosion_lim)
    all_diff_conn_choose_conv = nchoosek(1:size(all_diff_conn,1),converter_num);
    for j = 1:nchoosek(size(all_diff_conn,1),converter_num)
        linked_conn{end+1} =  all_diff_conn(all_diff_conn_choose_conv(j,:),:);  % if the search space is too large, do random sampling
    end    
    fprintf('%s%d','The search space',nchoosek(size(all_diff_conn,1),converter_num),'is complete');
    fprintf('\n');
else 
    for j = 1:trial_num
        linked_conn{end+1} = datasample(all_diff_conn,converter_num,'Replace',false);   % if the search space is too large, do random sampling
    end
    fprintf('%s%d%s%d','The search space ',nchoosek(size(all_diff_conn,1),converter_num),...
        ' is partially reached by ',trial_num);
    fprintf('\n');
end
%% Solve the average power flow by parallel computing
for i = 1:min(trial_num,exhaust_num)
    P_diff_mat_in = zeros(CL.Stat.Bat_num, CL.Stat.Bat_num);    % P_diff_mat_in connection should be reset in every cycle
    for j = 1:converter_num
         P_diff_mat_in(linked_conn{i}{j,1}.ind, linked_conn{i}{j,2}.ind) = ENV.Avg_Conv.p_lim_vec(j);
    end
    OptRes(i) = func_ac_powerflow_avg_layer(P_diff_mat_in + P_diff_mat_in', CL, ENV);
end
%% Pick the Optimal Data
max_u_p = [];
for i = 1:min(trial_num,exhaust_num)
    max_u_p = [max_u_p,OptRes(i).u_p];
end
OptAvg_power_collection = {};
for i = 1:min(trial_num,exhaust_num)
    if (OptRes(i).u_p == max(max_u_p))
        
        OptAvg_power_collection{end+1}.u_p = abs(OptRes(i).u_p);
        OptAvg_power_collection{end}.Avg_power_flow = OptRes(i).Avg_power_flow;
        
        OptAvg_power_collection{end}.diff_avg_mat = (OptRes(i).Avg_power_flow ~= 0);
        OptAvg_power_collection{end}.Conv_power_rating_partition_mat = OptRes(i).Conv_power_rating_partition_mat;
        OptAvg_power_collection{end}.Avg_conv_power_processed = abs(OptRes(i).Avg_power_flow);

        OptAvg_power_collection{end}.Bat_power_rating = Bat_power_rating;
        OptAvg_power_collection{end}.Bat_uratio = abs(OptRes(i).bus_plim)/sum(Bat_power_rating);
        OptAvg_power_collection{end}.Avg_conv_uratio = sum(sum(abs(OptRes(i).Avg_power_flow)))/sum(sum(OptRes(i).Conv_power_rating_partition_mat));
%        OptAvg.Conv_power_rate_partition = func_rate_partition(OptRes(i).Conv_power_rate,ENV.Avg_Conv.partition);         % sort the cnverter rating in ascending orders
%       fprintf('%s%d\n','Battery 1 to Battery 2: ',OptRes(i).sol((CL.Bat{2}.ind-1)*Node_num+CL.Bat{1}.ind));
%       fprintf('%s%d\n','Battery 2 to Battery 3: ',OptRes(i).sol((CL.Bat{3}.ind-1)*Node_num+CL.Bat{2}.ind));
%       fprintf('%s%d\n','Battery 3 to Capacitor 1: ',OptRes(i).sol((CL.Cap{1}.ind-1)*Node_num+CL.Bat{3}.ind));
%         fprintf('%s%d\n','Battery 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{1}.ind));
%         fprintf('%s%d\n','Battery 2 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{2}.ind));
%         fprintf('%s%d\n','Battery 3 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{3}.ind));
%         fprintf('%s%d\n','Capacitor 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Cap{1}.ind));
    end
end

if(size(OptAvg_power_collection,2) > 1)
    fprintf('Multiple Optimal Points, use lexicographic optimization\n');
    % optimize for the converter rating
    for i = 1:size(OptAvg_power_collection,2)
        converter_rating(i) = max(max(OptAvg_power_collection{i}.Avg_conv_power_processed));
    end
    OptAvg_power_rating_optimized = [];
    for i = 1:size(OptAvg_power_collection,2)
        if(max(max(OptAvg_power_collection{i}.Avg_conv_power_processed)) == min(converter_rating))
           OptAvg_power_rating_optimized{end+1} = OptAvg_power_collection{i};
        end
    end
    % optimize for the processed power
    for i = 1:size(OptAvg_power_rating_optimized,2)
        processed_power(i) = sum(sum(OptAvg_power_rating_optimized{i}.Avg_conv_power_processed))/2;
    end
    OptAvg_power_processed_optimized = [];
    for i = 1:size(OptAvg_power_rating_optimized,2)
        if((sum(sum(OptAvg_power_rating_optimized{i}.Avg_conv_power_processed))/2) == min(processed_power))
           OptAvg_power_processed_optimized{end+1} = OptAvg_power_rating_optimized{i};
        end
    end
     % if still multiple optimal, select by hands...
    if(size(OptAvg_power_processed_optimized,2)==1)
        OptAvg = OptAvg_power_processed_optimized{1};
    else
        save('PreDesign_power_collection', 'OptAvg_power_processed_optimized');
        fprintf('Still have multiple Optimal Points after lexicographic optimization\n');
        fprintf('Open PreDesign_power_collection.m \n');
        fprintf('The number of optimal points is: \n');
        size(OptAvg_power_processed_optimized,2)
        u_in = input('Assign xth entry of OptAvg_power_processed_optimized to OptAvg, Type a number: \n');
        OptAvg = OptAvg_power_processed_optimized{u_in};
    end
else
    OptAvg = OptAvg_power_collection{1};
end

function y = Flatten_battery_stat(CL)
    temp = norminv(linspace(1e-4,1-1e-4,CL.Stat.Bat_num+1),CL.Bat{1}.curlim_mu, CL.Bat{1}.curlim_var);   %Normal inverse cumulative distribution function
    for k = 1:CL.Stat.Bat_num
       x_int = linspace(temp(k),temp(k+1),100);
       y_int = normpdf(x_int,CL.Bat{1}.curlim_mu, CL.Bat{1}.curlim_var).*x_int;
       CL.Bat{k}.curlim = CL.Stat.Bat_num*trapz(x_int,y_int);
    end
    y = CL;

%     % temp = CL.Bat{1}.curinvcdf(linspace(1e-4,1-1e-4,CL.Stat.Bat_num+1));   %Normal inverse cumulative distribution function
%     N_int = 100;
%     for k = 1:CL.Stat.Bat_num
%        x_int = linspace(temp(k),temp(k+1),N_int);
%        y_int = CL.Bat{1}.curcdf(x_int);
%        % we do integral by part, int xp(x)dx = delta xc(x) - int c(x)dx
%        CL.Bat{k}.curlim = ((x_int(N_int)*CL.Bat{1}.curcdf(x_int(N_int)) - ...
%            x_int(1)*CL.Bat{1}.curcdf(x_int(1)) - trapz(x_int, y_int)))/(1/CL.Stat.Bat_num);
%        bat_vec(k) = CL.Bat{k}.curlim;
%     end
%     y = CL;

end

end

