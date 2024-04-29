function OptRes = func_ac_powerflow_avg_layer(P_diff_mat_in, CL, ENV)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
P_diff_mat = (P_diff_mat_in>0) | (P_diff_mat_in>0)';           % the interconnection matrix must be symmetric
Con_diff_num = sum(sum(P_diff_mat))/2;      % count the number of differential power delivery connections
eta = 1; % assuming the converter efficiency is 100% 
Conv_power_rating_partition_mat = zeros(CL.Stat.Bat_num,CL.Stat.Bat_num);
Node_num = CL.Stat.Bat_num;
T_num = CL.Stat.Delta_t_num;
Delta_t = CL.Stat.Delta_t;
Output_waveform = CL.Stat.Output;
%% Construct the optimization problems and solve the linear programming problems 
I_B = optimvar('I_B', Node_num,T_num); % battery charge
I_C_s = optimvar('I_C_s',Node_num,T_num); % sparse converter charge
% I_C_d = optimvar('I_C_d',Node_num-1,T_num); % dense converter charge

Q_B = optimvar('Q_B', Node_num,T_num); % battery charge
Q_C_s = optimvar('Q_C_s',Node_num,T_num); % sparse converter charge
% Q_C_d = optimvar('Q_C_d',Node_num-1,T_num); % dense converter charge
Q_L = optimvar('Q_L',Node_num,T_num);  % load charge

P = optimvar('P',Node_num, Node_num,T_num);    % dimensional power transfer variables
% bus_curlim = optimvar('bus_curlim',size(CL.Bus,2),1);     % bus current capability
% bus_volt = [];
% for i = 1:size(CL.Bus,2)
%     CL.Bus{i}.curlim = bus_curlim(i);
%     CL.Bus{i}.volt = 0;
%     for  j = 1:size(CL.Bat,2)
%        CL.Bus{i}.volt = CL.Bus{i}.volt + CL.Bat{j}.volt;
%     end
%     bus_volt(i) = CL.Bus{i}.volt;
% end
prob = optimproblem('Objective',sum(sum(Q_L)),'ObjectiveSense','min');  % minimize the summations

% Construct the constraint for charge matrix
topo_ql_con = optimconstr(0);
for i = 1:size(I_B,1)
    for j = 1:T_num
        topo_ql_con(end + 1) = Q_B(i,j) == I_B(i,j) * Delta_t(j);
        topo_ql_con(end + 1) = Q_C_s(i,j) == I_C_s(i,j) * Delta_t(j);
        topo_ql_con(end + 1) = Q_L(i,j) == Q_B(i,j) + Q_C_s(i,j);
    end
end
% topo_ql_con(end + 1) = Q_C_d == I_C_d .* Delta_t;

prob.Constraints.topo_ql_con = topo_ql_con;

% Construct the constriant for differential power processing links
topo_diff_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_mat(i,j) == 1)
            for k = 1 :T_num
                topo_diff_con(end+1) = P(i,j,k) == -P(j,i,k);
            end
        else
            for k = 1 :T_num
                topo_diff_con(end+1) = P(i,j,k) == 0;
                topo_diff_con(end+1) = P(j,i,k) == 0; 
            end
        end
    end
end
 prob.Constraints.topo_diff_con = topo_diff_con;

% Construct the Average Converter Power Limit
avg_conv_con = optimconstr(0);
temp_ct = 1;
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_mat(i,j) == 1)
            for k = 1 :T_num
    %             avg_conv_con(end+1) = P(i,j) <= ENV.Avg_Conv.e_lim(temp_ct);
    %             avg_conv_con(end+1) = -ENV.Avg_Conv.e_lim(temp_ct) <= P(i,j);
    %             Conv_power_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
                avg_conv_con(end+1) = P(i,j,k) <= P_diff_mat_in(i,j);
                avg_conv_con(end+1) = -P_diff_mat_in(i,j) <= P(i,j,k);
            end
            Conv_power_rating_partition_mat(i,j) = P_diff_mat_in(i,j);
            temp_ct = temp_ct + 1;
        end
    end
end
prob.Constraints.avg_conv_con = avg_conv_con;

% Construct the constriant for Qc
power_qc_con = optimconstr(0);
P_sum = reshape(sum(P,2),Node_num,T_num);
for i = 1:size(CL.Bus,2)
        if (sum(P_diff_mat(i,:)) ~= 0)
            for k = 1 :T_num
                power_qc_con(end+1) = I_C_s(i,k)* CL.Bat{i}.volt ==  P_sum(i,k);
            end
        else
            for k = 1 :T_num
                power_qc_con(end+1) = I_C_s(i,k) == 0;
            end
        end
end
prob.Constraints.power_qc_con = power_qc_con;

% Construct the constriant for qb_I
power_qb_con = optimconstr(0);
for i = 1: size(CL.Bat,2)
    for k = 1 :T_num
        power_qb_con(end+1) = CL.Bat{i}.curlim >= I_B(i,k) + I_C_s(i,k);
        power_qb_con(end+1) = I_B(i,k)+I_C_s(i,k) >= -CL.Bat{i}.curlim;
    end
end
prob.Constraints.power_qb_con = power_qb_con;

% Construct the constriant for U_e
ue_con = optimconstr(0);
for i = 2:size(CL.Bat,2)
   ue_con(end+1) =  sum(Q_B(i,:)) / (CL.Bat{i}.qlim/CL.Bat{i}.volt) == sum(Q_B(1,:)) / (CL.Bat{1}.qlim/CL.Bat{1}.volt);
end
prob.Constraints.ue_con = ue_con;

% % Construct constant power flow of converters
% constant_powerflow_con = optimconstr(0);
% for i = 1: Node_num
%     for j = 1 : T_num 
%         constant_powerflow_con(end+1) = I_C_s(i,j) == I_C_s(i,1);
%     end
% end
% prob.Constraints.constant_powerflow_con = constant_powerflow_con;

output_con = optimconstr(0);
for k = 1 :T_num
    sum_u = 0;
    for i = 1: Node_num
        sum_u = sum_u + Q_L(i,k)/(CL.Bat{i}.qlim/CL.Bat{i}.volt);
    end
    output_con(end + 1) = sum_u == abs(Output_waveform(k));
end
prob.Constraints.output_con = output_con;

problem = prob2struct(prob);
[sol,fval,exitflag,output] = linprog(problem);

if(exitflag == 1)
    OptRes.sol = sol;
    OptRes.u_p = fval;
    OptRes.Conv_power_rating_partition_mat = Conv_power_rating_partition_mat + Conv_power_rating_partition_mat';
    OptRes.P_diff = P_diff_mat;
%    OptRes.P_direct = P_direct_mat;

%% decomposition sol
% I_B = optimvar('I_B', Node_num,T_num); % battery charge
% I_C_s = optimvar('I_C_s',Node_num,T_num); % sparse converter charge
% % I_C_d = optimvar('I_C_d',Node_num-1,T_num); % dense converter charge
% 
% Q_B = optimvar('Q_B', Node_num,T_num); % battery charge
% Q_C_s = optimvar('Q_C_s',Node_num,T_num); % sparse converter charge
% % Q_C_d = optimvar('Q_C_d',Node_num-1,T_num); % dense converter charge
% Q_L = optimvar('Q_L',Node_num,T_num);  % load charge
% 
% P = optimvar('P',Node_num, Node_num,T_num);    % dimensional power transfer variables
OptRes.I_B = reshape(OptRes.sol(0*Node_num*T_num+1:1*Node_num*T_num),[Node_num,T_num]);
    OptRes.I_C_s = reshape(OptRes.sol(1*Node_num*T_num+1:2*Node_num*T_num),[Node_num,T_num]);
    OptRes.Q_B = reshape(OptRes.sol(2*Node_num*T_num+1:3*Node_num*T_num),[Node_num,T_num]);
    OptRes.Q_C_s = reshape(OptRes.sol(3*Node_num*T_num+1:4*Node_num*T_num),[Node_num,T_num]);
    OptRes.Q_L = reshape(OptRes.sol(4*Node_num*T_num+1:5*Node_num*T_num),[Node_num,T_num]);
%    OptRes.diff_mat_avg = P_diff_mat;
%    OptRes.Conv_power_processed = abs(sol((Node_num*Node_num+1):(Node_num*Node_num+Con_diff_num)));
    
else
    error('Error in optimizing dc power flow of the averaging converters\n');
end

function [Indirect_power_flow,Direct_power_flow] = func_power_flow_decomp(Total_power_flow,CL)
    Indirect_power_flow = zeros(size(Total_power_flow,1),size(Total_power_flow,2));
    Direct_power_flow = zeros(size(Total_power_flow,1),size(Total_power_flow,2));
    for ic = 1:size(Total_power_flow,1)
        for jc = 1:size(Total_power_flow,2)
           if ((ic < CL.Bus{1}.ind)&&(jc < CL.Bus{1}.ind))
               Indirect_power_flow(ic,jc) = Total_power_flow(ic,jc);
           else
               Direct_power_flow(ic,jc) = Total_power_flow(ic,jc);
           end
        end
    end 
end

end

