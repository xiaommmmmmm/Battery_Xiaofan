function y = func_plot_task_lsv5_1(CL,ENV,ExpRes,ENV_trad,ExpRes_trad)
        figure(13);
        for i = 1:size(ExpRes.bat_uratio_energy,1)
            subplot(4,2,i);
            pic(i) = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
                100*mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_trad(i) = plot(ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum, ...
            100*mean(reshape(ExpRes_trad.bat_uratio_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'s--','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            xlabel('Converter Energy Rating/Battery Energy Capability');
            ylabel('Battery Capacity Utilization (%)');
            str1 = strcat({'Battery Capacity Variation (sorted)'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.qlim_var(i)/ENV.Sweep.Bat{1}.qlim_mu(i))},{'%'});
 %           str2 = {};
            title(str1,'Fontsize',9);
            legend([pic(i),pic_trad(i)],'LSV-PPP','C-PPP','Fontsize',10);
            grid on;
            grid minor;
        end
    % xlabel('Sum of Converter energy ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
 
            %% Plot    
        figure(14);
        for i = 1:size(ExpRes.bat_uratio_energy,1)
            subplot(4,2,i);
            pic(i) = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            100-100*mean(reshape(ExpRes.total_energy_process(i,:,:)*0.15./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_trad(i) = plot(ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum, ...
            100-100*mean(reshape(ExpRes_trad.total_energy_process(i,:,:)*0.15./ExpRes_trad.max_output_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'s--','linewidth',2,'color',[i/size(ExpRes_trad.bat_uratio_energy,1),0,0]);
            xlabel('Converter Energy Rating/Battery Energy Capability');
            ylabel('System Energy Efficiency (%)');
            str1 = strcat({'Battery Capacity Variation (sorted)'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.qlim_var(i)/ENV.Sweep.Bat{1}.qlim_mu(i))},{'%'});
            title(str1,'Fontsize',9);
            legend([pic(i),pic_trad(i)],'LSV-PPP','C-PPP','Fontsize',10);
            grid on;
            grid minor;
        end
end