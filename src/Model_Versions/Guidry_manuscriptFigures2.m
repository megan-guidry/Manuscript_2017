%% Figure 4: Two disting ES curves

%Plotting Isometric Bars
filename = 'Isometric_D_SL';
%slValue = {'1.955', '2.024', '2.093', '2.162', '2.231', '2.3'};
slValue = {'1.955', '2.3'};
length_slValue = size(slValue);
len_SLValue = length_slValue(2);
line_color = [0.75, 0.75, 0.75];
lw = 1.25;
ES_force = [];
ES_length = [];

for i=1:len_SLValue
    SLv = slValue{i};
    data = strcat(filename,SLv,'.csv')
    
    [time, SL_norm, F_total_norm, Ca_i, dTropTot] = reading_Isometric_DATA(data);
    [M, I] = max(F_total_norm); %%%
    
    figure(4)
    h1 = plot(SL_norm, F_total_norm, 'color', line_color,'LineWidth',lw); hold on;
    ES_force(i) = F_total_norm(I)
    ES_length(i) = SL_norm(I)
end
h1_ES = plot(ES_length, ES_force, '-o', 'color', line_color, 'markers', 10, 'LineWidth', 2); %%%

%Plotting Work-loops
filename = 'WL_D_afterload';
%afterloadValue = {'0.135', '0.1911', '0.2576', '0.3359', '0.4312', '0.6'};
afterloadValue = {'0.135', '0.1911'};
length_afterloadValue = size(afterloadValue);
len_afterloads = length_afterloadValue(2);
line_color = [0, 0, 0];
lw = 1.25;
ES_force = [];
ES_length = [];

    for i=1:len_afterloads
        afterload = afterloadValue{i};
        data = strcat(filename,afterload,'.csv');
        
        [time, SL_norm, F_total_norm, Ca_i, dTropTot, ESmarker] = reading_WL_DATA(data);
        ES_point = find(ESmarker) %%%

        hold on
        figure(4)
        h2 = plot(SL_norm, F_total_norm, 'k-','LineWidth',lw); hold on;
        ES_force(i) = F_total_norm(ES_point);
        ES_length(i) = SL_norm(ES_point);
    end
h2_ES = plot(ES_length, ES_force, '-o', 'color', line_color, 'markers', 10, 'LineWidth', 2); %%%

    ax = gca;
    ax.XTick = [0.85 0.9 0.95 1];
    ax.YTick = [0.2 0.4 0.6 0.8 1];
    axis([0.84 1.01 0 1.2]);
    set(gca,'fontsize',14)
    xlabel('Normalised Sarcomere Length');
    ylabel('Normalised Total Force');
    box off;
    
legend([h1_ES, h2_ES],{'Isometric', 'Work-loop'}, 'Location', 'northwest', 'EdgeColor','w');

%% Figure 5: WL vs. Isometric Ca2+

%Plotting Isometric Bars
filename = 'Isometric_D_SL';
%slValue = {'1.955', '2.024', '2.093', '2.162', '2.231', '2.3'};
slValue = {'1.955', '2.3'};
length_slValue = size(slValue);
len_SLValue = length_slValue(2);
line_color = [0.75, 0.75, 0.75];
lw = 1.25;


for i=1:len_SLValue
    SLv = slValue{i};
    data = strcat(filename,SLv,'.csv')
    
    [time, SL_norm, F_total_norm, Ca_i, dTropTot] = reading_Isometric_DATA(data);
    
    figure(5)
    h1 = plot(time, Ca_i, 'color', line_color,'LineWidth',lw); hold on;
    
    if strcmp(SLv, '2.3') == 1
        figure(9) % figure 5 inset
        h3 = plot(time, Ca_i, 'color', line_color,'LineWidth',2.5); hold on;
    end
end

%Plotting Work-loops
filename = 'WL_D_afterload';
%afterloadValue = {'0.135', '0.1911', '0.2576', '0.3359', '0.4312', '0.6'};
afterloadValue = {'0.135', '0.1911'};
length_afterloadValue = size(afterloadValue);
len_afterloads = length_afterloadValue(2);
line_color = [0, 0, 0];
lw = 1.25;


    for i=1:len_afterloads
        afterload = afterloadValue{i};
        data = strcat(filename,afterload,'.csv');
        
        [time, SL_norm, F_total_norm, Ca_i, dTropTot, ESmarker] = reading_WL_DATA(data);

        hold on; figure(5)
        h2 = plot(time, Ca_i, 'k-','LineWidth',lw); hold on;
        
        if strcmp(afterload, '0.135') == 1
            hold on
            figure(9)
            h4 = plot(time, Ca_i, 'color', line_color,'LineWidth',2.5); hold on;
        end
    end

hold on; figure(5)
set(gca,'fontsize',14, 'XTick', [0 50 100 150 200 250 300], 'YTick', [0 0.2 0.4 0.6 0.8 1 1.2])
axis([0 300 0 1.4]);
xlabel('Time (ms)');
ylabel('[Ca^2^+]_i (\muM)');
box off;   
legend([h1, h2],{'Isometric', 'Work-loop'}, 'Location', 'north', 'EdgeColor','w');

hold on; figure(9)
set(gca,'fontsize',25, 'XTick', [150 300], 'YTick', [0.6 1.2])
axis([0 300 0 1.4]);
xlabel('Time (ms)');
ylabel('[Ca^2^+]_i (\muM)');
box off;   

%% Figure 6 

%  %Plotting Isometric bars. dynamic Cai. 
% 
% filename = 'Isometric_D_SL';
% slValue = {'2.093', '2.3'};
% length_afterloadValue = size(slValue);
% len_afterloads = length_afterloadValue(2);
% line_color = [0.75, 0.75, 0.75];
% 
% for i=1:len_afterloads
%     SLv = slValue{i};
%     data = strcat(filename,SLv,'.csv')  
%
%     [time, SL_norm, F_total_norm, Ca_i, dTropTot] = reading_Isometric_DATA(data);
%     
%     if strcmp(SLv, '2.093') == 1
%         lw = 1.5;
%         line_type = '--';
%     elseif strcmp(SLv, '2.3') == 1
%         lw = 2.5;
%         line_type = '-';
%     end
%     
%     figure(6)
%     
%     subplot(3,1,1)
%     hold on
%     h1 = plot(time, Ca_i, line_type, 'color', line_color,'LineWidth',lw); hold on;
%     axis([0 500 0 1.5])
%     box off;
% 
%     subplot(3,1,2)
%     plot(time, F_total_norm, line_type, 'color', line_color, 'LineWidth',lw); hold on; 
%     axis([0 500 0 1.1])
%     box off;
%     
%     subplot(3,1,3)
%     hold on
%     plot(time, dTropTot, line_type, 'color', line_color,'LineWidth',lw); hold on;
%     axis([0 500 -0.5 1]) 
%     box off;
% 
% 
% end

filename = 'WL_D_afterload';
afterloadValue = {'0.2576'};
length_afterloadValue = size(afterloadValue);
len_afterloads = length_afterloadValue(2);
line_color = [0, 0, 0];
lw = 1.32;

    for i=1:len_afterloads
        afterload = afterloadValue{i};
        data = strcat(filename,afterload,'.csv');
        
        [time, SL_norm, F_total_norm, Ca_i, dTropTot] = reading_WL_DATA(data);

        hold on
        figure(6)

        subplot(3,1,1)
        hold on
        plot(time, Ca_i, 'color', line_color,'LineWidth',lw); hold on;
        set(gca,'fontsize',14);
        axis([0 500 0 1.5])
        ylabel('[Ca^2^+]_i (\muM)')
        box off;

        subplot(3,1,2)
        plot(time, F_total_norm, 'color', line_color, 'LineWidth',lw); hold on; 
        set(gca,'fontsize',14);
        axis([0 500 0 1.1])
        ylabel('Normalised Total Force')
        box off;

        subplot(3,1,3)
        hold on
        h3 = plot(time, dTropTot, 'color', line_color,'LineWidth',lw); hold on;
        set(gca,'fontsize',14);
        axis([0 500 -0.5 1])
        xlabel('Time (ms)'); % 
        ylabel('Ca-Trop Flux (\muM/ms)') %
        box off;
    end
    
legend([h2, h1, h3],{'Isometric L_o', 'Isometric 0.91 L_o', 'Work-loop'}, 'Location', 'northeast', 'Fontsize', 12);
