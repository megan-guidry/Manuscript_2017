%% Figure 4: Two disting ES curves AND Figure 5: WL vs. Isometric Ca2+

%Plotting Isometric Bars
filename = 'Isometric_dynamicCai_SL';
slValue = {'1.955', '2.024', '2.093', '2.162', '2.231', '2.3'};
%slValue = {'2.3'};
length_afterloadValue = size(slValue);
len_afterloads = length_afterloadValue(2)
line_color = [0.75, 0.75, 0.75];

start = '1';
stop = '28572';

for i=1:len_afterloads
    SLv = slValue{i};
    data = strcat(filename,SLv,'.csv')  
    
    time= xlsread(data,strcat('A',start,':','A',stop))-19000;%time
    SL_norm= xlsread(data,strcat('B',start,':','B',stop))/2.3;%SL
    F_total_norm= xlsread(data,strcat('F',start,':','F',stop))./0.556;
    Ca_i = xlsread(data,strcat('D',start,':','D',stop));
    
    lw = 1.25;
    
    hold on
    figure(4)
    h1 = plot(SL_norm, F_total_norm, 'color', line_color,'LineWidth',lw); hold on;
    ax = gca
    ax.XTick = [0.85 0.9 0.95 1]
    axis([0.84 1.01 0 0.8])
    xlabel('Normalised Sarcomere Length');
    ylabel('Normalised Total Force')
    set(gca,'fontsize',14);
    box off;
    
    figure(5)
    hold on
    h5 = plot(time, Ca_i, 'color', line_color,'LineWidth',lw); hold on;
    axis([0 300 0 1.4])
    xlabel('Time (ms)');
    ylabel('[Ca^2^+]_i (\muM)')
    set(gca,'fontsize',14);
    box off;
end

%Plotting Work-loops
filename = 'WL_dynamicCai_afterload';
%afterloadValue = {'0.135'};
afterloadValue = {'0.135', '0.1911', '0.2576', '0.3359', '0.4312', '0.6'};
length_afterloadValue = size(afterloadValue);
len_afterloads = length_afterloadValue(2);
line_color = [0, 0, 0];


    for i=1:len_afterloads
        afterload = afterloadValue{i};
        data = strcat(filename,afterload,'.csv');

        start = '2';
        stop = '28572';
        
        time= xlsread(data,strcat('A',start,':','A',stop))-19000;%time
        SL_norm= xlsread(data,strcat('B',start,':','B',stop))/2.3;%SL
        F_total_norm= xlsread(data,strcat('D',start,':','D',stop))./0.556;
        Ca_i = xlsread(data,strcat('E',start,':','E',stop));


        lw = 1.25;

        hold on
        figure(4)
        h2 = plot(SL_norm, F_total_norm, 'k-','LineWidth',lw); hold on;
        ax = gca
        ax.XTick = [0.85 0.9 0.95 1]
        axis([0.84 1.01 0 1.2])
        xlabel('Normalised Sarcomere Length');
        ylabel('Normalised Total Force')
        set(gca,'fontsize',14);
        box off;
        
        figure(5)
        hold on
        h6 = plot(time, Ca_i, 'color', line_color,'LineWidth',lw); hold on;
        axis([0 300 0 1.4])
        xlabel('Time (ms)');
        ylabel('[Ca^2^+]_i (\muM)')
        set(gca,'fontsize',14);
        %set(gca,'xtick',[],'ytick',[])
        box off;
    end
    
%Plotting ES curves
SL_isometric = [1, 0.97, 0.94, 0.91, 0.88, 0.85];
Ftotal_isometric = [0.556, 0.4312, 0.3359, 0.2576, 0.1911, 0.135]./0.556

SL_WL_dynamicCai = [1, 0.9899, 0.9717, 0.9457, 0.9094, 0.8598];
Ftotal_WL_dynamicCai = [0.5561, 0.4305, 0.3351, 0.2567, 0.1901, 0.1339]./0.556;

Iso_line_color = [0.75, 0.75, 0.75];
WL_line_color = [0, 0, 0];


figure(4); hold on 
h3 = plot(SL_isometric, Ftotal_isometric, '-o', 'color', Iso_line_color, 'markers', 10, 'LineWidth', 2); hold on; %NO MARKERS ISO
h4 = plot(SL_WL_dynamicCai, Ftotal_WL_dynamicCai, '-o', 'color', WL_line_color, 'markers', 10, 'LineWidth', 2);
legend([h3, h4],{'Isometric', 'Work-loop'}, 'Location', 'northwest');

figure(5); hold on
legend([h5, h6],{'Isometric', 'Work-loop'}, 'Location', 'northoutside');


%% Figure 6 

 %Plotting Isometric bars. dynamic Cai. 

filename = 'Isometric_dynamicCai_SL';
slValue = {'2.093'};
length_afterloadValue = size(slValue);
len_afterloads = length_afterloadValue(2)
line_color = [0.75, 0.75, 0.75];

start = '1';
stop = '28572';

for i=1:len_afterloads
    SLv = slValue{i};
    data = strcat(filename,SLv,'_Tropreg.csv')  
  
    time= xlsread(data,strcat('A',start,':','A',stop))-19000;%time
    F_total_norm= xlsread(data,strcat('F',start,':','F',stop))./0.556;
    dTropTot = xlsread(data,strcat('W',start,':','W',stop));
    Ca_i = xlsread(data,strcat('D',start,':','D',stop));
    
    lw = 1.75;
    
    hold on
    figure(6)
    
    subplot(3,1,1)
    hold on
    h1 = plot(time, Ca_i, '--', 'color', line_color,'LineWidth',lw); hold on;
    axis([0 500 0 1.5])
    %xlabel('Time (ms)');
    ylabel('[Ca^2^+]_i (\muM)')
    set(gca,'fontsize',14);
    box off;

    subplot(3,1,2)
    plot(time, F_total_norm, '--', 'color', line_color, 'LineWidth',lw); hold on; 
    axis([0 500 0 1.1])
    ylabel('Normalised Total Force')
    %xlabel('Time (ms)')
    set(gca,'fontsize',14);
    box off;
    
    subplot(3,1,3)
    hold on
    plot(time, dTropTot, '--', 'color', line_color,'LineWidth',lw); hold on;
    axis([0 500 -0.5 1])
    xlabel('Time (ms)');
    ylabel('dTropTot (\muM/ms)')
    set(gca,'fontsize',14);
    box off;


end

 %Plotting Isometric bars. dynamic Cai. 

filename = 'Isometric_dynamicCai_SL';
slValue = {'2.3'};
length_afterloadValue = size(slValue);
len_afterloads = length_afterloadValue(2)
line_color = [0.75, 0.75, 0.75];

start = '1';
stop = '28572';

for i=1:len_afterloads
    SLv = slValue{i};
    data = strcat(filename,SLv,'_Tropreg.csv')  
  
    time= xlsread(data,strcat('A',start,':','A',stop))-19000;%time
    F_total_norm= xlsread(data,strcat('F',start,':','F',stop))./0.556;
    dTropTot = xlsread(data,strcat('W',start,':','W',stop));
    Ca_i = xlsread(data,strcat('D',start,':','D',stop));
    
    lw = 2.5;
    
    hold on
    figure(6)
    
    subplot(3,1,1)
    hold on
    plot(time, Ca_i, 'color', line_color,'LineWidth',lw); hold on;
    axis([0 500 0 1.5])
    %xlabel('Time (ms)');
    ylabel('[Ca^2^+]_i (\muM)')
    set(gca,'fontsize',14);
    box off;

    subplot(3,1,2)
    h2 = plot(time, F_total_norm, 'color', line_color, 'LineWidth',lw); hold on; 
    axis([0 500 0 1.1])
    ylabel('Normalised Total Force')
    %xlabel('Time (ms)')
    set(gca,'fontsize',14);
    box off;
    
    subplot(3,1,3)
    hold on
    plot(time, dTropTot, 'color', line_color,'LineWidth',lw); hold on;
    axis([0 500 -0.5 1])
    xlabel('Time (ms)');
    ylabel('dTropTot (\muM/ms)')
    set(gca,'fontsize',14);
    box off;

end

filename = 'WL_dynamicCai_afterload';
afterloadValue = {'0.2576'};
length_afterloadValue = size(afterloadValue);
len_afterloads = length_afterloadValue(2);
line_color = [0, 0, 0];


    for i=1:len_afterloads
        afterload = afterloadValue{i};
        data = strcat(filename,afterload,'_Tropreg.csv');

        start = '2';
        stop = '28572';

        time= xlsread(data,strcat('A',start,':','A',stop))-19000;%time
        F_total_norm= xlsread(data,strcat('D',start,':','D',stop))./0.556;
        dTropTot = xlsread(data,strcat('AA',start,':','AA',stop));
        Ca_i = xlsread(data,strcat('E',start,':','E',stop));

        lw = 1.25;

        hold on
        figure(6)

        subplot(3,1,1)
        hold on
        plot(time, Ca_i, 'color', line_color,'LineWidth',lw); hold on;
        axis([0 500 0 1.5])
        %xlabel('Time (ms)');
        ylabel('[Ca^2^+]_i (\muM)')
        set(gca,'fontsize',14);
        box off;

        subplot(3,1,2)
        plot(time, F_total_norm, 'color', line_color, 'LineWidth',lw); hold on; 
        axis([0 500 0 1.1])
        ylabel('Normalised Total Force')
        %xlabel('Time (ms)')
        set(gca,'fontsize',14);
        box off;

        subplot(3,1,3)
        hold on
        h3 = plot(time, dTropTot, 'color', line_color,'LineWidth',lw); hold on;
        axis([0 500 -0.5 1])
        xlabel('Time (ms)');
        ylabel('Ca-Trop Flux (\muM/ms)')
        set(gca,'fontsize',14);
        box off;
    end
    
legend([h2, h1, h3],{'Isometric L_o', 'Isometric 0.91 L_o', 'Work-loop'}, 'Location', 'northeast', 'Fontsize', 12);


%% Figure 7: ES shift caused by Wider Ca2+ transients

SL_WL_fixedCaiSL23 = [1, 0.9902, 0.9724, 0.947, 0.9119, 0.8647];
Ftotal_WL_fixedCaiSL23 = [0.5542, 0.4304, 0.3351, 0.2567, 0.1901, 0.1339]./0.556;

SL_WL_fixedCaiSL2231 = [1, 0.9852, 0.9652, 0.9376, 0.9003, 0.8506];
Ftotal_WL_fixedCaiSL2231 = [0.5984, 0.4305, 0.3351, 0.2567, 0.1901, 0.1339]./0.556;

SL_WL_fixedCaiSL2162 = [1, 0.9985, 0.98, 0.958, 0.9284, 0.889, 0.8371];
Ftotal_WL_fixedCaiSL2162 = [0.641, 0.599, 0.4304, 0.3351, 0.2567, 0.1901, 0.1339]./0.556;
 
SL_WL_fixedCaiSL2093 = [1, 0.9962, 0.9749, 0.9511, 0.9198, 0.8785, 0.8248];
Ftotal_WL_fixedCaiSL2093 = [0.6802, 0.5996, 0.4304, 0.335, 0.2567, 0.1901, 0.1338]./0.556;

SL_WL_fixedCaiSL2024 = [1, 0.9938, 0.9703, 0.9449, 0.9121, 0.8692, 0.8141];
Ftotal_WL_fixedCaiSL2024 = [0.7143, 0.5994, 0.4304, 0.335, 0.2566, 0.19, 0.1338]./0.556;

SL_WL_fixedCaiSL1955 = [1, 0.9915, 0.9662, 0.9395, 0.9054, 0.8613, 0.8051];
Ftotal_WL_fixedCaiSL1955 = [0.7431, 0.5993, 0.4304, 0.335, 0.2566, 0.19, 0.1338]./0.556;

line_color = [0, 0, 0];

figure(7)
plot(SL_WL_fixedCaiSL23, Ftotal_WL_fixedCaiSL23, '-o', 'color', line_color, 'markers', 2, 'LineWidth', 1.25); hold on;
plot(SL_WL_fixedCaiSL2231, Ftotal_WL_fixedCaiSL2231, '-o', 'color', line_color, 'markers', 2, 'LineWidth', 1.25); hold on;
plot(SL_WL_fixedCaiSL2162, Ftotal_WL_fixedCaiSL2162, '-o', 'color', line_color, 'markers', 2, 'LineWidth', 1.25); hold on;
plot(SL_WL_fixedCaiSL2093, Ftotal_WL_fixedCaiSL2093, '-o', 'color', line_color, 'markers', 2, 'LineWidth', 1.25); hold on;
plot(SL_WL_fixedCaiSL2024, Ftotal_WL_fixedCaiSL2024, '-o', 'color', line_color, 'markers', 2, 'LineWidth', 1.25); hold on;
plot(SL_WL_fixedCaiSL1955, Ftotal_WL_fixedCaiSL1955, '-o', 'color', line_color, 'markers', 2, 'LineWidth', 1.25); hold on;

set(gca,'fontsize',14);
xlabel('Normalised Sarcomere Length');
ylabel('Normalised Total Force')
ax = gca
box off;
ax.XTick = [0.8 0.85 0.9 0.95 1]
axis([0.79 1.01 0 1.5])

%% Figure 8: Uniting the Isometric and WL ES curves

filename = 'ES_merge_afterload';
afterloadValue = {'0.135', '0.1911', '0.2576', '0.3359', '0.4312', '0.6'};

%afterloadValue = {'0.135', '0.2576', '0.4312'};
length_afterloadValue = size(afterloadValue);
len_afterloads = length_afterloadValue(2);
line_color = [0, 0, 0];
white_out = [1, 1, 1];

 for i=1:len_afterloads
        afterload = afterloadValue{i};
        %CaiSL = CaiSL_value{i};
        data = strcat(filename,afterload,'.csv');

        start = '114288';
        stop = '142859';
        
        time= xlsread(data,strcat('A',start,':','A',stop))-4000;%time
        SL_norm= xlsread(data,strcat('B',start,':','B',stop))/2.3;%SL
        F_total_norm= xlsread(data,strcat('D',start,':','D',stop))./0.556;
        XBpostR = xlsread(data,strcat('L',start,':','L',stop));
        XBpreR = xlsread(data,strcat('P',start,':','P',stop));
        Force_producing_state = XBpostR + XBpreR;

        lw = 1.25;

        hold on
        figure(8)
        h7 = plot(SL_norm, F_total_norm, '-', 'color', line_color, 'LineWidth',lw); hold on;
        %h7 = scatter(SL_norm, F_total_norm, Force_producing_state*10, Force_producing_state); hold on;
        ax = gca
        ax.XTick = [0.85 0.9 0.95 1]
        axis([0.84 1.01 0 1.2])
        %axis([0.83 1.01 -0.01 65 ])
        set(gca,'fontsize',14);
        xlabel('Normalised Sarcomere Length');
        ylabel('Normalised Total Force')
        set(gca,'fontsize',14);
        box off;
        
 end
dummyX = [0.5, 0.6];
dummyY = [0.9, 1];
h10 = plot(dummyX, dummyY, '-o', 'color', line_color, 'markers', 10, 'LineWidth',lw); hold on;

SL_isometric = [1, 0.97, 0.94, 0.91, 0.88, 0.85];
Ftotal_isometric = [0.556, 0.4312, 0.3359, 0.2576, 0.1911, 0.135]./0.556

WL_with_fixed_Cai_EScurve_X = [1, 0.9703, 0.9395, 0.9121, 0.8785, 0.8506]
WL_with_fixed_Cai_EScurve_Y = [0.5542, 0.4304, 0.335, 0.2566, 0.1901, 0.1339]./0.556

Iso_line_color = [0.75, 0.75, 0.75];

hold on
h8 = plot(SL_isometric, Ftotal_isometric, '-', 'color', Iso_line_color, 'markers', 10, 'LineWidth', 2); hold on; %NO MARKERS ISO
h9 =plot(WL_with_fixed_Cai_EScurve_X, WL_with_fixed_Cai_EScurve_Y, 'o', 'color', line_color, 'markers', 10, 'LineWidth', 2); hold on; %NO MARKERS ISO

legend([h10, h8],{'Work-loop', 'Isometric'}, 'Location', 'northwest');
