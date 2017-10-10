function [time, SL_norm, F_total_norm, Ca_i, dTropTot, ESmarker] = reading_WL_DATA(input_data)
    time= xlsread(input_data,strcat('A:A'))-19000;%time
    SL_norm= xlsread(input_data,strcat('B:B'))/2.3;%SL
    F_total_norm= xlsread(input_data,strcat('D:D'))./0.556;
    Ca_i = xlsread(input_data,strcat('E:E'));
    dTropTot = xlsread(input_data,strcat('AA:AA'));
    ESmarker = xlsread(input_data,strcat('AB:AB'));
end