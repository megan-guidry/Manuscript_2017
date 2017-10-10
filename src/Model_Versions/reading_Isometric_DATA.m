function [time, SL_norm, F_total_norm, Ca_i, dTropTot] = reading_Isometric_DATA(input_data)
    time= xlsread(input_data,strcat('A:A'))-19000;%time
    SL_norm= xlsread(input_data,strcat('B:B'))/2.3;%SL
    F_total_norm= xlsread(input_data,strcat('F:F'))./0.556;
    Ca_i = xlsread(input_data,strcat('D:D'));
    dTropTot = xlsread(input_data,strcat('W:W'));
end