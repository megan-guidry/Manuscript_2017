import sys
import FixedCai_parameters

#preset_afterloads = [0.135, 0.1911, 0.2576, 0.3359, 0.4312, 0.6, 0.9]
preset_CaiSL = [2.024]
preset_afterloads = [0.135, 0.1911]
print range(len(preset_afterloads))
CaiSL_count = 0

while CaiSL_count < len(preset_CaiSL):
    Afterload_count = 0
    value_CaiSL = preset_CaiSL[CaiSL_count]

    if value_CaiSL == 1.955:
        parameters = [0.089002126, 0.283359682, 39.423358, 2.84495395, 1.6649028, 37.6686509, 9.995708, 49.5660573, 100.759014, 50, 70.3621179]
    elif value_CaiSL == 2.024:
        parameters = [0.08968451, 0.27926453, 37.66516463, 2.80900835, 1.67887931, 36.03943529, 5.42060899, 45.08829658, 82.71763373, 50, 72.90385783]
    elif value_CaiSL == 2.093:
        parameters = [0.08981397, 0.23605426, 37.19706502, 2.36169272, 1.76341758, 31.40873567, 9.99030452, 40.1918327, 74.94501223, 50 ,79.1466021]
    elif value_CaiSL == 2.162:
        parameters = [0.08982691, 0.23170048, 35.4328648, 2.32764163, 1.78087606, 29.92922172, 9.99824712, 36.10844273, 70.54464675, 50, 87.58238235]
    elif value_CaiSL == 2.231:
        parameters = [0.089699541, 0.225542405, 34.1045897, 2.28321127, 1.80277288, 28.3291654, 0.016845165, 32.1908114, 66.8007906, 50, 97.3577797]
    elif value_CaiSL == 2.3:
        parameters = [0.0894001, 0.21621596, 33.1456627, 2.21842617, 1.83302661, 26.491962, 9.97086786, 28.7470311, 63.8698117, 50, 108.854814]

    value_tstart = parameters[3]
    value_Tau1 = parameters[5]
    value_Tau2 = parameters[7]
    value_Tau3 = parameters[10]
    value_Ca_dia = parameters[0]
    value_a1 = parameters[1]
    value_b = parameters[2]
    value_k1 = parameters[4]
    value_c = parameters[8]
    
    CaiSL_count += 1
    
    while Afterload_count < len(preset_afterloads):
        
        value_afterload = preset_afterloads[Afterload_count]
        user_number = Afterload_count + 1
        execfile(r'Workloop_protocol_fixedCai.py')

        Afterload_count += 1
        print value_CaiSL
        print value_afterload

#from mlabwrap import mlab
#mlab.plotting_figures('wrkloop.csv')

sys.exit()
#Grab the csv workloops generated above (they will be in the same folder
#as main.py and Workloop_protocol_with_passive_cycling2.py) and plot
     
