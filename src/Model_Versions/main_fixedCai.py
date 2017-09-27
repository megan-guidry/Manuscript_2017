#the original executed file is 'Workloop_protocol_with_passive_cycling2.py'
import sys
#preset_afterloads = [0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
preset_afterloads = [0.2576]
print range(len(preset_afterloads))
count = 0

value_tstart = 2.80900835
value_Tau1 = 36.03943529
value_Tau2 = 45.08829658
value_Tau3 = 72.90385783
value_Ca_dia = 0.08968451
value_a1 = 0.27926453
value_b = 37.66516463
value_k1 = 1.67887931
value_c = 82.71763373


while count < len(preset_afterloads):
    
    value_afterload = preset_afterloads[count-1]
    user_number = count + 1
    execfile(r'Workloop_protocol_fixedCai.py')

    #for the above file, choose Workloop_protocol_with_passive_cycling2.py for fixed calcium.
    count += 1
    print count

#from mlabwrap import mlab
#mlab.plotting_figures('wrkloop.csv')

sys.exit()
#Grab the csv workloops generated above (they will be in the same folder
#as main.py and Workloop_protocol_with_passive_cycling2.py) and plot
     
