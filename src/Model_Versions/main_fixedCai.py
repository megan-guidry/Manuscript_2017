#the original executed file is 'Workloop_protocol_with_passive_cycling2.py'
import sys
#preset_afterloads = [0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
preset_afterloads = [0.6]
print range(len(preset_afterloads))
count = 0

value_tstart = 2.21842617
value_Tau1 = 26.491962
value_Tau2 = 28.7470311
value_Tau3 = 108.854814
value_Ca_dia = 0.0894001
value_a1 = 0.21621596
value_b = 33.1456627
value_k1 = 1.83302661
value_c = 63.8698117


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
     
