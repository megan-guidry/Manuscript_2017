#the original executed file is 'Workloop_protocol_with_passive_cycling2.py'
import sys
#preset_afterloads = [0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
preset_afterloads = [0.2576]
print range(len(preset_afterloads))
count = 0

value_tstart = 2.36169272
value_Tau1 = 31.40873567
value_Tau2 = 40.1918327
value_Tau3 = 79.1466021
value_Ca_dia = 0.08981397
value_a1 = 0.23605426
value_b = 37.19706502
value_k1 = 1.76341758
value_c = 74.94501223


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
     
