#the original executed file is 'Workloop_protocol_with_passive_cycling2.py'
import sys
#preset_afterloads = [0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
preset_afterloads = [0.6]
print range(len(preset_afterloads))
count = 0

value_tstart = 2.93884465
value_Tau1 = 41.4155557
value_Tau2 = 78.01759413
value_Tau3 = 52.84891992
value_Ca_dia = 0.09014398
value_a1 = 0.29336819
value_b = 39.99993005
value_k1 = 1.63406024
value_c = 47.32170775


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
     
