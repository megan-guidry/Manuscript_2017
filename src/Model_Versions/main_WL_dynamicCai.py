#the original executed file is 'Workloop_protocol_with_passive_cycling2.py'
import sys
#preset_afterloads = [0.135, 0.1911, 0.2576, 0.3359, 0.4312, 0.6, 0.9]
preset_afterloads = [0.2576, 0.3359, 0.4312, 0.6]
print range(len(preset_afterloads))
count = 0

while count < len(preset_afterloads):
    
    value_afterload = preset_afterloads[count-1]
    user_number = count + 1
    execfile(r'Workloop_protocol_dynamicCai.py')
    count += 1
    print count

#from mlabwrap import mlab
#mlab.plotting_figures('wrkloop.csv')

sys.exit()
#Grab the csv workloops generated above (they will be in the same folder
#as main.py and Workloop_protocol_with_passive_cycling2.py) and plot
     
