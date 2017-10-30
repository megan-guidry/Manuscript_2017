#the original executed file is 'Workloop_protocol_with_passive_cycling2.py'
import sys
preset_SL = [2.024, 2.093, 2.162, 2.231]
print range(len(preset_SL))
count = 0

while count < len(preset_SL):
    
    value_SL = preset_SL[count-1]
    user_number = count + 1
    execfile(r'Isometric_protocol_dynamicCai.py')
    count += 1
    print count

#from mlabwrap import mlab
#mlab.plotting_figures('wrkloop.csv')

sys.exit()
#Grab the csv workloops generated above (they will be in the same folder
#as main.py and Workloop_protocol_with_passive_cycling2.py) and plot
     
