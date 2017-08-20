#!/usr/bin/env python

""" Date: 13-7-2015.  This is the first attempt at getting the Hinch-Rice-Tran
model working in CellML/Python/OpenCMISS"""

#DOC-START imports
from __future__ import print_function
import sys, os, math
# Make sure $OPENCMISS_ROOT/cm/bindings/python is first in our PYTHONPATH.
sys.path.insert(1, os.path.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))
import main_dynamicCai

periodTime = 0
switch = 0
it = 1
WL_iteration = []

SL_top = []
active_tension_top = []
active_tension_WL = [] # WL stands for "whole loop"
F_total_WL = []
SL_bottom = []
SL_WL = [] # WL stands for "whole loop"
active_tension_bottom = []
step = 0
yes = 0
IT = []
iteration_value= 20

def setting_t_interval():
    global iterationTime
    change_currentTime_to_iterationTime = math.modf(currentTime/1000)
    iterationTime = change_currentTime_to_iterationTime[0]*1000

# get a specific loop from the multiple work-loop iterations that the python code generates:
def get_loop_info(): # SL used to be SL_T
    global it
    global loop_number
    global SEon
    global SL
    global T
    global SL_top
    global SL_bottom
    global SL_WL
    
    if it == iteration_value and SEon[-1] == 1:
        SL_top.append(SL) # this makes the top curve of the workloop (SL[-1] was originally SL_T)
        active_tension_top.append(T) # this makes the top curve of the workloop

    if it == iteration_value and SEon[-1] == 2:
        SL_bottom.append(SL_T) # this makes the bottom curve of the workloop
        active_tension_bottom.append(T) # this makes the bottom curve of the workloop

    if it == iteration_value:
        #print(SL_T)
        SL_WL.append(SL_T)
        F_total_WL.append(F)

# If Force is greater than afterload, switch SE on to 1 (isotonic)
def force_greater_than_afterload():
    global step
    if F_total[-1] > (value_afterload): # if total force is greater than the afterload, then switch to isotonic state (shortening)
        if cellmlNodeThisComputationalNode:
            cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 1)
            step = 1         
            
# If the Sarcomere length starts increasing then SEon will switch to 0 (isometric)
##def SL_increasing():
##    global yes
##    global YES
##    def SEon_Zero():
##        yes = 0
##        for element in [(int(SEon[-2]), int(SEon[-3]), int(SEon[-4]))]:
##            if element == int(0):
##                yes = 1
##                
##    global step
##    global element
##    global SEon
##    
##    SEon_Zero()
##    if len(SL) > 5: # This is just to make sure there are 10 elements in the list SL (10 is a random number)
##        if int(SEon[-1]) == 1 and yes == 1: #yes was originally SEon_Zero
##            if cellmlNodeThisComputationalNode:
##                cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 1)
##                SL_increasing = 1
##                return SL_increasing
##        else:
##            if SL[-1] > SL[-2]: # if the last value in SL ([-1]) is greater than the 2nd to last value in SL we know SL is increasing, so change to isometric (non shortening)
##                
##                if cellmlNodeThisComputationalNode:
##                    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 0)
##                    step = 2

def SL_increasing():
    global step
    if len(SL) > 5 and SL[-1] > SL[-2] > SL[-3]:
        if cellmlNodeThisComputationalNode:
                    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 0)
                    step = 2

##def SL_increasing():
##    global step
##    if time[-1] >= 251.7:
##        if cellmlNodeThisComputationalNode:
##                    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 0)
##                    step = 2

#If calcium has returned to diastolic levels then the sarcomere can return to its resting length (linearly)
def sarcomere_lengthen():
    global step
    if Ca_i[-1] <= 0.1 and iterationTime >= 700: # check just makes it so that the previous if statement must be achieved before this one can be achieved
        # this if statement is --> if the last value of intracellular calcium calculated is less than or equal to the initial value of calcium then allow the sarcomere to return to its starting length (SLset)
        if cellmlNodeThisComputationalNode:
            cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 2)
            step = 3

    
# Switch to isometric mode when sarcomere reaches its resting length
def stop_SL_recovery():
    global it
    global step
    global SL_T
    global SL
    
    if SL[-1] > value_SLset: # make sure the sarcomere stops relengthening when the length reaches SLset
        if cellmlNodeThisComputationalNode:
            cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, 0)
            periodTime = 0
            step = 0            
        if SL[-2] < value_SLset: # note this if statement is embedded in the previous if statement. This statement can be read as: If the current sarcomere length is greater than SLset but the previous value of
            # SL is less than SL rest then increase 'it' by 1.  This is the way I go about keeping track of each work-loop iteration. it will only increase by 1 once at the end of each work-loop.
            
            it += 1 #it correlates to the work-loop iteration (if it == 5 then its the 5th work-loop iteration)
            # this is done so that we can look at a work loop iteration that has achieved a steady state
            WL_iteration.append(it)

            


# Functions for the Workloop code##################################################################################
###################################################################################################################


# Intialise OpenCMISS
from opencmiss import CMISS
#DOC-END imports

# Set problem parameters
#DOC-START parameters

# get the model file
if len(sys.argv) > 1:
    modelFile = sys.argv[1]
else:
    modelFile = "../Hinch_Model3_quickrelease.xml"  

# 1D domain size
width = 1.0
numberOfXElements = 2


# Materials parameters
Am = 193.6
Cm = 0.014651
conductivity = 0.0


# Simulation parameters
StimPeriod = 1000
switch = 0
value_SEon = 0
value_mass = 50 #0.00666 #0.005

value_SLset = 2.3 #original value is 2.27
value_afterload = main_dynamicCai.value_afterload
timeStop = 20000 #ms
value_TmpC =  23 #Celsius
value_x_0 = 0.007 #micrometre
value_stimTime = 333 #ms
value_stepStart = 400 #ms
value_finalLambda = 1 #dimensionless
# this is used in the integration of the CellML model
odeTimeStep = 0.01 # [ms?] 0.01 originally
# this is used in the dummy monodomain problem/solver
pdeTimeStep = 0.01 # [ms]
# this is the step at which we grab output from the solver
outputTimeStep = 0.035 #WAS 0.1
# set this to 1 to get exfiles written out during solve
outputFrequency = 0
#DOC-END parameters




#Setup field number handles --> This is for setup... no change necessary
coordinateSystemUserNumber = main_dynamicCai.user_number
regionUserNumber = main_dynamicCai.user_number
basisUserNumber = main_dynamicCai.user_number
pressureBasisUserNumber = main_dynamicCai.user_number + 1
generatedMeshUserNumber = main_dynamicCai.user_number
meshUserNumber = main_dynamicCai.user_number
cellMLUserNumber = main_dynamicCai.user_number
decompositionUserNumber = main_dynamicCai.user_number
equationsSetUserNumber = main_dynamicCai.user_number
problemUserNumber = main_dynamicCai.user_number
#Mesh component numbers
linearMeshComponentNumber = 1
#Fields
geometricFieldUserNumber = main_dynamicCai.user_number
fibreFieldUserNumber = main_dynamicCai.user_number + 1
dependentFieldUserNumber = main_dynamicCai.user_number + 2
materialsFieldUserNumber = main_dynamicCai.user_number + 3
equationsSetFieldUserNumber = main_dynamicCai.user_number + 4
cellMLModelsFieldUserNumber = main_dynamicCai.user_number + 5
cellMLStateFieldUserNumber = main_dynamicCai.user_number + 6
cellMLParametersFieldUserNumber = main_dynamicCai.user_number + 7
cellMLIntermediateFieldUserNumber = main_dynamicCai.user_number + 8

#DOC-START parallel information
# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()
#DOC-END parallel information

#DOC-START initialisation
# Create a 2D rectangular cartesian coordinate system
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(1)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = CMISS.Region()
region.CreateStart(regionUserNumber,CMISS.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()
#DOC-END initialisation

#DOC-START basis
# Define a bilinear Lagrange basis
basis = CMISS.Basis()
basis.CreateStart(basisUserNumber)
basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 1
basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
basis.CreateFinish()
#DOC-END basis

#DOC-START generated mesh
# Create a generated mesh
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = CMISS.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width]
generatedMesh.numberOfElements = [numberOfXElements]

mesh = CMISS.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)
#DOC-END generated mesh

#DOC-START decomposition
# Create a decomposition for the mesh
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = CMISS.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()
#DOC-END decomposition

#DOC-START geometry
# Create a field for the geometry
geometricField = CMISS.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(CMISS.FieldVariableTypes.U, "coordinates")
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 1, linearMeshComponentNumber)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)
#DOC-END geometry

#DOC-START equations set
# Create the equations_set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
        CMISS.EquationsSetClasses.BIOELECTRICS,
        CMISS.EquationsSetTypes.MONODOMAIN_EQUATION,
        CMISS.EquationsSetSubtypes.NONE,
        equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()
#DOC-END equations set

#DOC-START equations set fields
# Create the dependent Field
dependentField = CMISS.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

# Create the materials Field
materialsField = CMISS.Field()
equationsSet.MaterialsCreateStart(materialsFieldUserNumber, materialsField)
equationsSet.MaterialsCreateFinish()

# Set the materials values
# Set Am
materialsField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,Am)
# Set Cm
materialsField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2,Cm)
# Set conductivity
materialsField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,3,conductivity)
#DOC-END equations set fields

#DOC-START create cellml environment
# Create the CellML environment
cellML = CMISS.CellML()
cellML.CreateStart(cellMLUserNumber, region)
# Import the Hinch-Rice-Tran cell model from a file
MeganModel = cellML.ModelImport(modelFile)
#DOC-END create cellml environment


#Flagging variables for MeganModel: known= variables we want to be able to set//
#cellML.VariableSetAsKnown(MeganModel, "parameters/stimTime")
cellML.VariableSetAsKnown(MeganModel, "parameters/TmpC")
cellML.VariableSetAsKnown(MeganModel, "parameters/stepStart")
cellML.VariableSetAsKnown(MeganModel, "parameters/finalLambda")
cellML.VariableSetAsKnown(MeganModel, "parameters/SLset")
cellML.VariableSetAsKnown(MeganModel, "parameters/afterload")
cellML.VariableSetAsKnown(MeganModel, "parameters/SEon")
cellML.VariableSetAsKnown(MeganModel, "parameters/mass")




#Below are the variables that we want to get from the CellML model 
cellML.VariableSetAsWanted(MeganModel, "parameters/active")
cellML.VariableSetAsWanted(MeganModel, "parameters/passive")
cellML.VariableSetAsWanted(MeganModel, "parameters/preload")
cellML.VariableSetAsWanted(MeganModel, "parameters/F_total")
cellML.VariableSetAsWanted(MeganModel, "parameters/dTropTot")
cellML.VariableSetAsWanted(MeganModel, "parameters/XB_cycling")

cellML.VariableSetAsWanted(MeganModel, "parameters/gxbT")
cellML.VariableSetAsWanted(MeganModel, "parameters/hfT")
cellML.VariableSetAsWanted(MeganModel, "parameters/SOVFThick")
cellML.VariableSetAsWanted(MeganModel, "parameters/SOVFThin")
cellML.VariableSetAsWanted(MeganModel, "parameters/fxbT")
cellML.VariableSetAsWanted(MeganModel, "parameters/hbT")
cellML.VariableSetAsWanted(MeganModel, "parameters/gappT")
cellML.VariableSetAsWanted(MeganModel, "parameters/fappT")


cellML.VariableSetAsWanted(MeganModel, "parameters/P")
cellML.VariableSetAsWanted(MeganModel, "parameters/kn_pT")
cellML.VariableSetAsWanted(MeganModel, "parameters/kp_nT")

cellML.VariableSetAsWanted(MeganModel, "parameters/I_LCC")
cellML.VariableSetAsWanted(MeganModel, "parameters/I_RyR")
cellML.VariableSetAsWanted(MeganModel, "parameters/I_SERCA")
cellML.VariableSetAsWanted(MeganModel, "parameters/I_SR")
cellML.VariableSetAsWanted(MeganModel, "parameters/I_NaCa")
cellML.VariableSetAsWanted(MeganModel, "parameters/I_pCa")
cellML.VariableSetAsWanted(MeganModel, "parameters/I_CaB")


#DOC-START create cellml finish
cellML.CreateFinish()
#DOC-END create cellml finish


#DOC-START map Vm components
# Start the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateStart()
#Now we can set up the field variable component <--> CellML model variable mappings.
#Map Vm (We only need to map Vm to make the monodomain happy)

cellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U,1, CMISS.FieldParameterSetTypes.VALUES,MeganModel,"parameters/Vm", CMISS.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(MeganModel,"parameters/Vm", CMISS.FieldParameterSetTypes.VALUES,dependentField,CMISS.FieldVariableTypes.U,1,CMISS.FieldParameterSetTypes.VALUES)

#Finish the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateFinish()

# Set the initial Vm values
dependentField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,-92.5)
#DOC-END map Vm components

#DOC-START define CellML models field
#Create the CellML models field
cellMLModelsField = CMISS.Field()
cellML.ModelsFieldCreateStart(cellMLModelsFieldUserNumber, cellMLModelsField)
cellML.ModelsFieldCreateFinish()
#DOC-END define CellML models field

#DOC-START define CellML state field
#Create the CellML state field 
cellMLStateField = CMISS.Field()
cellML.StateFieldCreateStart(cellMLStateFieldUserNumber, cellMLStateField)
cellML.StateFieldCreateFinish()
#DOC-END define CellML state field

#DOC-START define CellML parameters and intermediate fields
#Create the CellML parameters field 
cellMLParametersField = CMISS.Field()
cellML.ParametersFieldCreateStart(cellMLParametersFieldUserNumber, cellMLParametersField)
cellML.ParametersFieldCreateFinish()

#  Create the CellML intermediate field 
cellMLIntermediateField = CMISS.Field()
cellML.IntermediateFieldCreateStart(cellMLIntermediateFieldUserNumber, cellMLIntermediateField)
cellML.IntermediateFieldCreateFinish()
#DOC-END define CellML parameters and intermediate fields

# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
equations.outputType = CMISS.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Find the domains of the first and last nodes
firstNodeNumber = 1
lastNodeNumber = (numberOfXElements+1)
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber, 1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber, 1)

# set up the indices of the fields we want to grab
#stimTime_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/stimTime")
stepStart_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/stepStart")
finalLambda_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/finalLambda")
SLset_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/SLset")
afterload_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/afterload")
SEon_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/SEon")


mass_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/mass")

TmpC_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/TmpC")
active_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/active")
passive_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/passive")
preload_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/preload")
SL_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/SL")
integral_force_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/integral_force")
F_total_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/F_total")
dTropTot_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/dTropTot")
XB_cycling_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/XB_cycling")
#Vm_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/Vm")
Ca_i_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/Ca_i")
#dSL_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/dSL")

gxbT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/gxbT")
XBpostr_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/XBpostr")
hfT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/hfT")
SOVFThick_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/SOVFThick")
SOVFThin_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/SOVFThin")
xXBpostr_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/xXBpostr")
XBprer_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/XBprer")
xXBprer_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/xXBprer")
fxbT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/fxbT")
hbT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/hbT")
gappT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/gappT")
fappT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/fappT")

P_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/P")
N_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/N")
kn_pT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/kn_pT")
kp_nT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/kp_nT")

I_LCC_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_LCC")
I_RyR_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_RyR")
I_SERCA_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_SERCA")
I_SR_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_SR")
I_NaCa_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_NaCa")
I_pCa_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_pCa")
I_CaB_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/I_CaB")

# and the arrays to store them in

stepStart = []
finalLambda = []
stimTime = []
TmpC = []
Vm = []
active = [] 
Ca_i = []
dSL = []
SLset = []
afterload = []
integral_force = []
F_total = []
dTropTot = []
XB_cycling = []

gxbT = []
XBpostr = []
hfT = []
SOVFThick = []
SOVFThin = []
xXBpostr = []
XBprer = []
xXBprer = []

fxbT = []
hbT = []
gappT = []
fappT = []

SEon = []
SL = []
passive = []
preload = []

P = []
N = []
kn_pT = []
kp_nT = []

I_LCC = []
I_RyR = []
I_SERCA = []
I_SR = []
I_NaCa = []
I_pCa = []
I_CaB = []

# Do i need to put SLset and afterload here??

# We are using node 1 as the point in our dummy monodomain problem to integrate the CellML model
cellmlNode = 2
cellmlNodeDomain = decomposition.NodeDomainGet(cellmlNode, 1)
cellmlNodeThisComputationalNode = False
if cellmlNodeDomain == computationalNodeNumber:
    cellmlNodeThisComputationalNode = True

#if cellmlNodeThisComputationalNode:
    #cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stimTime_component, value_stimTime)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, TmpC_component, value_TmpC)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stepStart_component, value_stepStart)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, finalLambda_component, value_finalLambda)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SLset_component, value_SLset)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, afterload_component, value_afterload)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component, value_SEon)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, mass_component, value_mass)
    

#DOC-START define monodomain problem
#Define the problem
problem = CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.BIOELECTRICS,
    CMISS.ProblemTypes.MONODOMAIN_EQUATION,
    CMISS.ProblemSubTypes.MONODOMAIN_GUDUNOV_SPLIT)
problem.CreateFinish()
#DOC-END define monodomain problem



currentTime = 0.0
# grab initial results
#print "Time: " + str(currentTime)
time = []
value = 0.0
if cellmlNodeThisComputationalNode:
    time.append(currentTime)
    #stimTime.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stimTime_component))
    TmpC.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, TmpC_component))
    SEon.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component))
    stepStart.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stepStart_component))
    finalLambda.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, finalLambda_component))
    SLset.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SLset_component))
    afterload.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, afterload_component))
    active.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, active_component))
    passive.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, passive_component))
    preload.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, preload_component))
    SL.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SL_component))
    integral_force.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, integral_force_component))
    F_total.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, F_total_component))
    dTropTot.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, dTropTot_component))
    XB_cycling.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XB_cycling_component))

    gxbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, gxbT_component))
    XBpostr.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XBpostr_component))
    hfT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, hfT_component))
    SOVFThick.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SOVFThick_component))
    SOVFThin.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SOVFThin_component))
    xXBpostr.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, xXBpostr_component))
    XBprer.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XBprer_component))
    xXBprer.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, xXBprer_component))

    fxbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, fxbT_component))
    hbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, hbT_component))
    gappT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, gappT_component))
    fappT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, fappT_component))

    P.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, P_component))
    N.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, N_component))
    kn_pT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, kn_pT_component))
    kp_nT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, kp_nT_component))

    I_LCC.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_LCC_component))
    I_RyR.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_RyR_component))
    I_SERCA.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_SERCA_component))
    I_SR.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_SR_component))
    I_NaCa.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_NaCa_component))
    I_pCa.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_pCa_component))
    I_CaB.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_CaB_component))
    
    #Vm.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Vm_component))
    Ca_i.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Ca_i_component))
    #dSL.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, dSL_component))

#Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = CMISS.ControlLoop()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE],controlLoop)
if currentTime >= 0 and currentTime < 500:
    controlLoop.TimesSet(0.0,1.5,pdeTimeStep)
else:
    controlLoop.TimesSet(500,501.5,pdeTimeStep)

# controlLoop.OutputTypeSet(CMISS.ControlLoopOutputTypes.TIMING)
controlLoop.OutputTypeSet(CMISS.ControlLoopOutputTypes.NONE)

controlLoop.TimeOutputSet(outputFrequency)
problem.ControlLoopCreateFinish()

#Create the problem solvers
daeSolver = CMISS.Solver()
dynamicSolver = CMISS.Solver()
problem.SolversCreateStart()
# Get the first DAE solver
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,daeSolver)
daeSolver.DAETimeStepSet(odeTimeStep)
daeSolver.OutputTypeSet(CMISS.SolverOutputTypes.NONE)
# Get the second dynamic solver for the parabolic problem
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],2,dynamicSolver)
dynamicSolver.OutputTypeSet(CMISS.SolverOutputTypes.NONE)
problem.SolversCreateFinish()

#DOC-START define CellML solver
#Create the problem solver CellML equations
cellMLEquations = CMISS.CellMLEquations()
problem.CellMLEquationsCreateStart()
daeSolver.CellMLEquationsGet(cellMLEquations)
cellmlIndex = cellMLEquations.CellMLAdd(cellML)
problem.CellMLEquationsCreateFinish()
#DOC-END define CellML solver

#Create the problem solver PDE equations
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
dynamicSolver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe any boundary conditions 
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
solverEquations.BoundaryConditionsCreateFinish()

# Main time loop
# Here we are really doing something that Iron is not designed to handle, so this is not optimal
# but does allow up to test model integration in Iron.
F = F_total[0]
XBc = XB_cycling[0]
T = active[0]
SL_T = SL[0]
periodTime = 0
switch = 0
it = 1
WL_iteration = []

SL_top = []
active_tension_top = []
active_tension_WL = [] # WL stands for "whole loop"
SL_bottom = []
SL_WL = [] # WL stands for "whole loop"
active_tension_bottom = []
step = 0
IT = [1]
testt = [0]
YES = [0]

while currentTime < timeStop:
    #print(currentTime)
    setting_t_interval()
    
    if step <= 3:
        if step == 0:
            force_greater_than_afterload()
            # If Force is greater than afterload, switch SEon to 1 (isotonic)
        if step == 1:
            SL_increasing()
            # If the Sarcomere length starts increasing then SEon will switch to 0 (isometric)
            force_greater_than_afterload() #ADDED**
        if step == 2:
            sarcomere_lengthen()
            #If calcium has returned to diastolic levels then the sarcomere can return to its resting length (linearly)
            force_greater_than_afterload() #ADDED**
        if step == 3:
            stop_SL_recovery()
            

           

    # set the next solution interval
    nextTime = currentTime + outputTimeStep
    controlLoop.TimesSet(currentTime, nextTime, outputTimeStep)
    
    # integrate the model
    problem.Solve()

    currentTime = nextTime
    periodTime += 1
    
    # grab results
    #print "Time: " + str(currentTime)
    
    if cellmlNodeThisComputationalNode:
        time.append(currentTime)
        #stimTime.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stimTime_component))
        TmpC.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, TmpC_component))
        SEon.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SEon_component))
        stepStart.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stepStart_component))
        finalLambda.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, finalLambda_component))
        SLset.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SLset_component))
        afterload.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, afterload_component))
        T = cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, active_component)
        F = cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, F_total_component)
        
        XBc = cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XB_cycling_component)
        passive.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, passive_component))
        preload.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, preload_component))

        gxbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, gxbT_component))
        XBpostr.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XBpostr_component))
        hfT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, hfT_component))
        SOVFThick.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SOVFThick_component))
        SOVFThin.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SOVFThin_component))
        xXBpostr.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, xXBpostr_component))
        XBprer.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XBprer_component))
        xXBprer.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, xXBprer_component))

        fxbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, fxbT_component))
        hbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, hbT_component))
        gappT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, gappT_component))
        fappT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, fappT_component))

        P.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, P_component))
        N.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, N_component))
        kn_pT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, kn_pT_component))
        kp_nT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, kp_nT_component))
        
        active.append(T)
        SL_T = cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SL_component)
        SL.append(SL_T)
        IT.append(it)
        integral_force.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, integral_force_component))
        F_total.append(F)
        XB_cycling.append(XBc)
        #Vm.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Vm_component))
        Ca_i.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Ca_i_component))
        #dSL.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, dSL_component))
        testt.append(SEon[-1])
        YES.append(yes)

        dTropTot.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, dTropTot_component))
        I_LCC.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_LCC_component))
        I_RyR.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_RyR_component))
        I_SERCA.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_SERCA_component))
        I_SR.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_SR_component))
        I_NaCa.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_NaCa_component))
        I_pCa.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_pCa_component))
        I_CaB.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, I_CaB_component))
        
        
        
        print(it)
        print(step)
        
        
    get_loop_info() # This function generates the top and bottom curve in the selected workloop

    
   
        
"""
import numpy as np
from scipy.integrate import simps
from numpy import trapz

minSL_top = min(SL_top)
maxSL_top = max(SL_top)
SLstep_top = (maxSL_top[1] - minSL_top[1])/ len(SL_top)

minSL_bottom = min(SL_bottom)
maxSL_bottom = max(SL_bottom)
SLstep_bottom = (maxSL_bottom - minSL_bottom)/ len(SL_bottom)



#Separately finding the area under the top curve and the area under the bottom curve
area_top = trapz(active_tension_top, dx=SLstep_top)
area_bottom = trapz(active_tension_bottom, dx=SLstep_bottom)

print('area_top= ', area_top)
print('area_bottom= ', area_bottom)

workloop_area = area_top - area_bottom
"""


# save the results
import csv



## saving the area(work) and afterload results

"""
area_afterload_results = [value_afterload, workloop_area]

f = open('area_afterload.csv','a')
w = csv.writer(f, dialect='excel')
w.writerow(area_afterload_results)
f.close()
"""
with open('resultsALL.csv', "a") as csvfile:
    resultswriter = csv.writer(csvfile, dialect='excel')
    header_row = ["time", "SL", "active_tension", "F_total", "Ca_i", "integral_force", "value_afterload", "passive", "SEon", "XB_cycling", "gxbT", "XBpostr", "hfT", "SOVFThick", "xXBpostr", "XBprer", "xXBprer", "fxbT", "hbT", "gappT", "fappT", "P", "N", "kn_pT", "kp_nT", "SOVFThin", "dTropTot", "I_LCC", "I_RyR", "I_SERCA", "I_SR", "I_NaCa", "I_pCa", "I_CaB"]
    resultswriter.writerow(header_row)
    for i in range(0, len(time)):
        results_row = [time[i], SL[i], active[i], F_total[i], Ca_i[i], integral_force[i], value_afterload, passive[i], SEon[i], XB_cycling[i], gxbT[i], XBpostr[i], hfT[i], SOVFThick[i], xXBpostr[i], XBprer[i], xXBprer[i], fxbT[i], hbT[i], gappT[i], fappT[i], P[i], N[i], kn_pT[i], kp_nT[i], SOVFThin[i], dTropTot[i], I_LCC[i], I_RyR[i], I_SERCA[i], I_SR[i], I_NaCa[i], I_pCa[i], I_CaB[i]]
        resultswriter.writerow(results_row)


with open('20th_wrkloops.csv', "a") as csvfile: # the a means append ... this keeps the old data
    resultswriter = csv.writer(csvfile, dialect='excel')
    for i in range(542865, 571431):
        results_row = [time[i], SL_WL[i], active[i], F_total_WL[i], XB_cycling[i], Ca_i[i]]
        resultswriter.writerow(results_row)
    csvfile.close()
        
##end_point = iteration_value*28571.5 + 1
##print(end_point)
##with open('resultsSINGLE.csv', "a") as csvfile:
##    resultswriter = csv.writer(csvfile, dialect='excel')
##    for i in range(end_point-28571, end_point):
##        results_row = [time[i], active_tension[i], Ca_i[i], SL[i], F_total[i], integral_force[i], value_afterload, passive[i], SEon[i], XB_cycling[i], gxbT[i], XBpostr[i], hfT[i], SOVFThick[i], xXBpostr[i], XBprer[i], xXBprer[i], fxbT[i], hbT[i], gappT[i], fappT[i], P[i], N[i], kn_pT[i], kp_nT[i], SOVFThin[i]]
##        resultswriter.writerow(results_row)
##        i+= 0.1
        
# Export the results, here we export them as standard exnode, exelem files
if outputFrequency != 0:
    fields = CMISS.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("Monodomain","FORTRAN")
    fields.ElementsExport("Monodomain","FORTRAN")
    fields.Finalise()


##A: time[i]
##B: SL[i],
##C: active_tension[i],
##D: F_total[i],
##E: Ca_i[i],
##F: integral_force[i],
##G: value_afterload,
##H: passive[i],
##I: SEon[i],
##J: XB_cycling[i],
##K: gxbT[i],
##L: XBpostr[i],
##M: hfT[i],
##N: SOVFThick[i],
##O: xXBpostr[i],
##P: XBprer[i],
##Q: xXBprer[i],
##R: fxbT[i],
##S: hbT[i],
##T: gappT[i],
##U: fappT[i],
##V: P[i],
##W: N[i],
##X: kn_pT[i],
##Y: kp_nT[i],
##Z: SOVFThin[i]
##AA: dTropTot[i]
