#!/usr/bin/env python

""" Date: 13-7-2015.  This is the first attempt at getting the Hinch-Rice-Tran
model working in CellML/Python/OpenCMISS"""

#DOC-START imports
import sys, os, math
# Make sure $OPENCMISS_ROOT/cm/bindings/python is first in our PYTHONPATH.
sys.path.insert(1, os.path.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss import CMISS
#DOC-END imports

# Set problem parameters
#DOC-START parameters

# get the model file
if len(sys.argv) > 1:
    modelFile = sys.argv[1]
else:
    modelFile = "../Hinch_Model3_quickrelease.xml"  # changed from ../n98.xml

# 1D domain size
width = 1.0
numberOfXElements = 2



# Materials parameters and Simulation parameters below will need to be adjusted**************
# Materials parameters
Am = 193.6
Cm = 0.014651
conductivity = 0.0


# Simulation parameters
value_mass = 50
value_viscosity = 0.003
value_SLset = 2.3 #original value is 2.27
value_afterload = 0.1
value_stepStart = 400 #ms

stimstop = 1.5
timeStop = 5000 #ms
stimValue = 1.5
value_TmpC =  22.5 #Celsius
value_x_0 = 0.007 #micrometre
value_stimTime = 333 #ms
# this is used in the integration of the CellML model
odeTimeStep = 0.01 # [ms?]
# this is used in the dummy monodomain problem/solver
pdeTimeStep = 0.01 # [ms]
# this is the step at which we grab output from the solver
outputTimeStep = 0.035
# set this to 1 to get exfiles written out during solve
outputFrequency = 0
#DOC-END parameters




#Setup field number handles --> This is for setup... no change necessary
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
cellMLUserNumber = 1
decompositionUserNumber = 1
equationsSetUserNumber = 1
problemUserNumber = 1
#Mesh component numbers
linearMeshComponentNumber = 1
#Fields
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
materialsFieldUserNumber = 4
equationsSetFieldUserNumber = 5
cellMLModelsFieldUserNumber = 6
cellMLStateFieldUserNumber = 7
cellMLParametersFieldUserNumber = 8
cellMLIntermediateFieldUserNumber = 9

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


#Flagging variables for MeganModel: known= variables we want to be able to set// right now we just want to flag calcium and force, so we do
#not need to flag any variables as known...
#cellML.VariableSetAsKnown(MeganModel, "parameters/finalLambda")
cellML.VariableSetAsKnown(MeganModel, "parameters/mass")
cellML.VariableSetAsKnown(MeganModel, "parameters/SLset")
cellML.VariableSetAsKnown(MeganModel, "parameters/afterload")


#Below are the variables that we want to get from the CellML model (not entirely sure these are correct)
cellML.VariableSetAsWanted(MeganModel, "parameters/active_tension")
cellML.VariableSetAsWanted(MeganModel, "parameters/F_total")
cellML.VariableSetAsWanted(MeganModel, "parameters/XB_cycling")

cellML.VariableSetAsWanted(MeganModel, "parameters/gxbT")
cellML.VariableSetAsWanted(MeganModel, "parameters/hfT")
cellML.VariableSetAsWanted(MeganModel, "parameters/SOVFThick")
cellML.VariableSetAsWanted(MeganModel, "parameters/fxbT")
cellML.VariableSetAsWanted(MeganModel, "parameters/hbT")
cellML.VariableSetAsWanted(MeganModel, "parameters/gappT")
cellML.VariableSetAsWanted(MeganModel, "parameters/fappT")

cellML.VariableSetAsWanted(MeganModel, "parameters/P")
cellML.VariableSetAsWanted(MeganModel, "parameters/kn_pT")
cellML.VariableSetAsWanted(MeganModel, "parameters/kp_nT")


#DOC-START create cellml finish
cellML.CreateFinish()
#DOC-END create cellml finish


#DOC-START map Vm components
# Start the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateStart()
#Now we can set up the field variable component <--> CellML model variable mappings.
#Map Vm

cellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U,1, CMISS.FieldParameterSetTypes.VALUES,MeganModel,"parameters/Vm", CMISS.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(MeganModel,"parameters/Vm", CMISS.FieldParameterSetTypes.VALUES,dependentField,CMISS.FieldVariableTypes.U,1,CMISS.FieldParameterSetTypes.VALUES)

#_#_#cellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U,1, CMISS.FieldParameterSetTypes.VALUES,MeganModel,"parameters/SL", CMISS.FieldParameterSetTypes.VALUES)
#_#_#cellML.CreateCellMLToFieldMap(MeganModel,"parameters/SL", CMISS.FieldParameterSetTypes.VALUES,dependentField,CMISS.FieldVariableTypes.U,1,CMISS.FieldParameterSetTypes.VALUES)


#Finish the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateFinish()

# Set the initial Vm values
dependentField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,-92.5)
#_#_#dependentField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,2.27)
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
SLset_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/SLset")
afterload_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/afterload")
mass_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.PARAMETERS, "parameters/mass")
F_total_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/F_total")
XB_cycling_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/XB_cycling")

active_tension_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/active_tension")
SL_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/SL")
Vm_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/Vm")
Ca_i_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/Ca_i")

gxbT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/gxbT")
XBpostr_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.STATE, "parameters/XBpostr")
hfT_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/hfT")
SOVFThick_component = cellML.FieldComponentGet(MeganModel, CMISS.CellMLFieldTypes.INTERMEDIATE, "parameters/SOVFThick")
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
# and the arrays to store them in

Vm = []
active_tension = [] 
Ca_i = []
afterload = []
F_total = []
SLset = []
SL = []
XB_cycling = []

gxbT = []
XBpostr = []
hfT = []
SOVFThick = []
xXBpostr = []
XBprer = []
xXBprer = []

fxbT = []
hbT = []
gappT = []
fappT = []

P = []
N = []
kn_pT = []
kp_nT = []

# We are using node 1 as the point in our dummy monodomain problem to integrate the CellML model
cellmlNode = 2
cellmlNodeDomain = decomposition.NodeDomainGet(cellmlNode, 1)
cellmlNodeThisComputationalNode = False
if cellmlNodeDomain == computationalNodeNumber:
    cellmlNodeThisComputationalNode = True

#if cellmlNodeThisComputationalNode:
    #cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, stimTime_component, value_stimTime)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SLset_component, value_SLset)
    cellMLParametersField.ParameterSetUpdateNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, afterload_component, value_afterload)
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
    active_tension.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, active_tension_component))
    SL.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SL_component)) 
    Vm.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Vm_component))
    Ca_i.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Ca_i_component))

    afterload.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, afterload_component))
    F_total.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, F_total_component))
    SLset.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SLset_component))
    XB_cycling.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XB_cycling_component))

    gxbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, gxbT_component))
    XBpostr.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XBpostr_component))
    hfT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, hfT_component))
    SOVFThick.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SOVFThick_component))
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


### The below while loop comes from the original python script... Not sure if this is correct... look more into setting the stimulus correctly
XBc = XB_cycling[0]

while currentTime < timeStop:

    # set the next solution interval
    nextTime = currentTime + outputTimeStep
    controlLoop.TimesSet(currentTime, nextTime, outputTimeStep)
    
    # integrate the model
    problem.Solve()

    currentTime = nextTime
    
    # grab results
    #print "Time: " + str(currentTime)
    if cellmlNodeThisComputationalNode:
        time.append(currentTime)
        active_tension.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, active_tension_component))
        SL.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SL_component))
        Vm.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Vm_component))
        Ca_i.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, Ca_i_component))
        F = cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, F_total_component)
        F_total.append(F)
        XBc = cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XB_cycling_component)
        XB_cycling.append(XBc)
        SLset.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SLset_component))
        afterload.append(cellMLParametersField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, afterload_component))

        gxbT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, gxbT_component))
        XBpostr.append(cellMLStateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, XBpostr_component))
        hfT.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, hfT_component))
        SOVFThick.append(cellMLIntermediateField.ParameterSetGetNode(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1, 1, cellmlNode, SOVFThick_component))
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
        
        
# save the results
import csv
with open('results_isometric_total.csv', "a") as csvfile:
    resultswriter = csv.writer(csvfile, dialect='excel')
    for i in range(0, len(time)):
        results_row = [time[i], SL[i], afterload[i], Ca_i[i], active_tension[i],  F_total[i], XB_cycling[i], gxbT[i], XBpostr[i], hfT[i], SOVFThick[i], xXBpostr[i], XBprer[i], xXBprer[i], fxbT[i], hbT[i], gappT[i], fappT[i], P[i], N[i], kn_pT[i], kp_nT[i]]
        resultswriter.writerow(results_row)

selection= 15
end_point= selection*1000
with open('results_isometric.csv', "a") as csvfile:
    resultswriter = csv.writer(csvfile, dialect='excel')
    for i in range(end_point-1000, end_point):
        results_row = [time[i], SL[i], afterload[i], Ca_i[i], active_tension[i],  F_total[i], XB_cycling[i], gxbT[i], XBpostr[i], hfT[i], SOVFThick[i], xXBpostr[i], XBprer[i], xXBprer[i], fxbT[i], hbT[i], gappT[i], fappT[i], P[i], N[i], kn_pT[i], kp_nT[i]]
        resultswriter.writerow(results_row)

##with open('End_Systolic.csv', "a") as csvfile: # the a means append ... this keeps the old data
##    resultswriter = csv.writer(csvfile, dialect='excel')
##    for i in range(0, len(active_tension)):
##        results_row = [time[i], SLset[i], F_total[i], Ca_i[i]]
##        resultswriter.writerow(results_row)
##    csvfile.close()



# Export the results, here we export them as standard exnode, exelem files
if outputFrequency != 0:
    fields = CMISS.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("Monodomain","FORTRAN")
    fields.ElementsExport("Monodomain","FORTRAN")
    fields.Finalise()

