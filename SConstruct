import os
import matlab.engine
import glob
ml = matlab.engine.start_matlab()

# Global Parameter Values
pathToolbox = "../../software/matlab/toolbox-qst/"
trgData = "dataLeCroy.mat"
srcLoadLeCroy = pathToolbox+"loadLeCroy.m"
rawData = "C1flucwires*.txt"

trgPlotPointwiseVariance = "pointwiseVariance.png"
srcPlotPointwiseVariance = ["dataLeCroy.mat", pathToolbox+"pointwiseVariance.m", pathToolbox+"plotPointwiseVariance.m"]

trgPlotStackedWaveforms = "stackedWaveforms.png"
srcPlotStackedWaveforms = ["dataLeCroy.mat", pathToolbox+"plotStackedWaveforms.m"]

# Builder Functions
def bldDataLeCroy(target, source, env):
	ml.cd(env['ACTIVEDIR'])
	ml.loadLeCroy(str(target[0]), env['RAWDATA'])
	return None
	
def bldPlotPointwiseVariance(target, source, env):
	ml.cd(env['ACTIVEDIR'])
	ml.eval("data = load('"+str(source[0])+"','data');",nargout=0)
	ml.eval("data = data.data;",nargout=0)
	ml.eval("plotPointwiseVariance(data,'"+str(target[0])+"');",nargout=0)
	return None
	
def bldPlotStackedWaveforms(target, source, env):
	ml.cd(env['ACTIVEDIR'])
	ml.eval("data = load('"+str(source[0])+"','data');",nargout=0)
	ml.eval("data = data.data;",nargout=0)
	ml.eval("plotStackedWaveforms(data,'"+str(target[0])+"');",nargout=0)
	return None

# Builder Objects
builderDataLeCroy = Builder(action=bldDataLeCroy,suffix='.mat',src_suffix='.m')
builderPlotPointwiseVariance = Builder(action=bldPlotPointwiseVariance,suffix='.png',src_suffix='.m')
builderPlotStackedWaveforms = Builder(action=bldPlotStackedWaveforms,suffix='.png',src_suffix='.m')

# Environemnts
DefaultEnvironment(tools=[]) # Prevent SCons from searching for standard compilers
env = Environment(	BUILDERS = {'dataLeCroy' : builderDataLeCroy, 'plotPointwiseVariance' : builderPlotPointwiseVariance, 'plotStackedWaveforms' : builderPlotStackedWaveforms},
					ACTIVEDIR = os.getcwd(),
					RAWDATA = rawData,
					tools = [])

# Creating
env.plotPointwiseVariance(trgPlotPointwiseVariance, srcPlotPointwiseVariance, bldPlotPointwiseVariance)
env.dataLeCroy(trgData, srcLoadLeCroy, bldDataLeCroy)
env.plotStackedWaveforms(trgPlotStackedWaveforms, srcPlotStackedWaveforms, bldPlotStackedWaveforms)