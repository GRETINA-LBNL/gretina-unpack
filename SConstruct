import sys
import os

########################### You shouldn't need to change anything past this point..... ########

colors = {}
colors['cyan']   = '\033[36m'
colors['purple'] = '\033[95m'
colors['blue']   = '\033[94m'
colors['green']  = '\033[92m'
colors['yellow'] = '\033[93m'
colors['red']    = '\033[91m'
colors['end']    = '\033[0m'

## If the output is not a terminal, remove the colors
if not sys.stdout.isatty():
   for key, value in colors.iteritems():
      colors[key] = ''

compile_source_message = '%sCompiling %s===================> %s$SOURCE%s' % \
   (colors['blue'], colors['purple'], colors['cyan'], colors['end'])

compile_shared_source_message = '%sCompiling shared %s============> %s$SOURCE%s' % \
   (colors['blue'], colors['purple'], colors['cyan'], colors['end'])

link_program_message = '%sLinking Program %s=============> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

link_library_message = '%sLinking Static Library %s=====> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

ranlib_library_message = '%sRanlib Library %s===========> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

link_shared_library_message = '%sLinking Shared Library %s======> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

java_library_message = '%sCreating Java Archive %s======> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])

root_dictionary_message = '%sGenerating ROOT dictionary %s==> %s$TARGET%s' % \
   (colors['red'], colors['purple'], colors['cyan'], colors['end'])
 
def rootcint(target, source, env):
    """Executes the ROOT dictionary generator over a list of headers. """
    dictname = target[0]
    headers = ""
    for f in source:
    	headers += str(f) + " "

    command = "rootcint -f %s -c -p %s" %(dictname, headers)
    ok = os.system(command)
    return ok

def remember_pcm(target, source, env):
    new_target = os.path.splitext(str(target[0]))[0]+'_rdict.pcm'
    target.append(new_target)
    return target, source

## Create construction environment propagating the external environment
env = Environment(ENV=os.environ, 
      		  CXXCOMSTR = compile_source_message,
  		  CCCOMSTR = compile_source_message,
  		  SHCCCOMSTR = compile_shared_source_message,
  		  SHCXXCOMSTR = compile_shared_source_message,
  		  ARCOMSTR = link_library_message,
  		  RANLIBCOMSTR = ranlib_library_message,
  		  SHLINKCOMSTR = link_shared_library_message,
  		  LINKCOMSTR = link_program_message,
  		  JARCOMSTR = java_library_message,
  		  JAVACCOMSTR = compile_source_message) 

## Create a rootcint builder and attach it to the environment
bld = Builder(action=Action(rootcint,root_dictionary_message), emitter = remember_pcm)
env.Append(BUILDERS = {'RootCint':bld})

## Optimization flags ##################################################
env.Append(CCFLAGS = ['-O2', '-D_FILE_OFFSET_64', '-g'], LINKFLAGS=[])

## Link auxiliary detector system analysis #############################

## Finding dependencies (ROOT)
try:
    env.ParseConfig('root-config --cflags')
    env.ParseConfig('root-config --glibs')
except OSError:
    print "scons: ROOT not found!"
    exit(1)

env.Append(CPPPATH=['.', './include', './src', './src/hfc/'])
env.Append(LIBPATH='./lib')

envUnpack = env.Clone()
envDecompView = env.Clone()
envAddGH = env.Clone()
envClover = env.Clone()

## Building GRETINADict and libGRETINA #################################
gretinaDictTarget = 'lib/GRETINADict.cpp'
gretinaDictHeaders = ['include/GRETINA.h', 'include/SortingStructures.h', 
		      'include/GRETINAWavefunction.h',
		      'include/Histos.h', 'include/LinkDefGRETINA.h']
env.RootCint(gretinaDictTarget, gretinaDictHeaders)

gretinaLibTarget = 'lib/GRETINA'
gretinaLibSources = ['lib/GRETINADict.cpp', 'src/GRETINA.cpp', 
		     'src/SortingStructures.cpp', 'src/SuperPulse.cpp', 
		     'src/G3Waveform.cpp', 'src/Histos.cpp']
env.SharedLibrary(target = gretinaLibTarget, source = gretinaLibSources, 
                  SHLIBPREFIX='lib')

## Building Unpack executable ###########################################
unpackTarget = 'Unpack'
unpackSources = ['src/Unpack.cpp', 
	         'src/Globals.cpp', 'src/UnpackUtilities.cpp']

envUnpack.Append(LIBS=['GRETINA'], LIBPATH=['lib'])
envUnpack.Program(target = unpackTarget, source = unpackSources)

## Building GEB_HFC executable ###########################################
gebTarget = 'GEB_HFC'
gebSources = ['src/hfc/GEB_HFC.cpp', 'src/hfc/HFC.cpp']
env.Program(target = gebTarget, source = gebSources)

####################### OLD STUFF....but might need... #########################

## Building DecompViewerDict and libDecompViewer ########################
#decompViewerDictTarget = 'src/DecompViewerDict.cpp'
#decompViewerDictHeaders = ['include/decompViewer.h', 'include/decompViewerLinkDef.h'] 
#env.RootCint(decompViewerDictTarget, decompViewerDictHeaders)

#decompViewerLibTarget = 'DecompViewer'
#decompViewerLibSources = ['src/DecompViewerDict.cpp', 'src/decompViewer.cpp']
#env.SharedLibrary(target = decompViewerLibTarget, 
#	 	  source = decompViewerLibSources, SHLIBPREFIX='lib')

## Building AddGlobalHeaders executable ##################################
#addHeadersTarget = 'AddGlobalHeaders'
#addHeadersSources = ['src/AddGlobalHeaders.cpp', 
#		     'src/GlobalsAGH.cpp']
#envAddGH.Append(LIBS=['GRETINA'])
#envAddGH.Program(target = addHeadersTarget, source = addHeadersSources)

## Building DecompViewer executable ######################################
#viewerTarget = 'DecompViewer'
#viewerSources = ['src/decompViewer.cpp']
#envDecompView.Append(LIBS=['DecompViewer', 'GRETINA'])
#envDecompView.Program(target = viewerTarget, source = viewerSources)

## Building MergeArbFiles executable #####################################
#mergeTarget = 'MergeArbFiles'
#mergeSources = ['MergeSrc/MergeModes.cpp']
#env.Program(target = mergeTarget, source = mergeSources)

## Building StripClover executable #######################################
#cloverTarget = 'StripClover'
#cloverSources = ['src/stripCloverData4Coinc.cpp']
#envClover.Append(LIBS=['GRETINA'])
#envClover.Program(target = cloverTarget, source = cloverSources)

## Building StripDataType executable #######################################
#pullDataTarget = 'StripDataType'
#pullDataSources = ['src/stripDataType.cpp']
#envClover.Program(target = pullDataTarget, source = pullDataSources)
