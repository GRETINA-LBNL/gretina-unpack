import sys
import os

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
bld = Builder(action=Action(rootcint,root_dictionary_message))
env.Append(BUILDERS = {'RootCint':bld})

## Optimization flags ##################################################
env.Append(CCFLAGS = ['-O2', '-D_FILE_OFFSET_64', '-pg', '-g'], LINKFLAGS=['-pg'])
#env.Append(CCFLAGS = ['-O2', '-D_FILE_OFFSET_64'])

## Link auxiliary detector system analysis #############################
#env.Append(CPPDEFINES=['-DWITH_PWALL'])
#env.Append(CPPDEFINES=['-DWITH_S800'])
#env.Append(CPPDEFINES=['-DWITH_LENDA'])
env.Append(CPPDEFINES=['-DWITH_DFMA'])

## Finding dependencies (ROOT)
try:
    env.ParseConfig('root-config --cflags')
    env.ParseConfig('root-config --glibs')
except OSError:
    print "scons: ROOT not found!"
    exit(1)

env.Append(CPPPATH='#')
env.Append(LIBPATH='.')

envUnpack = env.Clone()
envDecompView = env.Clone()
envAddGH = env.Clone()
envClover = env.Clone()

## Building GRETINADict and libGRETINA #################################
gretinaDictTarget = 'src/GRETINADict.cpp'
gretinaDictHeaders = ['src/GRETINA.h', 'src/SortingStructures.h', 
		      'src/GRETINAWavefunction.h',
		      'src/INLCorrection.h', 'src/Histos.h', 
		      'src/Track.h', 'src/LinkDefGRETINA.h']
env.RootCint(gretinaDictTarget, gretinaDictHeaders)

gretinaLibTarget = 'GRETINA'
gretinaLibSources = ['src/GRETINADict.cpp', 'src/GRETINA.cpp', 
		     'src/SortingStructures.cpp', 'src/SuperPulse.cpp', 
		     'src/INLCorrection.cpp', 'src/G3Waveform.cpp',
		     'src/Histos.cpp', 'src/Track.cpp']
env.SharedLibrary(target = gretinaLibTarget, source = gretinaLibSources, 
                  SHLIBPREFIX='lib')

## Building TrackDict and libTrack #####################################
#trackDictTarget = 'src/TrackDict.cpp'
#trackDictHeaders = ['src/Track.h', 'src/LinkDefTrack.h']
#env.RootCint(trackDictTarget, trackDictHeaders)

#trackLibTarget = 'Track'
#trackLibSources = ['src/TrackDict.cpp', 'src/Track.cpp']
#env.SharedLibrary(target = trackLibTarget, source = trackLibSources,
#	          SHLIBPREFIX='lib')

## Building S800Dict and libS800 #######################################
s800DictTarget = 'src/S800Dict.cpp'
s800DictHeaders = ['src/S800Parameters.h', 'src/LinkDefS800.h'] 
env.RootCint(s800DictTarget, s800DictHeaders)

s800LibTarget = 'S800'
s800LibSources = ['src/S800Dict.cpp', 'src/S800Parameters.cpp']
env.SharedLibrary(target = s800LibTarget, source = s800LibSources, 
                  SHLIBPREFIX='lib')

## Building LendaDict and libLenda ######################################
lendaDictTarget = 'src/LendaDict.cpp'
lendaDictHeaders = ['src/LENDA-DDAS.h', 'src/LENDA-Controls.h', 'src/ddasChannel.h', 'src/LinkDefLenda.h'] 
env.RootCint(lendaDictTarget, lendaDictHeaders)

lendaLibTarget = 'Lenda'
lendaLibSources = ['src/LendaDict.cpp', 'src/LENDA-DDAS.cpp', 'src/LENDA-Controls.cpp', 'src/ddasChannel.cpp']
env.SharedLibrary(target = lendaLibTarget, source = lendaLibSources, 
                  SHLIBPREFIX='lib')

## Building CHICODict and libCHICO ######################################
chicoDictTarget = 'src/chicoDict.cpp'
chicoDictHeaders = ['src/CHICO.h', 'src/LinkDefCHICO.h'] 
env.RootCint(chicoDictTarget, chicoDictHeaders)

chicoLibTarget = 'chico'
chicoLibSources = ['src/chicoDict.cpp', 'src/CHICO.cpp']
env.SharedLibrary(target = chicoLibTarget, source = chicoLibSources, 
                  SHLIBPREFIX='lib')

## Building DFMADict and libDFMA ######################################
fmaDictTarget = 'src/fmaDict.cpp'
fmaDictHeaders = ['src/DFMA.h', 'src/LinkDefDFMA.h'] 
env.RootCint(fmaDictTarget, fmaDictHeaders)

fmaLibTarget = 'fma'
fmaLibSources = ['src/fmaDict.cpp', 'src/DFMA.cpp']
env.SharedLibrary(target = fmaLibTarget, source = fmaLibSources, 
                  SHLIBPREFIX='lib')

## Building PhosWallDict and libPhosWall ######################################
pwallDictTarget = 'src/phosWallDict.cpp'
pwallDictHeaders = ['src/PhosWall.h', 'src/LinkDefPhosWall.h'] 
env.RootCint(pwallDictTarget, pwallDictHeaders)

pwallLibTarget = 'phosWall'
pwallLibSources = ['src/phosWallDict.cpp', 'src/PhosWall.cpp']
env.SharedLibrary(target = pwallLibTarget, source = pwallLibSources, 
                  SHLIBPREFIX='lib')

## Building DecompViewerDict and libDecompViewer ########################
#decompViewerDictTarget = 'src/DecompViewerDict.cpp'
#decompViewerDictHeaders = ['src/decompViewer.h', 'src/decompViewerLinkDef.h'] 
#env.RootCint(decompViewerDictTarget, decompViewerDictHeaders)

#decompViewerLibTarget = 'DecompViewer'
#decompViewerLibSources = ['src/DecompViewerDict.cpp', 'src/decompViewer.cpp']
#env.SharedLibrary(target = decompViewerLibTarget, 
#	 	  source = decompViewerLibSources, SHLIBPREFIX='lib')

## Building Unpack executable ###########################################
unpackTarget = 'Unpack'
unpackSources = ['src/Unpack.cpp', 
	         'src/Globals.cpp', 'src/UnpackUtilities.cpp', 
		 'src/S800Functions.cpp']
envUnpack.Append(LIBS=['GRETINA', 'S800', 'chico', 'phosWall', 'Lenda'])
envUnpack.Program(target = unpackTarget, source = unpackSources)

## Building AddGlobalHeaders executable ##################################
#addHeadersTarget = 'AddGlobalHeaders'
#addHeadersSources = ['src/AddGlobalHeaders.cpp', 
#		     'src/GlobalsAGH.cpp']
#envAddGH.Append(LIBS=['GRETINA'])
#envAddGH.Program(target = addHeadersTarget, source = addHeadersSources)

## Building GEB_HFC executable ###########################################
gebTarget = 'GEB_HFC'
gebSources = ['src/hfc/GEB_HFC.cpp', 'src/hfc/HFC.cpp']
env.Program(target = gebTarget, source = gebSources)

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
