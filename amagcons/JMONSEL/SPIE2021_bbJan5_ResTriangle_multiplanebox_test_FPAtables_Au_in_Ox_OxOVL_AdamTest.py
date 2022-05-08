# This script determines the yield for electrons incident on 1 or more
# trapezoidal (with top corner radii) resist lines on a 3-layer substrate.								####Note great Jython syntax website:  http://www.jython.org/jythonbook/en/1.0/LangSyntax.html#operators
WinOrLin = 5	#=1 for WINDOWS/StandardSetup/homePC/workPC/Prolith1/Marcy/Spock, =2 for LINUX, =3 for WINDOWS/Metrosim group sim computer, =4 for WINDOWS/Prolith2 computer, =5 for WINDOWS/StandardSetup/homePC/workPC/Prolith1/Marcy/Spock with FPA tables
if (WinOrLin == 2):
	import sys															# add for LINUX, make comment if Windows
	files = [file.strip() for file in open( "./jar.list", "r" )]		# add for LINUX, make comment if Windows
	for file in files:													# add for LINUX, make comment if Windows
		sys.path.append(file)											# add for LINUX, make comment if Windows

import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as ept
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.nanoscalemetrology.JMONSEL as mon
import gov.nist.microanalysis.Utility as nmu
import java.io as jio
import java.util as jutil
import java.lang as jl
import jarray
import java.nio.charset as cs

# determine where to save the results
dest=DefaultOutput;
jio.File(dest).mkdirs()
filename = DefaultOutput+PathSep+"Results.txt"
file = open(filename,'w')
print "Output will be to file: ",filename

if (WinOrLin == 1): 
	tablePathOuter = "C:\NIST\JMONSEL\ScatteringTables"       			   # for home/workPC, "C:\NIST\JMONSEL\ScatteringTables"
if (WinOrLin == 2): 
	tablePathOuter = "/home/bundayb/NIST/JMONSEL/ScatteringTables"          # Only if LINUX on Silvaco1 or Silvaco2
if (WinOrLin == 3): 
	tablePathOuter = "D:\Program Files\NIST\JMONSEL\ScatteringTables"	   # for group sim computer, "D:\Program Files\NIST\JMONSEL\ScatteringTables"
if (WinOrLin == 4): 
	tablePathOuter = "C:\Users\prolith2\Desktop\NIST\JMONSEL"			   # for Prolith2, "C:\Users\prolith2\Desktop\NIST\JMONSEL"
if (WinOrLin == 5): 
	tablePathOuter = "C:\NIST\JMONSEL\FPAScatteringTables"       			   # for home/workPC, "C:\NIST\JMONSEL\FPAScatteringTables"

# Model parameters
nTrajectories = 10
mExp = 1		# # of experiments

# Shape parameters.
CTdiameter=20.
CTheight=20.
Lbox=20.
Wbox=20.
vertexZ=30.

pitchnm = 30. 	# Distance between line centers in nm
nlines = 12		# number of lines [for doe, 1 is iso, 7 is dense]
RotAngDeg = 0.	# Rotation Angle of Grating (in degrees)
ScanOrigX = -26.   # -240.
ScanOrigY = -10.
ScanStepX = 1.
ScanStepY = 10.
NumPixX = 31   #32
NumPixY = 103
linklengthnm = 45.
linkspacenm = 15.
#hnmvals = [50.] 	#resist line height in nm  [10,20,30,40,50,100]
wnmvals = [15.]	#resist line bottom width in nm
linelengthnm = 4000 # resist line length in nm
# Note that sidewall angles are specified with respect to vertical,
# so 0. is vertical, positive angles have bottom wider than top, and
# negative angles are the reverse (undercut).
SWAvals = [0.]  #can loop thru SWA with multiple values if SWAdecision=1, in this case these values overwrite values in both thetardeg and thetaldeg and make them equal.  
#Set to only one value and have SWAdecision=0 if using set thetardeg and thetaldeg.  
SWAdecision = 0   #If =1, then both right and left vals will loop thru SWAvals. if =0 then thetardeg and thetaldeg carry thru.
thetardeg = 3. #0.01	#resist line right sidewall angle in degrees
thetaldeg = 3. #0.01	#resist line left sidewall angle in degrees
radrnm = 0.000001		#resist line top right corner radius in nm
radlnm = 0.000001		#resist line top left corner radius in nm
layer1thicknessnm = 50.	# Thickness in nm of the 1st layer (immediately below the lines)
layer2thicknessnm = 500.# Thickness in nm of the 2nd layer
# We'll make the substrate infinitely thick.
beamEeVvals = [1000.] # Beam energies in eV   [100.,150.,200.,250.,300.,500.,800.]
beamsizenmvals = [0.5] # beam size in nm  [0.,0.1,0.2,0.3,0.5,1.]

DefectType = 0     #0=No Defect, 1=type A defect, 2=type By defect, 3=type Bx defect, 4=type Bx2 defect at link end, 5=line extension, 6=sidewall bump, 
#                        7=roughened link, 8=misaligned link, 9=CD variation&misalign, 10=missing link, 11=mid-link gap, 12=mouse-bite, 13=shortened link

CylDefDiam = 6.			#DefectType=1
CylDefHeight=CylDefDiam
BridgeDefWidth=4.		#DefectType=2/3/4
BridgeDefHeight=20.		#DefectType=2/3/4/5/6/7/8/9/12/13
LineExtDefLength=3.		#DefectType=5
LinkShorten=3.			#DefectType=13
SidewallBumpWidth=2.    #DefectType=6, in x
SidewallBumpLength=4.   #DefectType=6, in y

LERTabAmplitude=1.		#DefectType=7
LERTabLength=2.		    #DefectType=7
LERTabPitch=10.	        #DefectType=7

LinkMisAlign=2.		    #Center Link MisAlignment, goes with DefectType=8 and also DefectType=9
CDvar=9.			    #input width of center link if different from standard link width, goes with DefectType=9
MouseBiteWidth=2.       #in x, amplitude of mousebite, DefectType=12
MouseBiteLength=5.      #in y, how long the mousebite is along link, DefectType=12
MidLinkGap=8.		    #DefectType=11
TriangleLength=1000.
TriangleWidth=100.
TriangleHeight=100.
hnmvals=[0.,20.,50.,100.]  ###OverlayerThk!!!!!   0.,1.,2.,5.,10.,15.,

# Make a record of the random seed we use so we can exactly repeat this calculation
# (same random number sequence) if necessary.
seed = nmu.Math2.rgen.nextLong() # Pick a random seed

# To exactly repeat a previous calculation (e.g., for bug fix) uncomment the next line
# and replace the question marks with the seed that was recorded in the previous calculation's 
# output file.
#seed = -2769499132481846261L
print >>file,"Random number seed: ",seed
print "Random number seed: ",seed
nmu.Math2.initializeRandom(seed)
for i in range(0,10):
    r = nmu.Math2.rgen.nextDouble()
    print >>file, r
    print r

# The following parameter needs a bit of explaining. In the model we'll build below, the infinitely deep
# layer 3 will be artificially divided into two parts, a skin layer that is close to the surface and a 
# deep part that is farther from the surface. The definition of "deep" is set by 
# the deepnm parameter on the next line. Both the skin region and the deep region will contain the same
# material (Si), but in the deep region we'll make a model in which electrons with energies less than
# 50 eV are dropped from the simulation. This can save lots of time (particularly if beam energies are large)
# because there are lots of secondary electrons with energies < 50 eV, and lots of simulation time
# must be devoted to tracking them. We can't afford to drop them when they are generated near the surface,
# because they might escape and be detected. I.e., they're important there. However, low energy electrons
# that are deep inside the sample can't escape, so there is no harm done in not tracking them. Thus, the 
# parameter below should be set to several times the typical escape depth (so there's little chance of 
# dropping an electron that would have escaped). This is only important for high beam energies, because
# only then will the electrons have sufficient range to reach the deeper layer, but in such cases there
# can be a significant time savings.
deepnm = 15. # Depth below which to use the "deep model."

trajImg = 0
trajImgMaxTraj = 500     
trajImgSize = 150.e-9

VRML = 1
VRMLImgMaxTraj = 0 # Include no trajectories in VRML (Show sample only.) Leaving trajectories
# out makes a VRML that displays easily. It's good for checking the sample. I find that adding
# trajectories significantly slows down the display. If you want to try it, keep the number of
# displayed trajectories small (20 is a reasonable number) and turn off collision detection in
# your VRML viewer.


# Make materials

#	A Secondary Electron vaccum
vacuum = mon.SEmaterial()
vacuum.setName("SE vacuum")
vacuumBarrier = mon.ExpQMBarrierSM(vacuum)
vacuumMSM = mon.MONSEL_MaterialScatterModel(vacuum)
vacuumMSM.setBarrierSM(vacuumBarrier)

# PMMA
# Scattering tables to use the DFT model with PMMA are not yet available.
# Instead the code below specifies a back-up model based on the FittedInelSM
# and GanachaudMokraniPolaronTrapSM classes. Each of these classes has two
# free parameters. The parameters given below were chosen to provide a good
# fit to measured yield vs. energy on PMMA. There is no guarantee that these
# parameters will also provide a good fit to topographic yield (yield vs. angle)
# since the data employed for the fit were all acquired at normal incidence.

breakE = epq.ToSI.eV(45.) 
density = 1190.
workfun = 5.5
bandgap = 5. # width of band gap in eV, based on TPP. Svorcik, Lyutakov, Huttel get about 5.3
EFermi = -bandgap # This puts the Fermi level at the top of the valence band.
potU = -workfun-EFermi # Gives 0.5 for this material for now. This is based on Sayyah et al., 
# Int. J. Polymeric Mat. 54, p 505 (2005). It's about mid-range for the values they report for
# different forms of PMMA.

# Material defined in terms of its constituent elements and their weight fractions
# Elemental Constituents
C = epq.Element.C
Ox = epq.Element.O
H = epq.Element.H
PMMAcomp = epq.Composition()
PMMAcomp.defineByMoleFraction([C,Ox,H],[5,2,8])
PMMA = mon.SEmaterial(PMMAcomp,density)
PMMA.setName("PMMA")
PMMAWorkfunction=epq.ToSI.eV(workfun)
PMMA.setWorkfunction(PMMAWorkfunction)
PMMA.setBandgap(epq.ToSI.eV(bandgap))
PMMA.setEnergyCBbottom(epq.ToSI.eV(potU))

# Create scatter mechanisms
PMMANISTMott = mon.SelectableElasticSM(PMMA,mon.NISTMottRS.Factory)
PMMACSD = mon.JoyLuoNieminenCSD(PMMA,breakE)
PMMAfittedInel = mon.FittedInelSM(PMMA,epq.ToSI.eV(65.4),PMMACSD)
# Parameters in the next line are from my fits. (See PMMAOptimization.nb) 
PMMApolaron = mon.GanachaudMokraniPolaronTrapSM(2.e7,1./epq.ToSI.eV(4.))
#PMMAClassicalBarrier = mon.ExpQMBarrierSM(PMMA)
# Make a material scatter model
# to be used in thin layer
PMMAMSM = mon.MONSEL_MaterialScatterModel(PMMA)
PMMAMSM.addScatterMechanism(PMMANISTMott)
PMMAMSM.addScatterMechanism(PMMAfittedInel)
PMMAMSM.addScatterMechanism(PMMApolaron)
PMMAMSM.setCSD(PMMACSD)
#PMMAMSM.setBarrierSM(PMMAClassicalBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
PMMAMSMDeep = mon.MONSEL_MaterialScatterModel(PMMA)
PMMAMSMDeep.addScatterMechanism(PMMANISTMott)
PMMAMSMDeep.addScatterMechanism(PMMAfittedInel)
PMMAMSMDeep.addScatterMechanism(PMMApolaron)
PMMAMSMDeep.setCSD(PMMACSD)
#PMMAMSMDeep.setBarrierSM(PMMAClassicalBarrier)
PMMAMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

# TODO: Generate an ARC model.

# BEGIN TEMPORARY
# Replace the following lines with an ARC model when available.
# I'm replacing the ARC with PMMA during this test phase.
ARCMSM = PMMAMSM
# END TEMPORARY


# Si
phononE = 0.063 # I've seen the number reported as 510 cm^-1. this is conversion of that to eV.
phononStrength = 3. # Turner & Inkson dispersion curves appear to show 3 LO phonon modes converging to the same
# energy at the Gamma point. 
density = 2330.
workfun = 4.85
bandgap = 1.1 # width of band gap in eV
EFermi = -bandgap # This puts the Fermi level at the top of the valence band.
potU = -workfun-EFermi
Si = mon.SEmaterial([epq.Element.Si],[1.],density,"Silicon")
SiWorkfunction=epq.ToSI.eV(workfun)
Si.setWorkfunction(SiWorkfunction)
Si.setEnergyCBbottom(epq.ToSI.eV(potU))
Si.setBandgap(epq.ToSI.eV(bandgap))
Si.setCoreEnergy([epq.ToSI.eV(99.2),epq.ToSI.eV(99.8),epq.ToSI.eV(149.7),epq.ToSI.eV(1839.)])
# Edit the string below so it is the path to the folder where you have stored the silicon scattering tables
# that I provide
if (WinOrLin<>2): 
	tablePath = tablePathOuter+"\SiTables"+PathSep											#If Windows:      tablePath = tablePathOuter+"\SiTables"+PathSep            #If LINUX:      tablePath = tablePathOuter+PathSep+"SiTables"+PathSep                      #for LINUX, change all \ in all lines to "PathSep+"
if (WinOrLin==2): 
	tablePath = tablePathOuter+PathSep+"SiTables"+PathSep
if (WinOrLin<>5):
    SiTables = [tablePath +"IIMFPFullPennInterpSiSI.csv", 
    tablePath +"interpNUSimReducedDeltaEFullPennSiSI.csv", 
    tablePath +"interpNUThetaFullPennSiBGSI.csv", 
    tablePath +"interpSimESE0NUSiBGSI.csv"]
if (WinOrLin==5):
    SiTables = [tablePath +"IIMFP.csv", 
    tablePath +"dered.csv", 
    tablePath +"thscat.csv", 
    tablePath +"ese0.csv"]
# Create scatter mechanisms
SiNISTMott = mon.SelectableElasticSM(Si,mon.NISTMottRS.Factory)
SiDS = mon.TabulatedInelasticSM(Si,3,SiTables,epq.ToSI.eV(13.54))
# The eps0 value is n^2, where n=3.4155 is taken from Palik. epsInfinity is from Palik's 2000 eV value of
# n = 0.9999048
Siphonon = mon.GanachaudMokraniPhononInelasticSM(phononStrength,epq.ToSI.eV(phononE),300.,11.7,1.)
#SiAbruptBarrier = mon.ExpQMBarrierSM(Si,0.)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
SiMSM = mon.MONSEL_MaterialScatterModel(Si)
SiMSM.addScatterMechanism(SiNISTMott)
SiMSM.addScatterMechanism(SiDS)
SiMSM.addScatterMechanism(Siphonon)
#SiMSM.setBarrierSM(SiAbruptBarrier) # Omitting this line causes the barrier to default to gradual/classical
# MSM to be used deep inside (drops electrons with E<50 eV)
SiMSMDeep = mon.MONSEL_MaterialScatterModel(Si)
SiMSMDeep.addScatterMechanism(SiNISTMott)
SiMSMDeep.addScatterMechanism(SiDS)
SiMSMDeep.addScatterMechanism(Siphonon)
#SiMSMDeep.setBarrierSM(SiAbruptBarrier) # Omitting this line causes the barrier to default to gradual/classical
SiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))	#default


#	glassy carbon, set up for TabulatedInelasticSM mode 3 with energy levels as follows:
#	CB bottom at -25.4 eV relative to vacuum = 0 eV.
#	Interface barrier is gradual. 
density = 1800.
workfun = 5.0
bandgap = 0. # width of band gap in eV
EFermi = 20.4 # 
potU = -workfun-EFermi
glC = mon.SEmaterial([epq.Element.C],[1.],density,"glassy Carbon")
glCWorkfunction=epq.ToSI.eV(workfun)
glC.setWorkfunction(glCWorkfunction)
glC.setEnergyCBbottom(epq.ToSI.eV(potU))
glC.setBandgap(epq.ToSI.eV(bandgap))
glC.setCoreEnergy([epq.ToSI.eV(284.2)])
# Edit the string below so it is the path to the folder where you have stored the glassy carbon 
# scattering tables that I provide
if (WinOrLin<>2): 
	tablePath = tablePathOuter + "\glassyCTables"+PathSep								#If Windows:      tablePath = tablePathOuter + "\glassyCTables"+PathSep          #If LINUX:    tablePath = tablePathOuter + PathSep+"glassyCTables"+PathSep
if (WinOrLin==2): 
	tablePath = tablePathOuter + PathSep+"glassyCTables"+PathSep

if (WinOrLin<>5):
    glassyCTables = [tablePath +"IIMFPPennInterpglassyCSI.csv", tablePath +"interpNUSimReducedDeltaEglassyCSI.csv", tablePath +"interpsimTableThetaNUglassyCSI.csv", tablePath +"interpSimESE0NUglassyCSI.csv"]
if (WinOrLin==5):
    glassyCTables = [tablePath+"IIMFP.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
# Create scatter mechanisms
glCNISTMott = mon.SelectableElasticSM(glC,mon.NISTMottRS.Factory)
glCDS = mon.TabulatedInelasticSM(glC,3,glassyCTables)
#glCAbruptBarrier = mon.ExpQMBarrierSM(glC,0.)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
glCMSM = mon.MONSEL_MaterialScatterModel(glC)
glCMSM.addScatterMechanism(glCNISTMott)
glCMSM.addScatterMechanism(glCDS)
#glCMSM.setBarrierSM(glCAbruptBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
glCMSMDeep = mon.MONSEL_MaterialScatterModel(glC)
glCMSMDeep.addScatterMechanism(glCNISTMott)
glCMSMDeep.addScatterMechanism(glCDS)
#glCMSMDeep.setBarrierSM(glCAbruptBarrier)
glCMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

# Cu
#             Copper
density = 8933.
nve = 11
#plasmonE = 9.11
workfun = 4.65
EFermi = 8.7
potU = -workfun-EFermi # Assumes Cu Fermi energy is 8.7 eV
Cu = mon.SEmaterial([epq.Element.Cu],[1.],density,"Copper")
CuWorkfunction=epq.ToSI.eV(workfun)
Cu.setWorkfunction(CuWorkfunction)
Cu.setEnergyCBbottom(epq.ToSI.eV(potU))
Cu.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
if (WinOrLin<>2): 
	tablePath = tablePathOuter + "\CuTables"+PathSep									#If Windows:      tablePath = tablePathOuter + "\CuTables"+PathSep				#If LINUX:      tablePath = tablePathOuter + PathSep+"CuTables"+PathSep
if (WinOrLin==2): 
	tablePath = tablePathOuter + PathSep+"CuTables"+PathSep

if (WinOrLin<>5):
    CuTables = [tablePath+"IIMFPPennInterpCuSI.csv",tablePath+"interpNUSimReducedDeltaECuSI.csv",tablePath+"interpsimTableThetaNUCuSI.csv",tablePath+"interpSimESE0NUCuSI.csv"]
if (WinOrLin==5):
    CuTables = [tablePath+"IIMFP.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
Cu.setCoreEnergy([epq.ToSI.eV(75.1),epq.ToSI.eV(77.3),epq.ToSI.eV(122.5),epq.ToSI.eV(932.7),epq.ToSI.eV(1096.7),epq.ToSI.eV(8979.)])
#density1electron = Cu.getDensity()/epq.Element.Cu.getMass()
#Cu.addBindingEnergy(epq.ToSI.eV(0.)+CuWorkfunction,nve*density1electron)
# Create scatter mechanisms
CuNISTMott = mon.SelectableElasticSM(Cu,mon.NISTMottRS.Factory)
CuDS = mon.TabulatedInelasticSM(Cu,3,CuTables)
#CuMoller = mon.MollerInelasticSM(Cu)
#CuPlasmon = mon.KoteraPlasmonInelasticSM(Cu,1.)
CuBarrier = mon.ExpQMBarrierSM(Cu)
CuCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
CuMSM = mon.MONSEL_MaterialScatterModel(Cu)
CuMSM.addScatterMechanism(CuNISTMott)
CuMSM.addScatterMechanism(CuDS)
CuMSM.setCSD(CuCSD)
CuMSM.setBarrierSM(CuBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
CuMSMDeep = mon.MONSEL_MaterialScatterModel(Cu)
CuMSMDeep.addScatterMechanism(CuNISTMott)
CuMSMDeep.addScatterMechanism(CuDS)
CuMSMDeep.setCSD(CuCSD)
CuMSMDeep.setBarrierSM(CuBarrier)
CuMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

#SiO2
# SiO2 with care taken for TabulatedInelasticSM mode 3 with energy levels as follows:
#             CB bottom at -1.1 eV
#             VB top at -10 eV (i.e., 8.9 eV bandgap)
#             VB bottom at -30 eV (i.e., 28.9 eV offset between CB and VB bottoms)
#             All above energies are relative to vacuum = 0 eV.
# Two phonon modes near 0.145 eV and a gradual barrier are assumed.
density = 2200.
workfun = 10.
phononStrength = 2. # Number of phonon modes
phononE = 0.145 # Phonon mode energy in eV
bandgap = 8.9 # width of band gap in eV
EFermi = -bandgap # This puts the Fermi level at the top of the valence band.
potU = -workfun-EFermi
Si = epq.Element.Si
Ox = epq.Element.O
SiWeight = Si.getAtomicWeight()
OxWeight = 2*Ox.getAtomicWeight()
totalWeight = SiWeight+OxWeight
SiO2 = mon.SEmaterial([Si,Ox],[SiWeight/totalWeight,OxWeight/totalWeight],density,"Silicon Dioxide")
SiO2.setEpsr(3.9)
SiO2Workfunction=epq.ToSI.eV(workfun)
SiO2.setWorkfunction(SiO2Workfunction)
SiO2.setBandgap(epq.ToSI.eV(bandgap))
SiO2.setEnergyCBbottom(epq.ToSI.eV(potU))
SiO2.setCoreEnergy([epq.ToSI.eV(41.6),epq.ToSI.eV(99.2),epq.ToSI.eV(99.8),epq.ToSI.eV(149.7),epq.ToSI.eV(543.1),epq.ToSI.eV(1839.)])
if (WinOrLin<>2): 
	tablePath = tablePathOuter + "\SiO2Tables"+PathSep									#If Windows:      tablePath = tablePathOuter + "\SiO2Tables"+PathSep					#If LINUX:        tablePath = tablePathOuter + PathSep+"SiO2Tables"+PathSep
if (WinOrLin==2): 
	tablePath = tablePathOuter + PathSep+"SiO2Tables"+PathSep

if (WinOrLin<>5):
    SiO2Tables = [tablePath+"IIMFPPennInterpSiO2SI.csv",tablePath+"interpNUSimReducedDeltaESiO2SI.csv",tablePath+"interpsimTableThetaNUSiO2SI.csv",tablePath+"interpSimESE0NUSiO2SI.csv"]
if (WinOrLin==5):
    SiO2Tables = [tablePath+"IIMFP.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
# Create scatter mechanisms
SiO2NISTMott = mon.SelectableElasticSM(SiO2,mon.NISTMottRS.Factory)
SiO2DS = mon.TabulatedInelasticSM(SiO2,3,SiO2Tables,epq.ToSI.eV(20.+bandgap))
SiO2phonon = mon.GanachaudMokraniPhononInelasticSM(phononStrength,epq.ToSI.eV(phononE),300.,3.82,1.)
SiO2polaron = mon.GanachaudMokraniPolaronTrapSM(1.0e9,1./epq.ToSI.eV(1.))
SiO2Barrier = mon.ExpQMBarrierSM(SiO2,1.e-9)
#SiO2CSD = mon.ZeroCSD()                                        # The default. No need to actually execute this line.
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
SiO2MSM = mon.MONSEL_MaterialScatterModel(SiO2)
SiO2MSM.addScatterMechanism(SiO2NISTMott)
SiO2MSM.addScatterMechanism(SiO2DS)
SiO2MSM.addScatterMechanism(SiO2phonon)
# SiO2MSM.addScatterMechanism(SiO2polaron)
#SiO2MSM.setCSD(SiO2CSD)
SiO2MSM.setBarrierSM(SiO2Barrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
SiO2MSMDeep = mon.MONSEL_MaterialScatterModel(SiO2)
SiO2MSMDeep.addScatterMechanism(SiO2NISTMott)
SiO2MSMDeep.addScatterMechanism(SiO2DS)
SiO2MSMDeep.addScatterMechanism(SiO2phonon)
#SiO2MSMDeep.setCSD(SiO2CSD)
SiO2MSMDeep.setBarrierSM(SiO2Barrier)
SiO2MSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

#    Tungsten
density = 19300.
#plasmonE = 9.11
workfun = 4.55
EFermi = 10.1
potU = -workfun-EFermi 
W = mon.SEmaterial([epq.Element.W],[1.],density,"Tungsten")
WWorkfunction=epq.ToSI.eV(workfun)
W.setWorkfunction(WWorkfunction)
W.setEnergyCBbottom(epq.ToSI.eV(potU))
tablePath = "C:\NIST\JMONSEL\ScatteringTables\WTables"+PathSep
WTables = [tablePath+"IIMFPPennInterpWSI.csv",tablePath+"interpNUSimReducedDeltaEWSI.csv",tablePath+"interpsimTableThetaNUWSI.csv",tablePath+"interpSimESE0NUWSI.csv"]
coreEnergies = [31.4, 33.6, 36.8, 45.3, 75.6, 243.5, 255.9, 423.6, 490.4, 594.1,
1809., 1949., 2281., 2575., 2820., 10207., 11544., 12100., 69525.]
for i in range(len(coreEnergies)):
    coreEnergies[i] = epq.ToSI.eV(coreEnergies[i])
W.setCoreEnergy(coreEnergies)

#density1electron = W.getDensity()/epq.Element.W.getMass()
#W.addBindingEnergy(epq.ToSI.eV(0.)+WWorkfunction,nve*density1electron)
# Create scatter mechanisms
WNISTMott = mon.SelectableElasticSM(W,mon.NISTMottRS.Factory)
WDS = mon.TabulatedInelasticSM(W,3,WTables)
#WMoller = mon.MollerInelasticSM(W)
#WPlasmon = mon.KoteraPlasmonInelasticSM(W,1.)
#WBarrier = mon.ExpQMBarrierSM(W)
WCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
WMSM = mon.MONSEL_MaterialScatterModel(W)
WMSM.addScatterMechanism(WNISTMott)
WMSM.addScatterMechanism(WDS)
WMSM.setCSD(WCSD)
#WMSM.setBarrierSM(WBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
WMSMDeep = mon.MONSEL_MaterialScatterModel(W)
WMSMDeep.addScatterMechanism(WNISTMott)
WMSMDeep.addScatterMechanism(WDS)
WMSMDeep.setCSD(WCSD)
#WMSMDeep.setBarrierSM(WBarrier)
WMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

#	Gold
density = 19282.
#plasmonE = 9.11
workfun = 5.1
EFermi = 9.0
potU = -workfun-EFermi # 
Au = mon.SEmaterial([epq.Element.Au],[1.],density,"Gold")
AuWorkfunction=epq.ToSI.eV(workfun)
Au.setWorkfunction(AuWorkfunction)
Au.setEnergyCBbottom(epq.ToSI.eV(potU))
Au.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
if (WinOrLin<>2): 
	tablePath = tablePathOuter + "\AuTables"+PathSep							#If Windows:  tablePath = tablePathOuter + "\AuTables"+PathSep			If LINUX:      tablePath = tablePathOuter + PathSep+"AuTables"+PathSep
if (WinOrLin==2): 
	tablePath = tablePathOuter + PathSep+"AuTables"+PathSep
AuTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
if (WinOrLin<>5):
    AuTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
if (WinOrLin==5):
    AuTables = [tablePath+"IIMFP.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
for en in AuCoreEnergieseV:
    Au.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Au.getDensity()/epq.Element.Au.getMass()
#Au.addBindingEnergy(epq.ToSI.eV(0.)+AuWorkfunction,nve*density1electron)
# Create scatter mechanisms
AuNISTMott = mon.SelectableElasticSM(Au,mon.NISTMottRS.Factory)
AuDS = mon.TabulatedInelasticSM(Au,3,AuTables)
#AuMoller = mon.MollerInelasticSM(Au)
#AuPlasmon = mon.KoteraPlasmonInelasticSM(Au,1.)
AuBarrier = mon.ExpQMBarrierSM(Au,0.05e-9)
AuCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
AuMSM = mon.MONSEL_MaterialScatterModel(Au)
AuMSM.addScatterMechanism(AuNISTMott)
AuMSM.addScatterMechanism(AuDS)
AuMSM.setCSD(AuCSD)
AuMSM.setBarrierSM(AuBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
AuMSMDeep = mon.MONSEL_MaterialScatterModel(Au)
AuMSMDeep.addScatterMechanism(AuNISTMott)
AuMSMDeep.addScatterMechanism(AuDS)
AuMSMDeep.setCSD(AuCSD)
AuMSMDeep.setBarrierSM(AuBarrier)


#	graphite, set up for TabulatedInelasticSM mode 3 with energy levels as follows:
#	CB bottom at -25.4 eV relative to vacuum = 0 eV.
#	Interface barrier is gradual. 
density = 2250.
workfun = 5.0
bandgap = 0. # width of band gap in eV
EFermi = 20.4 # 
potU = -workfun-EFermi
Graphite = mon.SEmaterial([epq.Element.C],[1.],density,"Graphite")
GraphiteWorkfunction=epq.ToSI.eV(workfun)
Graphite.setWorkfunction(GraphiteWorkfunction)
Graphite.setEnergyCBbottom(epq.ToSI.eV(potU))
Graphite.setBandgap(epq.ToSI.eV(bandgap))
Graphite.setCoreEnergy([epq.ToSI.eV(284.2)])
if (WinOrLin<>2): 
	tablePath = tablePathOuter+"\GraphiteTables"+PathSep   #"D:\Program Files\NIST\JMONSEL\ScatteringTables\AuTables"          #If Windows:      tablePath = tablePathOuter+"\GraphiteTables"+PathSep         #If LINUX:      tablePath = tablePathOuter+PathSep+"GraphiteTables"+PathSep
if (WinOrLin==2): 
	tablePath = tablePathOuter+PathSep+"GraphiteTables"+PathSep

if (WinOrLin<>5):
    GraphiteTables = [tablePath +"IIMFPPennInterpgrCSI.csv", tablePath +"interpNUSimReducedDeltaEgrCSI.csv", tablePath +"interpsimTableThetaNUgrCSI.csv", tablePath +"interpSimESE0NUgrCSI.csv"]
if (WinOrLin==5):
    GraphiteTables = [tablePath+"IIMFP.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
# Create scatter mechanisms
GraphiteNISTMott = mon.SelectableElasticSM(Graphite,mon.NISTMottRS.Factory)
GraphiteDS = mon.TabulatedInelasticSM(Graphite,3,GraphiteTables)
#GraphiteDS.setBranchingRatios([0.0415246])
GraphiteClassicalBarrier = mon.ExpQMBarrierSM(Graphite)
GraphiteCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
GraphiteMSM = mon.MONSEL_MaterialScatterModel(Graphite)
GraphiteMSM.addScatterMechanism(GraphiteNISTMott)
GraphiteMSM.addScatterMechanism(GraphiteDS)
GraphiteMSM.setCSD(GraphiteCSD)
GraphiteMSM.setBarrierSM(GraphiteClassicalBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
GraphiteMSMDeep = mon.MONSEL_MaterialScatterModel(Graphite)
GraphiteMSMDeep.addScatterMechanism(GraphiteNISTMott)
GraphiteMSMDeep.addScatterMechanism(GraphiteDS)
GraphiteMSMDeep.setCSD(GraphiteCSD)
GraphiteMSMDeep.setBarrierSM(GraphiteClassicalBarrier)
GraphiteMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))



# Conversions of shape parameters to SI units.
# Shape parameters.
meterspernm = 1.e-9	# conversion from nanometers to meters
pitch = pitchnm*meterspernm

linelength = linelengthnm*meterspernm
# Note that sidewall angles are specified with respect to vertical,
# so 0. is vertical, positive angles have bottom wider than top, and
# negative angles are the reverse (undercut).
radperdeg = jl.Math.PI/180.	# conversion from degrees to radians
thetar = thetardeg*radperdeg
thetal = thetaldeg*radperdeg
RotAng = RotAngDeg*radperdeg
radr = radrnm*meterspernm
radl = radlnm*meterspernm
layer1thickness = layer1thicknessnm*meterspernm
layer2thickness = layer2thicknessnm*meterspernm

deep = deepnm*meterspernm 
RotPitch = pitch/(jl.Math.cos(RotAng))


linklengthnm = 36.
linkspacenm = 12.


# Print simulation parameters to window and output file.
print >>file, "# Trajectories at each landing position: ",nTrajectories
print "# Trajectories at each landing position: ",nTrajectories
print >>file, "# Experiments at each landing position: ",mExp
print "# Experiments at each landing position: ",mExp
print >>file, "# Defect Type: ",DefectType
print "# Defect Type: ",DefectType
print >>file, "# 0=No Defect, 1=type A defect, 2=type By defect, 3=type Bx defect, 4=type Bx2 defect at link end, 5=line extension, 6=sidewall bump"
print "# 0=No Defect, 1=type A defect, 2=type By defect, 3=type Bx defect, 4=type Bx2 defect at link end, 5=line extension, 6=sidewall bump"
print >>file, "# 7=roughened link, 8=misaligned link, 9=CD variation&misalign, 10=missing link, 11=mid-link gap, 12=mouse-bite, 13=shortened link"
print "# 7=roughened link, 8=misaligned link, 9=CD variation&misalign, 10=missing link, 11=mid-link gap, 12=mouse-bite, 13=shortened link"
print >>file, "# lines: ",nlines
print "# lines: ",nlines
print >>file, "Pitch of lines (nm): ",pitchnm
print "Pitch of lines (nm): ",pitchnm
print >>file, "Rotated Pitch of lines (nm): ",pitchnm/(jl.Math.cos(RotAng))
print "Rotated Pitch of lines (nm): ",pitchnm/(jl.Math.cos(RotAng))
print >>file, "Grating Rotation Angle (deg): ",RotAngDeg
print "Grating Rotation Angle (deg): ",RotAngDeg
print >>file, "Line height (nm): ",hnmvals
print "Line height: ",hnmvals
print >>file, "Line bottom width (nm): ",wnmvals
print "Line bottom width (nm): ",wnmvals
print >>file, "# LinkLength y[nm]: ",linklengthnm
print "# LinkLength y[nm]: ",linklengthnm
print >>file, "# LinkSpace y[nm]: ",linkspacenm
print "# LinkSpace y[nm]: ",linkspacenm
print >>file, "Line length (nm): ",linelengthnm
print "Line length (nm): ",linelengthnm
print >>file, "SWA range (deg): ",SWAvals
print "SWA range (deg): ",SWAvals
print >>file, "Left and right SWA (deg): ",thetaldeg,thetardeg
print "Left and right SWA (deg): ",thetaldeg,thetardeg
print >>file, "Left and right top corner radii (nm): ",radlnm,radrnm
print "Left and right top corner radii (nm): ",radlnm,radrnm
print >>file, "Thicknesses of 1st and second layers (nm): ",layer1thicknessnm,layer2thicknessnm
print "Thicknesses of 1st and second layers (nm): ",layer1thicknessnm,layer2thicknessnm
print >>file, "Beam landing energies (eV): ",beamEeVvals
print "Beam landing energies (eV): ",beamEeVvals
print >>file, "Beam size (standard deviation, in nm): ",beamsizenmvals
print "Beam size (standard deviation, in nm): ",beamsizenmvals
print >>file, "Scan Origin (x,y, in nm): ",ScanOrigX,ScanOrigY
print "Scan Origin (x,y, in nm): ",ScanOrigX,ScanOrigY
print >>file, "Scan Step (x,y, in nm): ",ScanStepX, ScanStepY
print "Scan Step (x,y, in nm): ",ScanStepX,ScanStepY
print >>file, "# Pixels (x,y): ",NumPixX, NumPixY
print "# Pixels (x,y): ",NumPixX, NumPixY
print >>file, "# CylDefDiam[nm]: ",CylDefDiam
print "# Def1 CylDefDiam[nm]: ",CylDefDiam
print >>file, "# Def1 CylDefHeight[nm]: ",CylDefHeight
print "# Def1 CylDefHeight[nm]: ",CylDefHeight
print >>file, "# Def2/3/4 BridgeDefWidth[nm]: ",BridgeDefWidth
print "# Def2/3/4 BridgeDefWidth[nm]: ",BridgeDefWidth
print >>file, "# Def2/3/4/5/6/7/8/9/12/13 BridgeDefHeight[nm]: ",BridgeDefHeight
print "# Def2/3/4/5/6/7/8/9/12/13 BridgeDefHeight[nm]: ",BridgeDefHeight
print >>file, "# Def5 LineExtDefLength[nm]: ",LineExtDefLength
print "# Def5 LineExtDefLength[nm]: ",LineExtDefLength
print >>file, "# Def13 LinkShorten[nm]: ",LinkShorten
print "# Def13 LinkShorten[nm]: ",LinkShorten
print >>file, "# Def6 SidewallBumpWidth x[nm]: ",SidewallBumpWidth
print "# Def6 SidewallBumpWidth x[nm]: ",SidewallBumpWidth
print >>file, "# Def6 SidewallBumpLength y[nm]: ",SidewallBumpLength
print "# Def6 SidewallBumpLength y[nm]: ",SidewallBumpLength
print >>file, "# Def7 LERTabAmplitude[nm]: ",LERTabAmplitude
print "# Def7 LERTabAmplitude[nm]: ",LERTabAmplitude
print >>file, "# Def7 LERTabLength[nm]: ",LERTabLength
print "# Def7 LERTabLength[nm]: ",LERTabLength
print >>file, "# Def7 LERTabPitch[nm]: ",LERTabPitch
print "# Def7 LERTabPitch[nm]: ",LERTabPitch
print >>file, "# Def8/9 LinkMisAlign[nm]: ",LinkMisAlign
print "# Def8/9 LinkMisAlign[nm]: ",LinkMisAlign
print >>file, "# Def9 CDvar[nm]: ",CDvar
print "# Def9 CDvar[nm]: ",CDvar
print >>file, "# Def12 MouseBiteWidth x[nm]: ",MouseBiteWidth
print "# Def12 MouseBiteWidth x[nm]: ",MouseBiteWidth
print >>file, "# Def12 MouseBiteLength y[nm]: ",MouseBiteLength
print "# Def12 MouseBiteLength y[nm]: ",MouseBiteLength
print >>file, "# Def11 MidLinkGap[nm]: ",MidLinkGap
print "# Def11 MidLinkGap[nm]: ",MidLinkGap

print >>file	# Blank line before start of calculation results
print 

print >>file, "    beamsznm     h(nm)     w(nm)  beamE(eV)   Exp#	x(nm)	 y(nm)	BSEyield	SEyield	  ElapsedTime"     #print >>file, "beamsize(nm) h(nm) w(nm) SWAl(deg) SWAr(deg) beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield	Elapsed Time"
print         "    beamsznm     h(nm)     w(nm)  beamE(eV)   Exp#	x(nm)	 y(nm)	BSEyield	SEyield	  ElapsedTime"     #print "beamsize(nm) h(nm) w(nm) SWAl(deg) SWAr(deg) beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield	Elapsed Time"

t0 = jl.System.currentTimeMillis()





#loop beamsizenmvals
for beamsizenm in beamsizenmvals:

	beamsize = beamsizenm*meterspernm

	# create an instance of the model
	monte=nm.MonteCarloSS() #creates an instance of the model with all default characteristics
	eg = nm.GaussianBeam(beamsize) # makes electron gun, Gaussian with standard deviation = beamsize
	monte.setElectronGun(eg) # This gun is "attached" to the model.

	# SAMPLE DESCRIPTION

	# NISTMonte provides us with a "chamber" in the form of a 0.1 m sphere inside of which we build
	# out sample. Replace the default vacuum in the chamber with SEvacuum. (SEmaterials define additional
	# properties, such as work function, that are needed by JMONSEL.)
	chamber = monte.getChamber()
	chamber.updateMaterial(chamber.getScatterModel(),vacuumMSM)

	# Generate the sample. The GaussianBeam electron gun has this pecularity: it defines the +z axis to
	# be in the direction of travel of the electrons. When we describe the sample in this coordinate
	# system, it is inverted along the z direction.

	# Make sample component shapes
	normalvector = [0.,0.,-1.]

	# First we make the layers. The simplest way to do this is to define each as a MultiPlaneShape with a single
	# plane, each nested inside the previous one.
	layer1 = mon.NormalMultiPlaneShape()
	layer1.addPlane(normalvector,[0.,0.,0.]) #layer 1 is now the half space of everything above the x-y plane
	# This region has shape defined by layer1, scattering properties defined for ARC, and is a subregion of (is 
	# wholly contained within) the chamber.
	layer1Region = monte.addSubRegion(chamber,SiO2MSM,layer1)  

	layer2 = mon.NormalMultiPlaneShape()
	layer2.addPlane(normalvector,[0.,0.,layer1thickness]) #layer 2 starts layer1thickness farther up.
	# We give it the properties of PMMA, and make it a subregion of layer1Region. At this point, layer1Region
	# extends only for 0<=z<=layer1thickness.
	layer2Region = monte.addSubRegion(layer1Region,SiMSMDeep,layer2)  

	layer3 = mon.NormalMultiPlaneShape()
	layer3.addPlane(normalvector,[0.,0.,layer1thickness+layer2thickness]) #layer 3 starts 
	# yet another layer2thickness farther up.
	# We give it the properties of Si, and make it a subregion of layer2Region. At this point, layer2Region
	# extends only for layer1thickness<=z<=layer1thickness+layer2thickness.
	layer3Region = monte.addSubRegion(layer2Region,SiMSMDeep,layer3)
	
	
	#loop hnm in hnmvals
	for hnm in hnmvals:
			
		#loop wnm in wnmvals
		for wnm in wnmvals:

			h = hnm*meterspernm 
			w = wnm*meterspernm
			#loop SWAtemp in SWAvals
			for SWAtemp in SWAvals:
				if SWAdecision: 
					thetardeg=SWAtemp
					thetaldeg=SWAtemp

				thetar = thetardeg*radperdeg
				thetal = thetaldeg*radperdeg
				
				# Make the array of lines. The integer divide in (nlines/2) truncates fractions.
				# The result always has one line centered at (x,y)=(0,0). If nlines is odd
				# the remaining lines are placed symmetrically left and right of this one.
				# If nlines is even, there will be one more line on the right side than on the
				# left side.
				
				#standard grating		###############comment out this IDA section or the grating section
				#PitchWalk = 2*1e-9
				#LERtabAmplitude = 2*1e-9
				#LERperiod = 50*1e-9
				#nLERtabs = linelength//LERperiod
				#leftmostLineCenterx = -pitch * (nlines/2)
				#for i in range(nlines):
				#	if (i % 2 == 0):
				#		xcenter = leftmostLineCenterx+i*pitch
				#	if (i % 2 == 1):
				#		xcenter = leftmostLineCenterx+i*pitch+PitchWalk
				#	line = mon.NShapes.createLine(-h,w,linelength,thetal,thetar,radl,radr)
				#	#line.rotate([0.,0.,0.],RotAng, 0., 0.)          #rotating chamber at end
				#	line.translate([xcenter,0.,0.])
				#	lineRegion = monte.addSubRegion(chamber,CuMSM,line)  # set material for the lines here!!!!!!!!
				######Add LER to center line
				#If (LERtabAmplitude > 1):
				#	bottommostLERtabcentery = -LERperiod * (nLERtabs/2)
				#	for k in range(nLERtabs):
				#		line = mon.NShapes.createLine(-h,LERtabAmplitude,LERperiod/2,thetal,thetar,radl,radr)
				#		line.translate([w/2+LERtabAmplitude/2,bottommostLERtabcentery+k*LERperiod,0.])
				#		lineRegion = monte.addSubRegion(chamber,CuMSM,line)    # set material for the lines here!!!!!!!!

				#IDA SRAM cell			###############comment out this IDA section or the grating section
                
                #########New SPIE2020 cell#########
				#######tall   h=50nm########
                #block = mon.NShapes.createLine(-CTheight*1E-9,2000*1E-9,2000*1E-9,0.,0.,0.,0.)  #	line = mon.NShapes.createLine(-h,w,linelength,thetal,thetar,radl,radr)
                #block = mon.NormalMultiPlaneShape()  #	line = mon.NShapes.createLine(-h,w,linelength,thetal,thetar,radl,radr)
                #block.addPlane([-1.,0.,0.],[-Wbox*1e-9/2,0,0])
                #block.addPlane([+1.,0.,0.],[+Wbox*1e-9/2,0,0])
                #block.addPlane([0.,-1.,0.],[0.,-Lbox*1e-9/2,0])
                #block.addPlane([0.,+1.,0.],[0.,+Lbox*1e-9/2,0])
                #block.addPlane([0.,0.,-1.],[0,0,-CTheight*1e-9])
                #block.addPlane([0.,0.,+1.],[0,0,0])
                #block.translate([-20.*1E-9,0.,0.])
                #blockregion= monte.addSubRegion(chamber,AuMSM,block)
                
                #link.translate([0*1E-9,0*1E-9,0.])
                #linkRegion = monte.addSubRegion(chamber,SiO2MSM,link)  # set material for the lines here!!!!!!!!
                #CylHole1 = mon.NormalCylindricalShape([0*1e-9,0*1e-9,0.],[0*1e-9,0*1e-9,-CTheight*1e-9],CTdiameter*1e-9/2)  #CylDefect1 = mon.NormalCylindricalShape([0.,0.,0.],[0.,0.,-h],w/2)  #SphDefect1 = mon.NormalSphereShape([0.,0.,-w/2],w/2)
                #CylDefectRegion = monte.addSubRegion(chamber,SiMSM,CylDefect1)     #CylDefect1Region = monte.addSubRegion(chamber,SiMSM,CylDefect1) 
                #pyramidhole=mon.NormalMultiPlaneShape()  
                #pyramidhole.addPlane([-CTdiameter*1E-9/2,-CTdiameter*1E-9/2,0.],[-CTdiameter*1E-9/2,+CTdiameter*1E-9/2,0.],[0.,0.,vertexZ*1E-9]) #SW corner-->NW corner-->nadir
                #pyramidhole.addPlane([+CTdiameter*1E-9/2,+CTdiameter*1E-9/2,0.],[+CTdiameter*1E-9/2,-CTdiameter*1E-9/2,0.],[0.,0.,vertexZ*1E-9]) #NE corner-->SE corner-->nadir
                #pyramidhole.addPlane([+CTdiameter*1E-9/2,-CTdiameter*1E-9/2,0.],[-CTdiameter*1E-9/2,-CTdiameter*1E-9/2,0.],[0.,0.,vertexZ*1E-9]) #SE corner-->SW corner-->nadir
                #pyramidhole.addPlane([-CTdiameter*1E-9/2,+CTdiameter*1E-9/2,0.],[+CTdiameter*1E-9/2,+CTdiameter*1E-9/2,0.],[0.,0.,vertexZ*1E-9]) #NW corner-->NE corner-->nadir
                #pyramidhole.addPlane([+CTdiameter*1E-9/2,+CTdiameter*1E-9/2,0.],[-CTdiameter*1E-9/2,+CTdiameter*1E-9/2,0.],[-CTdiameter*1E-9/2,-CTdiameter*1E-9/2,0.]) #NE corner-->NW corner-->SW corner (if inside solid you looking up at it)
                #pyramidhole.addPlane([-CTdiameter*1E-9/2,-CTdiameter*1E-9/2,CTheight*1E-9],[-CTdiameter*1E-9/2,+CTdiameter*1E-9/2,CTheight*1E-9],[+CTdiameter*1E-9/2,+CTdiameter*1E-9/2,CTheight*1E-9]) #SW corner-->NW corner-->NE corner (looking down at it)
                
                triangle=mon.NormalMultiPlaneShape()  #vertices N/S/E/W and will be 100nm tall
                triangle.addPlane([+TriangleWidth*1E-9/2,0.,0.],[-TriangleWidth*1E-9/2,0.,0.],[-TriangleWidth*1E-9/2,0.,-TriangleHeight*1e-9]) #N corner-->W corner-->up
                triangle.addPlane([-TriangleWidth*1E-9/2,0.,0.],[0.,TriangleLength*1e-9,0.],[-TriangleWidth*1E-9/2,0.,-TriangleHeight*1e-9]) #W corner-->S corner-->up
                triangle.addPlane([0.,TriangleLength*1e-9,0.],[+TriangleWidth*1E-9/2,0.,0.],[+TriangleWidth*1E-9/2,0.,-TriangleHeight*1e-9]) #S corner-->E corner-->up
                triangle.addPlane([0.,0.,+1.],[0.,0.,0.]) #plane pointing down thru origin for bottom of pyramid
                triangle.addPlane([0.,0.,-1.],[0.,0.,-TriangleHeight*1e-9]) #plane pointing down thru origin for bottom of pyramid
                #pyramidhole.addPlane([-CTdiameter*1E-9/2,-CTdiameter*1E-9/2,CTheight*1E-9],[-CTdiameter*1E-9/2,+CTdiameter*1E-9/2,CTheight*1E-9],[+CTdiameter*1E-9/2,+CTdiameter*1E-9/2,CTheight*1E-9]) #SW corner-->NW corner-->NE corner (looking down at it)
                triangle.translate([0.,0.,0.])
                
                film=mon.NormalMultiPlaneShape()
                film.addPlane([0.,0.,+1.],[0.,0.,0.]) #plane pointing down thru origin for bottom of pyramid
                film.addPlane([0.,0.,-1.],[0.,0.,-TriangleHeight*1e-9]) #plane pointing up thru origin for bottom of pyramid
                FilmMinusTriangle=mon.NormalDifferenceShape(film,triangle)
                FilmMinusTriangleRegion = monte.addSubRegion(chamber,SiO2MSM,FilmMinusTriangle)
                triangleregion= monte.addSubRegion(chamber,AuMSM,triangle)
                film2=mon.NormalMultiPlaneShape()
                film2.addPlane([0.,0.,+1.],[0.,0.,-TriangleHeight*1e-9]) #plane pointing down thru origin for bottom of pyramid
                film2.addPlane([0.,0.,-1.],[0.,0.,-(TriangleHeight+hnm)*1e-9]) #plane pointing up thru origin for bottom of pyramid
                Film2Region = monte.addSubRegion(chamber,SiO2MSM,film2)
                
                #link = mon.NShapes.createLine(-50*1E-9,30*1E-9,25*1E-9,0,0,radl,radr)
                #link.translate([70*1E-9,47.5*1E-9,0.])
                #linkRegion = monte.addSubRegion(chamber,SiMSM,link)  # set material for the lines here!!!!!!!!
                #CylDefect1 = mon.NormalCylindricalShape([-10*1e-9,55*1e-9,0.],[-10*1e-9,55*1e-9,-50*1e-9],10*1e-9/2)  #CylDefect1 = mon.NormalCylindricalShape([0.,0.,0.],[0.,0.,-h],w/2)  #SphDefect1 = mon.NormalSphereShape([0.,0.,-w/2],w/2)
                #CylDefectRegion = monte.addSubRegion(chamber,SiMSM,CylDefect1)     #CylDefect1Region = monte.addSubRegion(chamber,SiMSM,CylDefect1) 
                #CylDefect1 = mon.NormalCylindricalShape([-10*1e-9,35*1e-9,0.],[-10*1e-9,35*1e-9,-50*1e-9],8*1e-9/2)  #CylDefect1 = mon.NormalCylindricalShape([0.,0.,0.],[0.,0.,-h],w/2)  #SphDefect1 = mon.NormalSphereShape([0.,0.,-w/2],w/2)
                
                
                
                ##From IDA, no longr used but still here in case of orphaned variables.  The placement of the links is disabled.
                #nlinks = 6   #31
				#linklength = linklengthnm*1.e-9
				#linkspace = linkspacenm*1.e-9
				#linkpitch = linklength + linkspace
				#PitchWalk = 0 * 1.e-9
				#leftmostLineCenterx = -pitch * (nlines/2)
				#bottommostLinkCentery = -linkpitch * (nlinks/2)
				#for i in range(nlines):
					#if (i % 2 == 0):
						#xcenter = leftmostLineCenterx + i * pitch
					#if (i % 2 == 1):
						#xcenter = leftmostLineCenterx+i*pitch+PitchWalk
					#for k in range(nlinks):
						#ycenter = bottommostLinkCentery + k * linkpitch
						##if (i % 4 > 1):                                                                #this makes the links jog half a period in pairs; % means mod function.  If remainder =2 or 3 then no jog, if remainder = 0 or 1 then half linkpitch jog
							##ycenter = bottommostLinkCentery + k * linkpitch
						##if (i % 4 < 2):
							##ycenter = bottommostLinkCentery + (k + 0.5) * linkpitch	
						#if (DefectType < 8):    #note this makes normal links everywhere for DefectType<8
							#link = mon.NShapes.createLine(-h,w,linklength,thetal,thetar,radl,radr)
							#link.translate([xcenter,ycenter,0.])
							#linkRegion = monte.addSubRegion(chamber,SiMSM,link)  # set material for the lines here!!!!!!!!	
						#if ((DefectType >= 8) and ((xcenter<>0) or (ycenter<>0))):         #note If DefectType>=8 and at center link, will not place a link (leaves room for custom link for subtractor defects)
							#link = mon.NShapes.createLine(-h,w,linklength,thetal,thetar,radl,radr)
							#link.translate([xcenter,ycenter,0.])
							#linkRegion = monte.addSubRegion(chamber,SiMSM,link)  # set material for the lines here!!!!!!!!
							

				#if (DefectType == 0):     #No Defect
					#ntab=linklength % LERTabPitch
					#No Defect, "Perfect Cell", no further action necessary
				
				#Cyl Defects                                                                                RotPitch = pitch/(jl.Math.cos(RotAng))
				#if (DefectType == 1):      #Type A defect
					#CylDefect1 = mon.NormalCylindricalShape([-pitchnm/2*1e-9,0*1e-9,0.],[-12*1e-9,0*1e-9,-CylDefHeight*1e-9],CylDefDiam*1e-9/2)  #CylDefect1 = mon.NormalCylindricalShape([0.,0.,0.],[0.,0.,-h],w/2)  #SphDefect1 = mon.NormalSphereShape([0.,0.,-w/2],w/2)
					#CylDefectRegion = monte.addSubRegion(chamber,SiMSM,CylDefect1)     #CylDefect1Region = monte.addSubRegion(chamber,SiMSM,CylDefect1)  #SphDefect1Region = monte.addSubRegion(chamber,SiMSM,SphDefect1)    #This is where defect is defined.  Also where defect material is changed.
				
				#Line defects
				#if (DefectType == 2):     #By defect
					#LineDefect1 = mon.NShapes.createLine(-BridgeDefHeight*1e-9,BridgeDefWidth*1e-9,linkspace,thetal,thetar,radl,radr)	#bridge-x   (w.r.t. grating)
					#LineDefect1.translate([0.,-24e-9,0.])
					#LineDefectRegion = monte.addSubRegion(chamber,SiMSM,LineDefect1)
				
				#if (DefectType == 3):     #Bx defect at link center
					#LineDefect2 = mon.NShapes.createLine(-BridgeDefHeight*1e-9,BridgeDefWidth*1e-9,pitch-w,thetal,thetar,radl,radr)	#bridge-y   (w.r.t. grating)
					#LineDefect2.rotate([0.,0.,0.],90.*radperdeg,0.,0.)
					#LineDefect2.translate([pitch/2,0.,0.])
					#LineDefectRegion = monte.addSubRegion(chamber,SiMSM,LineDefect2)	
					
				#if (DefectType == 4):     #Bx2 defect at link end
					#LineDefect3 = mon.NShapes.createLine(-BridgeDefHeight*1e-9,BridgeDefWidth*1e-9,pitch-w,thetal,thetar,radl,radr)	#bridge-y   (w.r.t. grating)
					#LineDefect3.rotate([0.,0.,0.],90.*radperdeg,0.,0.)
					#LineDefect3.translate([pitch/2,linklength/2-BridgeDefWidth/2*1e-9,0.])
					#LineDefectRegion = monte.addSubRegion(chamber,SiMSM,LineDefect3)
					
				#if (DefectType == 5):     #Line Extension
					#LineDefect4 = mon.NShapes.createLine(-BridgeDefHeight*1.e-9,w,LineExtDefLength*1.e-9,thetal,thetar,radl,radr)	#bridge-y   (w.r.t. grating)
					#LineDefect4.rotate([0.,0.,0.],0.,0.,0.)
					#LineDefect4.translate([0.,linklength/2+LineExtDefLength/2*1.e-9,0.])
					#LineDefectRegion = monte.addSubRegion(chamber,SiMSM,LineDefect4)
					
				#if (DefectType == 6):     #Sidewall Bump
					#LineDefect5 = mon.NShapes.createLine(-BridgeDefHeight*1e-9,SidewallBumpWidth*1.e-9,SidewallBumpLength*1.e-9,thetal,thetar,radl,radr)	#bridge-y   (w.r.t. grating)
					#LineDefect5.rotate([0.,0.,0.],0.,0.,0.)
					#LineDefect5.translate([w/2+SidewallBumpWidth/2*1.e-9,0.,0.])
					#LineDefectRegion = monte.addSubRegion(chamber,SiMSM,LineDefect5)
					
				#if (DefectType == 7):     #Roughened Link
					#ntab = (linklengthnm / LERTabPitch - 1)
					#bottommostTabCentery = (-1*LERTabPitch*1.e-9)*(ntab/2) + (LERTabLength*1.e-9)
					#for iLER in range(ntab):
						#ycenterLER = bottommostTabCentery + iLER*LERTabPitch*1.e-9
						#linetab6 = mon.NShapes.createLine(-BridgeDefHeight*1.e-9,(LERTabAmplitude*1.e-9),(LERTabLength*1.e-9),0.,0.,0.,0.)
						#linetab6.translate([(-w/2-LERTabAmplitude*1.e-9/2),ycenterLER,0.])
						#linetab6Region = monte.addSubRegion(chamber,SiMSM,linetab6)  # set material for the lines here!!!!!!!!
						#linetab7 = mon.NShapes.createLine(-BridgeDefHeight*1.e-9,(LERTabAmplitude*1.e-9),(LERTabLength*1.e-9),0.,0.,0.,0.)
						#linetab7.translate([(+w/2+LERTabAmplitude*1.e-9/2),ycenterLER,0.])
						#linetab7Region = monte.addSubRegion(chamber,SiMSM,linetab7)  # set material for the lines here!!!!!!!!
				
				#if (DefectType == 8):     #Mis-Aligned Link              #note if DefectType=8 it adds center link, shifted over a little.  If DefectType>8 and at center link, no link (leaves room for custom link for subtractor defects)
						#link8 = mon.NShapes.createLine(-BridgeDefHeight*1.e-9,w,linklength,thetal,thetar,radl,radr)
						#link8.translate([0.+LinkMisAlign*1.e-9,0.,0.])
						#link8Region = monte.addSubRegion(chamber,SiMSM,link8)  # set material for the lines here!!!!!!!!
						
				#if (DefectType == 9):     #CD Variation & MisAlign              #note this builds an alternative link at the center link
						#link9 = mon.NShapes.createLine(-BridgeDefHeight*1.e-9,CDvar*1.e-9,linklength,thetal,thetar,radl,radr)
						#link9.translate([0.+LinkMisAlign*1.e-9,0.,0.])
						#link9Region = monte.addSubRegion(chamber,SiMSM,link9)  # set material for the lines here!!!!!!!!
					
				#if (DefectType == 10):     #Missing Link
					#ntab=linklength % LERTabPitch
					#Missing Link, Link at (0,0) skipped above in link logic, no further action necessary
				
				#if (DefectType == 11):     #Mid-Link Gap                 #note this builds an alternative link at the center link
					#LineDefect8 = mon.NShapes.createLine(-h,w,(linklength/2-MidLinkGap*1e-9/2),thetal,thetar,radl,radr)	#bridge-x   (w.r.t. grating)
					#LineDefect8.translate([0.,(-linklength/4-(MidLinkGap*1e-9/4)),0.])
					#LineDefect8Region = monte.addSubRegion(chamber,SiMSM,LineDefect8)
					#LineDefect9 = mon.NShapes.createLine(-h,w,(linklength/2-MidLinkGap*1e-9/2),thetal,thetar,radl,radr)	#bridge-x   (w.r.t. grating)
					#LineDefect9.translate([0.,(+linklength/4+(MidLinkGap*1e-9/4)),0.])
					#LineDefect9Region = monte.addSubRegion(chamber,SiMSM,LineDefect9)
				
				#if (DefectType == 12):     #MouseBite                    #note this builds an alternative link at the center link
					#LineDefect10 = mon.NShapes.createLine(-h,w,(linklength/2-(MouseBiteLength*1e-9)/2),thetal,thetar,radl,radr)	#bridge-x   (w.r.t. grating)
					#LineDefect10.translate([0.,(-linklength/4-(MouseBiteLength*1e-9/4)),0.])
					#LineDefect10Region = monte.addSubRegion(chamber,SiMSM,LineDefect10)
					#LineDefect11 = mon.NShapes.createLine(-h,w,(linklength/2-(MouseBiteLength*1e-9)/2),thetal,thetar,radl,radr)	#bridge-x   (w.r.t. grating)
					#LineDefect11.translate([0.,(+linklength/4+(MouseBiteLength*1e-9/4)),0.])
					#LineDefect11Region = monte.addSubRegion(chamber,SiMSM,LineDefect11)
					#LineDefect12 = mon.NShapes.createLine(-BridgeDefHeight*1.e-9,(w-MouseBiteWidth*1e-9),MouseBiteLength*1e-9,thetal,thetar,radl,radr)	#bridge-x   (w.r.t. grating)
					#LineDefect12.translate([(-w/2+(w-MouseBiteWidth*1e-9)/2),0.,0.])
					#LineDefect12Region = monte.addSubRegion(chamber,SiMSM,LineDefect12)
				
				#if (DefectType == 13):     #Shortened Link              #note if DefectType=13 it adds center link, shifted up a little.  If DefectType>8 and at center link, no link (leaves room for custom link for subtractor defects)
					#link13 = mon.NShapes.createLine(-BridgeDefHeight,w,linklength-LinkShorten*1.e-9,thetal,thetar,radl,radr)
					#link13.translate([0.,LinkShorten/2*1.e-9,0.])
					#link13Region = monte.addSubRegion(chamber,SiMSM,link13)  # set material for the lines here!!!!!!!!				

					
                chamber.rotate([0.,0.,0.],RotAng,0.,0.)         #Entire sample rotated!
				
				# Scan parameters
                m5=1
                deltam5=1
                m5vals = []
                while m5-1 < mExp:
                    m5vals.append(m5)
                    m5 += deltam5
				
                #If scan, use this block.  If manual, comment it and fill in yvals manually.
                yvals = [1000.,990.,980.,970.,960.,950.,940.,930.,920.,910.,900.,890.,880.,870.,860.,850.,840.,830.,820.,810.,800.,790.,780.,770.,760.,750.,740.,730.,720.,710.,700.,650.,600.,550.,500.]
                #yend = ScanOrigY + NumPixY * ScanStepY
                #y = ScanOrigY
                #while y < yend:
                    #yvals.append(y)
                    #y += ScanStepY
				
				#yvals = [0.,24.,48.,72.,96.,120.,144.,168.,192.,216.,240.] #A single line scan at y = 0 (center of lines)
				# The following parameters are set for a 201 nm scan centered on the position where the right
				# edge of the center line intersects the substrate. We could sample x at equal intervals (e.g., every nm)
				# but fine sampling is mainly important where the topography changes rapidly, so the code here samples
				# every 5 nm except within the zone that starts 25 nm to the left of the top corner and ends 25 nm to the
				# right of the bottom corner.

                #If scan, use this block.  If manual, comment it and fill in yvals manually.
				#xvals = [-480.,-456.,-432.,-408.,-384.,-360.,-336.,-312.,-288.,-264.,-240.,-216.,-192.,-168.,-144.,-120.,-96.,-72.,-48.,-24.,0.,24.,48.,72.,96.,120.,144.,168.,192.,216.,240.,264.,288.,312.,336.,360.,384.,408.,432.,456.,480.]
                xvals = []
                xend = ScanOrigX + NumPixX * ScanStepX
                x = ScanOrigX
                while x < xend:
                    xvals.append(x)
                    x += ScanStepX
					

                binSizeEV = 10.	# Width (in eV) of bins in energy histogram
                for beamEeV in beamEeVvals:
                    beamE = epq.ToSI.eV(beamEeV)
                    monte.setBeamEnergy(beamE) # sets this model's beam energy
                    for m5 in m5vals:
                        for xnm in xvals:          #####################To switch x and y, do it here and in 4 lines below
                            x = xnm*meterspernm
                            for ynm in yvals:
                                y = ynm*meterspernm
                                eg.setCenter([x,y,-1500.*meterspernm]) # Aims the gun at x,y.

								# Define our backscatter detector.
                                back=nm.BackscatterStats(monte)	
                                nbins = int(beamEeV/binSizeEV)
                                monte.addActionListener(back)
                                back.setEnergyBinCount(nbins)
									
								# Add a trajectory image
                                if trajImg:  # output the trajectory image
                                    img=nm.TrajectoryImage(512,512,trajImgSize)
                                    img.setMaxTrajectories(trajImgMaxTraj)
                                    img.setYRange(-h-trajImgSize/10,.9*trajImgSize-h)
                                    img.setXRange(x-trajImgSize/2.,x+trajImgSize/2.)
                                    monte.addActionListener(img)

								# Add a vrml
                                if VRML:  # output the trajectory image
									#fos=jio.FileOutputStream("%s/angle - %lg deg.wrl" % (dest, phi))
                                    fos=jio.FileOutputStream( "%s.wrl" % (dest))
                                    tw=jio.OutputStreamWriter(fos,cs.Charset.forName("UTF-8"))
                                    vrml=nm.TrajectoryVRML(monte,tw)
                                    vrml.setMaxTrajectories(VRMLImgMaxTraj)
                                    vrml.setTrajectoryWidth(0.25e-9)
                                    vrml.setDisplayBackscatter(1)
                                    vrml.addView("Gun",[0.0,0.0,-5.0e-7],[0.0,0.0,0.0])
                                    vrml.addView("X-Axis",[1.0e-6,0.0,0.0,0.0],[0.0,0.0,0.0])
                                    vrml.addView("Y-Axis",[0.0,1.0e-6,0.0,0.0],[0.0,0.0,0.0])
                                    vrml.addView("Close perspective",[-110.e-8,100.e-8,-100.e-8],[-100.e-9,0.,0.])
                                    vrml.renderSample()
                                    monte.addActionListener(vrml)

								# Run the simulation
                                monte.runMultipleTrajectories(nTrajectories)

                                hist = back.backscatterEnergyHistogram()
								#fhist = back.forwardscatterEnergyHistogram()
                                energyperbineV = beamEeV/hist.binCount()
                                maxSEbin = 50./energyperbineV	# bin number of the one with 50 eV
                                totalSE = 0
                                for j in range(0,int(maxSEbin)):
                                    totalSE = totalSE+hist.counts(j)
									#print j,hist.counts(j)
									#totalSE = totalSE+hist.counts(j)+fhist.counts(j)

                                SEf = float(totalSE)/nTrajectories
                                bsf = back.backscatterFraction()-SEf#+float(totalFwd)/nTrajectories
                                t = (jl.System.currentTimeMillis()-t0)/3600000.
                                print >>file, "%8.2f \t%8.2f \t%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%8.6f" % (beamsizenm,hnm,wnm,beamEeV,m5,xnm,ynm,bsf,SEf,t)
                                print "%8.2f \t%8.2f \t%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%8.6f" % (beamsizenm,hnm,wnm,beamEeV,m5,xnm,ynm,bsf,SEf,t)             #print "%8.2f \t%8.2f \t%8.2f \t%8.2f \t%8.2f \t%8.1f \t%8.2f \t%8.2f \t%3.4f \t%3.4f \t%8.4f" % (beamsizenm,hnm,wnm,thetaldeg,thetardeg,beamEeV,xnm,ynm,bsf,SEf,t)
                                rdirname = jl.Integer(int(beamEeV)).toString()+"eV_"+"x"+jl.Integer(int(xnm)).toString()+'y'+jl.Integer(int(ynm)).toString()

                                #back.dump(jio.FileOutputStream(dest+PathSep+rdirname+"backscatter.prn"))    ###Turn On/Off saving all the BSE hist files
                                monte.removeActionListener(back)
					
                                if trajImg:  # output the trajectory image
                                    img.dumpToFile(dest+PathSep+rdirname)
                                    monte.removeActionListener(img)

                                if VRML:
                                    tw.flush()
                                    fos.close()
                                    monte.removeActionListener(vrml)

#**************************************************************
#ROUND 2				
#**************************************************************
# This script determines the yield for electrons incident on 1 or more
# trapezoidal (with top corner radii) resist lines on a 3-layer substrate.
 # determine where to save the results
#dest=DefaultOutput;
#jio.File(dest).mkdirs()
#filename = DefaultOutput+PathSep+"Results2.txt"
#file = open(filename,'w')
#print "Output will be to file: ",filename

