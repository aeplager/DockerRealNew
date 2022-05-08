import gov.nist.microanalysis.EPQLibrary as epq
#import gov.nist.microanalysis.EPQTools as ept
#import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.nanoscalemetrology.JMONSEL as mon

import cmath

PathSep="/"

#tablePathOuter = "C:\\NIST\\JMONSEL\\ScatteringTables\\"
tablePathOuter = '/work/ScatteringTables/'
#tablePathOuterFPA = "C:\\NIST\\JMONSEL\\FPAScatteringTables\\"
tablePathOuterFPA ='/work/FPAScatteringTables/'
tablePathgl = tablePathOuter+"glassyCTables"+PathSep
tablePathSi = tablePathOuter+"SiTables"+PathSep
tablePathSiO2 = tablePathOuter+"SiO2Tables"+PathSep
tablePathCu = tablePathOuter+"CuTables"+PathSep
tablePathW = tablePathOuter+"WTables"+PathSep
tablePathAu = tablePathOuter+"AuTables"+PathSep
tablePathGraphite = tablePathOuter+"GraphiteTables"+PathSep
tablePathRu = tablePathOuter+"RuTables"+PathSep
tablePathMo = tablePathOuter+"MoTables"+PathSep
tablePathAg = tablePathOuter+"AgTables"+PathSep
tablePathAl = tablePathOuter+"AlTables"+PathSep
tablePathFe = tablePathOuter+"FeTables"+PathSep
tablePathPt = tablePathOuter+"PtTables"+PathSep
tablePathTa = tablePathOuter+"TaTables"+PathSep
tablePathTi = tablePathOuter+"TiTables"+PathSep
tablePathwhiteSn = tablePathOuter+"whiteSnTables"+PathSep
tablePathNi = tablePathOuter+"NiTables"+PathSep
tablePathTiN = tablePathOuter+"TiNTables"+PathSep
tablePathGaAs = tablePathOuter+"GaAsTables"+PathSep
tablePathDiamond = tablePathOuter+"DiamondTables"+PathSep
tablePathCo = tablePathOuter+"CoTables"+PathSep
tablePathGe = tablePathOuter+"GeTables"+PathSep
tablePathCr = tablePathOuter+"CrTables"+PathSep


# all shape parameters are in meters
# they are converted in the definesample.py script

def createBlock(xmin,xmax,ymin,ymax,zmin,zmax):
    block=mon.NormalMultiPlaneShape()
    block.addPlane([-1.,0.,0.],[xmin,0,0])
    block.addPlane([+1.,0.,0.],[xmax,0,0])
    block.addPlane([0.,-1.,0.],[0.,ymin,0])
    block.addPlane([0.,+1.,0.],[0.,ymax,0])
    block.addPlane([0.,0.,-1.],[0,0,zmin])
    block.addPlane([0.,0.,+1.],[0,0,zmax])
    return block

# For pyramid "above" aubstrate, vertexZ is positive, CTheight is positive, CTdiameter is a distance
# CTHeight is translated in -z, cone points to axes
def createPyramid(NDivisions,CTheight,CTdiameter,vertexZ):
    pyramid=mon.NormalMultiPlaneShape()
    pyramidPlanes = makePyramidPlanes(NDivisions,CTdiameter/2,[0.,0.,vertexZ])
    for plane in pyramidPlanes:
	pyramid.addPlane(*plane)                
    pyramid.addPlane([0.,0.,-1.],[0.,0.,0.]) #plane pointing up thru origin for base of pyramid
    pyramid.addPlane([0,0,1],[0,0,CTheight]) #cutoff planes
#    pyramid.translate([+0.,0.,-CTheight])
    return pyramid

def makePyramidPlanes(Ndivisions,radius,apex,min=1e-15):
#provides point combinations for wall planes of a multiplanar cone/N-gonal pyramid; points down = +z
#assuming: planes' normal vectors from .addPlane(p1,p2,p3) == (p2-p1) cross (p3-p1)
#preferably create planar shape inside
#(1) allowing apex = (x,y,z) will allow for skewed pyramids, but also inward pointing normal vectors; fix with apex = distance
#would like to use point transformations
    coords = []
    planes = []
    for n in range(Ndivisions):
        z = radius*cmath.exp(1j*2*cmath.pi*n/Ndivisions)
        coords.append([z.real if abs(z.real)>=min else 0.,z.imag if abs(z.imag)>=min else 0.,0.])
    
    for p in range(Ndivisions):
        planes.append([coords[p],coords[(p+1)%Ndivisions],apex])    
    return planes

def vacuumModel():
#	A Secondary Electron vaccum
	vacuum = mon.SEmaterial()
	vacuum.setName("SE vacuum")
	vacuumBarrier = mon.ExpQMBarrierSM(vacuum)
	vacuumMSM = mon.MONSEL_MaterialScatterModel(vacuum)
	vacuumMSM.setBarrierSM(vacuumBarrier)
	return(vacuumMSM)

# PMMA: Scattering tables for the DFT model of PMMA are not yet available.
# Instead the code below implements a backup model based on the FittedInelSM class,
# supplemented with a charge trapping model. Each of these has two free parameters.
# These are chosen to provide the best fit to measured SE yield vs. energy. 
# There is no guarantee that a model so constructed will match topographic yield
# (yield vs. angle of incidence) since all of the data used to determine the 
# parameters were measured at normal incidence. Nevertheless, it's the best we
# can do for now.

def PMMAModel():
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
	return (PMMAMSM, PMMAMSMDeep)

def PMMACSDModel():
# PMMA is inherently a CSD model
        PMMAMSM,PMMAMSMDeep = PMMAModel()
	return (PMMAMSM, PMMAMSMDeep)

# TODO: Generate an ARC model.

# BEGIN TEMPORARY
# Replace the following lines with an ARC model when available.
# I'm replacing the ARC with PMMA during this test phase.
#ARCMSM = PMMAMSM
# END TEMPORARY


def glCModel():
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
#tablePath = "C:\Program Files\NIST\JMONSEL\ScatteringTables\glassyCTables"+PathSep
#	tablePathgl = tablePath+"glassyCTables"+PathSep
	glCTables = [tablePathgl +"IIMFPPennInterpglassyCSI.csv", tablePathgl +"interpNUSimReducedDeltaEglassyCSI.csv", tablePathgl +"interpsimTableThetaNUglassyCSI.csv", tablePathgl +"interpSimESE0NUglassyCSI.csv"]
# Create scatter mechanisms
	glCNISTMott = mon.SelectableElasticSM(glC,mon.NISTMottRS.Factory)
	glCDS = mon.TabulatedInelasticSM(glC,3,glCTables)
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
	return (glCMSM, glCMSMDeep)

def SiModel():
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
#tablePath = "C:\Program Files\NIST\JMONSEL\ScatteringTables\SiTables"+PathSep
	SiTables = [tablePathSi +"IIMFPFullPennInterpSiSI.csv", 
	tablePathSi +"interpNUSimReducedDeltaEFullPennSiSI.csv", 
	tablePathSi +"interpNUThetaFullPennSiBGSI.csv", 
	tablePathSi +"interpSimESE0NUSiBGSI.csv"]
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
	return (SiMSM,SiMSMDeep)

def CuModel():
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
	tablePath=tablePathCu
	CuTables = [tablePath+"IIMFPPennInterpCuSI.csv",tablePath+"interpNUSimReducedDeltaECuSI.csv",tablePath+"interpsimTableThetaNUCuSI.csv",tablePath+"interpSimESE0NUCuSI.csv"]
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
	return(CuMSM,CuMSMDeep)

def SiO2Model():
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
	tablePath=tablePathSiO2
	SiO2Tables = [tablePath+"IIMFPPennInterpSiO2SI.csv",tablePath+"interpNUSimReducedDeltaESiO2SI.csv",tablePath+"interpsimTableThetaNUSiO2SI.csv",tablePath+"interpSimESE0NUSiO2SI.csv"]
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
	return (SiO2MSM,SiO2MSMDeep)

#    Tungsten
def WModel():
	density = 19300.
#plasmonE = 9.11
	workfun = 4.55
	EFermi = 10.1
	potU = -workfun-EFermi 
	W = mon.SEmaterial([epq.Element.W],[1.],density,"Tungsten")
	WWorkfunction=epq.ToSI.eV(workfun)
	W.setWorkfunction(WWorkfunction)
	W.setEnergyCBbottom(epq.ToSI.eV(potU))
	tablePath = tablePathW
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
	return(WMSM,WMSMDeep)

def AuModel():
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
	tablePath=tablePathAu
	AuTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
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
	return (AuMSM,AuMSMDeep)

def GraphiteModel():
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
	tablePath = tablePathGraphite
	GraphiteTables = [tablePath+"IIMFPPennInterpgrCSI.csv",tablePath+"interpNUSimReducedDeltaEgrCSI.csv",tablePath+"interpsimTableThetaNUgrCSI.csv",tablePath+"interpSimESE0NUgrCSI.csv"]
#        GraphiteTables = [tablePath +"iimfp.csv", tablePath +"dered.csv", tablePath +"thscat.csv", tablePath +"ese0.csv"]
# Create scatter mechanisms
	GraphiteNISTMott = mon.SelectableElasticSM(Graphite,mon.NISTMottRS.Factory)
	GraphiteDS = mon.TabulatedInelasticSM(Graphite,3,GraphiteTables)
	GraphiteDS.setBranchingRatios([0.0453593])
	#GraphiteDS.setE0fromDispersion(1);
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
	return(GraphiteMSM,GraphiteMSMDeep)

def AgModel():
#	Silver
    density = 10501.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 7.2
    potU = -workfun-EFermi # 
    Ag = mon.SEmaterial([epq.Element.Ag],[1.],density,"Silver")
    AgWorkfunction=epq.ToSI.eV(workfun)
    Ag.setWorkfunction(AgWorkfunction)
    Ag.setEnergyCBbottom(epq.ToSI.eV(potU))
    Ag.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
    tablePath = tablePathAg
    AgTables = [tablePath+"IIMFPPennInterpAgSI.csv",tablePath+"interpNUSimReducedDeltaEAgSI.csv",tablePath+"interpsimTableThetaNUAgSI.csv",tablePath+"interpSimESE0NUAgSI.csv"]
    AgCoreEnergieseV = [66., 368.3, 374, 3351, 3524, 3806, 25514]
    for en in AgCoreEnergieseV:
	Ag.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ag.getDensity()/epq.Element.Ag.getMass()
#Ag.addBindingEnergy(epq.ToSI.eV(0.)+AgWorkfunction,nve*density1electron)
# Create scatter mechanisms
    AgNISTMott = mon.SelectableElasticSM(Ag,mon.NU_ELSEPA_DCS.FactoryMT)
    AgDS = mon.TabulatedInelasticSM(Ag,3,AgTables)
    AgDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#AgMoller = mon.MollerInelasticSM(Ag)
#AgPlasmon = mon.KoteraPlasmonInelasticSM(Ag,1.)
    AgBarrier = mon.ExpQMBarrierSM(Ag)
    AgCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    AgMSM = mon.MONSEL_MaterialScatterModel(Ag)
    AgMSM.addScatterMechanism(AgNISTMott)
    AgMSM.addScatterMechanism(AgDS)
    AgMSM.setCSD(AgCSD)
    AgMSM.setBarrierSM(AgBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    AgMSMDeep = mon.MONSEL_MaterialScatterModel(Ag)
    AgMSMDeep.addScatterMechanism(AgNISTMott)
    AgMSMDeep.addScatterMechanism(AgDS)
    AgMSMDeep.setCSD(AgCSD)
    AgMSMDeep.setBarrierSM(AgBarrier)
    AgMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(AgMSM,AgMSMDeep)

def AlModel():
#	Aluminum
	density = 2700.
#plasmonE = 9.11
	workfun = 4.28
	EFermi = 11.2
	potU = -workfun-EFermi # 
	Al = mon.SEmaterial([epq.Element.Al],[1.],density,"Aluminum")
	AlWorkfunction=epq.ToSI.eV(workfun)
	Al.setWorkfunction(AlWorkfunction)
	Al.setEnergyCBbottom(epq.ToSI.eV(potU))
	Al.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
	tablePath = tablePathAl
	AlTables = [tablePath+"IIMFPPennInterpAlSI.csv",tablePath+"interpNUSimReducedDeltaEAlSI.csv",tablePath+"interpsimTableThetaNUAlSI.csv",tablePath+"interpSimESE0NUAlSI.csv"]
	AlCoreEnergieseV = [72.7, 1559.]
	for en in AlCoreEnergieseV:
		Al.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Al.getDensity()/epq.Element.Al.getMass()
#Al.addBindingEnergy(epq.ToSI.eV(0.)+AlWorkfunction,nve*density1electron)
# Create scatter mechanisms
	AlNISTMott = mon.SelectableElasticSM(Al,mon.NU_ELSEPA_DCS.FactoryMT)
	AlDS = mon.TabulatedInelasticSM(Al,3,AlTables)
	AlDS.setBranchingRatios([0.3243, 0.800729])
#AlMoller = mon.MollerInelasticSM(Al)
#AlPlasmon = mon.KoteraPlasmonInelasticSM(Al,1.)
	AlBarrier = mon.ExpQMBarrierSM(Al)
	AlCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
	AlMSM = mon.MONSEL_MaterialScatterModel(Al)
	AlMSM.addScatterMechanism(AlNISTMott)
	AlMSM.addScatterMechanism(AlDS)
	AlMSM.setCSD(AlCSD)
	AlMSM.setBarrierSM(AlBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
	AlMSMDeep = mon.MONSEL_MaterialScatterModel(Al)
	AlMSMDeep.addScatterMechanism(AlNISTMott)
	AlMSMDeep.addScatterMechanism(AlDS)
	AlMSMDeep.setCSD(AlCSD)
	AlMSMDeep.setBarrierSM(AlBarrier)
	AlMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
        return(AlMSM,AlMSMDeep)

def RuModel():
#	Ruthenium
	density = 12410.
#plasmonE = 28.54
	workfun = 4.7
	EFermi = 6.9
	potU = -workfun-EFermi # 
	Ru = mon.SEmaterial([epq.Element.Ru],[1.],density,"Gold")
	RuWorkfunction=epq.ToSI.eV(workfun)
	Ru.setWorkfunction(RuWorkfunction)
	Ru.setEnergyCBbottom(epq.ToSI.eV(potU))
	tablePath = tablePathRu
	RuTables = [tablePath+"IIMFPPennInterpRuSI.csv",tablePath+"interpNUSimReducedDeltaERuSI.csv",tablePath+"interpsimTableThetaNURuSI.csv",tablePath+"interpSimESE0NURuSI.csv"]
	RuCoreEnergieseV = [284.2, 2838, 22117]
	for en in RuCoreEnergieseV:
		Ru.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ru.getDensity()/epq.Element.Ru.getMass()
#Ru.addBindingEnergy(epq.ToSI.eV(0.)+RuWorkfunction,nve*density1electron)
# Create scatter mechanisms
	RuNISTMott = mon.SelectableElasticSM(Ru,mon.NISTMottRS.Factory)
	RuDS = mon.TabulatedInelasticSM(Ru,3,RuTables)
	RuDS.setBranchingRatios([0.66, 0.318851, 0.167881])
#RuMoller = mon.MollerInelasticSM(Ru)
#RuPlasmon = mon.KoteraPlasmonInelasticSM(Ru,1.)
	RuBarrier = mon.ExpQMBarrierSM(Ru)
	RuCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
	RuMSM = mon.MONSEL_MaterialScatterModel(Ru)
	RuMSM.addScatterMechanism(RuNISTMott)
	RuMSM.addScatterMechanism(RuDS)
	RuMSM.setCSD(RuCSD)
	RuMSM.setBarrierSM(RuBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
	RuMSMDeep = mon.MONSEL_MaterialScatterModel(Ru)
	RuMSMDeep.addScatterMechanism(RuNISTMott)
	RuMSMDeep.addScatterMechanism(RuDS)
	RuMSMDeep.setCSD(RuCSD)
	RuMSMDeep.setBarrierSM(RuBarrier)
	RuMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
        return(RuMSM,RuMSMDeep)

def MoModel():
#	Molybdenum
    density = 10280.
#plasmonE = 23.09
    workfun = 5.1
    EFermi = 6.5
    potU = -workfun-EFermi # 
    Mo = mon.SEmaterial([epq.Element.Mo],[1.],density,"Molybdenum")
    MoWorkfunction=epq.ToSI.eV(workfun)
    Mo.setWorkfunction(MoWorkfunction)
    Mo.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathMo
    MoTables = [tablePath+"IIMFPPennInterpMoSI.csv",tablePath+"interpNUSimReducedDeltaEMoSI.csv",tablePath+"interpsimTableThetaNUMoSI.csv",tablePath+"interpSimESE0NUMoSI.csv"]
    MoCoreEnergieseV = [35.5, 227.9, 2520., 20000.]
    for en in MoCoreEnergieseV:
	Mo.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Mo.getDensity()/epq.Element.Mo.getMass()
#Mo.addBindingEnergy(epq.ToSI.eV(0.)+MoWorkfunction,nve*density1electron)
# Create scatter mechanisms
    MoNISTMott = mon.SelectableElasticSM(Mo,mon.NU_ELSEPA_DCS.FactoryMT)
    MoDS = mon.TabulatedInelasticSM(Mo,3,MoTables)
    MoDS.setBranchingRatios([0.78, 0.62, 0.441262, 0.199441])
    #MoMoller = mon.MollerInelasticSM(Mo)
    #MoPlasmon = mon.KoteraPlasmonInelasticSM(Mo,1.)
    MoBarrier = mon.ExpQMBarrierSM(Mo)
    MoCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    MoMSM = mon.MONSEL_MaterialScatterModel(Mo)
    MoMSM.addScatterMechanism(MoNISTMott)
    MoMSM.addScatterMechanism(MoDS)
    MoMSM.setCSD(MoCSD)
    MoMSM.setBarrierSM(MoBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    MoMSMDeep = mon.MONSEL_MaterialScatterModel(Mo)
    MoMSMDeep.addScatterMechanism(MoNISTMott)
    MoMSMDeep.addScatterMechanism(MoDS)
    MoMSMDeep.setCSD(MoCSD)
    MoMSMDeep.setBarrierSM(MoBarrier)
    MoMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(MoMSM,MoMSMDeep)

def FeModel():
#	Iron
    density = 7874.
#plasmonE = 9.11
    workfun = 4.67
    EFermi = 8.9
    potU = -workfun-EFermi # 
    Fe = mon.SEmaterial([epq.Element.Fe],[1.],density,"Iron")
    FeWorkfunction=epq.ToSI.eV(workfun)
    Fe.setWorkfunction(FeWorkfunction)
    Fe.setEnergyCBbottom(epq.ToSI.eV(potU))
    Fe.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
    tablePath = tablePathFe
    FeTables = [tablePath+"IIMFPPennInterpFeSI.csv",tablePath+"interpNUSimReducedDeltaEFeSI.csv",tablePath+"interpsimTableThetaNUFeSI.csv",tablePath+"interpSimESE0NUFeSI.csv"]
    FeCoreEnergieseV = [710., 7112.]
    for en in FeCoreEnergieseV:
	Fe.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Fe.getDensity()/epq.Element.Fe.getMass()
#Fe.addBindingEnergy(epq.ToSI.eV(0.)+FeWorkfunction,nve*density1electron)
# Create scatter mechanisms
    FeNISTMott = mon.SelectableElasticSM(Fe,mon.NU_ELSEPA_DCS.FactoryMT)
    FeDS = mon.TabulatedInelasticSM(Fe,3,FeTables)
#FeDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#FeMoller = mon.MollerInelasticSM(Fe)
#FePlasmon = mon.KoteraPlasmonInelasticSM(Fe,1.)
    FeBarrier = mon.ExpQMBarrierSM(Fe)
    FeCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    FeMSM = mon.MONSEL_MaterialScatterModel(Fe)
    FeMSM.addScatterMechanism(FeNISTMott)
    FeMSM.addScatterMechanism(FeDS)
    FeMSM.setCSD(FeCSD)
    FeMSM.setBarrierSM(FeBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    FeMSMDeep = mon.MONSEL_MaterialScatterModel(Fe)
    FeMSMDeep.addScatterMechanism(FeNISTMott)
    FeMSMDeep.addScatterMechanism(FeDS)
    FeMSMDeep.setCSD(FeCSD)
    FeMSMDeep.setBarrierSM(FeBarrier)
    FeMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(FeMSM,FeMSMDeep)

def PtModel():
#	Platinum
    density = 21450.
#plasmonE = 9.11
    workfun = 5.64
    EFermi = 10.6
    potU = -workfun-EFermi # 
    Pt = mon.SEmaterial([epq.Element.Pt],[1.],density,"Platinum")
    PtWorkfunction=epq.ToSI.eV(workfun)
    Pt.setWorkfunction(PtWorkfunction)
    Pt.setEnergyCBbottom(epq.ToSI.eV(potU))
    Pt.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
    tablePath = tablePathPt
    PtTables = [tablePath+"IIMFPPennInterpPtSI.csv",tablePath+"interpNUSimReducedDeltaEPtSI.csv",tablePath+"interpsimTableThetaNUPtSI.csv",tablePath+"interpSimESE0NUPtSI.csv"]
    PtCoreEnergieseV = [2122.]
    for en in PtCoreEnergieseV:
	Pt.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Pt.getDensity()/epq.Element.Pt.getMass()
#Pt.addBindingEnergy(epq.ToSI.eV(0.)+PtWorkfunction,nve*density1electron)
# Create scatter mechanisms
    PtNISTMott = mon.SelectableElasticSM(Pt,mon.NU_ELSEPA_DCS.FactoryMT)
    PtDS = mon.TabulatedInelasticSM(Pt,3,PtTables)
#PtDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#PtMoller = mon.MollerInelasticSM(Pt)
#PtPlasmon = mon.KoteraPlasmonInelasticSM(Pt,1.)
    PtBarrier = mon.ExpQMBarrierSM(Pt)
    PtCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    PtMSM = mon.MONSEL_MaterialScatterModel(Pt)
    PtMSM.addScatterMechanism(PtNISTMott)
    PtMSM.addScatterMechanism(PtDS)
    PtMSM.setCSD(PtCSD)
    PtMSM.setBarrierSM(PtBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    PtMSMDeep = mon.MONSEL_MaterialScatterModel(Pt)
    PtMSMDeep.addScatterMechanism(PtNISTMott)
    PtMSMDeep.addScatterMechanism(PtDS)
    PtMSMDeep.setCSD(PtCSD)
    PtMSMDeep.setBarrierSM(PtBarrier)
    PtMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(PtMSM,PtMSMDeep)

def TaModel():
#	Tantalum
    density = 16690.
#plasmonE = 9.11
    workfun = 4.25
    EFermi = 8.4
    potU = -workfun-EFermi # 
    Ta = mon.SEmaterial([epq.Element.Ta],[1.],density,"Tantalum")
    TaWorkfunction=epq.ToSI.eV(workfun)
    Ta.setWorkfunction(TaWorkfunction)
    Ta.setEnergyCBbottom(epq.ToSI.eV(potU))
    Ta.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
    tablePath = tablePathTa
    TaTables = [tablePath+"IIMFPPennInterpTaSI.csv",tablePath+"interpNUSimReducedDeltaETaSI.csv",tablePath+"interpsimTableThetaNUTaSI.csv",tablePath+"interpSimESE0NUTaSI.csv"]
    TaCoreEnergieseV = [1735.]
    for en in TaCoreEnergieseV:
	Ta.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ta.getDensity()/epq.Element.Ta.getMass()
#Ta.addBindingEnergy(epq.ToSI.eV(0.)+TaWorkfunction,nve*density1electron)
# Create scatter mechanisms
    TaNISTMott = mon.SelectableElasticSM(Ta,mon.NU_ELSEPA_DCS.FactoryMT)
    TaDS = mon.TabulatedInelasticSM(Ta,3,TaTables)
#TaDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#TaMoller = mon.MollerInelasticSM(Ta)
#TaPlasmon = mon.KoteraPlasmonInelasticSM(Ta,1.)
    TaBarrier = mon.ExpQMBarrierSM(Ta)
    TaCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    TaMSM = mon.MONSEL_MaterialScatterModel(Ta)
    TaMSM.addScatterMechanism(TaNISTMott)
    TaMSM.addScatterMechanism(TaDS)
    TaMSM.setCSD(TaCSD)
    TaMSM.setBarrierSM(TaBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    TaMSMDeep = mon.MONSEL_MaterialScatterModel(Ta)
    TaMSMDeep.addScatterMechanism(TaNISTMott)
    TaMSMDeep.addScatterMechanism(TaDS)
    TaMSMDeep.setCSD(TaCSD)
    TaMSMDeep.setBarrierSM(TaBarrier)
    TaMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(TaMSM,TaMSMDeep)

def TiModel():
#	Titanium
    density = 4506.
#plasmonE = 9.11
    workfun = 4.33
    EFermi = 6.
    potU = -workfun-EFermi # 
    Ti = mon.SEmaterial([epq.Element.Ti],[1.],density,"Titanium")
    TiWorkfunction=epq.ToSI.eV(workfun)
    Ti.setWorkfunction(TiWorkfunction)
    Ti.setEnergyCBbottom(epq.ToSI.eV(potU))
    Ti.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
    tablePath = tablePathTi
    TiTables = [tablePath+"IIMFPPennInterpTiSI.csv",tablePath+"interpNUSimReducedDeltaETiSI.csv",tablePath+"interpsimTableThetaNUTiSI.csv",tablePath+"interpSimESE0NUTiSI.csv"]
    TiCoreEnergieseV = [460.2, 4966.]
    for en in TiCoreEnergieseV:
	Ti.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ti.getDensity()/epq.Element.Ti.getMass()
#Ti.addBindingEnergy(epq.ToSI.eV(0.)+TiWorkfunction,nve*density1electron)
# Create scatter mechanisms
    TiNISTMott = mon.SelectableElasticSM(Ti,mon.NU_ELSEPA_DCS.FactoryMT100Lin)
    TiDS = mon.TabulatedInelasticSM(Ti,3,TiTables)
#TiDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#AgMoller = mon.MollerInelasticSM(Ag)
#AgPlasmon = mon.KoteraPlasmonInelasticSM(Ag,1.)
    TiBarrier = mon.ExpQMBarrierSM(Ti)
    TiCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    TiMSM = mon.MONSEL_MaterialScatterModel(Ti)
    TiMSM.addScatterMechanism(TiNISTMott)
    TiMSM.addScatterMechanism(TiDS)
    TiMSM.setCSD(TiCSD)
    TiMSM.setBarrierSM(TiBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    TiMSMDeep = mon.MONSEL_MaterialScatterModel(Ti)
    TiMSMDeep.addScatterMechanism(TiNISTMott)
    TiMSMDeep.addScatterMechanism(TiDS)
    TiMSMDeep.setCSD(TiCSD)
    TiMSMDeep.setBarrierSM(TiBarrier)
    TiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return (TiMSM,TiMSMDeep)

def WhiteSnModel():
#	whiteSn
    density = 7265.
#plasmonE = 9.11
    workfun = 4.42
    EFermi = 10.2
    potU = -workfun-EFermi # 
    whiteSn = mon.SEmaterial([epq.Element.Sn],[1.],density,"white Tin")
    whiteSnWorkfunction=epq.ToSI.eV(workfun)
    whiteSn.setWorkfunction(whiteSnWorkfunction)
    whiteSn.setEnergyCBbottom(epq.ToSI.eV(potU))
    whiteSn.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
    tablePath =tablePathwhiteSn
    whiteSnTables = [tablePath+"IIMFPPennInterpwhiteSnSI.csv",tablePath+"interpNUSimReducedDeltaEwhiteSnSI.csv",tablePath+"interpsimTableThetaNUwhiteSnSI.csv",tablePath+"interpSimESE0NUwhiteSnSI.csv"]
    whiteSnCoreEnergieseV = [23.9,493.2,3929.,29150.]
    for en in whiteSnCoreEnergieseV:
	whiteSn.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = whiteSn.getDensity()/epq.Element.whiteSn.getMass()
#whiteSn.addBindingEnergy(epq.ToSI.eV(0.)+whiteSnWorkfunction,nve*density1electron)
# Create scatter mechanisms
    whiteSnNISTMott = mon.SelectableElasticSM(whiteSn,mon.NU_ELSEPA_DCS.FactoryMT)
    whiteSnDS = mon.TabulatedInelasticSM(whiteSn,3,whiteSnTables)
#whiteSnDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#whiteSnMoller = mon.MollerInelasticSM(whiteSn)
#whiteSnPlasmon = mon.KoteraPlasmonInelasticSM(whiteSn,1.)
    whiteSnBarrier = mon.ExpQMBarrierSM(whiteSn)
    whiteSnCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    whiteSnMSM = mon.MONSEL_MaterialScatterModel(whiteSn)
    whiteSnMSM.addScatterMechanism(whiteSnNISTMott)
    whiteSnMSM.addScatterMechanism(whiteSnDS)
    whiteSnMSM.setCSD(whiteSnCSD)
    whiteSnMSM.setBarrierSM(whiteSnBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    whiteSnMSMDeep = mon.MONSEL_MaterialScatterModel(whiteSn)
    whiteSnMSMDeep.addScatterMechanism(whiteSnNISTMott)
    whiteSnMSMDeep.addScatterMechanism(whiteSnDS)
    whiteSnMSMDeep.setCSD(whiteSnCSD)
    whiteSnMSMDeep.setBarrierSM(whiteSnBarrier)
    whiteSnMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(whiteSnMSM,whiteSnMSMDeep)

def SiCSDModel():
# Si CSD (from John V)
    phononE = 0.063 # I've seen the number reported as 510 cm^-1. this is conversion of that to eV.
    phononStrength = 3. # Turner & Inkson dispersion curves appear to show 3 LO phonon modes converging to the same
# energy at the Gamma point. 
    density = 2330.
    workfun = 4.85
    bandgap = 1.1 # width of band gap in eV
    EFermi = -bandgap # This puts the Fermi level at the top of the valence band.
    potU = -workfun-EFermi
    Si = mon.SEmaterial([epq.Element.Si],[1.],density,"Silicon CSD")
    SiWorkfunction=epq.ToSI.eV(workfun)
    Si.setWorkfunction(SiWorkfunction)
    Si.setEnergyCBbottom(epq.ToSI.eV(potU))
    Si.setBandgap(epq.ToSI.eV(bandgap))
    Si.setCoreEnergy([epq.ToSI.eV(99.9),epq.ToSI.eV(1839.)])

# Create scatter mechanisms
    SiNISTMott = mon.SelectableElasticSM(Si,mon.NU_ELSEPA_DCS.FactoryMT100Lin)
    SiCSD = mon.JoyLuoNieminenCSD(Si,epq.ToSI.eV(45.))                 
    SifittedInel= mon.FittedInelSM(Si,epq.ToSI.eV(50.),SiCSD)         

# MSM to be used in thin layer (includes SE generation)
    SiMSM = mon.MONSEL_MaterialScatterModel(Si)
    SiMSM.addScatterMechanism(SiNISTMott)
    SiMSM.setCSD(SiCSD)                             
    SiMSM.addScatterMechanism(SifittedInel)         

# MSM to be used deep inside (drops electrons with E<50 eV)
    SiMSMDeep = mon.MONSEL_MaterialScatterModel(Si)
    SiMSMDeep.addScatterMechanism(SiNISTMott)
    SiMSMDeep.setCSD(SiCSD)   # tells the material scatter model to use your defined CSD
    SiMSMDeep.addScatterMechanism(SifittedInel)      #Inelastic scattering
#SiMSMDeep.addScatterMechanism(SiDS)         #In-elastic scattering. Table
#SiMSMDeep.addScatterMechanism(Siphonon)      #Phonon scattering
#SiMSMDeep.setBarrierSM(SiAbruptBarrier)     #Omitting this line causes the barrier to default to gradual/classical

    SiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))	#default
    return(SiMSM,SiMSMDeep)

def CuCSDModel():
# Cu CSD
#             Copper CSD
    density = 8933.
    nve = 11
#plasmonE = 9.11
    workfun = 4.65
    EFermi = 8.7
    potU = -workfun-EFermi # Assumes Cu Fermi energy is 8.7 eV
    Cu = mon.SEmaterial([epq.Element.Cu],[1.],density,"Copper CSD")
    CuWorkfunction=epq.ToSI.eV(workfun)
    Cu.setWorkfunction(CuWorkfunction)
    Cu.setEnergyCBbottom(epq.ToSI.eV(potU))
    Cu.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#if (WinOrLin<>2): 
	#tablePath = tablePathOuter + "\CuTables"+PathSep									#If Windows:      tablePath = tablePathOuter + "\CuTables"+PathSep				#If LINUX:      tablePath = tablePathOuter + PathSep+"CuTables"+PathSep
#if (WinOrLin==2): 
	#tablePath = tablePathOuter + PathSep+"CuTables"+PathSep
#CuTables = [tablePath+"IIMFPPennInterpCuSI.csv",tablePath+"interpNUSimReducedDeltaECuSI.csv",tablePath+"interpsimTableThetaNUCuSI.csv",tablePath+"interpSimESE0NUCuSI.csv"]
    Cu.setCoreEnergy([epq.ToSI.eV(75.1),epq.ToSI.eV(77.3),epq.ToSI.eV(122.5),epq.ToSI.eV(932.7),epq.ToSI.eV(1096.7),epq.ToSI.eV(8979.)])
#density1electron = Cu.getDensity()/epq.Element.Cu.getMass()
#Cu.addBindingEnergy(epq.ToSI.eV(0.)+CuWorkfunction,nve*density1electron)
# Create scatter mechanisms
    CuNISTMott = mon.SelectableElasticSM(Cu,mon.NISTMottRS.Factory)
#CuDS = mon.TabulatedInelasticSM(Cu,3,CuTables)
#CuMoller = mon.MollerInelasticSM(Cu)
#CuPlasmon = mon.KoteraPlasmonInelasticSM(Cu,1.)
#CuBarrier = mon.ExpQMBarrierSM(Cu)
    CuCSD = mon.JoyLuoNieminenCSD(Cu,epq.ToSI.eV(45.))
    CufittedInel= mon.FittedInelSM(Cu,epq.ToSI.eV(50.),CuCSD)         
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    CuMSM = mon.MONSEL_MaterialScatterModel(Cu)
    CuMSM.addScatterMechanism(CuNISTMott)
    CuMSM.addScatterMechanism(CufittedInel)         
#CuMSM.addScatterMechanism(CuDS)
    CuMSM.setCSD(CuCSD)
#CuMSM.setBarrierSM(CuBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    CuMSMDeep = mon.MONSEL_MaterialScatterModel(Cu)
    CuMSMDeep.addScatterMechanism(CuNISTMott)
#CuMSMDeep.addScatterMechanism(CuDS)
    CuMSMDeep.setCSD(CuCSD)
#CuMSMDeep.setBarrierSM(CuBarrier)
    CuMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(CuMSM,CuMSMDeep)

def WCSDModel():
#    Tungsten CSD
    density = 19300.
#plasmonE = 9.11
    workfun = 4.55
    EFermi = 10.1
    potU = -workfun-EFermi 
    W = mon.SEmaterial([epq.Element.W],[1.],density,"Tungsten CSD")
    WWorkfunction=epq.ToSI.eV(workfun)
    W.setWorkfunction(WWorkfunction)
    W.setEnergyCBbottom(epq.ToSI.eV(potU))
#tablePath = "C:\NIST\JMONSEL\ScatteringTables\WTables"+PathSep
#WTables = [tablePath+"IIMFPPennInterpWSI.csv",tablePath+"interpNUSimReducedDeltaEWSI.csv",tablePath+"interpsimTableThetaNUWSI.csv",tablePath+"interpSimESE0NUWSI.csv"]
    coreEnergies = [31.4, 33.6, 36.8, 45.3, 75.6, 243.5, 255.9, 423.6, 490.4, 594.1, 1809., 1949., 2281., 2575., 2820., 10207., 11544., 12100., 69525.]
    for i in range(len(coreEnergies)):
        coreEnergies[i] = epq.ToSI.eV(coreEnergies[i])
    W.setCoreEnergy(coreEnergies)
    density1electron = W.getDensity()/epq.Element.W.getMass()
#W.addBindingEnergy(epq.ToSI.eV(0.)+WWorkfunction,nve*density1electron)
# Create scatter mechanisms
    WNISTMott = mon.SelectableElasticSM(W,mon.NISTMottRS.Factory)
#WDS = mon.TabulatedInelasticSM(W,3,WTables)
#WMoller = mon.MollerInelasticSM(W)
#WPlasmon = mon.KoteraPlasmonInelasticSM(W,1.)
#WBarrier = mon.ExpQMBarrierSM(W)
    WCSD = mon.JoyLuoNieminenCSD(W,epq.ToSI.eV(45.))
    WfittedInel= mon.FittedInelSM(W,epq.ToSI.eV(50.),WCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    WMSM = mon.MONSEL_MaterialScatterModel(W)
    WMSM.addScatterMechanism(WNISTMott)
#WMSM.addScatterMechanism(WDS)
    WMSM.setCSD(WCSD)
    WMSM.addScatterMechanism(WfittedInel)
#WMSM.setBarrierSM(WBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    WMSMDeep = mon.MONSEL_MaterialScatterModel(W)
    WMSMDeep.addScatterMechanism(WNISTMott)
    #WMSMDeep.addScatterMechanism(WDS)
    WMSMDeep.setCSD(WCSD)
    WMSMDeep.addScatterMechanism(WfittedInel)
    #WMSMDeep.setBarrierSM(WBarrier)
    WMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(WMSM,WMSMDeep)

def AuCSDModel():
#	Gold CSD
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
    Au = mon.SEmaterial([epq.Element.Au],[1.],density,"Gold CSD")
    AuWorkfunction=epq.ToSI.eV(workfun)
    Au.setWorkfunction(AuWorkfunction)
    Au.setEnergyCBbottom(epq.ToSI.eV(potU))
    Au.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#if (WinOrLin<>2): 
	#tablePath = tablePathOuter + "\AuTables"+PathSep							#If Windows:  tablePath = tablePathOuter + "\AuTables"+PathSep			If LINUX:      tablePath = tablePathOuter + PathSep+"AuTables"+PathSep
#if (WinOrLin==2): 
	#tablePath = tablePathOuter + PathSep+"AuTables"+PathSep
#AuTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
        762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    for en in AuCoreEnergieseV:
        Au.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Au.getDensity()/epq.Element.Au.getMass()
#Au.addBindingEnergy(epq.ToSI.eV(0.)+AuWorkfunction,nve*density1electron)
# Create scatter mechanisms
    AuNISTMott = mon.SelectableElasticSM(Au,mon.NISTMottRS.Factory)
#AuDS = mon.TabulatedInelasticSM(Au,3,AuTables)
#AuMoller = mon.MollerInelasticSM(Au)
#AuPlasmon = mon.KoteraPlasmonInelasticSM(Au,1.)
#AuBarrier = mon.ExpQMBarrierSM(Au,0.05e-9)
    AuCSD = mon.JoyLuoNieminenCSD(Au,epq.ToSI.eV(45.)) 
    AufittedInel= mon.FittedInelSM(Au,epq.ToSI.eV(50.),AuCSD)                         
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    AuMSM = mon.MONSEL_MaterialScatterModel(Au)
    AuMSM.addScatterMechanism(AuNISTMott)
#AuMSM.addScatterMechanism(AuDS)
    AuMSM.setCSD(AuCSD)
#AuMSM.setBarrierSM(AuBarrier)
    AuMSM.addScatterMechanism(AufittedInel)         
# MSM to be used deep inside (drops electrons with E<50 eV)
    AuMSMDeep = mon.MONSEL_MaterialScatterModel(Au)
    AuMSMDeep.addScatterMechanism(AuNISTMott)
#AuMSMDeep.addScatterMechanism(AuDS)
    AuMSMDeep.setCSD(AuCSD)
#AuMSMDeep.setBarrierSM(AuBarrier)
    AuMSMDeep.addScatterMechanism(AufittedInel)
    AuMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(AuMSM,AuMSMDeep)

def SiO2CSDModel():
#SiO2 CSD
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
    SiO2 = mon.SEmaterial([Si,Ox],[SiWeight/totalWeight,OxWeight/totalWeight],density,"Silicon Dioxide CSD")
    SiO2.setEpsr(3.9)
    SiO2Workfunction=epq.ToSI.eV(workfun)
    SiO2.setWorkfunction(SiO2Workfunction)
    SiO2.setBandgap(epq.ToSI.eV(bandgap))
    SiO2.setEnergyCBbottom(epq.ToSI.eV(potU))
    SiO2.setCoreEnergy([epq.ToSI.eV(41.6),epq.ToSI.eV(99.2),epq.ToSI.eV(99.8),epq.ToSI.eV(149.7),epq.ToSI.eV(543.1),epq.ToSI.eV(1839.)])
#if (WinOrLin<>2): 
	#tablePath = tablePathOuter + "\SiO2Tables"+PathSep									#If Windows:      tablePath = tablePathOuter + "\SiO2Tables"+PathSep					#If LINUX:        tablePath = tablePathOuter + PathSep+"SiO2Tables"+PathSep
#if (WinOrLin==2): 
	#tablePath = tablePathOuter + PathSep+"SiO2Tables"+PathSep
#SiO2Tables = [tablePath+"IIMFPPennInterpSiO2SI.csv",tablePath+"interpNUSimReducedDeltaESiO2SI.csv",tablePath+"interpsimTableThetaNUSiO2SI.csv",tablePath+"interpSimESE0NUSiO2SI.csv"]
# Create scatter mechanisms
    SiO2NISTMott = mon.SelectableElasticSM(SiO2,mon.NISTMottRS.Factory)
#SiO2DS = mon.TabulatedInelasticSM(SiO2,3,SiO2Tables,epq.ToSI.eV(20.+bandgap))
#SiO2phonon = mon.GanachaudMokraniPhononInelasticSM(phononStrength,epq.ToSI.eV(phononE),300.,3.82,1.)
#SiO2polaron = mon.GanachaudMokraniPolaronTrapSM(1.0e9,1./epq.ToSI.eV(1.))
#SiO2Barrier = mon.ExpQMBarrierSM(SiO2,1.e-9)
    SiO2CSD = mon.JoyLuoNieminenCSD(SiO2,epq.ToSI.eV(45.))                                        # The default. No need to actually execute this line.
    SiO2fittedInel= mon.FittedInelSM(SiO2,epq.ToSI.eV(50.),SiO2CSD)         
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    SiO2MSM = mon.MONSEL_MaterialScatterModel(SiO2)
    SiO2MSM.addScatterMechanism(SiO2NISTMott)
#SiO2MSM.addScatterMechanism(SiO2DS)
#SiO2MSM.addScatterMechanism(SiO2phonon)
# SiO2MSM.addScatterMechanism(SiO2polaron)
    SiO2MSM.setCSD(SiO2CSD)
    SiO2MSM.addScatterMechanism(SiO2fittedInel)         
#SiO2MSM.setBarrierSM(SiO2Barrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    SiO2MSMDeep = mon.MONSEL_MaterialScatterModel(SiO2)
    SiO2MSMDeep.addScatterMechanism(SiO2NISTMott)
#SiO2MSMDeep.addScatterMechanism(SiO2DS)
#SiO2MSMDeep.addScatterMechanism(SiO2phonon)
    SiO2MSMDeep.setCSD(SiO2CSD)
    SiO2MSMDeep.addScatterMechanism(SiO2fittedInel)
#SiO2MSMDeep.setBarrierSM(SiO2Barrier)
    SiO2MSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(SiO2MSM,SiO2MSMDeep)

def AlCSDModel():
#	Aluminum CSD
    density = 2700.
#plasmonE = 9.11
    workfun = 4.28
    EFermi = 11.2
    potU = -workfun-EFermi # 
    Al = mon.SEmaterial([epq.Element.Al],[1.],density,"Aluminum CSD")
    AlWorkfunction=epq.ToSI.eV(workfun)
    Al.setWorkfunction(AlWorkfunction)
    Al.setEnergyCBbottom(epq.ToSI.eV(potU))
    Al.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#tablePath = "C:\Program Files\NIST\JMONSEL\ScatteringTables\AlTables"+PathSep
#AlTables = [tablePath+"IIMFPPennInterpAlSI.csv",tablePath+"interpNUSimReducedDeltaEAlSI.csv",tablePath+"interpsimTableThetaNUAlSI.csv",tablePath+"interpSimESE0NUAlSI.csv"]
    AlCoreEnergieseV = [72.7, 1559.]
    for en in AlCoreEnergieseV:
	Al.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Al.getDensity()/epq.Element.Al.getMass()
#Al.addBindingEnergy(epq.ToSI.eV(0.)+AlWorkfunction,nve*density1electron)
# Create scatter mechanisms
    AlNISTMott = mon.SelectableElasticSM(Al,mon.NU_ELSEPA_DCS.FactoryMT)
#AlDS = mon.TabulatedInelasticSM(Al,3,AlTables)
#AlDS.setBranchingRatios([0.3243, 0.800729])
#AlMoller = mon.MollerInelasticSM(Al)
#AlPlasmon = mon.KoteraPlasmonInelasticSM(Al,1.)
#AlBarrier = mon.ExpQMBarrierSM(Al)
    AlCSD = mon.JoyLuoNieminenCSD(Al,epq.ToSI.eV(45.))
    AlfittedInel= mon.FittedInelSM(Al,epq.ToSI.eV(50.),AlCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    AlMSM = mon.MONSEL_MaterialScatterModel(Al)
    AlMSM.addScatterMechanism(AlNISTMott)
    AlMSM.setCSD(AlCSD)
    AlMSM.addScatterMechanism(AlfittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    AlMSMDeep = mon.MONSEL_MaterialScatterModel(Al)
    AlMSMDeep.addScatterMechanism(AlNISTMott)
    AlMSMDeep.setCSD(AlCSD)
    AlMSMDeep.addScatterMechanism(AlfittedInel)
    AlMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(AlMSM,AlMSMDeep)

def TiCSDModel():
#	Titanium CSD
    density = 4506.
#plasmonE = 9.11
    workfun = 4.33
    EFermi = 6.
    potU = -workfun-EFermi # 
    Ti = mon.SEmaterial([epq.Element.Ti],[1.],density,"Titanium CSD")
    TiWorkfunction=epq.ToSI.eV(workfun)
    Ti.setWorkfunction(TiWorkfunction)
    Ti.setEnergyCBbottom(epq.ToSI.eV(potU))
    Ti.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#tablePath = "C:\NIST\JMONSEL\ScatteringTables\TiTables"+PathSep
#TiTables = [tablePath+"IIMFPPennInterpTiSI.csv",tablePath+"interpNUSimReducedDeltaETiSI.csv",tablePath+"interpsimTableThetaNUTiSI.csv",tablePath+"interpSimESE0NUTiSI.csv"]
    TiCoreEnergieseV = [460.2, 4966.]
    for en in TiCoreEnergieseV:
	Ti.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ti.getDensity()/epq.Element.Ti.getMass()
#Ti.addBindingEnergy(epq.ToSI.eV(0.)+TiWorkfunction,nve*density1electron)
# Create scatter mechanisms
    TiNISTMott = mon.SelectableElasticSM(Ti,mon.NU_ELSEPA_DCS.FactoryMT100Lin)
#TiDS = mon.TabulatedInelasticSM(Ti,3,TiTables)
#TiDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#AgMoller = mon.MollerInelasticSM(Ag)
#AgPlasmon = mon.KoteraPlasmonInelasticSM(Ag,1.)
#TiBarrier = mon.ExpQMBarrierSM(Ti)
    TiCSD = mon.JoyLuoNieminenCSD(Ti,epq.ToSI.eV(45.))
    TifittedInel= mon.FittedInelSM(Ti,epq.ToSI.eV(50.),TiCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    TiMSM = mon.MONSEL_MaterialScatterModel(Ti)
    TiMSM.addScatterMechanism(TiNISTMott)
    TiMSM.setCSD(TiCSD)
    TiMSM.addScatterMechanism(TifittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    TiMSMDeep = mon.MONSEL_MaterialScatterModel(Ti)
    TiMSMDeep.addScatterMechanism(TiNISTMott)
    TiMSMDeep.setCSD(TiCSD)
    TiMSMDeep.addScatterMechanism(TifittedInel)
    TiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(TiMSM,TiMSMDeep)

def FeCSDModel():
#	Iron CSD
    density = 7874.
#plasmonE = 9.11
    workfun = 4.67
    EFermi = 8.9
    potU = -workfun-EFermi # 
    Fe = mon.SEmaterial([epq.Element.Fe],[1.],density,"Iron CSD")
    FeWorkfunction=epq.ToSI.eV(workfun)
    Fe.setWorkfunction(FeWorkfunction)
    Fe.setEnergyCBbottom(epq.ToSI.eV(potU))
    Fe.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#tablePath = "C:\NIST\JMONSEL\ScatteringTables\FeTables"+PathSep
#FeTables = [tablePath+"IIMFPPennInterpFeSI.csv",tablePath+"interpNUSimReducedDeltaEFeSI.csv",tablePath+"interpsimTableThetaNUFeSI.csv",tablePath+"interpSimESE0NUFeSI.csv"]
    FeCoreEnergieseV = [710., 7112.]
    for en in FeCoreEnergieseV:
	Fe.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Fe.getDensity()/epq.Element.Fe.getMass()
#Fe.addBindingEnergy(epq.ToSI.eV(0.)+FeWorkfunction,nve*density1electron)
# Create scatter mechanisms
    FeNISTMott = mon.SelectableElasticSM(Fe,mon.NU_ELSEPA_DCS.FactoryMT)
#FeDS = mon.TabulatedInelasticSM(Fe,3,FeTables)
#FeDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#FeMoller = mon.MollerInelasticSM(Fe)
#FePlasmon = mon.KoteraPlasmonInelasticSM(Fe,1.)
#FeBarrier = mon.ExpQMBarrierSM(Fe)
    FeCSD = mon.JoyLuoNieminenCSD(Fe,epq.ToSI.eV(45.))
    FefittedInel= mon.FittedInelSM(Fe,epq.ToSI.eV(50.),FeCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    FeMSM = mon.MONSEL_MaterialScatterModel(Fe)
    FeMSM.addScatterMechanism(FeNISTMott)
    FeMSM.setCSD(FeCSD)
    FeMSM.addScatterMechanism(FefittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    FeMSMDeep = mon.MONSEL_MaterialScatterModel(Fe)
    FeMSMDeep.addScatterMechanism(FeNISTMott)
    FeMSMDeep.setCSD(FeCSD)
    FeMSMDeep.addScatterMechanism(FefittedInel)
    FeMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(FeMSM,FeMSMDeep)

def MoCSDModel():
#	Molybdenum CSD
    density = 10280.
#plasmonE = 23.09
    workfun = 5.1
    EFermi = 6.5
    potU = -workfun-EFermi # 
    Mo = mon.SEmaterial([epq.Element.Mo],[1.],density,"Molybdenum CSD")
    MoWorkfunction=epq.ToSI.eV(workfun)
    Mo.setWorkfunction(MoWorkfunction)
    Mo.setEnergyCBbottom(epq.ToSI.eV(potU))
#tablePath = "C:\Program Files\NIST\JMONSEL\FPAScatteringTables\MoTables"+PathSep
#MoTables = [tablePath+"iimfp.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
    MoCoreEnergieseV = [227.9, 2520., 20000.]
    for en in MoCoreEnergieseV:
	Mo.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Mo.getDensity()/epq.Element.Mo.getMass()
#Mo.addBindingEnergy(epq.ToSI.eV(0.)+MoWorkfunction,nve*density1electron)
# Create scatter mechanisms
    MoNISTMott = mon.SelectableElasticSM(Mo,mon.NU_ELSEPA_DCS.FactoryMT)
#MoDS = mon.TabulatedInelasticSM(Mo,3,MoTables)
#MoDS.setBranchingRatios([0.78, 0.62, 0.441262, 0.199441])
#MoMoller = mon.MollerInelasticSM(Mo)
#MoPlasmon = mon.KoteraPlasmonInelasticSM(Mo,1.)
#MoBarrier = mon.ExpQMBarrierSM(Mo)
    MoCSD = mon.JoyLuoNieminenCSD(Mo,epq.ToSI.eV(45.))
    MofittedInel= mon.FittedInelSM(Mo,epq.ToSI.eV(50.),MoCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    MoMSM = mon.MONSEL_MaterialScatterModel(Mo)
    MoMSM.addScatterMechanism(MoNISTMott)
    MoMSM.setCSD(MoCSD)
    MoMSM.addScatterMechanism(MofittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    MoMSMDeep = mon.MONSEL_MaterialScatterModel(Mo)
    MoMSMDeep.addScatterMechanism(MoNISTMott)
    MoMSMDeep.setCSD(MoCSD)
    MoMSMDeep.addScatterMechanism(MofittedInel)
    MoMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(MoMSM,MoMSMDeep)

def PtCSDModel():
#	Platinum CSD
    density = 21450.
#plasmonE = 9.11
    workfun = 5.64
    EFermi = 10.6
    potU = -workfun-EFermi # 
    Pt = mon.SEmaterial([epq.Element.Pt],[1.],density,"Platinum CSD")
    PtWorkfunction=epq.ToSI.eV(workfun)
    Pt.setWorkfunction(PtWorkfunction)
    Pt.setEnergyCBbottom(epq.ToSI.eV(potU))
    Pt.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#tablePath = "C:\NIST\JMONSEL\ScatteringTables\PtTables"+PathSep
#PtTables = [tablePath+"IIMFPPennInterpPtSI.csv",tablePath+"interpNUSimReducedDeltaEPtSI.csv",tablePath+"interpsimTableThetaNUPtSI.csv",tablePath+"interpSimESE0NUPtSI.csv"]
    PtCoreEnergieseV = [2122.]
    for en in PtCoreEnergieseV:
	Pt.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Pt.getDensity()/epq.Element.Pt.getMass()
#Pt.addBindingEnergy(epq.ToSI.eV(0.)+PtWorkfunction,nve*density1electron)
# Create scatter mechanisms
    PtNISTMott = mon.SelectableElasticSM(Pt,mon.NU_ELSEPA_DCS.FactoryMT)
#PtDS = mon.TabulatedInelasticSM(Pt,3,PtTables)
#PtDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#PtMoller = mon.MollerInelasticSM(Pt)
#PtPlasmon = mon.KoteraPlasmonInelasticSM(Pt,1.)
#PtBarrier = mon.ExpQMBarrierSM(Pt)
    PtCSD = mon.JoyLuoNieminenCSD(Pt,epq.ToSI.eV(45.))
    PtfittedInel= mon.FittedInelSM(Pt,epq.ToSI.eV(50.),PtCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    PtMSM = mon.MONSEL_MaterialScatterModel(Pt)
    PtMSM.addScatterMechanism(PtNISTMott)
    PtMSM.setCSD(PtCSD)
    PtMSM.addScatterMechanism(PtfittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    PtMSMDeep = mon.MONSEL_MaterialScatterModel(Pt)
    PtMSMDeep.addScatterMechanism(PtNISTMott)
    PtMSMDeep.setCSD(PtCSD)
    PtMSMDeep.addScatterMechanism(PtfittedInel)
    PtMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(PtMSM,PtMSMDeep)

def RuCSDModel():
#	Ruthenium CSD
    density = 12410.
#plasmonE = 28.54
    workfun = 4.7
    EFermi = 6.9
    potU = -workfun-EFermi # 
    Ru = mon.SEmaterial([epq.Element.Ru],[1.],density,"Ruthenium CSD")
    RuWorkfunction=epq.ToSI.eV(workfun)
    Ru.setWorkfunction(RuWorkfunction)
    Ru.setEnergyCBbottom(epq.ToSI.eV(potU))
#tablePath = "C:\Program Files\NIST\JMONSEL\FPAScatteringTables\RuTables"+PathSep
#RuTables = [tablePath+"iimfp.csv",tablePath+"dered.csv",tablePath+"thscat.csv",tablePath+"ese0.csv"]
    RuCoreEnergieseV = [284.2, 2838., 22117.]
    for en in RuCoreEnergieseV:
	Ru.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ru.getDensity()/epq.Element.Ru.getMass()
#Ru.addBindingEnergy(epq.ToSI.eV(0.)+RuWorkfunction,nve*density1electron)
# Create scatter mechanisms
    RuNISTMott = mon.SelectableElasticSM(Ru,mon.NU_ELSEPA_DCS.FactoryMT)
#RuDS = mon.TabulatedInelasticSM(Ru,3,RuTables)
#RuDS.setBranchingRatios([0.66, 0.318851, 0.167881])
#RuMoller = mon.MollerInelasticSM(Ru)
#RuPlasmon = mon.KoteraPlasmonInelasticSM(Ru,1.)
#RuBarrier = mon.ExpQMBarrierSM(Ru)
    RuCSD = mon.JoyLuoNieminenCSD(Ru,epq.ToSI.eV(45.))
    RufittedInel= mon.FittedInelSM(Ru,epq.ToSI.eV(50.),RuCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    RuMSM = mon.MONSEL_MaterialScatterModel(Ru)
    RuMSM.addScatterMechanism(RuNISTMott)
    RuMSM.setCSD(RuCSD)
    RuMSM.addScatterMechanism(RufittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    RuMSMDeep = mon.MONSEL_MaterialScatterModel(Ru)
    RuMSMDeep.addScatterMechanism(RuNISTMott)
    RuMSMDeep.setCSD(RuCSD)
    RuMSMDeep.addScatterMechanism(RufittedInel)
    RuMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(RuMSM,RuMSMDeep)

def TaCSDModel():
#	Tantalum CSD
    density = 16690.
#plasmonE = 9.11
    workfun = 4.25
    EFermi = 8.4
    potU = -workfun-EFermi # 
    Ta = mon.SEmaterial([epq.Element.Ta],[1.],density,"Tantalum CSD")
    TaWorkfunction=epq.ToSI.eV(workfun)
    Ta.setWorkfunction(TaWorkfunction)
    Ta.setEnergyCBbottom(epq.ToSI.eV(potU))
    Ta.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#tablePath = "C:\NIST\JMONSEL\ScatteringTables\TaTables"+PathSep
#TaTables = [tablePath+"IIMFPPennInterpTaSI.csv",tablePath+"interpNUSimReducedDeltaETaSI.csv",tablePath+"interpsimTableThetaNUTaSI.csv",tablePath+"interpSimESE0NUTaSI.csv"]
    TaCoreEnergieseV = [1735.]
    for en in TaCoreEnergieseV:
	Ta.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ta.getDensity()/epq.Element.Ta.getMass()
#Ta.addBindingEnergy(epq.ToSI.eV(0.)+TaWorkfunction,nve*density1electron)
# Create scatter mechanisms
    TaNISTMott = mon.SelectableElasticSM(Ta,mon.NU_ELSEPA_DCS.FactoryMT)
#TaDS = mon.TabulatedInelasticSM(Ta,3,TaTables)
#TaDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#TaMoller = mon.MollerInelasticSM(Ta)
#TaPlasmon = mon.KoteraPlasmonInelasticSM(Ta,1.)
#TaBarrier = mon.ExpQMBarrierSM(Ta)
    TaCSD = mon.JoyLuoNieminenCSD(Ta,epq.ToSI.eV(45.))
    TafittedInel= mon.FittedInelSM(Ta,epq.ToSI.eV(50.),TaCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    TaMSM = mon.MONSEL_MaterialScatterModel(Ta)
    TaMSM.addScatterMechanism(TaNISTMott)
    TaMSM.setCSD(TaCSD)
    TaMSM.addScatterMechanism(TafittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    TaMSMDeep = mon.MONSEL_MaterialScatterModel(Ta)
    TaMSMDeep.addScatterMechanism(TaNISTMott)
    TaMSMDeep.setCSD(TaCSD)
    TaMSMDeep.addScatterMechanism(TafittedInel)
    TaMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(TaMSM,TaMSMDeep)

def AgCSDModel():
#	Silver CSD
    density = 10501.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 7.2
    potU = -workfun-EFermi # 
    Ag = mon.SEmaterial([epq.Element.Ag],[1.],density,"Silver CSD")
    AgWorkfunction=epq.ToSI.eV(workfun)
    Ag.setWorkfunction(AgWorkfunction)
    Ag.setEnergyCBbottom(epq.ToSI.eV(potU))
    Ag.setEpsr(50.0) # a large number, to mimic a metal by using large dielectric constant
#tablePath = "C:\Program Files\NIST\JMONSEL\ScatteringTables\AgTables"+PathSep
#AgTables = [tablePath+"IIMFPPennInterpAgSI.csv",tablePath+"interpNUSimReducedDeltaEAgSI.csv",tablePath+"interpsimTableThetaNUAgSI.csv",tablePath+"interpSimESE0NUAgSI.csv"]
    AgCoreEnergieseV = [368.3, 25514.]
    for en in AgCoreEnergieseV:
	Ag.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ag.getDensity()/epq.Element.Ag.getMass()
#Ag.addBindingEnergy(epq.ToSI.eV(0.)+AgWorkfunction,nve*density1electron)
# Create scatter mechanisms
    AgNISTMott = mon.SelectableElasticSM(Ag,mon.NU_ELSEPA_DCS.FactoryMT)
#AgDS = mon.TabulatedInelasticSM(Ag,3,AgTables)
#AgDS.setBranchingRatios([0.964539, 0.207453, 0.909605, 0.838235, 0.528771, 0.909477, 0.172381])
#AgMoller = mon.MollerInelasticSM(Ag)
#AgPlasmon = mon.KoteraPlasmonInelasticSM(Ag,1.)
#AgBarrier = mon.ExpQMBarrierSM(Ag)
    AgCSD = mon.JoyLuoNieminenCSD(Ag,epq.ToSI.eV(45.))
    AgfittedInel= mon.FittedInelSM(Ag,epq.ToSI.eV(50.),AgCSD)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    AgMSM = mon.MONSEL_MaterialScatterModel(Ag)
    AgMSM.addScatterMechanism(AgNISTMott)
    AgMSM.setCSD(AgCSD)
    AgMSM.addScatterMechanism(AgfittedInel)
# MSM to be used deep inside (drops electrons with E<50 eV)
    AgMSMDeep = mon.MONSEL_MaterialScatterModel(Ag)
    AgMSMDeep.addScatterMechanism(AgNISTMott)
    AgMSMDeep.setCSD(AgCSD)
    AgMSMDeep.addScatterMechanism(AgfittedInel)
    AgMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(AgMSM,AgMSMDeep)

###### still need to check properties of these materials #######
def CrModel():
#	Chromium
    print("Warning. Chromium is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
    Cr = mon.SEmaterial([epq.Element.Cr],[1.],density,"Chromium")
    CrWorkfunction=epq.ToSI.eV(workfun)
    Cr.setWorkfunction(CrWorkfunction)
    Cr.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathCr
#    CrTables = [tablePath+"IIMFPPennInterpCrSI.csv",tablePath+"interpNUSimReducedDeltaECrSI.csv",tablePath+"interpsimTableThetaNUCrSI.csv",tablePath+"interpSimESE0NUCrSI.csv"]
    tablePath=tablePathAu
    CrTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]

    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    CrCoreEnergieseV=AuCoreEnergieseV
    for en in CrCoreEnergieseV:
	Cr.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Cr.getDensity()/epq.Element.Cr.getMass()
#Cr.addBindingEnergy(epq.ToSI.eV(0.)+CrWorkfunction,nve*density1electron)
# Create scatter mechanisms
    CrNISTMott = mon.SelectableElasticSM(Cr,mon.NISTMottRS.Factory)
    CrDS = mon.TabulatedInelasticSM(Cr,3,CrTables)
#    CrDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441]) 
    #CrMoller = mon.MollerInelasticSM(Cr)
    #CrPlasmon = mon.KoteraPlasmonInelasticSM(Cr,1.)
    CrBarrier = mon.ExpQMBarrierSM(Cr,0.05e-9)
    CrCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    CrMSM = mon.MONSEL_MaterialScatterModel(Cr)
    CrMSM.addScatterMechanism(CrNISTMott)
    CrMSM.addScatterMechanism(CrDS)
    CrMSM.setCSD(CrCSD)
    CrMSM.setBarrierSM(CrBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    CrMSMDeep = mon.MONSEL_MaterialScatterModel(Cr)
    CrMSMDeep.addScatterMechanism(CrNISTMott)
    CrMSMDeep.addScatterMechanism(CrDS)
    CrMSMDeep.setCSD(CrCSD)
    CrMSMDeep.setBarrierSM(CrBarrier)
    CrMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(CrMSM,CrMSMDeep)

def NiModel():
#	Nickel
    print("Warning. Nickel is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
    Ni = mon.SEmaterial([epq.Element.Ni],[1.],density,"Nickel")
    NiWorkfunction=epq.ToSI.eV(workfun)
    Ni.setWorkfunction(NiWorkfunction)
    Ni.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathNi
#    NiTables = [tablePath+"IIMFPPennInterpNiSI.csv",tablePath+"interpNUSimReducedDeltaENiSI.csv",tablePath+"interpsimTableThetaNUNiSI.csv",tablePath+"interpSimESE0NUNiSI.csv"]
    tablePath=tablePathAu
    NiTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    NiCoreEnergieseV = AuCoreEnergieseV
    for en in NiCoreEnergieseV:
	Ni.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ni.getDensity()/epq.Element.Ni.getMass()
#Ni.addBindingEnergy(epq.ToSI.eV(0.)+NiWorkfunction,nve*density1electron)
# Create scatter mechanisms
    NiNISTMott = mon.SelectableElasticSM(Ni,mon.NISTMottRS.Factory)
    NiDS = mon.TabulatedInelasticSM(Ni,3,NiTables)
#    NiDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441]) 
    #NiMoller = mon.MollerInelasticSM(Ni)
    #NiPlasmon = mon.KoteraPlasmonInelasticSM(Ni,1.)
    NiBarrier = mon.ExpQMBarrierSM(Ni,0.05e-9)
    NiCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    NiMSM = mon.MONSEL_MaterialScatterModel(Ni)
    NiMSM.addScatterMechanism(NiNISTMott)
    NiMSM.addScatterMechanism(NiDS)
    NiMSM.setCSD(NiCSD)
    NiMSM.setBarrierSM(NiBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    NiMSMDeep = mon.MONSEL_MaterialScatterModel(Ni)
    NiMSMDeep.addScatterMechanism(NiNISTMott)
    NiMSMDeep.addScatterMechanism(NiDS)
    NiMSMDeep.setCSD(NiCSD)
    NiMSMDeep.setBarrierSM(NiBarrier)
    NiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(NiMSM,NiMSMDeep)

def GaAsModel():
#	Gallium Aresenide
    print("Warning. Gallium Arsenide is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
#    GaAs = mon.SEmaterial([epq.Element.GaAs],[1.],density,"GalliumArsenide")
    GaAs = mon.SEmaterial([epq.Element.Au],[1.],density,"GalliumArsenide")
    GaAsWorkfunction=epq.ToSI.eV(workfun)
    GaAs.setWorkfunction(GaAsWorkfunction)
    GaAs.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathGaAs
#    GaAsTables = [tablePath+"IIMFPPenGaAsnterpGaAsSI.csv",tablePath+"interpNUSimReducedDeltaEGaAsSI.csv",tablePath+"interpsimTableThetaNUGaAsSI.csv",tablePath+"interpSimESE0NUGaAsSI.csv"]
    tablePath=tablePathAu
    GaAsTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    GaAsCoreEnergieseV = AuCoreEnergieseV
    for en in GaAsCoreEnergieseV:
	GaAs.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = GaAs.getDensity()/epq.Element.GaAs.getMass()
#GaAs.addBindingEnergy(epq.ToSI.eV(0.)+GaAsWorkfunction,nve*density1electron)
# Create scatter mechanisms
    GaAsNISTMott = mon.SelectableElasticSM(GaAs,mon.NISTMottRS.Factory)
    GaAsDS = mon.TabulatedInelasticSM(GaAs,3,GaAsTables)
#    GaAsDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441])
    #GaAsMoller = mon.MollerInelasticSM(GaAs)
    #GaAsPlasmon = mon.KoteraPlasmonInelasticSM(GaAs,1.)
    GaAsBarrier = mon.ExpQMBarrierSM(GaAs,0.05e-9)
    GaAsCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    GaAsMSM = mon.MONSEL_MaterialScatterModel(GaAs)
    GaAsMSM.addScatterMechanism(GaAsNISTMott)
    GaAsMSM.addScatterMechanism(GaAsDS)
    GaAsMSM.setCSD(GaAsCSD)
    GaAsMSM.setBarrierSM(GaAsBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    GaAsMSMDeep = mon.MONSEL_MaterialScatterModel(GaAs)
    GaAsMSMDeep.addScatterMechanism(GaAsNISTMott)
    GaAsMSMDeep.addScatterMechanism(GaAsDS)
    GaAsMSMDeep.setCSD(GaAsCSD)
    GaAsMSMDeep.setBarrierSM(GaAsBarrier)
    GaAsMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(GaAsMSM,GaAsMSMDeep)

def GeModel():
#	Germanium
    print("Warning. Germanium is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
    Ge = mon.SEmaterial([epq.Element.Ge],[1.],density,"Geremanium")
    GeWorkfunction=epq.ToSI.eV(workfun)
    Ge.setWorkfunction(GeWorkfunction)
    Ge.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathGe
#    GeTables = [tablePath+"IIMFPPennInterpGeSI.csv",tablePath+"interpNUSimReducedDeltaEGeSI.csv",tablePath+"interpsimTableThetaNUGeSI.csv",tablePath+"interpSimESE0NUGeSI.csv"]
    tablePath=tablePathAu
    GeTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    GeCoreEnergieseV = AuCoreEnergieseV
    for en in GeCoreEnergieseV:
	Ge.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Ge.getDensity()/epq.Element.Ge.getMass()
#Ge.addBindingEnergy(epq.ToSI.eV(0.)+GeWorkfunction,nve*density1electron)
# Create scatter mechanisms
    GeNISTMott = mon.SelectableElasticSM(Ge,mon.NISTMottRS.Factory)
    GeDS = mon.TabulatedInelasticSM(Ge,3,GeTables)
#    GeDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441]) 
    #GeMoller = mon.MollerInelasticSM(Ge)
    #GePlasmon = mon.KoteraPlasmonInelasticSM(Ge,1.)
    GeBarrier = mon.ExpQMBarrierSM(Ge,0.05e-9)
    GeCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    GeMSM = mon.MONSEL_MaterialScatterModel(Ge)
    GeMSM.addScatterMechanism(GeNISTMott)
    GeMSM.addScatterMechanism(GeDS)
    GeMSM.setCSD(GeCSD)
    GeMSM.setBarrierSM(GeBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    GeMSMDeep = mon.MONSEL_MaterialScatterModel(Ge)
    GeMSMDeep.addScatterMechanism(GeNISTMott)
    GeMSMDeep.addScatterMechanism(GeDS)
    GeMSMDeep.setCSD(GeCSD)
    GeMSMDeep.setBarrierSM(GeBarrier)
    GeMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(GeMSM,GeMSMDeep)


def TiNModel():
#	Titanium Nitride
    print("Warning. Titanium Nitride is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
#    TiN = mon.SEmaterial([epq.Element.TiN],[1.],density,"Titanium Nitride")
    TiN = mon.SEmaterial([epq.Element.Au],[1.],density,"Titanium Nitride")
    TiNWorkfunction=epq.ToSI.eV(workfun)
    TiN.setWorkfunction(TiNWorkfunction)
    TiN.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathTiN
    TiNTables = [tablePath+"IIMFPPennInterpTiNSI.csv",tablePath+"interpNUSimReducedDeltaETiNSI.csv",tablePath+"interpsimTableThetaNUTiNSI.csv",tablePath+"interpSimESE0NUTiNSI.csv"]
#    tablePath=tablePathAu
#    TiNTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    TiNCoreEnergieseV = AuCoreEnergieseV
    for en in TiNCoreEnergieseV:
	TiN.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = TiN.getDensity()/epq.Element.TiN.getMass()
#TiN.addBindingEnergy(epq.ToSI.eV(0.)+TiNWorkfunction,nve*density1electron)
# Create scatter mechanisms
    TiNNISTMott = mon.SelectableElasticSM(TiN,mon.NISTMottRS.Factory)
    TiNDS = mon.TabulatedInelasticSM(TiN,3,TiNTables)
#    TiNDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441])
    #TiNMoller = mon.MollerInelasticSM(TiN)
    #TiNPlasmon = mon.KoteraPlasmonInelasticSM(TiN,1.)
    TiNBarrier = mon.ExpQMBarrierSM(TiN,0.05e-9)
    TiNCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    TiNMSM = mon.MONSEL_MaterialScatterModel(TiN)
    TiNMSM.addScatterMechanism(TiNNISTMott)
    TiNMSM.addScatterMechanism(TiNDS)
    TiNMSM.setCSD(TiNCSD)
    TiNMSM.setBarrierSM(TiNBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    TiNMSMDeep = mon.MONSEL_MaterialScatterModel(TiN)
    TiNMSMDeep.addScatterMechanism(TiNNISTMott)
    TiNMSMDeep.addScatterMechanism(TiNDS)
    TiNMSMDeep.setCSD(TiNCSD)
    TiNMSMDeep.setBarrierSM(TiNBarrier)
    TiNMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(TiNMSM,TiNMSMDeep)

def CoModel():
#	Cobalt
    print("Warning. Cobalt is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
    Co = mon.SEmaterial([epq.Element.Co],[1.],density,"Cobalt")
    CoWorkfunction=epq.ToSI.eV(workfun)
    Co.setWorkfunction(CoWorkfunction)
    Co.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathCo
#    CoTables = [tablePath+"IIMFPPennInterpCoSI.csv",tablePath+"interpNUSimReducedDeltaECoSI.csv",tablePath+"interpsimTableThetaNUCoSI.csv",tablePath+"interpSimESE0NUCoSI.csv"]
    tablePath=tablePathAu
    CoTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    CoCoreEnergieseV = AuCoreEnergieseV
    for en in CoCoreEnergieseV:
	Co.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Co.getDensity()/epq.Element.Co.getMass()
#Co.addBindingEnergy(epq.ToSI.eV(0.)+CoWorkfunction,nve*density1electron)
# Create scatter mechanisms
    CoNISTMott = mon.SelectableElasticSM(Co,mon.NISTMottRS.Factory)
    CoDS = mon.TabulatedInelasticSM(Co,3,CoTables)
#    CoDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441])
    #CoMoller = mon.MollerInelasticSM(Co)
    #CoPlasmon = mon.KoteraPlasmonInelasticSM(Co,1.)
    CoBarrier = mon.ExpQMBarrierSM(Co,0.05e-9)
    CoCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    CoMSM = mon.MONSEL_MaterialScatterModel(Co)
    CoMSM.addScatterMechanism(CoNISTMott)
    CoMSM.addScatterMechanism(CoDS)
    CoMSM.setCSD(CoCSD)
    CoMSM.setBarrierSM(CoBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    CoMSMDeep = mon.MONSEL_MaterialScatterModel(Co)
    CoMSMDeep.addScatterMechanism(CoNISTMott)
    CoMSMDeep.addScatterMechanism(CoDS)
    CoMSMDeep.setCSD(CoCSD)
    CoMSMDeep.setBarrierSM(CoBarrier)
    CoMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(CoMSM,CoMSMDeep)

def DiamondModel():
#	Diamond
    print("Warning. Diamond is not yet available.  Defaulting to Gold.  Expect nonphysical results.")
    density = 19282.
#plasmonE = 9.11
    workfun = 5.1
    EFermi = 9.0
    potU = -workfun-EFermi # 
#    Diamond = mon.SEmaterial([epq.Element.Diamond],[1.],density,"Diamond")
    Diamond = mon.SEmaterial([epq.Element.Au],[1.],density,"Diamond")
    DiamondWorkfunction=epq.ToSI.eV(workfun)
    Diamond.setWorkfunction(DiamondWorkfunction)
    Diamond.setEnergyCBbottom(epq.ToSI.eV(potU))
    tablePath = tablePathDiamond
    DiamondTables = [tablePath+"IIMFPPennInterpDiamondSI.csv",tablePath+"interpNUSimReducedDeltaEDiamondSI.csv",tablePath+"interpsimTableThetaNUDiamondSI.csv",tablePath+"interpSimESE0NUDiamondSI.csv"]
#    tablePath=tablePathAu
#    DiamondTables = [tablePath+"IIMFPPennInterpAuSI.csv",tablePath+"interpNUSimReducedDeltaEAuSI.csv",tablePath+"interpsimTableThetaNUAuSI.csv",tablePath+"interpSimESE0NUAuSI.csv"]
    AuCoreEnergieseV = [57.2, 74.2, 83.9, 87.6, 107.2, 335.1, 353.2, 546.3, 642.7, \
	762.1, 2206, 2291, 2743, 3148, 3425, 11919, 13734, 14353, 80725]
    DiamondCoreEnergieseV = AuCoreEnergieseV
    for en in DiamondCoreEnergieseV:
	Diamond.addCoreEnergy(epq.ToSI.eV(en))
#density1electron = Diamond.getDensity()/epq.Element.Diamond.getMass()
#Diamond.addBindingEnergy(epq.ToSI.eV(0.)+DiamondWorkfunction,nve*density1electron)
# Create scatter mechanisms
    DiamondNISTMott = mon.SelectableElasticSM(Diamond,mon.NISTMottRS.Factory)
    DiamondDS = mon.TabulatedInelasticSM(Diamond,3,DiamondTables)
#    DiamondDS.setBranchingRatios([0.78, 0.62, 0.441262]) #, 0.199441])
    #DiamondMoller = mon.MollerInelasticSM(Diamond)
    #DiamondPlasmon = mon.KoteraPlasmonInelasticSM(Diamond,1.)
    DiamondBarrier = mon.ExpQMBarrierSM(Diamond,0.05e-9)
    DiamondCSD = mon.ZeroCSD()
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
    DiamondMSM = mon.MONSEL_MaterialScatterModel(Diamond)
    DiamondMSM.addScatterMechanism(DiamondNISTMott)
    DiamondMSM.addScatterMechanism(DiamondDS)
    DiamondMSM.setCSD(DiamondCSD)
    DiamondMSM.setBarrierSM(DiamondBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
    DiamondMSMDeep = mon.MONSEL_MaterialScatterModel(Diamond)
    DiamondMSMDeep.addScatterMechanism(DiamondNISTMott)
    DiamondMSMDeep.addScatterMechanism(DiamondDS)
    DiamondMSMDeep.setCSD(DiamondCSD)
    DiamondMSMDeep.setBarrierSM(DiamondBarrier)
    DiamondMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))
    return(DiamondMSM,DiamondMSMDeep)


