# Copyright: Pursuant to title 17 Section 105 of the United States Code this
# software is not subject to copyright protection and is in the public domain.

# Institution: National Institute of Standards and Technology

# Author: John Villarrubia, July 23, 2010

# This script determines the yield for electrons incident on 1 or more
# trapezoidal (with top corner radii) lines on a 3-layer substrate with 
# an on-line defect in the central line.
 
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
filename = DefaultOutput+PathSep+"OnLineDefect_results.txt"
file = open(filename,'w')
print "Output will be to file: ",filename

# Make a record of the random seed we use so we can exactly repeat this calculation
# (same random number sequence) if necessary.
seed = nmu.Math2.rgen.nextLong() # Pick a random seed


# To exactly repeat a previous calculation (e.g., for bug fix) uncomment the next line
# and replace the question marks with the seed that was recorded in the previous calculation's 
# output file.

#seed = -5991603611489619985L
print >>file,"Random number seed: ",seed
print "Random number seed: ",seed
nmu.Math2.initializeRandom(seed)
#for i in range(0,10):
#    r = nmu.Math2.rgen.nextDouble()
#    print >>file, r
#    print r

# Model parameters
nTrajectories = 5000

# Shape parameters.
pitchnm = 50. 	# Distance between line centers in nm
nlines = 11		# number of lines
hnm = 100.  # line height in nm
wnm = 25.	# line bottom width in nm
linelengthnm = 1000 #  line length in nm
# Note that sidewall angles are specified with respect to vertical,
# so 0. is vertical, positive angles have bottom wider than top, and
# negative angles are the reverse (undercut).
thetardeg = 1.	# line right sidewall angle in degrees
thetaldeg = 5.	# line left sidewall angle in degrees
radrnm = 3.		# line top right corner radius in nm
radlnm = 3.		# line top left corner radius in nm
layer1thicknessnm = 15.	# Thickness in nm of the 1st layer (immediately below the lines)
layer2thicknessnm = 200.# Thickness in nm of the 2nd layer
# We'll make the substrate infinitely thick.
beamEeVvals = [1500.] # Beam energies in eV
beamsizenm = 10.0	# beam size in nm. This is the standard deviation of a Gaussian distribution.
# This is much larger than typical for a high resolution SEM. I've indicated the lower resolution
# of a defect detection tool by choosing this larger value, but it's probably still not large enough.
# Adjust as desired.

defectDepthFrac = 0.25 # Defect depth as a fraction of line height
defectLengthnm = 15. # Defect length in nm

trajImg = 1
trajImgMaxTraj = 50 
trajImgSize = 150.e-9

VRML = 1
VRMLImgMaxTraj = 0	# Include no trajectories in VRML (Show sample only.) Leaving trajectories
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
tablePath = "C:\Program Files\NIST\JMONSEL\ScatteringTables\SiTables"+PathSep
SiTables = [tablePath +"IIMFPFullPennInterpSiSI.csv", 
tablePath +"interpNUSimReducedDeltaEFullPennSiSI.csv", 
tablePath +"interpNUThetaFullPennSiBGSI.csv", 
tablePath +"interpSimESE0NUSiBGSI.csv"]
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
SiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

# Conversions of shape parameters to SI units.
# Shape parameters.
meterspernm = 1.e-9	# conversion from nanometers to meters
pitch = pitchnm*meterspernm
h = hnm*meterspernm 
w = wnm*meterspernm
linelength = linelengthnm*meterspernm
# Note that sidewall angles are specified with respect to vertical,
# so 0. is vertical, positive angles have bottom wider than top, and
# negative angles are the reverse (undercut).
radperdeg = jl.Math.PI/180.	# conversion from degrees to radians
thetar = thetardeg*radperdeg
thetal = thetaldeg*radperdeg
radr = radrnm*meterspernm
radl = radlnm*meterspernm
layer1thickness = layer1thicknessnm*meterspernm
layer2thickness = layer2thicknessnm*meterspernm
beamsize = beamsizenm*meterspernm
defectDepth = h*defectDepthFrac # Defect depth in meters
defectLength = defectLengthnm*meterspernm # Defect length in meters

# create an instance of the model
monte=nm.MonteCarloSS() #creates an instance of the model with all default characteristics
eg = nm.GaussianBeam(beamsize) # makes electron gun, Gaussian with standard deviation = beamsize
monte.setElectronGun(eg) # This gun is "attached" to the model.

# SAMPLE DESCRIPTION

# NISTMonte provides us with a "chamber" in the form of a 0.1 m sphere inside of which we build
# our sample. Replace the default vacuum in the chamber with SEvacuum. (SEmaterials define additional
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
# This region has shape defined by layer1, scattering properties defined for Si, and is a subregion of 
# the chamber.
#layer1Region = monte.addSubRegion(chamber,SiMSM,layer1)  
layer1Region = monte.addSubRegion(chamber,SiMSM,layer1)  

# The problem we have set ourselves calls for only a single substrate material. I.e., it's not layered.
# However, layer2 below illustrates a trick that speeds up the code excecution. Almost all secondary
# electrons have energies below 50 eV. Such electrons have short escape depths--on the order of a 
# few nanometers--in any solid material. The ones generated at greater depths than this have little
# to no chance of escaping the sample, so any time we spend computing and tracking their trajectories
# is wasted time. The waste of time is especially noticeable at larger landing energies, since high energy
# electrons spend more of their time at greater depth and they generate more secondary electrons.
# In this script, we avoid wasting as much time by making layer1 and layer2 out of the same material.
# We do full modeling, including all secondaries, in layer1, which is close enough to the surface that
# such electrons are relevant for our yield calculation. We make the thickness of layer1 several times
# the escape depth, just to be safe. Then in layer2 we use a model that drops electrons with 
# energies <50 eV, that is, all but the highest energy secondaries.

layer2 = mon.NormalMultiPlaneShape()
layer2.addPlane(normalvector,[0.,0.,layer1thickness]) #layer 2 starts layer1thickness farther up.
# We make it a subregion of layer1Region. At this point, layer1Region
# extends only for 0<=z<=layer1thickness.
layer2Region = monte.addSubRegion(layer1Region,SiMSMDeep,layer2)  

# The code for layer3 is commented out because we don't need it. To include such a layer, uncomment 
# the code and change the material scatter model parameter (the 2nd argument of addSubRegion) to be a 
# model for the material you wish to assign to the layer. You'll have to include a code block  
# analogous to the Si block above to define that new material. 
layer3 = mon.NormalMultiPlaneShape()
layer3.addPlane(normalvector,[0.,0.,layer1thickness+layer2thickness]) # Notice that the position of  
# a layer surface is equal to the sum of thicknesses of all the layers above it.
# We give it the properties of ???, and make it a subregion of layer2Region. At this point, layer2Region
# extends only for layer1thickness<=z<=layer1thickness+layer2thickness.
layer3Region = monte.addSubRegion(layer2Region,SiMSM,layer3)

# THE LINES

# We always place at least one line, at the center. This line has a more complicated structure
# because it's the one where we place our defect. The defect we make here is a cutout--a missing portion--
# at the top of the line. There is more than one way to generate this. Here, we choose to do it by 
# constructing the line in 4 parts: (1) The undefective positive y part of the line, (2) the undefective 
# negative y part, (3) The "bridge" that connects parts 1 and 2 below the defect and (4) the defect itself.
# Part 4 is not strictly necessary for our present simulation. Since the defect is empty, and it is a 
# subregion of the chamber, which is also empty, we could just leave it out. Instead, we include an explicit
# region for it and we "fill" this region with vacuum. This is equivalent, but it has the advantage that
# we can easily edit this script to fill it with something else should we so desire.

# An alternative method would be to define the entire undefected line, separately define the defect, 
# and then subtract the defect from the line using NormalDifferenceShape(). This second method is in 
# some ways more general, so keep it in mind, but a disadvantage is that NormalDifferenceShape does not
# implement the TrajectoryVRML.IRender interface. This means VRML renderings would not work.

# The next code line generates a line at the origin with length (linelength-defectLength)/2.
# in the y direction.
upperline =  mon.NShapes.createLine(-h,w,(linelength-defectLength)/2.,thetal,thetar,radl,radr)
# Move it up the y axis so it extends from defectLength/2. < y < linelength/2.
upperline.translate([0.,(linelength+defectLength)/4.,0.])
# Now do the same thing for the negative y part
lowerline =  mon.NShapes.createLine(-h,w,(linelength-defectLength)/2.,thetal,thetar,radl,radr)
lowerline.translate([0.,-(linelength+defectLength)/4.,0.])
# Now make the "bridge" between these two parts. We'll assume the bridge height does not extend into
# the rounded corners at the top of the line.
bridge = mon.NShapes.createLine(-(h-defectDepth),w,defectLength,thetal,thetar,0.,0.)
# Now make the defect...
# If the edges are sloped, the defect bottom corner positions have to be adjusted for height
xDefectLeft = -w/2+(h-defectDepth)*jl.Math.tan(thetal)
xDefectRight = w/2-(h-defectDepth)*jl.Math.tan(thetar)
defectBottWidth = xDefectRight-xDefectLeft
# The defect bottom width must be correspondingly adjusted
defect = mon.NShapes.createLine(-defectDepth,defectBottWidth,defectLength,thetal,thetar,radl,radr)
# ...and translate it into position
# As generated its bottom right edge is at defectBottWidth/2. We need to put it at xDefectRight
defect.translate([xDefectRight-defectBottWidth/2.,0.,-(h-defectDepth)])

# Now make regions (shape + material scatter model) for all these shapes, and put them in the chamber
upperlineRegion = monte.addSubRegion(chamber,SiMSM,upperline)
lowerlineRegion = monte.addSubRegion(chamber,SiMSM,lowerline)
bridgeRegion = monte.addSubRegion(chamber,SiMSM,bridge)
defectRegion = monte.addSubRegion(chamber,vacuumMSM,defect)

# If the simulation parameters above specify more than one line, we add the extra ones here.
# We alternate them first to the right and then to the left of the center line.
for i in range(nlines-1):
	# Note that in the integer arithmetic below, fractions are truncated to the next lowest integer.
	distanceIndex = 1+(i-1)/2 # distanceIndex = 1 for the 1st pair of lines, 2 for the 2nd pair, etc.
	if i - 2*(i/2) == 0:	# This means "if i is even"
		distanceIndex *= -1	# Put even-numbered lines on the left
	xcenter = distanceIndex*pitch
	line = mon.NShapes.createLine(-h,w,linelength,thetal,thetar,radl,radr)
	line.translate([xcenter,0.,0.])
	lineRegion = monte.addSubRegion(chamber,SiMSM,line)

# SCAN PARAMETERS

# The following block of code determines the landing positions we'll simulate. We'll put
# x coordinates in one list (called xvals) and the corresponding ycoordinates in another (yvals)
# of equal length. We'll include two linescans, one from left to right across the defect and
# one from below it to above. 

# We'll do a nonuniform scan, with closely spaced landing positions near points of particular interest
# (the defect in this case) and larger spacings farther away.

# First we generate the left to right scan, with 1 nm pixel spacing within about 20 nm of the
# right edge of the defect. (Note that 1 nm is probably overkill for large beam sizes.) Farther
# from the edge we'll use 5 nm spacing.

xbottom = wnm/2. # Coordinate of right edge of defect
xtop = wnm/2.-hnm*jl.Math.tan(thetar)
xstart = xbottom - 100.5
xstop = xbottom + 100.5
xfinestart = xtop-20.5
if thetar<0.:	# undercut line
	xfinestop = xtop+20.5
else: # normal line
	xfinestop = wnm/2.+20.5	# 20.5 nm to the right of the bottom corner

yconstantval = 0.	# y=0 for our horizontal scan
xvals = []
yvals = []
deltax = 5.
x = xstart
while x<xfinestart:
	xvals.append(x)
	yvals.append(yconstantval)
	x += deltax
x = xfinestart
deltax = 1.
while x<xfinestop:
	xvals.append(x)
	yvals.append(yconstantval)
	x += deltax
x = xfinestop
deltax = 5.
while x<xstop:
	xvals.append(x)
	yvals.append(yconstantval)
	x += deltax
xvals.append(xstop)
yvals.append(yconstantval)

# Now generate the bottom to top scan analogously. 
yEdge = defectLengthnm/2.	# y coordinate of the top edge of the defect
start = -yEdge - 50.5
stop = yEdge + 50.5
finestart = yEdge-20.5
finestop = yEdge+20.5

xconstantval = 0.	# x=0 for our vertical scan
delta = 5.
y = start
while y<finestart:
	xvals.append(xconstantval)
	yvals.append(y)
	y += delta
x = xfinestart
delta = 1.
while y<finestop:
	xvals.append(xconstantval)
	yvals.append(y)
	y += delta
y = finestop
delta = 5.
while y<stop:
	xvals.append(xconstantval)
	yvals.append(y)
	y += delta
yvals.append(stop)
xvals.append(xconstantval)

# Print simulation parameters to window and output file.
print >>file, "# Trajectories at each landing position: ",nTrajectories

print "# Trajectories at each landing position: ",nTrajectories
print >>file, "# Pitch of lines (nm): ",pitchnm
print "# Pitch of lines (nm): ",pitchnm
print >>file, "# lines: ",nlines
print "# lines: ",nlines
print >>file, "Line height (nm): ",hnm
print "Line height: ",hnm
print >>file, "Line bottom width (nm): ",wnm
print "Line bottom width (nm): ",wnm
print >>file, "Line length (nm): ",linelengthnm
print "Line length (nm): ",linelengthnm
print >>file, "Left and right sidewall angles (deg): ",thetaldeg,thetardeg
print "Left and right sidewall angles (deg): ",thetaldeg,thetardeg
print >>file, "Left and right top corner radii (nm): ",radlnm,radrnm
print "Left and right top corner radii (nm): ",radlnm,radrnm
print >>file, "Thicknesses of 1st and second layers (nm): ",layer1thicknessnm,layer2thicknessnm
print "Thicknesses of 1st and second layers (nm): ",layer1thicknessnm,layer2thicknessnm
print >>file, "Beam landing energies (eV): ",beamEeVvals
print "Beam landing energies (eV): ",beamEeVvals
print >>file, "Beam size (standard deviation, in nm): ",beamsizenm
print "Beam size (standard deviation, in nm): ",beamsizenm
print >>file, "Defect depth as a fraction of total line height: ",defectDepthFrac
print "Defect depth as a fraction of total line height: ",defectDepthFrac
print >>file, "Defect length: ",defectLengthnm
print "Defect length: ",defectLengthnm

print >>file	# Blank line before start of calculation results
print 

print >>file, "beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield"
print "beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield"

binSizeEV = 10.	# Width (in eV) of bins in energy histogram

for beamEeV in beamEeVvals:
	beamE = epq.ToSI.eV(beamEeV) # convert beam energy from eV to SI units (Joules)
	monte.setBeamEnergy(beamE) # sets this model's beam energy
	for lpindex in range(len(yvals)):	# lpindex is landing position index
		ynm = yvals[lpindex]
		xnm = xvals[lpindex]
		y = ynm*meterspernm
		x = xnm*meterspernm
		eg.setCenter([x,y,-h-20.*meterspernm]) # Aims the gun at x,y.

		# Define our backscatter detector.
		back=nm.BackscatterStats(monte)	
		nbins = int(beamEeV/binSizeEV)
		monte.addActionListener(back)
		back.setEnergyBinCount(nbins)
			
		# Add a trajectory image
		if trajImg:  # output the trajectory image
			img=nm.TrajectoryImage(2048,2048,trajImgSize)
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

		print >>file, "%8.1f \t%8.1f \t%8.1f \t%3.3f \t%3.3f" % (beamEeV,xnm,ynm,bsf,SEf)
		print "%8.1f \t%8.1f \t%8.1f \t%3.3f \t%3.3f" % (beamEeV,xnm,ynm,bsf,SEf)

		back.dump(jio.FileOutputStream(dest+PathSep+"backscatter.prn"))
		monte.removeActionListener(back)

		if trajImg:  # output the trajectory image
			img.dumpToFile(dest)
			monte.removeActionListener(img)

		if VRML:
			tw.flush()
			fos.close()
			monte.removeActionListener(vrml)
	
file.close()
