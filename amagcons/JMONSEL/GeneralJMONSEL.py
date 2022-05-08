####Note great Jython syntax website:  http://www.jython.org/jythonbook/en/1.0/LangSyntax.html#operators
WinOrLin = 1	#=1 for WINDOWS/StandardSetup/homePC/workPC/Prolith1/Marcy/Spock, =2 for LINUX, =3 for WINDOWS/Metrosim group sim computer, =4 for WINDOWS/Prolith2 computer
print("Position 0")
import sys
import os

if (WinOrLin == 2):
        print("Position 0.1adfaf")
	import sys
	print("Position 0.2")
	# add for LINUX, make comment if Windows
	print("Position 0.3")
	sys.path.append('/work/')
	#files = [file.strip() for file in open( "./jar.list", "r" )]		# add for LINUX, make comment if Windows
	print("Position 0.4")
	for file in files:
                # add for LINUX, make comment if Windows
                print("Position 0.5")
		sys.path.append(file)
		print("Position 0.6")
		# add for LINUX, make comment if Windows
		print("Position 0.7")
print("Position 0.1")
basepath = '/work/'
files = os.listdir(basepath)
for i in range(0,len(files),1):
        path = basepath + "/" + files[i]
        print(path)
        sys.path.append(path)
print('Added Path')
print('Position 0.2')
import gov.nist.microanalysis.EPQLibrary as epq
print('Position 0.3')
import gov.nist.microanalysis.EPQTools as ept
print('Position 0.4')
import gov.nist.microanalysis.NISTMonte as nm
print('Position 0.5')
import gov.nist.nanoscalemetrology.JMONSEL as mon
print('Position 0.6')
import gov.nist.microanalysis.Utility as nmu
print('Position 0.7')
import java.io as jio
print('Position 0.8')
import java.util as jutil
import java.lang as jl
import jarray
import java.nio.charset as cs
print("Position 1")
import sys
sys.modules.clear()
print("Position 2")
#sys.path.append('C:\\Users\\Shari\\Documents\\JMONSELMySimProject\\Examples')
#sys.path.append('C:\\NIST\\JMONSEL\\Examples')
#sys.path.append('C:\\NIST\\JMONSEL')
#sys.path.append('C:\\NIST\\JMONSEL\\')
sys.path.append('/work/')
print("Position 3")
import definesample as geometry
print("Position 4")
PathSep = "/"
DefaultOutput = basepath + PathSep + 'output_qkss'
# determine where to save the results
dest=DefaultOutput;
jio.File(dest).mkdirs()
PathSep = "/"
print("Position 64")
filename = DefaultOutput + PathSep + "Results.txt"
file = open(filename,'w')
print "Output will be to file: ",filename

# Conversions of shape parameters to SI units.
meterspernm = 1.e-9	# conversion from nanometers to meters

trajImg = 0
trajImgMaxTraj = 500     
trajImgSize = 150.e-9

VRML = 1
VRMLImgMaxTraj = 0 # Include no trajectories in VRML (Show sample only.) Leaving trajectories
# out makes a VRML that displays easily. It's good for checking the sample. I find that adding
# trajectories significantly slows down the display. If you want to try it, keep the number of
# displayed trajectories small (20 is a reasonable number) and turn off collision detection in
# your VRML viewer.

binSizeEV = 10.	# Width (in eV) of bins in energy histogram

# Model parameters
nTrajectories = 10
mExp = 1		# # of experiments
RotAngDeg = 0.	# Rotation Angle of Grating (in degrees)

def xygettrajectories(xnm,ynm,monte,beamEeV):
	x = xnm*meterspernm
	y = ynm*meterspernm
	eg.setCenter([x,y,-1050.*meterspernm]) # Aims the gun at x,y.

# Define our backscatter detector.
        back=nm.BackscatterStats(monte)	
        nbins = int(beamEeV/binSizeEV)
        monte.addActionListener(back)
        back.setEnergyBinCount(nbins)
									
	# Add a trajectory image
	if trajImg:  # output the trajectory image
		img=nm.TrajectoryImage(512,512,trajImgSize)
		img.setMaxTrajectories(trajImgMaxTraj)
		img.setYRange(-y-trajImgSize/10,.9*trajImgSize-y)
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
	print >>file, "%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%8.4f" % (beamsizenm,beamEeV,m5,xnm,ynm,bsf,SEf,t)
	print "%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%8.4f" % (beamsizenm,beamEeV,m5,xnm,ynm,bsf,SEf,t)
	rdirname = jl.Integer(int(beamEeV)).toString()+"eV_"+"x"+jl.Integer(int(xnm)).toString()+'y'+jl.Integer(int(ynm)).toString()+"traj.prn"

	back.dump(jio.FileOutputStream(dest+PathSep+rdirname+"backscatter.prn"))
	monte.removeActionListener(back)
#	print(dest+PathSep+rdirname)	
        if trajImg:  # output the trajectory image
#       	img.dumpToFile(dest+rdirname)
		img.dumpToFile(dest+PathSep+rdirname)
		monte.removeActionListener(img)

	if VRML:
		tw.flush()
		fos.close()
		monte.removeActionListener(vrml)
#end of xygettrajectories

def looptrajectories(m5vals,xvals,yvals,monte,beamEeV):
	for m5 in m5vals:
		for xnm in xvals:          #####################To switch x and y, do it here and in 4 lines below
			for ynm in yvals:
				xygettrajectories(xnm,ynm,monte,beamEeV)
# end of looptrajectories

# Beginning of main program

# We'll make the substrate infinitely thick.
beamEeVvals = [500.] # Beam energies in eV   [100.,150.,200.,250.,300.,500.,800.]
beamsizenmvals = [0.5] # beam size in nm  [0.,0.1,0.2,0.3,0.5,1.]


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

radperdeg = jl.Math.PI/180.	# conversion from degrees to radians

# Print simulation parameters to window and output file.
print >>file, "# Trajectories at each landing position: ",nTrajectories
print "# Trajectories at each landing position: ",nTrajectories
print >>file, "# Experiments at each landing position: ",mExp
print "# Experiments at each landing position: ",mExp

print >>file, "Beam landing energies (eV): ",beamEeVvals
print "Beam landing energies (eV): ",beamEeVvals
print >>file, "Beam size (standard deviation, in nm): ",beamsizenmvals
print "Beam size (standard deviation, in nm): ",beamsizenmvals

print 

print >>file, "beamsize(nm)	beamE(eV)   Experiment#	x(nm)	 y(nm)	BSEyield	SEyield	  ElapsedTime"     #print >>file, "beamsize(nm) h(nm) w(nm) SWAl(deg) SWAr(deg) beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield	Elapsed Time"
print         "beamsize(nm)	beamE(eV)   Experiment#	x(nm)	 y(nm)	BSEyield	SEyield	  ElapsedTime"     #print "beamsize(nm) h(nm) w(nm) SWAl(deg) SWAr(deg) beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield	Elapsed Time"

t0 = jl.System.currentTimeMillis()


#loop beamsizenmvals
for beamsizenm in beamsizenmvals:

	beamsize = beamsizenm*meterspernm

# define multiple experiments??
	m5=1
	deltam5=1
	m5vals = []
	while m5-1 < mExp:
		m5vals.append(m5)
		m5 += deltam5

	# SAMPLE DESCRIPTION
	monte,eg, xvals,yvals,nTrajectories=geometry.definesample(beamsize)                
	print "nTrajectories=",nTrajectories					
#       chamber.rotate([0.,0.,0.],RotAng,0.,0.)         #Entire sample rotated!

# Run the simulation
	for beamEeV in beamEeVvals:
		beamE = epq.ToSI.eV(beamEeV)
		monte.setBeamEnergy(beamE) # sets this model's beam energy
# reduce the number of loops in the main program to make it easier to add loops over other parameters
		looptrajectories(m5vals,xvals,yvals,monte,beamEeV)
