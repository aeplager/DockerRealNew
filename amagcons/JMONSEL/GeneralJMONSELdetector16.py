####Note great Jython syntax website:  http://www.jython.org/jythonbook/en/1.0/LangSyntax.html#operators
WinOrLin = 1	#=1 for WINDOWS/StandardSetup/homePC/workPC/Prolith1/Marcy/Spock, =2 for LINUX, =3 for WINDOWS/Metrosim group sim computer, =4 for WINDOWS/Prolith2 computer
print("Position 0")
import sys
import os
import shutil
from datetime import datetime

if (WinOrLin == 2):
        print("Position 0.1")
	import sys
	print("Position 0.2")
	# add for LINUX, make comment if Windows
	print("Position 0.3")
	sys.path.append('/work/')
	files = [file.strip() for file in open( "./jar.list", "r" )]		# add for LINUX, make comment if Windows
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
DefaultOutput = '/work/output_qkss/'
# Establish Default Data
DefaultOutput = '/output_data/'

print("Clearing Directory and Files")
i = 0
for root, dirs, files in os.walk(DefaultOutput):
    for f in files:
        i=i+1
        os.unlink(os.path.join(root, f))
    for d in dirs:
        print("Deleted Directory")
        i=i+1
        shutil.rmtree(os.path.join(root, d))
print("Cleared Directory and Files:  " + str(i))        
        
files = os.listdir(basepath)
for i in range(0,len(files),1):
        path = basepath + "/" + files[i]
        print(path)
        sys.path.append(path)										# add for LINUX, make comment if Windows
# Opening Status of File Opening

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

import math, sys
sys.modules.clear()
print("Position 2")
sys.path.append('/work/')
print("Position 3")
# determine where to save the results
dest=DefaultOutput;
jio.File(dest).mkdirs()
PathSep = "/"
print("Position 64")
filename = DefaultOutput + PathSep + "Results.txt"
file = open(filename,'w')
print "Output will be to file: ",filename

#open a file for each scan region, so move it from here to scan region loop

# add the direction in which this file is location to the path
# so that it will find definesample.py and makematerials.py
# adding the directory to the beginning of the path
# so that it does not find files in the working directory first
# makematerials.py could be in the working directory
currentpath=dest.rsplit(PathSep,1)[0]
sys.path.insert(0,currentpath)  #is prepending to path a good idea??
#sys.path.append(currentpath)

import definesample as geometry
print("Position 4")
PathSep = "/"
DefaultOutput = basepath + PathSep + 'output_qkss'
# AEP Changed Path Here
DefaultOutput = '/output_data/'
import savedata

# Starting Status
print("Started Writing Status")
statusfilename = DefaultOutput+"Status.csv"
statusfileid = open(statusfilename,'a')
now = datetime.now() # current date and time
date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
savedata.writeStatusFile(statusfileid, "Starting", date_time)
statusfileid.close()
print("Finished First Line")


# Conversions of shape parameters to SI units.
meterspernm = 1.e-9	# conversion from nanometers to meters

# Trajectory image parameters
# Set trajImg =1 to save trajectory file(s). Set it to 0 to neither record nor save them.
# Set overwriteTF=1 to keep only a single (the latest) trajectory file.
#trajImg = 0
trajImgMaxTraj = 500     
trajImgSize = 150.e-9
trajImgFile=""

# Trajectory log: To turn on, set trajLog on the next line to 1.
#trajLog = 1
trajLogMaxTraj = 1000
trajlogfile=""

histfile=""

VRML = 1
VRMLImgMaxTraj = 0 # Include no trajectories in VRML (Show sample only.) Leaving trajectories
# out makes a VRML that displays easily. It's good for checking the sample. I find that adding
# trajectories significantly slows down the display. If you want to try it, keep the number of
# displayed trajectories small (20 is a reasonable number) and turn off collision detection in
# your VRML viewer.

# Model parameters
tprevsec=0
#mExp = 1		# # of experiments
RotAngDeg = 0.	# Rotation Angle of Grating (in degrees)
thicknessvals = [0]   # if in nm, make sure it is converted to meters.  most parameters converted in definesample.py
paramvals = [0]

# Beam Parameters--most of these are defined in definesample.py
binSizeEV = 10.	# Width (in eV) of bins in energy histogram

# Scanning parameters
deltax = 10.	# size of fine x pixel
deltay = 10.	# size of fine y pixel
nx = 11L		# number of x pixels per linescan
ny = 1L		# number of linscans per image
xmin = -50.
ymin = 0.
pixelsPerFrame = nx*ny
#gunZnm = 1050.

def defineBackscatterChamber(monte,nbins):
        back=nm.BackscatterStats(monte)	
        monte.addActionListener(back)
        back.setEnergyBinCount(nbins)
        return(back)

def regionBackscatter(monte,detectorRegion,nbins,maxE,minE,destructive):
        back_regdet=mon.RegionDetector(monte,[detectorRegion],epq.ToSI.eV(minE),epq.ToSI.eV(maxE),destructive)	#regiondetectorTEST=RegionDetector(monte,[linedefect9region,reg2,reg3],minE,maxE,boolean0or1)
# conversion from eV to J (?)
        back_regdet.setLogDetected(0) #toggles on the saving the log info from detector, electronID, step#, energy, positionxyz, directionthetaphi
        monte.addActionListener(back_regdet)
        back_regdet.setEnergyBinCount(nbins)
        return(back_regdet)

def XfChamberBins(back,beamEeV,nTrajectories,minbin,maxbin):
	hist = back.backscatterEnergyHistogram()
	energyperbineV = beamEeV/hist.binCount()
	totalSE = 0
	for j in range(minbin/energyperbineV,maxbin/energyperbineV):
		totalSE = totalSE+hist.counts(j)
	SEf = float(totalSE)/nTrajectories
        return(SEf)

def getSfandBsf(back,beamEeV,nTrajectories):
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
# backscatterFraction() on the next line returns total electron fraction scattered back toward
# the electron gun regardless of energy. We divide these into SE (E<=50 eV) and BSE (E>50eV).
# The fraction of SE were counted above. To get the fraction of BSE we have to subtract SEf
# from the total backscatterFraction().
	bsf = back.backscatterFraction()-SEf#+float(totalFwd)/nTrajectories
        return (SEf,bsf)

def getSfandBsfRegion(back,beamEeV,nTrajectories,minbin,maxbin):
	hist = back.energyHistogram()
	energyperbineV = beamEeV/hist.binCount()
	maxSEbin = 50./energyperbineV	# bin number of the one with 50 eV
	totalSE = 0
 	for j in range(0,int(maxSEbin)):
		totalSE = totalSE+hist.counts(j)
#                print(j,totalSE,hist.counts(j))
		#totalSE = totalSE+hist.counts(j)+fhist.counts(j)
	SEf = float(totalSE)/nTrajectories
	bsf = back.detectedFraction()-SEf#+float(totalFwd)/nTrajectories
        totalXf = 0
 	for j in range(int(minbin/energyperbineV),int(maxbin/energyperbineV)):
		totalXf = totalXf+hist.counts(j)
        Xf=totalXf/nTrajectories
        return (SEf,bsf,Xf)

# The times in the following ScanGenerator are not used.
#change gunZnm to beamStartZnm
# this call should be in the scanregion loop, and the xvals and yvals are defined differently (untested)
#scanGen = mon.SimpleRasterScanGenerator(xmin, ymin, -beamStartZnm, deltax, deltay, nx, ny, 0., 1.e-6, 5.e-6, 10.e-6)

#class computeTrajectories(object):  #thinking about changing to object oriented program
#    def __init__(self):
#        pass

def xygettrajectories(mexp,xnm,ynm,pixelnum,monte,beamEeV,beamsize,paramlist,chamberdict,detectorList,fileid):
	x = xnm*meterspernm
	y = ynm*meterspernm
	eg.setCenter([x,y,-chamberdict["gunZnm"]*meterspernm]) # Aims the gun at x,y.
        nbins = int(beamEeV/binSizeEV)

# Define our backscatter detector.
	back=defineBackscatterChamber(monte,nbins)
        back_regdetList=[]
        for detector in detectorList:
             	back_regdet=regionBackscatter(monte,detector["Name"],nbins,detector["maxE"],detector["minE"],int(detector["destruct"]))
		detector["back_regdet"]=back_regdet
									
	# Add a trajectory image
        trajImgFile=""
	if trajImg:  # output the trajectory image
		img=nm.TrajectoryImage(512,512,trajImgSize)
		img.setMaxTrajectories(trajImgMaxTraj)
		img.setYRange(-y-trajImgSize/10,.9*trajImgSize-y)
		img.setXRange(x-trajImgSize/2.,x+trajImgSize/2.)
		monte.addActionListener(img)
                trajImgFile="trajectory.png"

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
        trajlogfile=dest+PathSep+"no_file"
	if trajLog:  # add a trajectory log listener
                trajlogfile=dest+PathSep+"trajectory"+str(xnm)+"_"+str(ynm)+".log"
	        log=mon.TrajectoryLoggingListener(monte,trajlogfile)
	        log.setMaxTrajectories(trajLogMaxTraj)
	        monte.addActionListener(log)
 
        tbefore = (jl.System.currentTimeMillis()-t0)/1000.
	monte.runMultipleTrajectories(nTrajectories)
        SEf,bsf=getSfandBsf(back,beamEeV,nTrajectories)
        Xfchamber=XfChamberBins(back,beamEeV,nTrajectories,chamberdict["XfminE"],chamberdict["XfmaxE"])

	t = (jl.System.currentTimeMillis()-t0)/3600000.
	rdirname = jl.Integer(int(beamEeV)).toString()+"eV_"+"x"+jl.Integer(int(xnm)).toString()+'y'+jl.Integer(int(ynm)).toString()
        histfile=dest+PathSep+rdirname+"backscatter.prn"
        if not histogram:
            histfile=dest+PathSep+"no_histogram"
        rotbeam="x0_y0_z0"
        savedata.writeChamberData(fileid,pixelnum,t*3600,t,t*3600-tbefore,chamberdict["scanregname"],xnm,ynm,mexp,nTrajectories,beamEeV,beamsize,beamsize,rotbeam,chamberdict['rotAng'],histfile,trajImgFile,trajlogfile,bsf,SEf,Xfchamber,paramlist)
        tprevsec=t*3600

        detBsf=0
        detBSf2 = 0
        for detector in detectorList:
            detSEf,detBsf,detXf=getSfandBsfRegion(detector["back_regdet"],beamEeV,nTrajectories,detector["minbin"],detector["maxbin"])
            detBSf2=detector["back_regdet"].detectedFraction()
 #           histfile is the same as for chamber data
            savedata.writeDetectorData(fileid,detBsf,detSEf,detXf,histfile,trajlogfile)
        savedata.endline(fileid)

#	print >>file, "%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%8.4f \t%8.2f" % (beamsizenm,beamEeV,mexp,xnm,ynm,bsf,SEf,t,paramlist[0])
	print "%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%3.5f \t%3.5f \t%8.4f \t%8.2f" % (beamsizenm,beamEeV,mexp,xnm,ynm,bsf,SEf,detBsf,detBSf2,t,paramlist[0])

        if histogram:  #check before saving histogram file
	    back.dump(jio.FileOutputStream(dest+PathSep+rdirname+"backscatter.prn"))
	    monte.removeActionListener(back)
        if trajImg:  # output the trajectory image
#       	img.dumpToFile(dest+rdirname)
		img.dumpToFile(dest+PathSep+rdirname)
		monte.removeActionListener(img)

	if VRML:
		tw.flush()
		fos.close()
		monte.removeActionListener(vrml)
	if trajLog:
		monte.removeActionListener(log)
		log.close()
#end of xygettrajectories

#class loopPixels(computeTrajectories): 

# loopPixels uses scanGen but is still untested
def loopPixels(m5vals,xvals,yvals,monte,beamEeV,beamsize,paramlist,chamberdict,detectorRegionList,fileid):
# used for E and B fields
        for m5 in m5vals:
            pixelCount = 0
            while pixelCount < pixelsPerFrame:
                pixelcount += 1				
                [xnm,ynm,znm,t] = scanGen.next() # Get the next beam start position (in nm)
                xygettrajectories(m5,xnm,ynm,pixelcount,monte,beamEeV,beamsize,paramlist,chamberdict,detectorRegionList,fileid)
# end of loopPixels

#class loopPoints(computeTrajectories): 
# iterate over a list of coordinate
def loopPoints(m5vals,xyvallist,monte,beamEeV,beamsize,paramlist,chamberdict,detectorRegionList,fileid):
        for m5 in m5vals:
           for xynm in xyvallist:
                xnm=xynm[1]
                ynm=xynm[2]
                pixelcount=xynm[0]        
                xygettrajectories(m5,xnm,ynm,pixelcount,monte,beamEeV,beamsize,paramlist,chamberdict,detectorRegionList,fileid)
# end of looppoints

def pointsfromxyvals(xvals,yvals,xyvallist):
        pixelcount=0
        for xnm in xvals:
            for ynm in yvals:
                xyvallist.append([pixelcount,xvals,yvals])
                pixelcount += 1
# end of points from xyvals

#class loopXYvals(computeTrajectories):

# iterate over a list of x values and y values
def loopXYvals(m5vals,xvals,yvals,monte,beamEeV,beamsize,paramlist,chamberdict,detectorRegionList,fileid):       
      for m5 in m5vals:
        pixelcount=0
        for xnm in xvals:          #####################To switch x and y, do it here and in 4 lines below
            for ynm in yvals:
                pixelcount +=1
                xygettrajectories(m5,xnm,ynm,pixelcount,monte,beamEeV,beamsize,paramlist,chamberdict,detectorRegionList,fileid)
# end of looptrajectories

# Beginning of main program

# We'll make the substrate infinitely thick.

# Make a record of the random seed we use so we can exactly repeat this calculation
# (same random number sequence) if necessary.
seed = nmu.Math2.rgen.nextLong() # Pick a random seed

# To exactly repeat a previous calculation (e.g., for bug fix) uncomment the next line
# and replace the question marks with the seed that was recorded in the previous calculation's 
# output file.
#seed = -1920631692402242952L
#print >>file,"Random number seed: ",seed
print "Random number seed: ",seed
nmu.Math2.initializeRandom(seed)
for i in range(0,10):
    r = nmu.Math2.rgen.nextDouble()
#    print >>file, r
    print r

radperdeg = jl.Math.PI/180.	# conversion from degrees to radians

t0 = jl.System.currentTimeMillis()

# open large file with all the results
filename = DefaultOutput+PathSep+"ResultsAll.csv"


fileid = open(filename,'w')

#headerfilename = DefaultOutput+PathSep+"Header.txt"
#headerfileid = open(filename,'w')

#get scan region parameters from list of scan region dictionaries
# one dictionary for each scan region
for scanregion in geometry.scanregionlist:
   beamsizenmvals=scanregion["beamsizenmvals"]
   beamEeVvals=scanregion["beamEeVvals"]
   nTrajectories=scanregion["nTrajectories"]
   mExp=scanregion["mExp"]
   trajImg=scanregion["trajImg"]
   trajLog=scanregion["trajLog"]
   histogram=scanregion["histogram"]
# create separate results file for each scan region
   filename = DefaultOutput+PathSep+"Results"+scanregion["Name"]+".txt"
#   file = open(filename,'w')
#   print "Output will be to file: ",filename
#   print >>file, "beamsize(nm)	beamE(eV)   Experiment#	x(nm)	 y(nm)	BSEyield	SEyield	det BSf	detectedFraction	ElapsedTime	variable 1"  
   print         "beamsize(nm)	beamE(eV)   Experiment#	x(nm)	 y(nm)	BSEyield	SEyield	det BSf	detectedFraction	ElapsedTime	variable 1"
   datetime=jl.System.currentTimeMillis()
   starttime=jl.System.currentTimeMillis()
   name="Project 1"
   jobfolder=currentpath
   numdetectors=12
   looplabellist=["thickness","None","None"]
   if "xyvals" in scanregion.keys():
       xyvals=scanregion["xyvals"]
       totalpixels=len(xyvals)
   else:
       xvals=scanregion["xvals"]
       yvals=scanregion["yvals"]
       totalpixels=len(xvals)*len(yvals) #scanregion["npixelsx"]*scanregion["npixelsy"]
   totaltrajs=totalpixels*nTrajectories*mExp
#   savedata.writeheader(fileid,name,datetime,jobfolder,materials_file,definesample_file,blend_file,geometry.GUI_Addon_Version,docker_config,userid,affiliation,starttime,seed,totalpixels,totaltrajs)
#  much of the header info does not need to be in the mail part of the program
   geometry.detectorInfoList[0]["gunZnm"]=scanregion["gunZnm"]
   geometry.detectorInfoList[0]["rotAng"]='x'+str(scanregion["rotAng"][0]*180/math.pi)+'_y'+str(scanregion["rotAng"][1]*180/math.pi)+'_z'+str(scanregion["rotAng"][2]*180/math.pi)+' at '+str(scanregion["pivotloc"]).replace(',',' ')
   # AEP Changed Data Here
   savedata.writeheader(fileid,datetime,geometry.GUI_Addon_Version,seed,totalpixels,totaltrajs,geometry.detectorInfoList,looplabellist)
   chamberdict=geometry.detectorInfoList[0]


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

	for parameter in paramvals:            
	    for thickness in geometry.thicknessvals:                
	# SAMPLE DESCRIPTION
		monte,eg,detectorRegionList=geometry.definesample(beamsize,thickness)
# The list of detectors have to be returned here of the program won't recognize the names of the regions with detectors
#  add thickness and geometry to definesample method call
                chamber = monte.getChamber()
       	        chamber.rotate(scanregion["pivotloc"],math.pi/2,scanregion["rotAng"][0],-math.pi/2)         #Entire sample rotated!
       	        chamber.rotate(scanregion["pivotloc"],0,scanregion["rotAng"][1],scanregion["rotAng"][2])         #Entire sample rotated!                
		for beamEeV in beamEeVvals:                    
		    beamE = epq.ToSI.eV(beamEeV)
		    monte.setBeamEnergy(beamE) # sets this model's beam energy
# Run the simulation
# reduce the number of loops in the main program to make it easier to add loops over other parameters
                    chamberdict["scanregname"]=scanregion["Name"]
                    if "xyvals" in scanregion.keys():
                        loopPoints(m5vals,xyvals,monte,beamEeV,beamsize,[thickness,0,0],chamberdict,detectorRegionList,fileid)
                    else:
		        loopXYvals(m5vals,xvals,yvals,monte,beamEeV,beamsizenm,[thickness,0,0],chamberdict,detectorRegionList,fileid)

print("Ended Writing Status")
statusfilename = DefaultOutput+"Status.csv"
statusfileid = open(statusfilename,'a')
date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
savedata.writeStatusFile(statusfileid, "Ended", date_time)
statusfileid.close()
print("Closed Status File")
fileid.close()
# AEP Command Copy all files to temp volume
# Copy Command Here For All PRN Files
#DestinationPath = '/work/output_qkss/'
#cmd = "cp " + DefaultOutput + ' ' + DestinationPath
cmd = "cp -r /work/output_qkss /output_data"
os.system(cmd)
# Write to File Finalized

#qkss_end_date = datetime.now()
#print("The start time is " + qkss_start_date.strftime(format))
#print("The end time is " + qkss_end_date.strftime(format))
