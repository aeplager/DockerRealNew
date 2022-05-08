####Note great Jython syntax website:  http://www.jython.org/jythonbook/en/1.0/LangSyntax.html#operators
WinOrLin = 1	#=1 for WINDOWS/StandardSetup/homePC/workPC/Prolith1/Marcy/Spock, =2 for LINUX, =3 for WINDOWS/Metrosim group sim computer, =4 for WINDOWS/Prolith2 computer
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

import math, sys
from datetime import datetime
from time import ctime
sys.modules.clear()

max_detectors=12
# determine where to save the results
dest=DefaultOutput;
jio.File(dest).mkdirs()
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
import savedata

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
trajLogMaxTraj = 20
trajlogfile=""

histfile=""

VRML = 0
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

def createBlock(xmin,xmax,ymin,ymax,zmin,zmax):
    block=mon.NormalMultiPlaneShape()
    block.addPlane([-1.,0.,0.],[xmin,0,0])
    block.addPlane([+1.,0.,0.],[xmax,0,0])
    block.addPlane([0.,-1.,0.],[0.,ymin,0])
    block.addPlane([0.,+1.,0.],[0.,ymax,0])
    block.addPlane([0.,0.,-1.],[0,0,zmin])
    block.addPlane([0.,0.,+1.],[0,0,zmax])
    return block

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

#class computeTrajectories(object):  #thinking about changing to object oriented program
#    def __init__(self):
#        pass

def xygettrajectories(mexp,xnm,ynm,pixelnum,imagenum,eg,monte,beamEeV,beamsize,paramlist,chamberdict,detectorList,fileid):
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
                trajImgFile="traj_"+str(pixelnum)+".png"

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
                trajlogfile=dest+PathSep+"trajlog_"+str(pixelnum)+".log"
	        log=mon.TrajectoryLoggingListener(monte,trajlogfile)
	        log.setMaxTrajectories(trajLogMaxTraj)
	        monte.addActionListener(log)
 
        tbefore = (jl.System.currentTimeMillis()-t0)/1000.
	monte.runMultipleTrajectories(nTrajectories)
        SEf,bsf=getSfandBsf(back,beamEeV,nTrajectories)
        Xfchamber=XfChamberBins(back,beamEeV,nTrajectories,chamberdict["XfminE"],chamberdict["XfmaxE"])

	t = (jl.System.currentTimeMillis()-t0)/3600000.
	rdirname = jl.Integer(int(beamEeV)).toString()+"eV_"+"x"+jl.Integer(int(xnm)).toString()+'y'+jl.Integer(int(ynm)).toString()
        histfile=dest+PathSep+"histogram_"+str(pixelnum)+".prn"
        if not histogram:
            histfile=dest+PathSep+"no_histogram"
        rotbeam="x0_y0_z0"
        energydep="NA"
        savedata.writeChamberData(fileid,pixelnum,imagenum,t*3600,t*3600-tbefore,xnm,ynm,mexp,nTrajectories,beamEeV,beamsize,beamsize,rotbeam,energydep,histfile,trajImgFile,trajlogfile,bsf,SEf,Xfchamber,paramlist)
        tprevsec=t*3600

        detBsf=0
        detBSf2 = 0
        for detector in detectorList:
            detSEf,detBsf,detXf=getSfandBsfRegion(detector["back_regdet"],beamEeV,nTrajectories,detector["minbin"],detector["maxbin"])
            detBSf2=detector["back_regdet"].detectedFraction()
 #           histfile is the same as for chamber data
            savedata.writeDetectorData(fileid,detBsf,detSEf,detXf,histfile,trajlogfile)
        for nodetector in range(len(detectorList),max_detectors):
            savedata.writeDetectorData(fileid,0,0,0,"NA","NA")
          
        savedata.endline(fileid)

#	print >>file, "%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%8.4f \t%8.2f" % (beamsizenm,beamEeV,mexp,xnm,ynm,bsf,SEf,t,paramlist[0])
	print "%8.2f \t%8.1f \t%8.1f \t%8.2f \t%8.2f \t%3.5f \t%3.5f \t%3.5f \t%3.5f \t%12.6f \t%8.2f" % (beamsizenm,beamEeV,mexp,xnm,ynm,bsf,SEf,detBsf,detBSf2,t,paramlist[0])

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
def loopPixels(m5vals,xvals,yvals,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid,pixelcount,imagenum):
        print("loopPixels is not working yet")
# used for E and B fields
        for m5 in m5vals:
            while pixelcount < pixelsPerFrame:
                pixelcount += 1				
                [xnm,ynm,znm,t] = scanGen.next() # Get the next beam start position (in nm)
                xygettrajectories(m5,xnm,ynm,pixelcount,imagenum,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid)
        return(pixelcount)
# end of loopPixels

#class loopPoints(computeTrajectories): 
# iterate over a list of coordinate
def loopPoints(m5vals,xyvallist,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid,pixelcount,imagenum):
        for m5 in m5vals:
           for xynm in xyvallist:
                xnm=xynm[1]
                ynm=xynm[2]
                pixelcount += 1 #=xynm[0]        
                xygettrajectories(m5,xnm,ynm,pixelcount,imagenum,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid)
        return(pixelcount)
# end of looppoints

def pointsfromxyvals(m5vals,xvals,yvals,xyvallist,pixelcount):
    for m5 in m5vals:
        for xnm in xvals:
            for ynm in yvals:
                xyvallist.append([pixelcount,xvals,yvals])
                pixelcount+=1
# end of points from xyvals

#class loopXYvals(computeTrajectories):

# iterate over a list of x values and y values
def loopXYvals(m5vals,xvals,yvals,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid,pixelcount,imagenum):
      for m5 in m5vals:
        for xnm in xvals:          #####################To switch x and y, do it here and in 4 lines below
            monte,eg,detectorRegionList=geometry.definesample(beamsizenm*meterspernm,thickness)
	    beamE = epq.ToSI.eV(beamEeV)
	    monte.setBeamEnergy(beamE) # sets this model's beam energy
            chamber = monte.getChamber()
#Entire sample rotated!
       	    chamber.rotate(chamberdict["pivotloc"],0,chamberdict["rotAng"][1],chamberdict["rotAng"][2])  
       	    chamber.rotate(chamberdict["pivotloc"],math.pi/2,chamberdict["rotAng"][0],-math.pi/2) 

            for ynm in yvals:
                pixelcount +=1
                xygettrajectories(m5,xnm,ynm,pixelcount,imagenum,eg,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid)
      return (pixelcount)
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
nmu.Math2.initializeRandom(seed)
print "Random number seed: ",seed
for i in range(0,10):
    r = nmu.Math2.rgen.nextDouble()
#    print >>file, r
    print r

radperdeg = jl.Math.PI/180.	# conversion from degrees to radians

t0 = jl.System.currentTimeMillis()

# open large file with all the results--have one file per scan region
#filename = DefaultOutput+PathSep+"ResultsAll.csv"
#fileid = open(filename,'w')

#date and time of start of project
dateTimeObj = datetime.now()
datestr = str(dateTimeObj.year)+ '/'+str(dateTimeObj.month)+'/'+str( dateTimeObj.day)
timestr=str(dateTimeObj.hour)+ ':'+str( dateTimeObj.minute)+ ':'+str( dateTimeObj.second)+ '.'+str(dateTimeObj.microsecond)
projectstart=datestr+'. '+timestr

#pixelnumber is unique for each pixel, frame and scan region in the simulation
pixelCount = 0

# compute total number of pixels in each scan region
totalimages=0
imagenum=0
totalpixels = 0
for scanregion in geometry.scanregionlist:
    mExp=scanregion["mExp"]
    beamsizenmvals=scanregion["beamsizenmvals"]
    beamEeVvals=scanregion["beamEeVvals"]
    images_region=mExp*len(beamsizenmvals)*len(beamEeVvals) #*len(geometry.thicknessvals)
    for loop_key in geometry.loopdict:
        images_region=images_region*len(geometry.loopdict[loop_key])
    if "xyvals" in scanregion.keys():
        totalpixels=len(scanregion["xyvals"])*mExp*len(geometry.loopdict.values()[0])*len(beamsizenmvals)*len(beamEeVvals)
    else:
        pixels_region=len(scanregion["xvals"])*len(scanregion["yvals"])*mExp*len(geometry.loopdict.values()[0])*len(beamsizenmvals)*len(beamEeVvals) #scanregion["npixelsx"]*scanregion["npixelsy"]
    scanregion["pixels_region"]=pixels_region
    scanregion["images_region"]=images_region
    totalimages=totalimages+images_region
    totalpixels = totalpixels+pixels_region

#get scan region parameters from list of scan region dictionaries
# one dictionary for each scan region
for scanregion in geometry.scanregionlist:
# for now reset pixelCount for each scan region, but fix in near future
#   pixelCount=0
#   startpixel=0
   scanregion["startpixel"]=pixelCount+1
   beamsizenmvals=scanregion["beamsizenmvals"]
   beamEeVvals=scanregion["beamEeVvals"]
   nTrajectories=scanregion["nTrajectories"]
   mExp=scanregion["mExp"]
   trajImg=scanregion["trajImg"]
   trajLog=scanregion["trajLog"]
   histogram=scanregion["histogram"]
# create separate results file for each scan region
   filename = DefaultOutput+PathSep+"Results_"+scanregion["Name"]+"_header"+".csv"
   fileid = open(filename,'w')
#   print >>file, "beamsize(nm)	beamE(eV)   Experiment#	x(nm)	 y(nm)	BSEyield	SEyield	det BSf	detectedFraction	ElapsedTime	variable 1"  
   starttime=jl.System.currentTimeMillis()
   dateTimeObj = datetime.now()
   datestr = str(dateTimeObj.year)+ '/'+str(dateTimeObj.month)+'/'+str( dateTimeObj.day)
   timestr=str(dateTimeObj.hour)+ ':'+str( dateTimeObj.minute)+ ':'+str( dateTimeObj.second)+ '.'+str(dateTimeObj.microsecond)
   regionstart=datestr+'. '+timestr
   name="Project 1"
   jobfolder=currentpath
   numdetectors=12
   looplabellist=["Loop1_thickness","Loop2_NA","Loop3_NA"]
   scanregioninfo={'Scan_Region_Name':scanregion["Name"],'pixel_size(x_y)':str(scanregion["deltax"])+', '+str(scanregion["deltay"]),\
       'images/region':scanregion["images_region"],'pixels/region':scanregion["pixels_region"],'pixels/image': 0,'frames': mExp, \
       'trajectories':scanregion["nTrajectories"],'center(x_y)':'0,0','start(x_y)':'0,0','end(x_y)':'0,0','pixels(ny_ny)':'0,0','filename':'none'}
   if "xyvals" in scanregion.keys():
       xyvals=scanregion["xyvals"]
       scanregioninfo["pixels/image"]=len(xyvals)
       scanregioninfo['pixels(nx_ny)'] = ' '+str(len(xyvals))
       scanregioninfo['start(x_y)'] = ' '+str(xvals[0])+', '+str(yvals[0])
       scanregioninfo['end(x_y)']= ' '+str(xyvals[-1][1])+','+str(xyvals[-1][2])
   else:
       xvals=scanregion["xvals"]
       yvals=scanregion["yvals"]
       scanregioninfo["pixels/image"]=len(xvals)*len(yvals)
       scanregioninfo['pixels(nx_ny)'] = ' '+str(len(xvals))+', '+str(len(yvals))
       scanregioninfo['start(x_y)'] = ' '+str(xvals[0])+', '+str(yvals[0])
       scanregioninfo['end(x_y)']= ' '+str(xvals[-1])+', '+str(yvals[-1])
   if len(xvals)==1 and "nx" in scanregion:
       xmin = xvals[0]
       ymin = yvals[0]
       deltax = scanregion["deltax"]
       deltay = scanregion["deltay"]
       nx = scanregion["nx"]
       ny = scanregion["ny"]
       scanregioninfo["pixels/image"]=nx*ny
# The times in the following ScanGenerator are not used.
#pixelDwell - time between pixels along the x direction
#retraceTime - time between completion of the last pixel of a line and start of the first pixel of the next line
#frameSettleTime - time between completion of the last line of a frame and the start of the first line of the next frame
#The time of the ith point is t = t0 + ix*pixelDwell + iy*(nx*pixelDwell+retraceTime) + if*(frameSettleTime + ny*(retraceTime + nx*dwellTime))
       t0=0
       pixelDwell = 1.e-6
       retraceTime = 5.e-6
       frameSettleTime = 10.e-6
       xmin=xvals[0]
       ymin=yvals[0]
       scanGen = mon.SimpleRasterScanGenerator(xmin, ymin, -gunZnm, deltax, deltay, nx, ny, t0, pixelDwell, retraceTime, frameSettleTime)
       scanregioninfo['pixels(ny_ny)'] = ' '+str(nx)+', '+str(ny)
       scanregioninfo['start(x_y)'] = ' '+str(xmin)+', '+str(ymin)

# create scanregioninfo dictionary with info for header
#   scanregioninfo['center(x_y)']= ' '+str(scanregion["centerx"])+', '+str(scanregion["centery"])
   scanregioninfo['Beam eV']=str(beamEeVvals).replace(',',' ')
   scanregioninfo['beamsize']=str(beamsizenmvals).replace(',',' ')

   loopdict={looplabellist[2]: [0], looplabellist[1]: [0]}
   loopdict[looplabellist[0]]=geometry.loopdict['thicknessvals']
   thicknessvals=geometry.loopdict['thicknessvals']

   if 'filename' in scanregion.keys():
       scanregioninfo['filename']=scanregion["filename"]
   print "total pixels = ",totalpixels
   print "scan region: ",scanregioninfo
   print         "beamsize(nm)	beamE(eV)   Experiment#	x(nm)	 y(nm)	BSEyield	SEyield	det BSf	detectedFraction	ElapsedTime	variable 1"
#  much of the header info does not need to be in the main part of the program
   geometry.detectorInfoList[0]["gunZnm"]=scanregion["gunZnm"]
   filedict=geometry.file_dict
   filedict["blend_modify"] = ctime(filedict["blend_epoch"])
   geometry.detectorInfoList[0]["rotAng"]='x'+str(scanregion["rotAng"][0]*180/math.pi)+'_y'+str(scanregion["rotAng"][1]*180/math.pi)+'_z'+str(scanregion["rotAng"][2]*180/math.pi)+' at '+str(scanregion["pivotloc"]).replace(',',' ')
   savedata.writeheader(fileid,[projectstart,regionstart],geometry.GUI_Addon_Version,seed,[totalpixels,totalimages],scanregion["startpixel"],geometry.detectorInfoList,scanregioninfo,filedict,geometry.commentList,loopdict)
   chamberdict=geometry.detectorInfoList[0]
   fileid.close()  #finished writing header and start new file for data

   filename = DefaultOutput+PathSep+"Results_"+scanregion["Name"]+"_data"+".csv"
   fileid = open(filename,'w')
   savedata.write_column_titles(fileid,looplabellist)
   chamberdict["rotAng"]=scanregion["rotAng"]
   chamberdict["pivotloc"]=scanregion["pivotloc"]
# define multiple experiments??
   m5=1
   deltam5=1
   m5vals = []
   while m5-1 < mExp:
       m5vals.append(m5)
       m5 += deltam5

#loop beamsizenmvals
   for beamsizenm in beamsizenmvals:
	beamsize = beamsizenm*meterspernm

	for parameter in paramvals:
	    for thickness in thicknessvals:
                imagenum += 1

	# SAMPLE DESCRIPTION
		monte,eg,detectorRegionList=geometry.definesample(beamsize,thickness)

# The list of detectors have to be returned here of the program won't recognize the names of the regions with detectors
#  add thickness and geometry to definesample method call
                chamber = monte.getChamber()
#Entire sample rotated!
       	        chamber.rotate(scanregion["pivotloc"],0,scanregion["rotAng"][1],scanregion["rotAng"][2])  
       	        chamber.rotate(scanregion["pivotloc"],math.pi/2,scanregion["rotAng"][0],-math.pi/2) 
		for beamEeV in beamEeVvals:
		    beamE = epq.ToSI.eV(beamEeV)
		    monte.setBeamEnergy(beamE) # sets this model's beam energy
# Run the simulation
# reduce the number of loops in the main program to make it easier to add loops over other parameters
                    chamberdict["scanregname"]=scanregion["Name"]
                    if "xyvals" in scanregion.keys():
                        pixelCount=loopPoints(m5vals,xyvals,monte,beamEeV,beamsizenm,[thickness,0,0],chamberdict,detectorRegionList,fileid,pixelCount,imagenum)
                    elif "nx" in scanregion.keys():
                        pixelCount=loopPixels(m5vals,xvals,yvals,monte,beamEeV,beamsizenm,[thickness,0,0],chamberdict,detectorRegionList,fileid,pixelCount,imagenum)                      
                    else:
#                        for m5 in m5vals:
#                            for xnm in xvals:          #####################To switch x and y, do it here and in 4 lines below
#                                monte,eg,detectorRegionList=geometry.definesample(beamsizenm*meterspernm,thickness)
#	                        beamE = epq.ToSI.eV(beamEeV)
#	                        monte.setBeamEnergy(beamE) # sets this model's beam energy
#                                for ynm in yvals:
#                                    pixelCount +=1
#                                    paramlist=[thickness,0,0]
#                                    xygettrajectories(m5,xnm,ynm,pixelCount,imagenum,monte,beamEeV,beamsizenm,paramlist,chamberdict,detectorRegionList,fileid)
		        pixelCount=loopXYvals(m5vals,xvals,yvals,monte,beamEeV,beamsizenm,[thickness,0,0],chamberdict,detectorRegionList,fileid,pixelCount,imagenum)

   fileid.close()
print ('Successful finish!')

