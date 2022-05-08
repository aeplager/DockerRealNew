# This script determines the yield for electrons incident on 1 or more
# trapezoidal (with top corner radii) resist lines on a 3-layer substrate.
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.nanoscalemetrology.JMONSEL as mon
import java.lang as jl
import math

import sys

# conversion factors
radperdeg = jl.Math.PI/180.	# conversion from degrees to radians
meterspernm = 1.e-9	# conversion from nanometers to meters


def definesample(beamsize,thickness):
	meterspernm = 1.e-9	# conversion from nanometers to meters
	radperdeg = jl.Math.PI/180.	# conversion from degrees to radians
	normalVector = [0.,0.,-1.]

# SAMPLE DESCRIPTION

# NISTMonte provides us with a "chamber" in the form of a 0.1 m sphere inside of which we build
# out sample. Replace the default vacuum in the chamber with SEvacuum. (SEmaterials define additional
# properties, such as work function, that are needed by JMONSEL.)

	eg = nm.GaussianBeam(beamsize) # makes electron gun, Gaussian with standard deviation = beamsize

# create an instance of the model
	monte=nm.MonteCarloSS() #creates an instance of the model with all default characteristics
	monte.setElectronGun(eg) # This gun is "attached" to the model.
	chamber = monte.getChamber()

#  Define materials


	import makematerials as mat
	vacuumMSM=mat.vacuumModel()
	AuMSM, AuMSMDeep = mat.AuModel()
	SiMSM, SiMSMDeep = mat.SiModel()

	chamber.updateMaterial(chamber.getScatterModel(),vacuumMSM)

# Pyramid
	location = [  0.000*meterspernm,  0.000*meterspernm,-15.000*meterspernm]
	rotation = [  1.000*math.pi,-  0.000*math.pi,  0.250*math.pi]
	dimensions = [ 40.000*meterspernm, 40.000*meterspernm, 30.000*meterspernm]
	pyramidsides=4
	height=dimensions[2]
	radius=20.0*meterspernm
	vertexZ=30.000001907348633*meterspernm
	pyramid=mat.createPyramid(pyramidsides,height,2*radius,vertexZ)
	pyramid.translate([location[0],location[1],location[2]-dimensions[2]/2])
	pyramid.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	pyramid.rotate(location,0,rotation[1],rotation[2])
	pyramidregion = monte.addSubRegion(chamber,AuMSM,pyramid)

# Layer.s1
	location = [  0.000*meterspernm,  0.000*meterspernm,250.000*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm,500.000*meterspernm]
	thickness1=location[2]+dimensions[2]/2
	layers1=mon.NormalMultiPlaneShape()
	layers1.addPlane(normalVector,[0,0,0])
	layers1.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	layers1.rotate(location,0,rotation[1],rotation[2])
	layers1region = monte.addSubRegion(chamber,SiMSM,layers1)

# Layer.s2
	location = [  0.000*meterspernm,  0.000*meterspernm,525.000*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm, 50.000*meterspernm]
	thickness2=location[2]+dimensions[2]/2
	layers2=mon.NormalMultiPlaneShape()
	layers2.addPlane(normalVector,[0,0,thickness1])
	layers2.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	layers2.rotate(location,0,rotation[1],rotation[2])
	layers2region = monte.addSubRegion(layers1region,SiMSMDeep,layers2)

# Layer.s3
	location = [  0.000*meterspernm,  0.000*meterspernm,557.500*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm, 15.000*meterspernm]
	thickness3=location[2]+dimensions[2]/2
	layers3=mon.NormalMultiPlaneShape()
	layers3.addPlane(normalVector,[0,0,thickness2])
	layers3.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	layers3.rotate(location,0,rotation[1],rotation[2])
	layers3region = monte.addSubRegion(layers2region,SiMSMDeep,layers3)

	detectorList=[]
	return(monte,eg,detectorList)

detectorInfoList=[]

detectorInfoList=[{'Name': 'chamber', 'XfminE': 0.0, 'XfmaxE': 50.0, 'rotAng': 'x0_y0_z0'}]

#save scan region and beam parameters
scanregionlist=[]
scanregiondict={
	"Name":"Scan_Region", 	"trajImg":0,	"trajLog":0,	"histogram":0,	"gunZnm":50,	"nTrajectories":100,	"mExp":1,
	"rotAng":[  0.00*radperdeg,  0.00*radperdeg,  0.00*radperdeg],	"pivotloc":[  0.00*meterspernm,  0.00*meterspernm,  0.00*meterspernm],
# create x and y vectors for simulation
	"xvals": [-23.00,-21.00,-19.00,-17.00,-15.00,-13.00,-11.00,-9.00,-7.00,-5.00,-3.00,-1.00,1.00,3.00,\
	5.00,7.00,9.00,11.00,13.00,15.00,17.00,19.00,21.00,],
	"yvals": [-23.00,-21.00,-19.00,-17.00,-15.00,-13.00,-11.00,-9.00,-7.00,-5.00,-3.00,-1.00,1.00,3.00,\
	5.00,7.00,9.00,11.00,13.00,15.00,17.00,19.00,21.00,],
 # Beam energies in eV
	"beamEeVvals" : [  500.0,],
 # beam size in nm
	"beamsizenmvals" : [  0.5,],
}
scanregionlist.append(scanregiondict)
thicknessvals=[0.0]
GUI_Addon_Version=(0, 9)
