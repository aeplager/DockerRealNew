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
	SiMSM, SiMSMDeep = mat.SiModel()
	AuMSM, AuMSMDeep = mat.AuModel()

	chamber.updateMaterial(chamber.getScatterModel(),vacuumMSM)

# Block
	location = [ 18.618*meterspernm, 50*meterspernm, -7.500*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [ 31.550*meterspernm, 30.000*meterspernm, 15.000*meterspernm]
	block=mon.NormalMultiPlaneShape()
	block.addPlane([-1.,0.,0.],[-dimensions[0]/2,0,0])
	block.addPlane([+1.,0.,0.],[dimensions[0]/2,0,0])
	block.addPlane([0.,-1.,0.],[0.,-dimensions[1]/2,0])
	block.addPlane([0.,+1.,0.],[0.,dimensions[1]/2,0])
	block.addPlane([0.,0.,-1.],[0,0,-dimensions[2]])
	block.addPlane([0.,0.,+1.],[0,0,0])
	block.translate([location[0],location[1],location[2]+dimensions[2]/2])
	block.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	block.rotate(location,0,rotation[1],rotation[2])
	blockregion = monte.addSubRegion(chamber,AuMSM,block)

# Cylinder
	location = [ 18.287*meterspernm,  12*meterspernm, -5.000*meterspernm]
	rotation = [  0.000*math.pi,  0.500*math.pi,  0.000*math.pi]
	dimensions = [ 10.000*meterspernm, 10.000*meterspernm, 30.000*meterspernm]
	cylinderheight=dimensions[2]
	bottomcenter=[0,0,0-cylinderheight/2.0]
	topcenter=[0,0,0+cylinderheight/2.0]
	radius=5.0*meterspernm
	cylinder= mon.NormalCylindricalShape(bottomcenter,topcenter,radius)
	cylinder.translate([location[0],location[1],location[2]])
	cylinder.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	cylinder.rotate(location,0,rotation[1],rotation[2])
	cylinderregion = monte.addSubRegion(chamber,AuMSM,cylinder)

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
	"rotAng":[  0.00*radperdeg, 20.00*radperdeg,  0.00*radperdeg],	"pivotloc":[  2.00*meterspernm,  2.00*meterspernm,  0.00*meterspernm],
# create x and y vectors for simulation
	"xvals": [-1.00,3.00,7.00,11.00,15.00,19.00,23.00,27.00,31.00,35.00,39.00,43.00,47.00,51.00,\
	55.00,59.00,63.00,67.00,],
	"yvals": [-1.00,3.00,7.00,11.00,15.00,19.00,23.00,27.00,31.00,35.00,39.00,43.00,47.00,51.00,\
	55.00,59.00,63.00,67.00,],
 # Beam energies in eV
	"beamEeVvals" : [  500.0,],
 # beam size in nm
	"beamsizenmvals" : [  0.5,],
}
scanregionlist.append(scanregiondict)
thicknessvals=[0.0]
GUI_Addon_Version=(0, 9)
