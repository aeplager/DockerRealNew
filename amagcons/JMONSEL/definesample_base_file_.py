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
	WMSM, WMSMDeep = mat.WModel()
	FeMSM, FeMSMDeep = mat.FeModel()
	AgMSM, AgMSMDeep = mat.AgModel()
	AuMSM, AuMSMDeep = mat.AuModel()
	SiMSM, SiMSMDeep = mat.SiModel()

	chamber.updateMaterial(chamber.getScatterModel(),vacuumMSM)

# Ring
	location = [  0.000*meterspernm,-190.000*meterspernm,-300.500*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm,  1.000*meterspernm]
	cylinderheight=dimensions[2]
	bottomcenter=[0,0,0-cylinderheight/2.0]
	topcenter=[0,0,0+cylinderheight/2.0]
	radinner=50.0*meterspernm
	radius=250.0*meterspernm
	angle=360*radperdeg
	cylinderbig= mon.NormalCylindricalShape(bottomcenter,topcenter,radius)
	cylindersmall= mon.NormalCylindricalShape(bottomcenter,topcenter,radinner)
	ring=mon.NormalDifferenceShape(cylinderbig,cylindersmall)
	ring.translate([location[0],location[1],location[2]])
	ring.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	ring.rotate(location,0,rotation[1],rotation[2])
	ringregion = monte.addSubRegion(chamber,vacuumMSM,ring)

# Line
	location = [ 11.068*meterspernm,-137.032*meterspernm, -5.000*meterspernm]
	rotation = [  0.500*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [ 30.000*meterspernm, 10.000*meterspernm, 20.000*meterspernm]
	h=-dimensions[2]
	w=30.0*meterspernm
	linelength=dimensions[1]
	rangle=15.0*radperdeg
	radr=2.0*meterspernm
	langle=10.0*radperdeg
	radl=2.0*meterspernm
	line=mon.NShapes.createLine(h,w,linelength,langle,rangle,radl,radr)
	line.translate([location[0],location[1],location[2]-location[2]/abs(location[2])*dimensions[2]/2])
	line.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	line.rotate(location,0,rotation[1],rotation[2])
	lineregion = monte.addSubRegion(chamber,FeMSM,line)

# Pyramid.001
	location = [-19.011*meterspernm,-137.886*meterspernm,-10.000*meterspernm]
	rotation = [  1.000*math.pi,-  0.000*math.pi,  0.000*math.pi]
	dimensions = [ 20.000*meterspernm, 20.000*meterspernm, 20.000*meterspernm]
	pyramidsides=4
	height=dimensions[2]
	radius=10.0*meterspernm
	vertexZ=20.0*meterspernm
	pyramid001=mat.createPyramid(pyramidsides,height,2*radius,vertexZ)
	pyramid001.translate([location[0],location[1],location[2]-dimensions[2]/2])
	pyramid001.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	pyramid001.rotate(location,0,rotation[1],rotation[2])
	pyramid001region = monte.addSubRegion(chamber,AgMSM,pyramid001)

# Sphere
	location = [ -4.088*meterspernm,-167.962*meterspernm, 27.111*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [ 20.000*meterspernm, 20.000*meterspernm, 20.000*meterspernm]
	radius=10.0*meterspernm
	sphere= mon.NormalSphereShape([0,0,0],radius)
	sphere.translate([location[0],location[1],location[2]])
	sphere.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	sphere.rotate(location,0,rotation[1],rotation[2])
	sphereregion = monte.addSubRegion(chamber,AuMSM,sphere)

# Cylinder
	location = [ -0.216*meterspernm,-207.394*meterspernm, 74.932*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [  6.000*meterspernm,  6.000*meterspernm, 10.000*meterspernm]
	cylinderheight=dimensions[2]
	bottomcenter=[0,0,0-cylinderheight/2.0]
	topcenter=[0,0,0+cylinderheight/2.0]
	radius=3.0*meterspernm
	cylinder= mon.NormalCylindricalShape(bottomcenter,topcenter,radius)
	cylinder.translate([location[0],location[1],location[2]])
	cylinder.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	cylinder.rotate(location,0,rotation[1],rotation[2])
	cylinderregion = monte.addSubRegion(chamber,WMSM,cylinder)

# Layer.s1
	location = [  0.000*meterspernm,  0.000*meterspernm,250.000*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm,500.000*meterspernm]
	thickness1=location[2]+dimensions[2]/2
	layers1=mon.NormalMultiPlaneShape()
	layers1.addPlane(normalVector,[0,0,0])
	layers1.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	layers1.rotate(location,0,rotation[1],rotation[2])
# Block
	location = [  0.000*meterspernm,-203.303*meterspernm, 17.890*meterspernm]
	rotation = [  0.000*math.pi,  0.000*math.pi,  0.000*math.pi]
	dimensions = [ 40.000*meterspernm,105.000*meterspernm, 42.000*meterspernm]
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
	layers1=mon.NormalDifferenceShape(layers1,block)
# Pyramid
	location = [  0.000*meterspernm,-207.564*meterspernm, 59.015*meterspernm]
	rotation = [  0.000*math.pi,-  0.000*math.pi,  0.000*math.pi]
	dimensions = [ 21.000*meterspernm, 21.000*meterspernm, 42.000*meterspernm]
	pyramidsides=12
	height=dimensions[2]
	radius=10.0*meterspernm
	vertexZ=80.0*meterspernm
	pyramid=mat.createPyramid(pyramidsides,height,2*radius,vertexZ)
	pyramid.translate([location[0],location[1],location[2]-dimensions[2]/2])
	pyramid.rotate(location,-math.pi/2,rotation[0],math.pi/2)
	pyramid.rotate(location,0,rotation[1],rotation[2])
	layers1=mon.NormalDifferenceShape(layers1,pyramid)
# add subregion with a boolean modifier(s)
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
	detectordict={'Name': ringregion, 'minE': 50.0, 'maxE': 100.0, 'destruct': '0', 'minbin': 0.0, 'maxbin': 50.0}
	detectorList.append(detectordict)
	return(monte,eg,detectorList)

detectorInfoList=[]

detectorInfoList=[{'Name': 'chamber', 'XfminE': 50.0, 'XfmaxE': 100.0, 'rotAng': 'x0_y0_z0'}]
detectorInfo={'Name': 'ring', 'minE': 50.0, 'maxE': 100.0, 'destruct': '0', 'minbin': 0.0, 'maxbin': 50.0}
detectorInfoList.append(detectorInfo)
#save scan region and beam parameters
scanregionlist=[]
scanregiondict={
	"Name":"Scan_Region", 	"trajImg":0,	"trajLog":0,	"histogram":0,	"gunZnm":50,	"nTrajectories":100,	"mExp":1,
	"rotAng":[  0.00*radperdeg,  0.00*radperdeg,  0.00*radperdeg],	"pivotloc":[  0.00*meterspernm,  0.00*meterspernm,  0.00*meterspernm],
# create x and y vectors for simulation
	"xvals": [-31.00,-29.00,-27.00,-25.00,-23.00,-21.00,-19.00,-17.00,-15.00,-13.00,-11.00,-9.00,-7.00,-5.00,\
	-3.00,-1.00,1.00,3.00,5.00,7.00,9.00,11.00,13.00,15.00,17.00,19.00,21.00,23.00,\
	25.00,27.00,29.00,31.00,],
	"yvals": [-256.00,-254.00,-252.00,-250.00,-248.00,-246.00,-244.00,-242.00,-240.00,-238.00,-236.00,-234.00,-232.00,-230.00,\
	-228.00,-226.00,-224.00,-222.00,-220.00,-218.00,-216.00,-214.00,-212.00,-210.00,-208.00,-206.00,-204.00,-202.00,\
	-200.00,-198.00,-196.00,-194.00,-192.00,-190.00,-188.00,-186.00,-184.00,-182.00,-180.00,-178.00,-176.00,-174.00,\
	-172.00,-170.00,-168.00,-166.00,-164.00,-162.00,-160.00,-158.00,-156.00,-154.00,-152.00,-150.00,-148.00,-146.00,\
	-144.00,-142.00,-140.00,-138.00,-136.00,-134.00,-132.00,-130.00,-128.00,-126.00,-124.00,],
 # Beam energies in eV
	"beamEeVvals" : [  500.0,],
 # beam size in nm
	"beamsizenmvals" : [  0.5,],
}
scanregionlist.append(scanregiondict)
thicknessvals=[0.0]
GUI_Addon_Version=(0, 9)
