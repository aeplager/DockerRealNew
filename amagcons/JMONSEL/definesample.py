# This script determines the yield for electrons incident a user created
# sample and substrate
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.nanoscalemetrology.JMONSEL as mon
import java.lang as jl
import java.util as jutil
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


# This addon was last updated on 1/6/2022
	import makematerials as mat
	vacuumMSM=mat.vacuumModel()
	AuMSM, AuMSMDeep = mat.AuModel()
	SiMSM, SiMSMDeep = mat.SiModel()

	chamber.updateMaterial(chamber.getScatterModel(),vacuumMSM)

# Sphere.011
	location = [-50.000*meterspernm,-25.000*meterspernm, -3.000*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.000*math.pi, -0.000*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere011= mon.NormalSphereShape([0,0,0],radius)
	sphere011= mon.AffinizedNormalShape(sphere011)
	sphere011.scale(stretch[0],stretch[1],stretch[2])
	sphere011.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere011.rotate(pivot,0,rotation[1],rotation[2])
	sphere011.translate([location[0],location[1],location[2]])
	sphere011region = monte.addSubRegion(chamber,AuMSM,sphere011)

# Sphere.010
	location = [-50.000*meterspernm, 25.000*meterspernm, -6.000*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [ -0.000*math.pi,  0.500*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere010= mon.NormalSphereShape([0,0,0],radius)
	sphere010= mon.AffinizedNormalShape(sphere010)
	sphere010.scale(stretch[0],stretch[1],stretch[2])
	sphere010.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere010.rotate(pivot,0,rotation[1],rotation[2])
	sphere010.translate([location[0],location[1],location[2]])
	sphere010region = monte.addSubRegion(chamber,AuMSM,sphere010)

# Sphere.009
	location = [-50.000*meterspernm,  0.000*meterspernm, -4.731*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.000*math.pi,  0.250*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere009= mon.NormalSphereShape([0,0,0],radius)
	sphere009= mon.AffinizedNormalShape(sphere009)
	sphere009.scale(stretch[0],stretch[1],stretch[2])
	sphere009.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere009.rotate(pivot,0,rotation[1],rotation[2])
	sphere009.translate([location[0],location[1],location[2]])
	sphere009region = monte.addSubRegion(chamber,AuMSM,sphere009)

# Sphere.008
	location = [-25.000*meterspernm,-25.000*meterspernm, -6.691*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.250*math.pi, -0.000*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere008= mon.NormalSphereShape([0,0,0],radius)
	sphere008= mon.AffinizedNormalShape(sphere008)
	sphere008.scale(stretch[0],stretch[1],stretch[2])
	sphere008.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere008.rotate(pivot,0,rotation[1],rotation[2])
	sphere008.translate([location[0],location[1],location[2]])
	sphere008region = monte.addSubRegion(chamber,AuMSM,sphere008)

# Sphere.007
	location = [-25.000*meterspernm, 25.000*meterspernm, -6.000*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.250*math.pi,  0.500*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere007= mon.NormalSphereShape([0,0,0],radius)
	sphere007= mon.AffinizedNormalShape(sphere007)
	sphere007.scale(stretch[0],stretch[1],stretch[2])
	sphere007.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere007.rotate(pivot,0,rotation[1],rotation[2])
	sphere007.translate([location[0],location[1],location[2]])
	sphere007region = monte.addSubRegion(chamber,AuMSM,sphere007)

# Sphere.006
	location = [-25.000*meterspernm,  0.000*meterspernm, -6.356*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.250*math.pi,  0.250*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere006= mon.NormalSphereShape([0,0,0],radius)
	sphere006= mon.AffinizedNormalShape(sphere006)
	sphere006.scale(stretch[0],stretch[1],stretch[2])
	sphere006.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere006.rotate(pivot,0,rotation[1],rotation[2])
	sphere006.translate([location[0],location[1],location[2]])
	sphere006region = monte.addSubRegion(chamber,AuMSM,sphere006)

# Sphere.005
	location = [ 25.000*meterspernm,-25.000*meterspernm, -3.000*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.000*math.pi, -0.000*math.pi,  0.250*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere005= mon.NormalSphereShape([0,0,0],radius)
	sphere005= mon.AffinizedNormalShape(sphere005)
	sphere005.scale(stretch[0],stretch[1],stretch[2])
	sphere005.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere005.rotate(pivot,0,rotation[1],rotation[2])
	sphere005.translate([location[0],location[1],location[2]])
	sphere005region = monte.addSubRegion(chamber,AuMSM,sphere005)

# Sphere.004
	location = [ 25.000*meterspernm, 25.000*meterspernm, -4.731*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.000*math.pi,  0.250*math.pi,  0.250*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere004= mon.NormalSphereShape([0,0,0],radius)
	sphere004= mon.AffinizedNormalShape(sphere004)
	sphere004.scale(stretch[0],stretch[1],stretch[2])
	sphere004.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere004.rotate(pivot,0,rotation[1],rotation[2])
	sphere004.translate([location[0],location[1],location[2]])
	sphere004region = monte.addSubRegion(chamber,AuMSM,sphere004)

# Sphere.003
	location = [ 25.000*meterspernm,  0.000*meterspernm, -6.691*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.250*math.pi, -0.000*math.pi,  0.250*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere003= mon.NormalSphereShape([0,0,0],radius)
	sphere003= mon.AffinizedNormalShape(sphere003)
	sphere003.scale(stretch[0],stretch[1],stretch[2])
	sphere003.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere003.rotate(pivot,0,rotation[1],rotation[2])
	sphere003.translate([location[0],location[1],location[2]])
	sphere003region = monte.addSubRegion(chamber,AuMSM,sphere003)

# Sphere.002
	location = [  0.000*meterspernm,-25.000*meterspernm, -9.000*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.500*math.pi, -0.000*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere002= mon.NormalSphereShape([0,0,0],radius)
	sphere002= mon.AffinizedNormalShape(sphere002)
	sphere002.scale(stretch[0],stretch[1],stretch[2])
	sphere002.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere002.rotate(pivot,0,rotation[1],rotation[2])
	sphere002.translate([location[0],location[1],location[2]])
	sphere002region = monte.addSubRegion(chamber,AuMSM,sphere002)

# Sphere.001
	location = [  0.000*meterspernm, 25.000*meterspernm, -6.000*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.500*math.pi,  0.500*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere001= mon.NormalSphereShape([0,0,0],radius)
	sphere001= mon.AffinizedNormalShape(sphere001)
	sphere001.scale(stretch[0],stretch[1],stretch[2])
	sphere001.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere001.rotate(pivot,0,rotation[1],rotation[2])
	sphere001.translate([location[0],location[1],location[2]])
	sphere001region = monte.addSubRegion(chamber,AuMSM,sphere001)

# Sphere
	location = [  0.000*meterspernm,  0.000*meterspernm, -7.649*meterspernm]
	stretch = [  1.000,  1.500,  0.500]
	rotbase = [  0.500*math.pi,  0.250*math.pi,  0.000*math.pi]
	dimensions = [ 12.000*meterspernm, 12.000*meterspernm, 12.000*meterspernm]
	rotation=[-rotbase[0],-rotbase[1],rotbase[2]]
	pivot=[0,0,0]
	radius=6.000*meterspernm
	sphere= mon.NormalSphereShape([0,0,0],radius)
	sphere= mon.AffinizedNormalShape(sphere)
	sphere.scale(stretch[0],stretch[1],stretch[2])
	sphere.rotate(pivot,math.pi/2,rotation[0],-math.pi/2)
	sphere.rotate(pivot,0,rotation[1],rotation[2])
	sphere.translate([location[0],location[1],location[2]])
	sphereregion = monte.addSubRegion(chamber,AuMSM,sphere)

# Layer.s1
	location = [  0.000*meterspernm,  0.000*meterspernm,250.000*meterspernm]
	rotbase = [  0.000*math.pi, -0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm,500.000*meterspernm]
	thickness1=location[2]+dimensions[2]/2
	layers1=mon.NormalMultiPlaneShape()
	layers1.addPlane(normalVector,[0,0,0])
	layers1region = monte.addSubRegion(chamber,SiMSM,layers1)

# Layer.s2
	location = [  0.000*meterspernm,  0.000*meterspernm,525.000*meterspernm]
	rotbase = [  0.000*math.pi, -0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm, 50.000*meterspernm]
	thickness2=location[2]+dimensions[2]/2
	layers2=mon.NormalMultiPlaneShape()
	layers2.addPlane(normalVector,[0,0,thickness1])
	layers2region = monte.addSubRegion(layers1region,SiMSMDeep,layers2)

# Layer.s3
	location = [  0.000*meterspernm,  0.000*meterspernm,557.500*meterspernm]
	rotbase = [  0.000*math.pi, -0.000*math.pi,  0.000*math.pi]
	dimensions = [500.000*meterspernm,500.000*meterspernm, 15.000*meterspernm]
	thickness3=location[2]+dimensions[2]/2
	layers3=mon.NormalMultiPlaneShape()
	layers3.addPlane(normalVector,[0,0,thickness2])
	layers3region = monte.addSubRegion(layers2region,SiMSMDeep,layers3)

	detectorList=[]
	return(monte,eg,detectorList)

detectorInfoList=[]

detectorInfoList=[{'Name': 'chamber', 'XfminE': 0.0, 'XfmaxE': 50.0, 'rotAng': 'x0_y0_z0'}]

#save scan region and beam parameters
scanregionlist=[]
scanregiondict={
	"Name":"Scan_Region.001", 	"trajImg":0,	"trajLog":0,	"histogram":0,	"gunZnm":50,	"nTrajectories":100,	"mExp":1,
	"centerx":-12.00,  	"centery":  0.00, 	"deltax":  2.00,  	"deltay":  2.00,
	"rotAng":[  0.00*radperdeg,  0.00*radperdeg,  0.00*radperdeg],	"pivotloc":[  0.00*meterspernm,  0.00*meterspernm,  0.00*meterspernm],
# create x and y vectors for simulation
	"xvals": [-58.00,-56.00,-54.00,-52.00,-50.00,-48.00,-46.00,-44.00,-42.00,-40.00,-38.00,-36.00,\
	-34.00,-32.00,-30.00,-28.00,-26.00,-24.00,-22.00,-20.00,-18.00,-16.00,-14.00,-12.00,],
	"yvals": [-38.00,-36.00,-34.00,-32.00,-30.00,-28.00,-26.00,-24.00,-22.00,-20.00,-18.00,-16.00,-14.00,\
	-12.00,-10.00,-8.00,-6.00,-4.00,-2.00,0.00,2.00,4.00,6.00,8.00,10.00,12.00,14.00,\
	16.00,18.00,20.00,22.00,24.00,26.00,28.00,30.00,32.00,34.00,36.00,38.00,],
 # Beam energies in eV
	"beamEeVvals" : [  500.0,],
 # beam size in nm
	"beamsizenmvals" : [  0.5,],
}
scanregionlist.append(scanregiondict)
loopdict={"thicknessvals": [0]}
GUI_Addon_Version=(1, 7)
file_dict={"blend": r'D:\JMONSEL_annex\Testing_June4\Jan16_Addon1_7_Phase2\sphere4.blend',"blend_epoch": 1642315841.9520166}
commentList=['Project Name','userid','affiliation',\
'','','',\
'','','']
