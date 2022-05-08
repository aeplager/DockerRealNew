import os
max_detectors=12
# only pass parameters that are in the main program.
def writeheader(fileid,starttime,GUI_Addon_Version,randseed,totals,startpixel,detectorInfoList,scanregioninfo,filedict,commentList,loopdict):
        
        jobfolder=os.path.dirname(filedict["blend"])
        materials_file="makematerials.py"
        definesample_file="definesample.py"
        docker_config= ["none","none","none"]
        userid= commentList[1]
        affiliation=commentList[2]
        numdetectors=len(detectorInfoList)

	print >> fileid, "Project_name, ",commentList[0]
	print >> fileid, "Project_date_time, ",starttime[0]
	print >> fileid, "Filepath_job_folder, ",jobfolder
	print >> fileid, "Filepath_materials, ",materials_file
	print >> fileid, "Filepath_definesample, ",definesample_file
	print >> fileid, "Filepath_BlenderProject, ",filedict["blend"]
	print >> fileid, "Blender_Save_date_time, ",filedict["blend_modify"]
	print >> fileid, "Filepath_Cont1, ","N/A"
	print >> fileid, "Filepath_Cont2, ","N/A"
	print >> fileid, "Filepath_Cont3, ","N/A"
	print >> fileid, "Filepath_Cont4, ","N/A"
	print >> fileid, "Filepath_Cont5, ","N/A"
	print >> fileid, "JMONSEL_GUI_Addon_version, ",str(GUI_Addon_Version).replace(',','.')
	print >> fileid, "BlenderVersionNum, ",2.92
        print >> fileid, "JMONSEL Sourcecode version number, ","Nov2017"
	print >> fileid, "Docker_config1, ",docker_config[0]
	print >> fileid, "Docker_config2, ",docker_config[1]
	print >> fileid, "Docker_config3, ",docker_config[2]
	print >> fileid, "UserID, ",userid
	print >> fileid, "UserAffiliation, ",affiliation
	print >> fileid, "scan_region_start, ",starttime[1]
        fileid.write("Rnd_seed %d\n" %randseed)
	print >> fileid, "TotalNumPixels, ",totals[0]
        print >> fileid, "TotalImages,",totals[1]
	print >> fileid, "StartPixels, ",startpixel # row 23
        print >> fileid, "Chamber ",",",str(detectorInfoList[0]).replace('{','').replace('}','')
        for detnum in range(1,len(detectorInfoList)):  # write for loop so column title is written for each detector
            print >> fileid, "Detector",detnum,", ",str(detectorInfoList[detnum]).replace('{','').replace('}','')
        for detnum in range(len(detectorInfoList),max_detectors+1):
            print >> fileid,"Detector",detnum,", None" 
 	print >> fileid, "JobComment1,",commentList[3]  # row 37
	print >> fileid, "JobComment2,",commentList[4]
	print >> fileid, "SampleDesc1,",commentList[5]
	print >> fileid, "SampleDesc2,",commentList[6]
	print >> fileid, "RunComment1,",commentList[7]
	print >> fileid, "RunComment2,",commentList[8]
	print >> fileid, "Scan Region info," #,str(scanregioninfo).replace('{','').replace('}','') # row 43
        for scan_region_key in scanregioninfo:
            print >> fileid,scan_region_key,",",scanregioninfo[scan_region_key]
        for loop_key in loopdict:
            print >> fileid,loop_key,",",str(loopdict[loop_key]).replace('[','').replace(']','')
	print >> fileid, "header_contintingency2,","N/A"
	print >> fileid, "header_contintingency3,","N/A"

def write_column_titles(fileid,looplabellist):
        fileid.write("PixelNum, ImageNum, TimeSeconds, TimePixel, Contingency, X_nm, Y_nm, FrameNum, N,")
        fileid.write("BeamV, BeamSizeX, BeamSizeY, BeamRotAng, %s, %s, %s, " %(looplabellist[0],looplabellist[1],looplabellist[2]))
#  Add loop labels here
        fileid.write("energyHist, traj_image, traj_log, BSEf, SEf, Xf, Energy Depo,")
#        fileid.write("contingency2, contingency3, contingency4, contingency5, ")
        for detnum in range(1,max_detectors+1):  # write for loop so column title is written for each detector
            fileid.write("Det%d_BSEf, Det%d_SEf, Det%d_Xf, energyHist, traj_log, "%(detnum,detnum,detnum))
        fileid.write("\n")
        
def writeChamberData(fileid,pixelnum,imagenum,timestampsec,timeelapsed,xnm,ynm,frame,nTrajectories,beamEeV,beamsizex,beamsizey,rotbeam,energydep,histfile,trajImgFile,trajlogfile,bsf,SEf,Xfchamber,paramList):
        fileid.write("%d, %d, %8.6f, %8.6f, %s, %7.3f, %7.3f, %d, %d, "%(pixelnum,imagenum,timestampsec,timeelapsed,"NA",xnm,ynm,frame,nTrajectories)) 
	fileid.write("%7.2f, %7.2f, %7.2f, %s,"%(beamEeV,beamsizex,beamsizey,rotbeam))
	fileid.write("%7.2f, %7.2f, %7.2f,"%(paramList[0],paramList[1],paramList[2]))
        fileid.write("%s, %s, %s, %8.6g, %8.6g, %8.6g, %s, "%(histfile,trajImgFile,trajlogfile,bsf,SEf,Xfchamber,energydep))
#        fileid.write("contingency2, contingency3, contingency4, contingency5, ")

def writeDetectorData(fileid,detBsf,detSEf,detXf,histfile,trajfile):
	fileid.write("%8.4g, %8.4g, %8.4g, %s, %s, "%(detBsf,detSEf,detXf,histfile,trajfile))

def writeNoDetectorData(fileid):
	fileid.write("%8.4g, %8.4g, %8.4g, %s, %s, "%(0,0,0,"NA","NA"))

def endline(fileid):
        print >> fileid