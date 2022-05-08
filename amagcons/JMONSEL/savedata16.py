# only pass parameters that are in the main program.
#def writeheader(fileid,name,datetime,jobfolder,materials_file,definesample_file,blend_file,GUI_Addon_version,docker_config,userid,affiliation,starttime,randseed,totalpixels,totaltrajs):
#import sqldata 
#def writeheader(fileid,GUI_Addon_Version,starttime,randseed,totalpixels,totaltrajs,detectorInfoList,looplabellist):
def writeheader(fileid,GUI_Addon_Version,starttime,randseed,totalpixels,totaltrajs,detectorInfoList,looplabellist):
        
        name="name"
        datetime="datetime"
        jobfolder="jobfolder"
        materials_file="makematerials.py"
        definesample_file="definesample.py"
        blend_file = "example.blend"
        docker_config= ["none","none","none"]
        userid= "usesrid"
        affiliation="affiliation"
        commentList=["job1","job2","desc1","desc2","run1","run2"]
        numdetectors=len(detectorInfoList)

	# print >> fileid, "Project_name, ",name
	# print >> fileid, "Project_date_time, ",datetime
	# print >> fileid, "Filepath_job_folder, ",jobfolder
	# print >> fileid, "Filepath_materials, ",materials_file
	# print >> fileid, "Filepath_definesample, ",definesample_file
	# print >> fileid, "Filepath_BlenderProject, ",blend_file
	# print >> fileid, "Filepath_Cont1, ","N/A"
	# print >> fileid, "Filepath_Cont2, ","N/A"
	# print >> fileid, "Filepath_Cont3, ","N/A"
	# print >> fileid, "Filepath_Cont4, ","N/A"
	# print >> fileid, "Filepath_Cont5, ","N/A"
	# print >> fileid, "JMONSEL_GUI_Addon_version, ",GUI_Addon_Version
	# print >> fileid, "BlenderVersionNum, ",2.92
    #     print >> fileid, "JMONSEL Sourcecode version number, ","Nov2017"
	# print >> fileid, "Docker_config1, ",docker_config[0]
	# print >> fileid, "Docker_config2, ",docker_config[1]
	# print >> fileid, "Docker_config3, ",docker_config[2]
	# print >> fileid, "UserID, ",userid
	# print >> fileid, "UserAffiliation, ",affiliation
	# print >> fileid, "Date_time_start, ",starttime
    #     fileid.write("Rnd_seed %d\n" %randseed)
	# print >> fileid, "TotalNumPixels, ",totalpixels
    #     for detnum in range(0,len(detectorInfoList)):  # write for loop so column title is written for each detector
    #         print >> fileid, "Detector",detnum,", ",str(detectorInfoList[detnum]).replace('{','').replace('}','')
            
 	# print >> fileid, "JobComment1",commentList[0]
	# print >> fileid, "JobComment2",commentList[1]
	# print >> fileid, "SampleDesc1",commentList[2]
	# print >> fileid, "SampleDesc2",commentList[3]
	# print >> fileid, "RunComment1",commentList[4]
	# print >> fileid, "RunComment2",commentList[5]
	# print >> fileid, "header_contintingency1","N/A"
	# print >> fileid, "header_contintingency2","N/A"
	# print >> fileid, "header_contintingency3","N/A"
 
 # Beginning Header Here
        fileid.write("PixelNum, TimeSeconds, TimeHours, TimePixel, ScanRegionName, X_nm, Y_nm, FrameNum, N,")
        fileid.write("BeamV, BeamSizeX, BeamSizeY, BeamRotAng, %s, %s, %s, " %(looplabellist[0],looplabellist[1],looplabellist[2]))
#  Add loop labels here
        fileid.write("energyHist, traj_image, traj_log, BSEf, SEf, Xf, rotChamber,")
#        fileid.write("contingency2, contingency3, contingency4, contingency5, ")
        for detnum in range(1,numdetectors+1):  # write for loop so column title is written for each detector
            fileid.write("Det%d_BSEf, Det%d_SEf, Det%d_Xf, energyHist, traj_log, "%(detnum,detnum,detnum))
        fileid.write("\n")
        
def writeChamberData(fileid,pixelnum,timestampsec,timestamphr,timeelapsed,scanregname,xnm,ynm,frame,nTrajectories,beamEeV,beamsizex,beamsizey,rotbeam,rotchamber,histfile,trajImgFile,trajlogfile,bsf,SEf,Xfchamber,paramList):
	PathSep='\\'
        fileid.write("%d, %7.2f, %8.4f, %7.2f, %s, %7.2f, %7.2f, %d, %d, "%(pixelnum,timestampsec,timestamphr,timeelapsed,scanregname,xnm,ynm,frame,nTrajectories)) 
	fileid.write("%7.2f, %7.2f, %7.2f, %s,"%(beamEeV,beamsizex,beamsizey,rotbeam))
	fileid.write("%7.2f, %7.2f, %7.2f,"%(paramList[0],paramList[1],paramList[2]))
		# AEP Commented it Out
        #fileid.write("%s, %s, %s, %8.4g, %8.4g, %8.4g, %s, "%(histfile.rsplit(PathSep,1)[1],trajImgFile,trajlogfile.rsplit(PathSep,1)[1],bsf,SEf,Xfchamber,rotchamber))
        fileid.write("%s, %s, %s, %8.4g, %8.4g, %8.4g, %s, "%(histfile,trajImgFile,trajlogfile,bsf,SEf,Xfchamber,rotchamber))
#        fileid.write("contingency2, contingency3, contingency4, contingency5, ")

def writeDetectorData(fileid,detBsf,detSEf,detXf,histfile,trajfile):
	fileid.write("%8.4g, %8.4g, %8.4g, %s, %s, "%(detBsf,detSEf,detXf,histfile,trajfile))

def endline(fileid):
        print >> fileid

def writeStatusFile(fileid, first_field, second_field):    
	fileid.write("%s, %s, "%(first_field,second_field))
	fileid.write("\n")