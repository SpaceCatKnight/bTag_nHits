#! /usr/bin/env python

import os, multiprocessing, math, sys
import ROOT as rt

useCondorBatch = False
import optparse

'''
#pT and PV on complex functional model on both samples and 4TeV with and without PU
job_dict = {0:'noPU_both', 1:'noPU_both_withPT', 2:'noPU_4TeV', 3:'noPU_4TeV_withPT', 4:'withPU_both', 5:'withPU_both_withPT', 6:'withPU_both_withPV', 7:'withPU_4TeV', 8:'withPU_4TeV_withPT', 9:'withPU_4TeV_withPV'}
Epoch_dict = {0:160, 1:180, 2:120, 3:135, 4:160, 5:180, 6:180, 7:130, 8:170, 9:100}
pathstring = "Submitted_Models/data/"

def submitJobs(cmd,queue,global_jobname):
    path = os.getcwd()
    joblist = []
    for i in range(0,10):
	local_jobname = job_dict[i]
	nEpochs = Epoch_dict[i]

	if i==0 or i==2 or i==4 or i==7:
		addFeature = "No"
		#continue
	elif i==1 or i==3 or i==5 or i==8:
		addFeature = "pT"
		#continue
	elif i==6 or i== 9:
		addFeature = "PV"
		#continue

       	workdir = "tmp"+global_jobname+"/job_{}".format(job_dict[i])
	#os.system("mkdir {}".format(workdir))
        os.makedirs(workdir)
	os.chdir(workdir)
	   
	with open('job_{}.sh'.format(job_dict[i]), 'w') as fout:
	    fout.write("#!/bin/sh\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	    fout.write("echo 'START---------------'\n")
	    fout.write("echo 'WORKDIR ' ${PWD}\n")
	    fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
	    fout.write("cd "+str(path)+"\n")
	    fout.write("cmsenv\n")
	    string = cmd+job_dict[i]+" "+pathstring+job_dict[i]+" "+addFeature+" "+str(nEpochs)+" "+str(i)+"\n"
	    fout.write(string)
	    print string
	       
	    fout.write("echo 'STOP---------------'\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	       
	    os.system("chmod 755 job_{}.sh".format(job_dict[i]))
	    os.system("bsub -q "+queue+" -o logs job_{}.sh -J {}".format(job_dict[i],local_jobname))
	    print "job nr {}: {} being submitted".format(i,job_dict[i])
	    os.chdir("../..")
'''
'''
#simple approaches on both samples without PU
job_dict = {0:'Li004', 1:'Li01', 2:'twoCones', 3:'fiveCones', 4:'fiveConesConv', 5:'fiveConesFunc'}
model_dict = {0:"simple", 1:"simple", 2:"simple", 3:"simple", 4:"conv", 5:"functional"}
Epoch_dict = {0:260, 1:145, 2:260, 3:260, 4:90, 5:250}

pathstring = "Submitted_Models/data/noPU_both/"

def submitJobs(cmd,queue,global_jobname):
    path = os.getcwd()
    joblist = []
    for i in range(0,6):
	if i == 1 or i == 4: continue
	local_jobname = job_dict[i]
	nEpochs = Epoch_dict[i]
	model = model_dict[i]
	addFeature = "No"

       	workdir = "tmp"+global_jobname+"/job_{}".format(job_dict[i])
	#os.system("mkdir {}".format(workdir))
        os.makedirs(workdir)
	os.chdir(workdir)
	   
	with open('job_{}.sh'.format(job_dict[i]), 'w') as fout:
	    fout.write("#!/bin/sh\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	    fout.write("echo 'START---------------'\n")
	    fout.write("echo 'WORKDIR ' ${PWD}\n")
	    fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
	    fout.write("cd "+str(path)+"\n")
	    fout.write("cmsenv\n")
	    string = cmd+job_dict[i]+" "+pathstring+" "+addFeature+" "+str(nEpochs)+" "+str(i)+" "+model+"\n"
	    fout.write(string)
	    print string
	       
	    fout.write("echo 'STOP---------------'\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	       
	    os.system("chmod 755 job_{}.sh".format(job_dict[i]))
	    os.system("bsub -q "+queue+" -o logs job_{}.sh -J {}".format(job_dict[i],local_jobname))
	    print "job nr {}: {} being submitted".format(i,job_dict[i])
	    os.chdir("../..")
'''
#pT and PV with functional model on both samples with and without PU
job_dict = {0:'noPU_functional', 1:'noPU_functional_withPT', 2:'withPU_functional', 3:'withPU_functional_withPT', 4:'withPU_functional_withPV'}
Epoch_dict = {0:250, 1:300, 2:300, 3:300, 4:300}
pathstring = "Submitted_Models/data/"

def submitJobs(cmd,queue,global_jobname):
    path = os.getcwd()
    joblist = []
    for i in range(0,5):
	if i == 0 or i==2: continue
	
	local_jobname = job_dict[i]
	nEpochs = Epoch_dict[i]
	

	if i==0 or i==2:
		addFeature = "No"
		if i==0:
			local_path = 'noPU_both'
		if i==2:
			local_path = 'withPU_both'
	elif i==1 or i==3:
		addFeature = "pT"
		if i==1:
			local_path = 'noPU_both_withPT'
		if i==3:
			local_path = 'withPU_both_withPT'

	elif i==4:
		addFeature = "PV"
		local_path = 'withPU_both_withPV'

       	workdir = "tmp"+global_jobname+"/job_{}".format(job_dict[i])
	#os.system("mkdir {}".format(workdir))
        os.makedirs(workdir)
	os.chdir(workdir)
	   
	with open('job_{}.sh'.format(job_dict[i]), 'w') as fout:
	    fout.write("#!/bin/sh\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	    fout.write("echo 'START---------------'\n")
	    fout.write("echo 'WORKDIR ' ${PWD}\n")
	    fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
	    fout.write("cd "+str(path)+"\n")
	    fout.write("cmsenv\n")
	    string = cmd+job_dict[i]+" "+pathstring+local_path+" "+addFeature+" "+str(nEpochs)+" "+str(i)+" "+"functional"+"\n"
	    fout.write(string)
	    print string
	       
	    fout.write("echo 'STOP---------------'\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	       
	    os.system("chmod 755 job_{}.sh".format(job_dict[i]))
	    os.system("bsub -q "+queue+" -o logs job_{}.sh -J {}".format(job_dict[i],local_jobname))
	    print "job nr {}: {} being submitted".format(i,job_dict[i])
	    os.chdir("../..")

        

if __name__ == "__main__":
  
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-i', '--input', action='store', type='string', dest='origin', default='')

  (options, args) = parser.parse_args()
  
  origin        = options.origin
  
  global_jobname = "ANN"
  
  cmd='python ANN4Grid.py '
  queue="1nd"
  
  submitJobs(cmd,queue,global_jobname)
 
  print
  print "your jobs:"
  os.system("bjobs")
  userName=os.environ['USER']
  
  print
  print 'Done submitting jobs!'
  print
  
