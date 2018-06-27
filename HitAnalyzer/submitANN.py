#! /usr/bin/env python

import os, multiprocessing, math, sys
import ROOT as rt

useCondorBatch = False
import optparse

job_dict = {0:'noPU_both', 1:'noPU_both_withPT', 2:'noPU_4TeV', 3:'noPU_4TeV_withPT', 4:'withPU_both', 5:'withPU_both_withPT', 6:'withPU_both_withPV', 7:'withPU_4TeV', 8:'withPU_4TeV_withPT', 9:'withPU_4TeV_withPV'}
pathstring = "Submitted_Models/data/"

def submitJobs(cmd,queue,jobname,nEpochs):
    path = os.getcwd()
    joblist = []
    for i in range(0,10):

	if i==0 or i==2 or i==4 or i==7:
		addFeature = "No"
		continue
	elif i==1 or i==3 or i==5 or i==8:
		addFeature = "pT"
	elif i==6 or i== 9:
		addFeature = "PV"

       	workdir = "tmp"+jobname+"/job_{}".format(job_dict[i])
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
	    os.system("bsub -q "+queue+" -o logs job_{}.sh -J {}".format(job_dict[i],jobname))
	    print "job nr {}: {} being submitted".format(i,job_dict[i])
	    os.chdir("../..")


        

if __name__ == "__main__":
  
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-i', '--input', action='store', type='string', dest='origin', default='QCD_PU_samples.txt')

  (options, args) = parser.parse_args()
  
  origin        = options.origin
  
  jobname = "ANN"
  nEpochs = 200
  
  cmd='python ANN4Grid.py '
  queue="1nd"
  
  submitJobs(cmd,queue,jobname,nEpochs)
 
  print
  print "your jobs:"
  os.system("bjobs")
  userName=os.environ['USER']
  
  print
  print 'Done submitting jobs!'
  print
  
