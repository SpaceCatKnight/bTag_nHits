#! /usr/bin/env python

import os, multiprocessing, math, sys
import ROOT as rt

useCondorBatch = False
import optparse


def submitJobs(cmd,queue,jobname):
    path = os.getcwd()
    joblist = []
    for i in range(1,48):
        #workdir= "/afs/cern.ch/user/m/msommerh/CMSSW_9_4_0_patch1/src/bTag_nHits/HitAnalyzer/tmp"+jobname+"/job_{}".format(i)
       	workdir = "tmp"+jobname+"/job_{}".format(i)
	#os.system("mkdir {}".format(workdir))
        os.makedirs(workdir)
	os.chdir(workdir)
	   
	with open('job_{}.sh'.format(i), 'w') as fout:
	    fout.write("#!/bin/sh\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	    fout.write("echo 'START---------------'\n")
	    fout.write("echo 'WORKDIR ' ${PWD}\n")
	    fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
	    fout.write("cd "+str(path)+"\n")
	    fout.write("cmsenv\n")
	    string = cmd+str(i)+"\n"
	    fout.write(string)
	    print string
	       
	    fout.write("echo 'STOP---------------'\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	       
	    os.system("chmod 755 job_%i.sh"%i )
	    os.system("bsub -q "+queue+" -o logs job_%i.sh -J %s"%(i,jobname))
	    print "job nr " + str(i) + " being submitted"
	    os.chdir("../..")


        

if __name__ == "__main__":
  
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-i', '--input', action='store', type='string', dest='origin', default='')

  (options, args) = parser.parse_args()
  
  origin        = options.origin
  
  jobname = "dR_dist"
  
  cmd='python dR_dist.py '
  queue="1nd"
  """
  ##### Creating and sending jobs #####
  filelist = []
  with open(origin) as f: 
      for line in f: 
          print "Adding file to list: " ,line
          fileslist.append(f)
  ###### loop for creating and sending jobs #####
  """

  submitJobs(cmd,queue,jobname)
 
  print
  print "your jobs:"
  os.system("bjobs")
  userName=os.environ['USER']
  
  print
  print 'Done submitting jobs!'
  print
 
