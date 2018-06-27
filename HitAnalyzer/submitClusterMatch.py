#! /usr/bin/env python

import os, multiprocessing, math, sys
import ROOT as rt

useCondorBatch = False
import optparse


def submitJobs(cmd,queue,jobname,filestring):
    path = os.getcwd()
    joblist = []
    for i in range(1,17):
        #workdir= "/afs/cern.ch/user/m/msommerh/CMSSW_9_4_0_patch1/src/bTag_nHits/HitAnalyzer/tmp"+jobname+"/job_{}".format(i)
       	workdir = "tmp"+jobname+"/job_4TeV_{}".format(i)
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
	    string = cmd+filestring.format(i)+" "+"4TeV-Signal_PU"+str(i)+"\n"
	    fout.write(string)
	    print string
	       
	    fout.write("echo 'STOP---------------'\n")
	    fout.write("echo\n")
	    fout.write("echo\n")
	       
	    os.system("chmod 755 job_%i.sh"%i )
	    os.system("bsub -q "+queue+" -o logs job_%i.sh -J %s"%(i,jobname))
	    print "job nr " + str(i) + " file " + filestring.format(i) + " being submitted"
	    os.chdir("../..")


        

if __name__ == "__main__":
  
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-i', '--input', action='store', type='string', dest='origin', default='QCD_PU_samples.txt')

  (options, args) = parser.parse_args()
  
  origin        = options.origin
  
  jobname = "CM"
  
  cmd='python FinalClusterMatcher2.py '
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
  #filestring = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_wPU-v2/180528_101050/0000/flatTuple_{}.root'	
  filestring1 = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_2000/btagHits_wPU/180604_113337/0000/flatTuple_{}.root'	
  filestring2 = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_4000/btagHits_wPU/180604_071651/0000/flatTuple_{}.root'


  submitJobs(cmd,queue,jobname,filestring2)
 
  print
  print "your jobs:"
  os.system("bjobs")
  userName=os.environ['USER']
  
  print
  print 'Done submitting jobs!'
  print
  
