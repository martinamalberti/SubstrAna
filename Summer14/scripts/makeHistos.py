import os
import sys
import numpy
import glob
import ROOT
from ROOT import TFile, TTree

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputdir"   , dest="inputdir"   , type="string", help="Path where ntuples are located. Example: /store/group/phys_jetmet/ntran/PUPPI/miniSamples/62x/qcd300-470_62x_PU40BX50")
parser.add_option("-w","--workdir"    , dest="workdir"    , type="string", default="mydir",help="Name of the directory for jobs")
parser.add_option("-o","--outputname" , dest="outputname" , type="string", default="histograms",help="Name of the output file. Default is: histograms")
parser.add_option("-q","--queue"      , dest="queue"      , type="string", default="1nh",help="Name of the queue on lxbatch")
parser.add_option(""  ,"--checkJobs"  , dest="checkJobs"  , action="store_true", default=False,help="Checks job status")
parser.add_option(""  ,"--resubmit"   , dest="resubmit"   , action="store_true", default=False,help="Resubmit job ")
parser.add_option("-j","--job"        , dest="job"        , type="int", help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("-n","--njobs"      , dest="njobs"      , type="int"   , help="Number of jobs")
parser.add_option("","--ptMin"        , dest="ptmin"      , type="float", default=25., help="Min jet pT")
parser.add_option("","--ptMax"        , dest="ptmax"      , type="float", default=1000., help="Max jet pT")
parser.add_option("","--etaMin"       , dest="etamin"     , type="float", default=0., help="Min jet eta")
parser.add_option("","--etaMax"       , dest="etamax"     , type="float", default=5., help="Max jet eta")
parser.add_option("","--useMatchedWJets", dest="useMatchedWJets", type="int", default=0, help="Use matched W jets, for WW sample only")

parser.add_option(""  ,"--dryRun"     , dest="dryRun"     , default=False, action="store_true",help="Do not submit jobs")


(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'

def makeFilesList(indir,wdir):
    list = []
    command = ('%s find -f %s | grep root > %s/list.txt' % (eos,indir,wdir))
    #print command
    os.system(command)
    file = open('%s/list.txt'%wdir, 'r')
    list =[line.replace('/eos/cms/','root://eoscms.cern.ch//').replace('\n','') for line in file ]
    print 'Found %d files' %len(list)
    #print list             
    return list


def writeJobs(wdir, indir, outputname, ptmin, ptmax, etamin, etamax, matchedwjets, njobs):
    #---------------------------------------------
    # --- prepare the list of files to be analyzed
    #---------------------------------------------
    listoffiles = []
    listoffiles = makeFilesList(indir,wdir) ## make the file list for the input directory on eos
    print listoffiles
    print listoffiles[0]

    #---------------------------------------------
    # --- now split the jobs
    #---------------------------------------------
    #njobs = len(listoffiles)
    
    for jobid in range(njobs):         
        jobdir = '%s/JOB_%d'%(wdir,jobid)
        os.system("mkdir -p "+jobdir)        

        #--- prepare the list of files for each job
        f = open('%s/input_%d.txt'%(jobdir,jobid), 'w')
        sublist = [file for i,file in enumerate(listoffiles) if (i%njobs==jobid)]
        for fname in sublist:
            f.write('%s \n'%fname)

        #--- prepare the jobs scripts
        jobscript = open('%s/sub_%d.sh'%(jobdir,jobid),'w')
        jobscript.write('cd %s \n'%jobdir)
        #jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc472 \n')
        jobscript.write('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
        jobscript.write('eval ` scramv1 runtime -sh ` \n')
        jobscript.write('cd - \n')
        jobscript.write('if ( \n')
        jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,jobid))
        #jobscript.write('\t TestJetTreeAnalysis %s %s_%d.root %f %f %f %f 0'%(listoffiles[jobid], outputname, jobid, ptmin, ptmax, etamin, etamax))
        jobscript.write('\t TestJetTreeAnalysis %s/input_%d.txt %s_%d.root %f %f %f %f %d'%(jobdir, jobid, outputname, jobid, ptmin, ptmax, etamin, etamax, matchedwjets))
        jobscript.write(') then \n')
        jobscript.write('\t cp ./%s_%d.root %s \n'%(outputname,jobid,jobdir))
        jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,jobid))
        jobscript.write('else \n')
        jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,jobid))
        jobscript.write('fi \n')
        os.system('chmod a+x %s/sub_%d.sh'%(jobdir,jobid))
        

def submitJobs(wdir, queue):
    jobs = glob.glob( '%s/JOB_*/sub*.sh'% (wdir) )
    njobs = len(jobs)
    print "for job in range(njobs): ",njobs
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        jobname = '%s/sub_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/sub_%d.log %s'%(queue,jobdir,job,jobname )
        os.system('bsub -q %s -o %s/sub_%d.log %s'%(queue,jobdir,job,jobname ))
        


def checkJobs(wdir, output, queue):
    jobs = glob.glob( '%s/JOB_*/sub*.sh'% (wdir) )
    print 'Total number of jobs: %d' %len(jobs)
    
    listdone = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.done' % (wdir,j,j))]
    print 'Total number of DONE jobs: %s ' % len(listdone)
    print '  %s' %listdone
    for j in listdone:
        f = '%s/JOB_%d/sub_%d.run'%(wdir,j,j)
        #print 'rm %s/JOB_%d/sub_%d.run'%(wdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
            
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.run' % (wdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
    print '   %s' %listrun
            
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.fail' % (wdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed
                        
    if (len(listdone) == len(jobs)):
        print "All jobs successful! Merging output files..."
        os.chdir(wdir)
        #        os.system('ls JOB*/*.root')
        # add automathically only if saved locally
        os.system('hadd -f %s.root JOB_*/*.root'%output)
    else:
        for j in listfailed:
            print 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(queue,wdir,j,j,wdir,j,j )
                    


#-----------------------------------------
#--- MAIN
#-----------------------------------------

path = os.getcwd()
workingdir = path+'/'+options.workdir

if not options.checkJobs and not options.resubmit:

    # -- write jobs scripts
    if os.path.isdir(workingdir):
     os.system("rm -r "+workingdir)

    os.system("mkdir -p "+workingdir)
    writeJobs(workingdir, options.inputdir, options.outputname, options.ptmin, options.ptmax, options.etamin, options.etamax, options.useMatchedWJets,options.njobs )
    
    # -- submit jobs 
    if not options.dryRun:
        submitJobs(workingdir, options.queue)

elif options.resubmit and options.job >-1 :
    print 'Resubmitting job %d ' %options.job
    resubcmd = 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(options.queue,workingdir,options.job,options.job,workingdir,options.job,options.job )
    #print resubcmd
    os.system(resubcmd)

elif options.checkJobs:
    checkJobs(workingdir,options.outputname, options.queue)
    
