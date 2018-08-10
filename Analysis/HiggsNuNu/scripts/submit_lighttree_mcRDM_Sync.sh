#!/bin/sh
DOCERN=0
DOSUBMIT=0
JETTYPE="pfJetsPFlow"
MYEXEC=LightTreeMakerFromMiniAODRDM
PRODUCTION=170201
PRODUSER=rdimaria
JPTCUTVAL=0

DATE=170704

## Try and take the JOBWRAPPER and JOBSUBMIT commands
## from the environment if set, otherwise use these defaults
: ${JOBWRAPPER:="./scripts/generate_job.sh $DOCERN $MYEXEC $PRODUCTION"}
: ${JOBSUBMIT:="eval"}

GRIDSETUP=1

if [ "$DOCERN" = "0" ]
  then
  JOBSCRIPT="./scripts/submit_ic_batch_job.sh"
else
  JOBSCRIPT="./scripts/submit_cern_batch_job.sh"
  GRIDSETUP=0
fi
export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE


echo "Using job-wrapper: " $JOBWRAPPER
echo "Using job-submission: " $JOBSUBMIT

INPUTPARAMS="filelists/$PRODUCTION/Params${PRODUCTION}_Sync.dat"
CONFIG=scripts/DefaultLightTreeConfig_mc_Sync.cfg



for SYST in central #JESUP JESDOWN JERBETTER JERWORSE UESUP UESDOWN #NOTE TO RUN JER DOSMEAR MUST BE SET TO TRUE IN THE CONFIG
  do
  SYSTOPTIONS="--dojessyst=false --dojersyst=false"

  JOBDIRPREFIX=./jobs_lighttree_${DATE}_Sync
  JOBDIR=$JOBDIRPREFIX/
  OUTPUTPREFIX=./output_lighttree_${DATE}_Sync

  OUTPUTDIR=$OUTPUTPREFIX/

  if [ "$SYST" != "central" ]
    then
    JOBDIR=$JOBDIRPREFIX/$SYST/
    OUTPUTDIR=$OUTPUTPREFIX/$SYST/
  fi


  echo "Config file: $CONFIG"
  mkdir -p $JOBDIR
  mkdir -p $OUTPUTDIR

  cp $CONFIG $OUTPUTDIR

  for QUEUEDIR in Sync #short medium
    do
    if [ "$DOCERN" = "0" ]
      then
      if [ "$QUEUEDIR" = "medium" ]
        then
        JOBQUEUE="5:59:0"
      elif [ "$QUEUEDIR" = "long" ]
        then
        JOBQUEUE="47:59:0"
      else
        JOBQUEUE="2:59:0"
      fi
    else
      if [ "$QUEUEDIR" = "medium" ]
        then
        JOBQUEUE="1nd"
      elif [ "$QUEUEDIR" = "long" ]
        then
        JOBQUEUE="2nd"
      else
        JOBQUEUE="1nh"
      fi
    fi
    export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE
    echo "Using job-submission: " $JOBSUBMIT

    #PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/${PRODUSER}/${PRODUCTION}_MC
    #PREFIX=root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/${PRODUSER}/${PRODUCTION}_MC
    PREFIX=/home/hep/rd1715/CMSSW_8_0_28/src/UserCode/ICHiggsTauTau

    #for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_Powheg-VBF*125.dat`
    #for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_WJetsToLNu*.dat`
    #for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_*`
    for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_*`
      do
      echo "Processing files in "$FILELIST

      echo $FILELIST
      echo $FILELIST > tmp.txt


      MCOPTION="--mc=summer16_80X"
      #PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/${PRODUSER}/${PRODUCTION}_MC
      #PREFIX=root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/${PRODUSER}/${PRODUCTION}_MC
      PREFIX=/home/hep/rd1715/CMSSW_8_0_28/src/UserCode/ICHiggsTauTau
      sed "s/filelists\/${PRODUCTION}\/$QUEUEDIR\/${PRODUCTION}_MC_//" tmp.txt > tmp2.txt
      INPUTPARAMS="filelists/$PRODUCTION/Params${PRODUCTION}_Sync.dat"

      echo $PREFIX

      JOB=MC_`sed "s/\.dat//" tmp2.txt`

      echo "JOB name = $JOB"

      JPTCUT=$JPTCUTVAL
      grep "Htoinv" tmp.txt
      if (( "$?" == 0 )); then
        JPTCUT=0
        MCOPTION="--mc=summer16_80X --donoskim=true"
      fi

      NEEDSSTREAM=0
      grep "EWKW" tmp.txt
      if (( "$?" == 0 )); then
        NEEDSSTREAM=1
      fi
      grep  "JetsToLNu" tmp.txt
      if (( "$?" == 0 )); then
        NEEDSSTREAM=1
      fi
      if (( "$NEEDSSTREAM" == 1 )); then
        for FLAVOUR in enu munu taunu
          do
          WJOB=$JOB"_"$FLAVOUR

          $JOBWRAPPER $JOBDIR $OUTPUTDIR "./bin/$MYEXEC --cfg=$CONFIG --prod="$PRODUCTION" $MCOPTION --filelist="$FILELIST" --input_prefix=$PREFIX --output_name=$JOB.root --output_folder=$OUTPUTDIR $SYSTOPTIONS --inputparams=$INPUTPARAMS --wstream=$FLAVOUR --jet1ptcut="$JPTCUT" --jet2ptcut="$JPTCUT" --jettype=$JETTYPE &> $JOBDIR/$WJOB.log" $JOBDIR/$WJOB.sh $GRIDSETUP
          if [ "$DOSUBMIT" = "1" ]; then
            $JOBSUBMIT $JOBDIR/$WJOB.sh
          else
            echo "$JOBSUBMIT $JOBDIR/$WJOB.sh"
          fi
        done

      else
        $JOBWRAPPER $JOBDIR $OUTPUTDIR "./bin/$MYEXEC --cfg=$CONFIG --prod="$PRODUCTION" $MCOPTION --filelist="$FILELIST" --input_prefix=$PREFIX --output_name=$JOB.root --output_folder=$OUTPUTDIR $SYSTOPTIONS --inputparams=$INPUTPARAMS --jet1ptcut="$JPTCUT" --jet2ptcut="$JPTCUT" --jettype=$JETTYPE &> $JOBDIR/$JOB.log" $JOBDIR/$JOB.sh  $GRIDSETUP
        if [ "$DOSUBMIT" = "1" ]; then
          $JOBSUBMIT $JOBDIR/$JOB.sh
        else
          echo "$JOBSUBMIT $JOBDIR/$JOB.sh"
        fi
      fi
      rm tmp.txt tmp2.txt

    done

  done

done
