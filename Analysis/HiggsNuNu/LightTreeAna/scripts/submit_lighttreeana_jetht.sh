#!/bin/sh
DOCERN=0
DOSUBMIT=1

DATE=170811

## Try and take the JOBWRAPPER and JOBSUBMIT commands
## from the environment if set, otherwise use these defaults
: ${JOBWRAPPER:="./scripts/generate_job.sh"}
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

QUEUEDIR=long #short #medium long

JOBDIRPREFIX=jobs_run2ana_${DATE}_CC
JOBDIR=$JOBDIRPREFIX/
OUTPUTPREFIX=output_run2ana_${DATE}_CC
OUTPUTDIR=$OUTPUTPREFIX/

OUTPUTNAME="output.root"

mkdir -p $JOBDIR
mkdir -p $OUTPUTDIR

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

echo "JOB name = $JOB"
for syst in ""
do
  mkdir -p $JOBDIR$syst
  mkdir -p $OUTPUTDIR$syst
  for channels in qcdA qcdC
    do
    JOB=$channels
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/weights.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/weights.hists`
    #executable expect strings separated by "!"
    ## To produce all of the hist
    HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}.hists`
    SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}.hists`
    ## To produce all of the hist for datacard
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}_datacard.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}_datacard.hists`
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}_sig.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}_sig.hists`
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}_approval.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}_approval.hists`
    ## To test for one hist
    #HISTSTRING=";(calo-pf)/recoil;Events!;(calo-pf)/recoil;Events"
    #SHAPESTRING="TMath::Abs(calomet-met)/metnomuons(40,0,2)!TMath::Abs(calomet-met)/metnoelectrons(40,0,2)"
    #HISTSTRING=";p_{T}^{j2} (GeV);Events"
    #SHAPESTRING="jet2_pt(12,40.,250.)"
    #HISTSTRING=";#Delta#phi_{jj};Events"
    #SHAPESTRING="dijet_dphi(50,0.,3.1416)"
    #HISTSTRING=";E_{T,no-#mu}^{miss} (GeV);Events"
    #SHAPESTRING="metnomuons(25,200.,600.)"
    #HISTSTRING=";E_{T,no-#mu}^{miss} (GeV);Events!;Forward tag jet #eta;Events"
    #SHAPESTRING="metnomuons(25,200.,600.)!forward_tag_eta(25,-5.,5.)"
    #HISTSTRING=";#Delta#phi(E_{T,no-#mu}^{miss},j);Events"
    #SHAPESTRING="fourjetsmetnomu_mindphi(14,2.3,3.1416)"
    echo "Making histograms: " $SHAPESTRING
    OUTPUTNAME="$channels.root"
    MINDPHICUT="fourjetsmetnomu_mindphi\>=0.5"
    DATASET="MET"
    if [ "$channels" = "qcdA" ]; then
	CONFIG=scripts/DefaultRun2ConfigQCDAC.cfg
	MINDPHICUT="fourjetsmetnomu_mindphi\<0.5"
	DATASET="JetHT"
    fi
    if [ "$channels" = "qcdB" ]; then
	CONFIG=scripts/DefaultRun2ConfigQCDBD.cfg
	MINDPHICUT="fourjetsmetnomu_mindphi\<0.5"
	DATASET="MET"
    fi
    if [ "$channels" = "qcdC" ]; then
	CONFIG=scripts/DefaultRun2ConfigQCDAC.cfg
	MINDPHICUT="fourjetsmetnomu_mindphi\>=0.5"
	DATASET="JetHT"
    fi
    if [ "$channels" = "qcdD" ]; then
	CONFIG=scripts/DefaultRun2ConfigQCDBD.cfg
	MINDPHICUT="fourjetsmetnomu_mindphi\>=0.5"
	DATASET="MET"
    fi
    if [ "$syst" = "" ]
      then
      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --dataset=$DATASET --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP
    else
	#if [ "$channels" = "nunu" ]
	   # then
	    #for other systs only run nunu
    #executable expect strings separated by "!"
	#HISTSTRING="Forward tag jet #eta;Events"
	#SHAPESTRING="forward_tag_eta(25,-5.,5.)"
	#echo "Making histograms: " $SHAPESTRING

      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --dataset=$DATASET --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' --syst=$syst -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP
	#fi
    fi
    if [ "$DOSUBMIT" = "1" ]; then 
      $JOBSUBMIT $JOBDIR$syst/$JOB.sh
    else 
      echo "$JOBSUBMIT $JOBDIR$syst/$JOB.sh"
    fi
  done
done
