# KinFitTester

# Development has moved to [GitLab](https://gitlab.cern.ch/cms-bph/KinFitTester)

## Install and run

Tested in CMSSW_10_6_10

In `CMSSW_10_6_10/src`:

* `git clone git@github.com:tpmccauley/KinFitTester.git KinFitTester/KinFitTester`
* `scram b`
* You may have to init your GRID proxy `voms-proxy-init --rfc --voms cms`
* `cmsRun KinFitTester/KinFitTester/test/KinFitTest_cfg.py`

## to-do

* Sensible fits
* Documentation
* Comparisons

## Links
*  [SWGuideKinematicFit](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit)
*  [SWGuideTransientTracks](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTransientTracks)
*  [KineExample](https://gitlab.cern.ch/cms-sw/cmssw/-/tree/009e1bdc5f6bdf33bbc367f3b6fd60df26d204bc/RecoVertex/KinematicFit/plugins)

