# Short recipe for CMGTools 

For the general recipe, [follow these instructions](https://twiki.cern.ch/twiki/bin/view/CMS/CMGToolsReleasesExperimental).

--------------

#### Set up CMSSW and the base git

```
cmsrel CMSSW_9_4_6_patch1
cd CMSSW_9_4_6_patch1/src 
cmsenv
git cms-init
```

# add the central cmg-cmssw repository to get the Heppy 94X_dev branch

```
git remote add cmg-central https://github.com/CERN-PH-CMG/cmg-cmssw.git  -f  -t heppy_94X_dev
```

#### Configure the sparse checkout, and get the base heppy packages

```
cp /afs/cern.ch/user/c/cmgtools/public/sparse-checkout_94X_heppy .git/info/sparse-checkout
git checkout -b heppy_94X_dev cmg-central/heppy_94X_dev
```

#### Add your mirror, and push the 80X branch to it

```
YOUR_GITHUB_REPOSITORY=$(git config user.github) # or set it manually if this doesn't work for you
git remote add origin git@github.com:$YOUR_GITHUB_REPOSITORY/cmg-cmssw.git
git push -u origin heppy_94X_dev
```

#### Now get the CMGTools subsystem from the cmgtools-lite repository

```
git clone -o mdunser https://github.com/mdunser/cmgtools-lite.git -b hinTTbar_marc_dev CMGTools
cd CMGTools 
```

#### Add your fork, and push the 80X branch to it

```
git remote add origin  git@github.com:$YOUR_GITHUB_REPOSITORY/cmgtools-lite.git 
git push -u origin hinTTbar_marc_dev
```

#### Compile

```
cd $CMSSW_BASE/src
scram b -j 8
```
