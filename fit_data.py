#Assumes data in form of excitation energy, counts
import ROOT
import sys
import numpy
import math
import json

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

################
#Load JSON file#
################
jsonFile = open(sys.argv[1],"r")
data = json.load(jsonFile)

###########################
#Check for JSON parameters#
###########################
headings=["data_source","qfc","lowerPeaks","midPeaks","ias","gtr","upperPeaks","backgrounds","normalization"]
subHeadings=[
  ["name","emin","emax","nbins"],
  ["eproj","ethresh","sp","eejectilefree","bcoulomb","exn","height"],
  ["nPeaks","means","sigma","height"],
  ["nPeaks","means","sigma","height"],
  ["mean","sigma","height"],
  ["mean","sigma","height"],
  ["nPeaks","means","sigma","height"],
  ["nPeaks","means","sigmas","heights"],
  ["A","Z","ratio"]
]
for i,heading in enumerate(headings):
  if not heading in data:
    print("No "+heading+" parameter found in json file, exiting")
    sys.exit()
  for subHeading in subHeadings[i]:
    if not subHeading in data[heading]:
      print("No "+subHeading+" parameter found in parameter "+heading+",exiting")
      sys.exit()
      
####################
##LOAD DATA POINTS##
####################
energies=[]
counts=[]
graph=ROOT.TGraph()
graph.SetPoint(0,0,0)
pointNum=1
#Load data sets
inpFile=open(data["data_source"]["name"],"r")
for line in inpFile:
  line=line.strip("\n")
  lineParts=line.split(",")
  E=float(lineParts[0])
  contents=float(lineParts[1])
  graph.SetPoint(pointNum,E,contents)
  pointNum+=1

###############
##Interpolate##
###############
eMin=float(data["data_source"]["emin"])
eMax=float(data["data_source"]["emax"])
nBins=int(data["data_source"]["nbins"])
hist=ROOT.TH1D("hist",";Energy (MeV);Counts",nBins,eMin,eMax)
for i in range(0,nBins):
  val=graph.Eval(hist.GetBinCenter(i+1))
  hist.SetBinContent(i+1,val)
  hist.SetBinError(i+1,0.001)

##########################
##Make RooFit Observable##
##########################
eVar=ROOT.RooRealVar("eVar","Excitation Energy (MeV)",eMin,eMax)
argSet=ROOT.RooArgSet(eVar)
argList=ROOT.RooArgList(eVar)

###################
##Make PDF to Fit##
###################
dataHist=ROOT.RooDataHist("dataHist","dataHist",argList,hist)
histPdf=ROOT.RooHistPdf("histPdf","histPdf",argSet,dataHist,1)

###################################
##Quasi-free continuum background##
###################################
#Read pars
E_projectile=float(data["qfc"]["eproj"])
E_thresh=float(data["qfc"]["ethresh"])
E_gs=E_projectile-E_thresh
S_p=float(data["qfc"]["sp"])
E_ejectile_free = float(data["qfc"]["eejectilefree"])
B_coulomb=float(data["qfc"]["bcoulomb"])
E_xn = float(data["qfc"]["exn"])

#Make PDF
E_shift=ROOT.RooRealVar("E_shift","E_shift",E_gs)
EQF = ROOT.RooRealVar("EQF","EQF",E_ejectile_free-(S_p+E_xn+B_coulomb))
E0 = ROOT.RooRealVar("E0","E0",E_gs-S_p)
T=ROOT.RooRealVar("T","T",100) #Recommeneded value
W=ROOT.RooRealVar("W","W",22) #Recommended value
argList=ROOT.RooArgList(eVar,E_shift,E0,EQF,T,W)
qfc=ROOT.RooGenericPdf("qfc","qfc","(1-exp(((E_shift - eVar) - E0)/T)) / (1+(((E_shift - eVar) - EQF)/W)^2)" ,argList)

qfc_amp=float(data["qfc"]["height"])
qfcAmp=ROOT.RooRealVar("qfcAmp","qfcAmp",qfc_amp,0.5*qfc_amp,2*qfc_amp)

###############
##Lower peaks##
###############
lowerPeak_sigma=float(data["lowerPeaks"]["sigma"])
lowerPeakSigma=ROOT.RooRealVar("lowerPeakSigma","lowerPeakSigma",lowerPeak_sigma,0.1*lowerPeak_sigma,2*lowerPeak_sigma)
lowerPeak_means=data["lowerPeaks"]["means"]
lowerPeakMeans=[]
lowerPeak_amp=float(data["lowerPeaks"]["height"])
lowerPeakAmplitudes=[]
lowerPeakGaussians=[]
nLowerPeaks=int(data["lowerPeaks"]["nPeaks"])
for i in range(0,nLowerPeaks):
  lowerPeakMeans.append(ROOT.RooRealVar("lowerPeakMean"+str(i),"lowerPeakMean"+str(i),lowerPeak_means[i],0.9*lowerPeak_means[i],1.1*lowerPeak_means[i]))
  lowerPeakAmplitudes.append(ROOT.RooRealVar("lowerPeakAmp"+str(i),"lowerPeakAmp"+str(i),lowerPeak_amp,0.1*lowerPeak_amp,2*lowerPeak_amp))
  lowerPeakGaussians.append(ROOT.RooGaussian("lowerPeakGauss"+str(i),"lowerPeakGauss"+str(i),eVar,lowerPeakMeans[i],lowerPeakSigma))
  
#############
##Mid peaks##
#############
midPeak_sigma=float(data["midPeaks"]["sigma"])
midPeakSigma=ROOT.RooRealVar("midPeakSigma","midPeakSigma",midPeak_sigma,0.1*midPeak_sigma,2*midPeak_sigma)
midPeak_means=data["midPeaks"]["means"]
midPeakMeans=[]
midPeak_amp=float(data["midPeaks"]["height"])
midPeakAmplitudes=[]
midPeakGaussians=[]
nMidPeaks=int(data["midPeaks"]["nPeaks"])
for i in range(0,nMidPeaks):
  midPeakMeans.append(ROOT.RooRealVar("midPeakMean"+str(i),"midPeakMean"+str(i),midPeak_means[i],0.9*midPeak_means[i],1.1*midPeak_means[i]))
  midPeakAmplitudes.append(ROOT.RooRealVar("midPeakAmp"+str(i),"midPeakAmp"+str(i),midPeak_amp,0.1*midPeak_amp,2*midPeak_amp))
  midPeakGaussians.append(ROOT.RooGaussian("midPeakGauss"+str(i),"midPeakGauss"+str(i),eVar,midPeakMeans[i],midPeakSigma))

#######
##IAS##
#######
ias_mean=float(data["ias"]["mean"])
ias_sigma=float(data["ias"]["sigma"])
ias_amp=float(data["ias"]["height"])
iasMean=ROOT.RooRealVar("iasMean","IAS",ias_mean,0.99*ias_mean,1.01*ias_mean)
iasAmp=ROOT.RooRealVar("iasAmp","iasAmp",ias_amp,0.1*ias_amp,2*ias_amp)
iasSigma=ROOT.RooRealVar("iasSigma","iasSigma",ias_sigma,0.1*ias_sigma,2*ias_sigma)
iasGauss=ROOT.RooGaussian("iasGauss","iasGauss",eVar,iasMean,iasSigma)

#######
##GTR##
#######
gtr_mean=float(data["gtr"]["mean"])
gtr_sigma=float(data["gtr"]["sigma"])
gtr_amp=float(data["gtr"]["height"])
gtrMean=ROOT.RooRealVar("gtrMean","gtr",gtr_mean,0.95*gtr_mean,1.05*gtr_mean)
gtrAmp=ROOT.RooRealVar("gtrAmp","gtrAmp",gtr_amp,0.1*gtr_amp,2*gtr_amp)
gtrSigma=ROOT.RooRealVar("gtrSigma","gtrSigma",gtr_sigma,0.1*gtr_sigma,2*gtr_sigma)
gtrGauss=ROOT.RooGaussian("gtrGauss","gtrGauss",eVar,gtrMean,gtrSigma)

##############
##upperPeaks##
##############
upperPeak_sigma=float(data["upperPeaks"]["sigma"])
upperPeakSigma=ROOT.RooRealVar("upperPeakSigma","upperPeakSigma",upperPeak_sigma,0.1*upperPeak_sigma,2*upperPeak_sigma)
upperPeak_means=data["upperPeaks"]["means"]
upperPeakMeans=[]
upperPeak_amp=float(data["upperPeaks"]["height"])
upperPeakAmplitudes=[]
upperPeakGaussians=[]
nUpperPeaks=int(data["upperPeaks"]["nPeaks"])
for i in range(0,nUpperPeaks):
  upperPeakMeans.append(ROOT.RooRealVar("upperPeakMean"+str(i),"upperPeakMean"+str(i),upperPeak_means[i],0.9*upperPeak_means[i],1.1*upperPeak_means[i]))
  upperPeakAmplitudes.append(ROOT.RooRealVar("upperPeakAmp"+str(i),"upperPeakAmp"+str(i),upperPeak_amp,0.1*upperPeak_amp,2*upperPeak_amp))
  upperPeakGaussians.append(ROOT.RooGaussian("upperPeakGauss"+str(i),"upperPeakGauss"+str(i),eVar,upperPeakMeans[i],upperPeakSigma))
  
###############
##backgrounds##
###############
nBackgroundPeaks=data["backgrounds"]["nPeaks"]
background_means=data["backgrounds"]["means"]
background_sigmas=data["backgrounds"]["sigmas"]
background_heights=data["backgrounds"]["heights"]
backgroundMeans=[]
backgroundSigmas=[]
backgroundAmplitudes=[]
backgroundGaussians=[]
for i in range(0,int(nBackgroundPeaks)):
  backgroundSigmas.append(ROOT.RooRealVar("backgroundSigma"+str(i),"backgroundSigma"+str(i),background_sigmas[i],0.1*background_sigmas[i],2*background_sigmas[i]))
  backgroundMeans.append(ROOT.RooRealVar("backgroundMean"+str(i),"backgroundMean"+str(i),background_means[i],0.9*background_means[i],1.1*background_means[i]))
  backgroundAmplitudes.append(ROOT.RooRealVar("backgroundAmp"+str(i),"backgroundAmp"+str(i),background_heights[i],0.1*background_heights[i],2*background_heights[i]))
  backgroundGaussians.append(ROOT.RooGaussian("backgroundGauss"+str(i),"backgroundGauss"+str(i),eVar,backgroundMeans[i],backgroundSigmas[i]))

##############
##Make model##
##############
pdfList=ROOT.RooArgList()
ampList=ROOT.RooArgList()
pdfList.add(qfc)
ampList.add(qfcAmp)
for i in range(0,nLowerPeaks):
  pdfList.add(lowerPeakGaussians[i])
  ampList.add(lowerPeakAmplitudes[i])
for i in range(0,nMidPeaks):
  pdfList.add(midPeakGaussians[i])
  ampList.add(midPeakAmplitudes[i])
pdfList.add(iasGauss)
ampList.add(iasAmp)
pdfList.add(gtrGauss)
ampList.add(gtrAmp)
for i in range(0,nUpperPeaks):
  pdfList.add(upperPeakGaussians[i])
  ampList.add(upperPeakAmplitudes[i])
for i in range(0,nBackgroundPeaks):
  pdfList.add(backgroundGaussians[i])
  ampList.add(backgroundAmplitudes[i])

model=ROOT.RooAddPdf("model","model",pdfList,ampList)

#######
##Fit##
#######
res = model.fitTo(dataHist,ROOT.RooFit.Save(1),ROOT.RooFit.Extended(1))

############
##Plotting##
############
c1=ROOT.TCanvas("c1","c1",1000,550)
pad1 = ROOT.TPad("pad1","pad1",0.0,0.0,0.5,1.0)
pad1.Draw()
pad1.cd()
frame=eVar.frame()
dataHist.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Binning(nBins),ROOT.RooFit.MarkerSize(0.5))
model.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlack))
model.plotOn(frame,ROOT.RooFit.Components("qfc"),ROOT.RooFit.LineColor(ROOT.kGray),ROOT.RooFit.Name("qfc"))
for i in range(0,nLowerPeaks):
  norm=lowerPeakAmplitudes[i].getVal()
  lowerPeakGaussians[i].plotOn(frame,ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.AddTo("qfc"),ROOT.RooFit.Normalization(norm,ROOT.RooAbsReal.NumEvent))
for i in range(0,nMidPeaks):
  model.plotOn(frame,ROOT.RooFit.Components("midPeakGauss"+str(i)),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.AddTo("qfc"))
model.plotOn(frame,ROOT.RooFit.Components("iasGauss"),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.AddTo("qfc"))
model.plotOn(frame,ROOT.RooFit.Components("gtrGauss"),ROOT.RooFit.LineColor(ROOT.kViolet),ROOT.RooFit.AddTo("qfc"))
for i in range(0,nUpperPeaks):
  model.plotOn(frame,ROOT.RooFit.Components("upperPeakGauss"+str(i)),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen),ROOT.RooFit.AddTo("qfc"))
for i in range(0,nBackgroundPeaks):
  model.plotOn(frame,ROOT.RooFit.Components("backgroundsGauss"+str(i)),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.AddTo("qfc"))
frame.Draw()


#############
##Normalize##
#############
#Use user-supplied ratio of cross sections at 0-degree or q-min?
ratio=float(data["normalization"]["ratio"])

#Calculate B(F)
A=data["normalization"]["A"]
Z=data["normalization"]["Z"]
N=A-Z
BF=N-Z
BF_amp = iasAmp.getVal()

#Calculate IAS wave number
kf_ias = math.sqrt(E_projectile-E_thresh-iasMean.getVal())

scaleFactor = BF/(ratio*BF_amp)

sum=0
scaledPdfList=[]
scaledPdfs=ROOT.RooArgList()
scaledAmplitudesList=[]
scaledAmplitudes=ROOT.RooArgList()
print("\n\nEnergy B(GT)")
for i in range(0,nLowerPeaks):

  #Wave number for this GT peak
  kf_gt = math.sqrt(E_projectile-E_thresh-lowerPeakMeans[i].getVal())
  ratio=kf_ias/kf_gt
  
  #Compute B(GT)
  BGT = lowerPeakAmplitudes[i].getVal()*scaleFactor*ratio
  
  #Scale
  scaledAmplitudesList.append(ROOT.RooRealVar("scaledLowerPeakAmp"+str(i),"scaledLowerPeakAmp"+str(i),BGT))
  scaledAmplitudes.add(scaledAmplitudesList[i])
  scaledPdfList.append(lowerPeakGaussians[i])
  scaledPdfs.add(scaledPdfList[i])
  
  print("{:.3f}".format(lowerPeakMeans[i].getVal())," {:.3f}".format(BGT))
  
  sum+=BGT
for i in range(0,nMidPeaks):
  #Wave number for this GT peak
  kf_gt = math.sqrt(E_projectile-E_thresh-midPeakMeans[i].getVal())
  ratio=kf_ias/kf_gt
  
  #Compute B(GT)
  BGT = midPeakAmplitudes[i].getVal()*scaleFactor*ratio

  #Scale
  scaledAmplitudesList.append(ROOT.RooRealVar("scaledMidPeakAmp"+str(i),"scaledMidPeakAmp"+str(i),BGT))
  scaledAmplitudes.add(scaledAmplitudesList[i+nLowerPeaks])
  scaledPdfList.append(midPeakGaussians[i])
  scaledPdfs.add(scaledPdfList[i+nLowerPeaks])

  print("{:.3f}".format(midPeakMeans[i].getVal())," {:.3f}".format(BGT))

  sum+=BGT
  

#Wave number for GTR
kf_gt = math.sqrt(E_projectile-E_thresh-gtrMean.getVal())
ratio=kf_ias/kf_gt

#Compute B(GT)
BGT = gtrAmp.getVal()*scaleFactor*ratio

#Scale
scaledAmplitudesList.append(ROOT.RooRealVar("scaledGtrAmp","scaledGtrAmp",BGT))
scaledAmplitudes.add(scaledAmplitudesList[-1])
scaledPdfList.append(gtrGauss)
scaledPdfs.add(scaledPdfList[-1])
sum+=BGT

print("{:.3f}".format(iasMean.getVal())+" (IAS)")
print("{:.3f}".format(gtrMean.getVal())," {:.3f}".format(BGT))

for i in range(0,nUpperPeaks):
  print("{:.3f}".format(upperPeakMeans[i].getVal())+" (non L=0)")

print("\n\nTotal B(GT) strength: "+str(sum)+"\n")

#######################
##PLot Normalized BGT##
#######################
c1.cd()
pad2 = ROOT.TPad("pad2","pad2",0.5,0.5,1.0,1.0)
pad2.Draw()
pad2.cd()
frame2=eVar.frame()
normalizedModel=ROOT.RooAddPdf("normalizedModel","normalizedModel",scaledPdfs,scaledAmplitudes)
normalizedModel.plotOn(frame2,ROOT.RooFit.Normalization(sum,ROOT.RooAbsReal.NumEvent))
frame2.Draw()

############
##Plot CDF##
############
c1.cd()
pad3 = ROOT.TPad("pad3","pad3",0.5,0.0,1.0,0.5)
pad3.Draw()
pad3.cd()
cfdFrame = eVar.frame();
cdf = normalizedModel.createCdf(ROOT.RooArgSet(eVar))
cdf.plotOn(cfdFrame,ROOT.RooFit.Normalization(sum,ROOT.RooAbsReal.NumEvent));
cfdFrame.Draw()
c1.Modified()
c1.Update()
try:
  input("Press enter to continue")
except SyntaxError:
  pass
