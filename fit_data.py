import ROOT
import sys
import numpy
import math

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

energies=[]
counts=[]
graph=ROOT.TGraph()
graph.SetPoint(0,0,0)
pointNum=1
#Load data sets
inpFile=open(sys.argv[1],"r")
for line in inpFile:
  line=line.strip("\n")
  lineParts=line.split(",")
  E=float(lineParts[0])
  contents=float(lineParts[1])
  graph.SetPoint(pointNum,E,contents)
  pointNum+=1

if not "208" in sys.argv[1]:
  nBins=600
else:
  nBins=300
hist=ROOT.TH1D("hist",";Energy (MeV);Counts",nBins,0,30)
for i in range(0,nBins):
  val=graph.Eval(hist.GetBinCenter(i+1))
  hist.SetBinContent(i+1,val)
  hist.SetBinError(i+1,0.001)

#Draw
c1=ROOT.TCanvas("c1","c1")
#hist.Draw()
c1.Modified()
c1.Update()

#Make RooFit Observable
eVar=ROOT.RooRealVar("eVar","Excitation Energy (MeV)",0,30)
argSet=ROOT.RooArgSet(eVar)
argList=ROOT.RooArgList(eVar)

#Make RooDataHist
dataHist=ROOT.RooDataHist("dataHist","dataHist",argList,hist)
histPdf=ROOT.RooHistPdf("histPdf","histPdf",argSet,dataHist,1)

#Quasi-free continuum background
if "208" in sys.argv[1]:
  E_projectile=200
  #Mass (208Bi + 3H) - (Mass 208Pb + 3He) in MeV/c^2
  E_thresh=2.897
  E_gs=E_projectile-E_thresh
  #Proton binding energy: Mass (208Bi [- electron if atomic mass]) - (207Pb + proton) in MeV/c^2
  S_p=3.707
  #Not sure how to calculate, but think it's not dependent on the nucleus
  E_ejectile_free = 200.761
  #Coulomb barrier for the proton: http://hyperphysics.phy-astr.gsu.edu/hbase/NucEne/coubar.html
  #interacting with the final state nucleus
  B_coulomb=14.38
  #E_xn = neutron emission excitation energy: Mass (207Bi + n) - (208Bi)
  E_xn = 6.887
else:
  E_projectile=200
  #Mass (208Bi + 3H) - (Mass 208Pb + 3He) in MeV/c^2
  E_thresh=0.518
  E_gs=E_projectile-E_thresh
  #Proton binding energy: Mass (208Bi [- electron if atomic mass]) - (207Pb + proton) in MeV/c^2
  S_p=5.158
  #Not sure how to calculate, but think it's not dependent on the nucleus
  E_ejectile_free = 200.761
  #Coulomb barrier for the proton: http://hyperphysics.phy-astr.gsu.edu/hbase/NucEne/coubar.html
  #interacting with the final state nucleus
  B_coulomb=15.2838
  #E_xn = neutron emission excitation energy: Mass (207Bi + n) - (208Bi)
  E_xn = 5.549
  

#Definitions
#E_t = E_shift - eVar
E_shift=ROOT.RooRealVar("E_shift","E_shift",E_gs)
EQF = ROOT.RooRealVar("EQF","EQF",E_ejectile_free-(S_p+E_xn+B_coulomb))
E0 = ROOT.RooRealVar("E0","E0",E_gs-S_p)
T=ROOT.RooRealVar("T","T",100) #Recommeneded value
W=ROOT.RooRealVar("W","W",22) #Recommended value
argList=ROOT.RooArgList(eVar,E_shift,E0,EQF,T,W)
qfc=ROOT.RooGenericPdf("qfc","qfc","(1-exp(((E_shift - eVar) - E0)/T)) / (1+(((E_shift - eVar) - EQF)/W)^2)" ,argList)

#
if "208" in sys.argv[1]:
  if not "Krasz" in sys.argv[1]:
    qfcAmp=ROOT.RooRealVar("qfcAmp","qfcAmp",10000,1000,20000)
  else:
    qfcAmp=ROOT.RooRealVar("qfcAmp","qfcAmp",1000,300,12000)
  
  if not "Krasz" in sys.argv[1]:
    iasMean=ROOT.RooRealVar("iasMean","IAS",16.2,16.0,16.5)
    iasAmp=ROOT.RooRealVar("iasAmp","iasAmp",700,50,900)
  else:
    iasMean=ROOT.RooRealVar("iasMean","IAS",15.0,13.8,15.0)
    iasAmp=ROOT.RooRealVar("iasAmp","iasAmp",30,30,400)
  iasSigma=ROOT.RooRealVar("iasSigma","iasSigma",0.16,0.14,0.20)
  iasGauss=ROOT.RooGaussian("iasGauss","iasGauss",eVar,iasMean,iasSigma)

  gtrMean=ROOT.RooRealVar("gtrMean","GTR",15,14,18)
  if not "Krasz" in sys.argv[1]:
    gtrSigma=ROOT.RooRealVar("gtrSigma","gtrSigma",2,1,2.4)
    gtrAmp=ROOT.RooRealVar("gtrAmp","gtrAmp",2000,500,33000)
  else:
    gtrSigma=ROOT.RooRealVar("gtrSigma","gtrSigma",1.5,0.8,1.6)
    gtrAmp=ROOT.RooRealVar("gtrAmp","gtrAmp",700,150,700)
  gtrGauss=ROOT.RooGaussian("gtrGauss","gtrGauss",eVar,gtrMean,gtrSigma)
  
  sdrMean=ROOT.RooRealVar("sdrMean","SDR",22,20,29)
  sdrSigma=ROOT.RooRealVar("sdrSigma","sdrSigma",3,1.4,8)
  sdrAmp=ROOT.RooRealVar("sdrAmp","sdrAmp",100,100,3000)
  sdrGauss=ROOT.RooGaussian("sdrGauss","sdrGauss",eVar,sdrMean,sdrSigma)
  
  gtr1Mean=ROOT.RooRealVar("gtr1Mean","GTR1",8,5.0,8.5)
  gtr1Sigma=ROOT.RooRealVar("gtr1Sigma","gtr1Sigma",1.0,.3,5.5)
  gtr1Amp=ROOT.RooRealVar("gtr1Amp","gtr1Amp",20,10,2000)
  gtr1Gauss=ROOT.RooGaussian("gtr1Gauss","gtr1Gauss",eVar,gtr1Mean,gtr1Sigma)
  
  gtr2Mean=ROOT.RooRealVar("gtr2Mean","GTR2",9.5,9.0,10)
  gtr2Sigma=ROOT.RooRealVar("gtr2Sigma","gtr2Sigma",1.0,.5,3.0)
  gtr2Amp=ROOT.RooRealVar("gtr2Amp","gtr2Amp",20,5,2000)
  gtr2Gauss=ROOT.RooGaussian("gtr2Gauss","gtr2Gauss",eVar,gtr2Mean,gtr2Sigma)
  
  if not "Krasz" in sys.argv[1]:
    n12Mean=ROOT.RooRealVar("n12Mean","12N (bgnd)",15.5,15.4,15.6)
    n12Sigma=ROOT.RooRealVar("n12Sigma","n12Sigma",0.05,0.05,0.10 )
    n12Amp=ROOT.RooRealVar("n12Amp","n12Amp",20,20,2000)
    n12Gauss=ROOT.RooGaussian("n12Gauss","n12Gauss",eVar,n12Mean,n12Sigma)
  
  ###########
  #Gaussians#
  ###########
  g1Mean=ROOT.RooRealVar("g1Mean","g1Mean",0.6,0.4,0.8)
  g1Sigma=ROOT.RooRealVar("g1Sigma","g1Sigma",0.1,.08,0.12)
  g1Amp=ROOT.RooRealVar("g1Amp","g1Amp",15,2,100)
  g1Gauss=ROOT.RooGaussian("g1Gauss","g1Gauss",eVar,g1Mean,g1Sigma)
  
  g2Mean=ROOT.RooRealVar("g2Mean","g2Mean",1.5,1.3,1.7)
  g2Sigma=ROOT.RooRealVar("g2Sigma","g2Sigma",0.1,0.08,0.15)
  g2Amp=ROOT.RooRealVar("g2Amp","g2Amp",20,2,300)
  g2Gauss=ROOT.RooGaussian("g2Gauss","g2Gauss",eVar,g2Mean,g2Sigma)
  
  g3Mean=ROOT.RooRealVar("g3Mean","g3Mean",2.75,2.6,2.85)
  g3Sigma=ROOT.RooRealVar("g3Sigma","g3Sigma",0.12,0.1,0.16)
  g3Amp=ROOT.RooRealVar("g3Amp","g3Amp",20,2,300)
  g3Gauss=ROOT.RooGaussian("g3Gauss","g3Gauss",eVar,g3Mean,g3Sigma)
  
  g4Mean=ROOT.RooRealVar("g4Mean","g4Mean",3.15,3.05,3.2)
  g4Sigma=ROOT.RooRealVar("g4Sigma","g4Sigma",0.05,0.05,0.06)
  g4Amp=ROOT.RooRealVar("g4Amp","g4Amp",10,1,200)
  g4Gauss=ROOT.RooGaussian("g4Gauss","g4Gauss",eVar,g4Mean,g4Sigma)
  
  g5Mean=ROOT.RooRealVar("g5Mean","g5Mean",4.0,3.8,4.2)
  g5Sigma=ROOT.RooRealVar("g5Sigma","g5Sigma",0.1,0.09,0.13)
  g5Amp=ROOT.RooRealVar("g5Amp","g5Amp",10,1,300)
  g5Gauss=ROOT.RooGaussian("g5Gauss","g5Gauss",eVar,g5Mean,g5Sigma)
  
  g6Mean=ROOT.RooRealVar("g6Mean","g6Mean",4.8,4.6,5.0)
  g6Sigma=ROOT.RooRealVar("g6Sigma","g6Sigma",0.08,0.05,0.09)
  g6Amp=ROOT.RooRealVar("g6Amp","g6Amp",10,1,200)
  g6Gauss=ROOT.RooGaussian("g6Gauss","g6Gauss",eVar,g6Mean,g6Sigma)
  
  g7Mean=ROOT.RooRealVar("g7Mean","g7Mean",5.5,5.0,7.0)
  g7Sigma=ROOT.RooRealVar("g7Sigma","g7Sigma",0.28,0.20,0.31)
  g7Amp=ROOT.RooRealVar("g7Amp","g7Amp",20,10,500)
  g7Gauss=ROOT.RooGaussian("g7Gauss","g7Gauss",eVar,g7Mean,g7Sigma)
  
  if not "Krasz" in sys.argv[1]:
    pdfList=ROOT.RooArgList(qfc,gtrGauss,iasGauss,sdrGauss,gtr1Gauss,gtr2Gauss,n12Gauss)
  else:
    pdfList=ROOT.RooArgList(qfc,gtrGauss,iasGauss,sdrGauss,gtr1Gauss,gtr2Gauss)
  gaussPdfList=ROOT.RooArgList(g1Gauss,g2Gauss,g3Gauss,g4Gauss,g5Gauss,g6Gauss,g7Gauss)
  pdfList.add(gaussPdfList)
  
  if not "Krasz" in sys.argv[1]:
    ampList=ROOT.RooArgList(qfcAmp,gtrAmp,iasAmp,sdrAmp,gtr1Amp,gtr2Amp,n12Amp)#,g5Amp,g6Amp)
  else:
    ampList=ROOT.RooArgList(qfcAmp,gtrAmp,iasAmp,sdrAmp,gtr1Amp,gtr2Amp)#,g5Amp,g6Amp)
  gaussAmpList=ROOT.RooArgList(g1Amp,g2Amp,g3Amp,g4Amp,g5Amp,g6Amp,g7Amp)
  ampList.add(gaussAmpList)
  
  model=ROOT.RooAddPdf("model","model",pdfList,ampList)
else:
  qfcAmp=ROOT.RooRealVar("qfcAmp","qfcAmp",500000,450000,600000)
  
  iasMean=ROOT.RooRealVar("iasMean","IAS",17.5,17,18)
  iasSigma=ROOT.RooRealVar("iasSigma","iasSigma",0.18,0.15,0.20)
  iasAmp=ROOT.RooRealVar("iasAmp","iasAmp",100000,10000,250000)
  iasGauss=ROOT.RooGaussian("iasGauss","iasGauss",eVar,iasMean,iasSigma)

  gtrMean=ROOT.RooRealVar("gtrMean","GTR",17,16.5,18.5)
  gtrSigma=ROOT.RooRealVar("gtrSigma","gtrSigma",1.0,1.0,2)
  gtrAmp=ROOT.RooRealVar("gtrAmp","gtrAmp",35000,30000,200000)
  gtrGauss=ROOT.RooGaussian("gtrGauss","gtrGauss",eVar,gtrMean,gtrSigma)
  
  gtr1Mean=ROOT.RooRealVar("gtr1Mean","GTR1",6,3,8.5)
  gtr1Sigma=ROOT.RooRealVar("gtr1Sigma","gtr1Sigma",1.0,1.0,4)
  gtr1Amp=ROOT.RooRealVar("gtr1Amp","gtr1Amp",5000,4000,70000)
  gtr1Gauss=ROOT.RooGaussian("gtr1Gauss","gtr1Gauss",eVar,gtr1Mean,gtr1Sigma)
  '''
  gtr2Mean=ROOT.RooRealVar("gtr2Mean","GTR2",9.5,9.0,10)
  gtr2Sigma=ROOT.RooRealVar("gtr2Sigma","gtr2Sigma",2.0,1.5,5.0)
  gtr2Amp=ROOT.RooRealVar("gtr2Amp","gtr2Amp",3000,2000,50000)
  gtr2Gauss=ROOT.RooGaussian("gtr2Gauss","gtr2Gauss",eVar,gtr2Mean,gtr2Sigma)
  '''
  sdrMean=ROOT.RooRealVar("sdrMean","SDR",24,23,25)
  sdrSigma=ROOT.RooRealVar("sdrSigma","sdrSigma",3.5,2.0,5)
  sdrAmp=ROOT.RooRealVar("sdrAmp","sdrAmp",250000,150000,350000)
  sdrGauss=ROOT.RooGaussian("sdrGauss","sdrGauss",eVar,sdrMean,sdrSigma)
  
  ############
  #Background#
  ############
  n12Mean=ROOT.RooRealVar("n12Mean","12N (bgnd)",15.5,15.4,15.95)
  n12Sigma=ROOT.RooRealVar("n12Sigma","n12Sigma",0.18,0.02,0.19)
  n12Amp=ROOT.RooRealVar("n12Amp","n12Amp",20,20,2000)
  n12Gauss=ROOT.RooGaussian("n12Gauss","n12Gauss",eVar,n12Mean,n12Sigma)
  
  n12Mean2=ROOT.RooRealVar("n12Mean2","12N (bgnd)",16.8,16.6,17.)
  n12Sigma2=ROOT.RooRealVar("n12Sigma2","n12Sigma2",0.10,0.01,0.15)
  n12Amp2=ROOT.RooRealVar("n12Amp2","n12Amp2",1000,30,2000)
  n12Gauss2=ROOT.RooGaussian("n12Gauss2","n12Gauss2",eVar,n12Mean2,n12Sigma2)
  
  f16Mean=ROOT.RooRealVar("f16Mean","16F (bgnd)",20.6,20,21)
  f16Sigma=ROOT.RooRealVar("f16Sigma","f16Sigma",0.18,0.05,0.19)
  f16Amp=ROOT.RooRealVar("f16Amp","f16Amp",20,20,8000)
  f16Gauss=ROOT.RooGaussian("f16Gauss","f16Gauss",eVar,f16Mean,f16Sigma)
  
  f16Mean2=ROOT.RooRealVar("f16Mean2","16F (bgnd)",18.75,18.4,19.)
  f16Sigma2=ROOT.RooRealVar("f16Sigma2","f16Sigma2",0.18,0.10,0.19)
  f16Amp2=ROOT.RooRealVar("f16Amp2","f16Amp2",20,20,8000)
  f16Gauss2=ROOT.RooGaussian("f16Gauss2","f16Gauss2",eVar,f16Mean2,f16Sigma2)
  
  f16Mean3=ROOT.RooRealVar("f16Mean3","16F (bgnd)",14.0,13.5,14.5)
  f16Sigma3=ROOT.RooRealVar("f16Sigma3","f16Sigma3",0.18,0.05,0.19)
  f16Amp3=ROOT.RooRealVar("f16Amp3","f16Amp3",20,20,8000)
  f16Gauss3=ROOT.RooGaussian("f16Gauss3","f16Gauss3",eVar,f16Mean3,f16Sigma3)
  
  pdfList=ROOT.RooArgList(qfc,iasGauss,gtrGauss,gtr1Gauss,sdrGauss)
  backgroundPdfList=ROOT.RooArgList(n12Gauss,n12Gauss2,f16Gauss,f16Gauss2,f16Gauss3)
  pdfList.add(backgroundPdfList)
  
  ampList=ROOT.RooArgList(qfcAmp,iasAmp,gtrAmp,gtr1Amp,sdrAmp)
  backgroundAmpList=ROOT.RooArgList(n12Amp,n12Amp2,f16Amp,f16Amp2,f16Amp3)
  ampList.add(backgroundAmpList)
  model=ROOT.RooAddPdf("model","model",pdfList,ampList)

res = model.fitTo(dataHist,ROOT.RooFit.Save(1),ROOT.RooFit.Extended(1))

c2=ROOT.TCanvas("c2","c2")
frame=eVar.frame()
dataHist.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Binning(300))
model.plotOn(frame)
model.plotOn(frame,ROOT.RooFit.Components("qfc"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGray))
model.plotOn(frame,ROOT.RooFit.Components("iasGauss"),ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(frame,ROOT.RooFit.Components("gtr1Gauss"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kViolet))
model.plotOn(frame,ROOT.RooFit.Components("sdrGauss"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(frame,ROOT.RooFit.Components("gtrGauss"),ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.kMagenta))
#model.plotOn(frame,ROOT.RooFit.Components("gtr2Gauss"),ROOT.RooFit.LineStyle(ROOT.kDotted),ROOT.RooFit.LineColor(ROOT.kViolet))
if not "Krasz" in sys.argv[1]:
  model.plotOn(frame,ROOT.RooFit.Components("n12Gauss"),ROOT.RooFit.LineStyle(ROOT.kDotted),ROOT.RooFit.LineColor(ROOT.kRed))
if "232" in sys.argv[1]:
  model.plotOn(frame,ROOT.RooFit.Components("n12Gauss2"),ROOT.RooFit.LineStyle(ROOT.kDotted),ROOT.RooFit.LineColor(ROOT.kRed))
  model.plotOn(frame,ROOT.RooFit.Components("f16Gauss"),ROOT.RooFit.LineStyle(ROOT.kDotted),ROOT.RooFit.LineColor(ROOT.kGreen))
  model.plotOn(frame,ROOT.RooFit.Components("f16Gauss2"),ROOT.RooFit.LineStyle(ROOT.kDotted),ROOT.RooFit.LineColor(ROOT.kGreen))
  model.plotOn(frame,ROOT.RooFit.Components("f16Gauss3"),ROOT.RooFit.LineStyle(ROOT.kDotted),ROOT.RooFit.LineColor(ROOT.kGreen))

frame.Draw()
c2.Modified()
c2.Update()

if "208" in sys.argv[1]:
  A=208.
  Z=82.
  N=A-Z
  
  '''
  #R^2 = unit GT x-section / unit Fermi x-section
  unitGT=109./math.pow(A,0.65) #From fit
  unitF = 72./math.pow(A,1.06) #From fit
  Rsquared=unitGT/unitF
  
  #Alternate approach
  Rsquared = math.pow(E_projectile/55.0,2)
  
  print("R2 = "+str(Rsquared))
  
  BF=N-Z
  
  BF_amp = iasAmp.getVal()
  
  scaleFactor = BF / (BF_amp * Rsquared)
  '''
  BF=N-Z
  BF_amp = iasAmp.getVal()
  if not "Krasz" in sys.argv[1]:
    scaleFactor = BF/(1.47*BF_amp) #multiply by ratio GT/F 0-degree cross sections
  else:
    scaleFactor=1
    #TODO
  
  GT1_amp = ROOT.RooRealVar("GT1_amp","GT1_amp",g1Amp.getVal()*scaleFactor)
  GT2_amp = ROOT.RooRealVar("GT2_amp","GT2_amp",g2Amp.getVal()*scaleFactor)
  GT3_amp = ROOT.RooRealVar("GT3_amp","GT3_amp",g3Amp.getVal()*scaleFactor)
  GT4_amp = ROOT.RooRealVar("GT4_amp","GT4_amp",g4Amp.getVal()*scaleFactor)
  GT5_amp = ROOT.RooRealVar("GT5_amp","GT5_amp",g5Amp.getVal()*scaleFactor)
  GT6_amp = ROOT.RooRealVar("GT6_amp","GT6_amp",g6Amp.getVal()*scaleFactor)
  GT7_amp = ROOT.RooRealVar("GT7_amp","GT7_amp",g7Amp.getVal()*scaleFactor)
  GT8_amp = ROOT.RooRealVar("GT8_amp","GT8_amp",gtr1Amp.getVal()*scaleFactor)
  GT9_amp = ROOT.RooRealVar("GT9_amp","GT9_amp",gtr2Amp.getVal()*scaleFactor)
  GT10_amp = ROOT.RooRealVar("GT10_amp","GT10_amp",gtrAmp.getVal()*scaleFactor)
  sum=GT1_amp.getVal()+GT2_amp.getVal()+GT3_amp.getVal()+GT4_amp.getVal()+GT5_amp.getVal()+GT6_amp.getVal()+GT7_amp.getVal()+GT8_amp.getVal()+GT9_amp.getVal()+GT10_amp.getVal()
  print("Total GT strength: "+str(sum))
  
  normalizedGaussPdfList=ROOT.RooArgList(g1Gauss,g2Gauss,g3Gauss,g4Gauss,g5Gauss,g6Gauss,g7Gauss)
  normalizedPdfList=ROOT.RooArgList(gtr1Gauss,gtr2Gauss,gtrGauss)
  normalizedPdfList.add(normalizedGaussPdfList)
  
  normalizedGaussAmpList=ROOT.RooArgList(GT1_amp,GT2_amp,GT3_amp,GT4_amp,GT5_amp,GT6_amp,GT7_amp)
  normalizedAmpList=ROOT.RooArgList(GT8_amp,GT9_amp,GT10_amp)
  normalizedAmpList.add(normalizedGaussAmpList)
  
  normalizedModel=ROOT.RooAddPdf("normalizedModel","normalizedModel",normalizedPdfList,normalizedAmpList)
  
  c1.cd()
  frame2=eVar.frame(100)
  normalizedModel.plotOn(frame2,ROOT.RooFit.Normalization(sum,ROOT.RooAbsReal.NumEvent))
  frame2.Draw()
  c1.Modified()
  c1.Update()
  
  c3=ROOT.TCanvas("c3","c3")
  cfdFrame = eVar.frame();
  cdf = normalizedModel.createCdf(ROOT.RooArgSet(eVar))
  cdf.plotOn(cfdFrame,ROOT.RooFit.Normalization(sum,ROOT.RooAbsReal.NumEvent));
  cfdFrame.Draw()
  c3.Modified()
  c3.Update()
  
else:
  A=232.
  Z=90.
  N=A-Z
  
  BF=N-Z
  BF_amp = iasAmp.getVal()
  scaleFactor = BF/(1.47*BF_amp) #multiply by ratio GT/F 0-degree cross sections
  
  GT1_amp = ROOT.RooRealVar("GT1_amp","GT1_amp",gtr1Amp.getVal()*scaleFactor)
  #GT2_amp = ROOT.RooRealVar("GT2_amp","GT2_amp",gtr2Amp.getVal()*scaleFactor)
  GT3_amp = ROOT.RooRealVar("GT3_amp","GT3_amp",gtrAmp.getVal()*scaleFactor)
  sum=GT1_amp.getVal()+GT3_amp.getVal()
  print("Total GT strength: "+str(sum))

  normalizedPdfList=ROOT.RooArgList(gtr1Gauss,gtrGauss)

  normalizedAmpList=ROOT.RooArgList(GT1_amp,GT3_amp)

  normalizedModel=ROOT.RooAddPdf("normalizedModel","normalizedModel",normalizedPdfList,normalizedAmpList)

  c1.cd()
  frame2=eVar.frame(30)
  normalizedModel.plotOn(frame2,ROOT.RooFit.Normalization(sum,ROOT.RooAbsReal.NumEvent))
  frame2.Draw()
  c1.Modified()
  c1.Update()
  
  c3=ROOT.TCanvas("c3","c3")
  cfdFrame = eVar.frame();
  cdf = normalizedModel.createCdf(ROOT.RooArgSet(eVar))
  cdf.plotOn(cfdFrame,ROOT.RooFit.Normalization(sum,ROOT.RooAbsReal.NumEvent));
  cfdFrame.Draw()
  c3.Modified()
  c3.Update()
  
try:
  input("Press enter to continue")
except SyntaxError:
  pass
