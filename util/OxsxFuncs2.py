import re
import os
import subprocess
import string
import time
import numpy as np
import argparse

import sensitivity_util as sutil
import ROOT
ROOT.gROOT.SetBatch(True) # ROOT in batch mode

import scipy.stats as st

#ENVIRONMENT_PATH = "/home/blakei/env_rat-6.3.6.sh"
ENVIRONMENT_PATH = "/home/blakei/Env_rat-KL.sh"

def write_shell_script(commands, shell_script_name):
    with open(shell_script_name, "w") as f:
        f.write(commands)
    f.close()


# Pruning down root files to an ntuple with entries needed for coincidence tagging
def PruneforCoincTag():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh prune_for_coinctag.cpp
./prune_for_coinctag '/data/snoplus/blakei/antinu/mc/rootfiles/Bruce_LAB_5day_flux1000_s*' /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                           input root                                                                       output pruned ntuple

# produces a fake data set to be used as data input for a fit, a single entry _oxsx.ntuple.root file containing oscillated (or unoscillated) EPrompt, produced by applying coincidence tagging on the pruned ntuple made with above func.
def FakeDataEPrompt():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh prune_to_KE_EPrompt.cpp
./prune_to_KE_EPrompt /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root /data/snoplus/blakei/antinu/test/Bruce5yr1000flux 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 E1""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                    input pruned ntp                                     output _oxsx fake data ntuple(no .root)   nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 s13 E1orKE

def FakeDataKE():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh prune_to_KE_EPrompt.cpp
./prune_to_KE_EPrompt /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root /data/snoplus/blakei/antinu/test/Bruce5yr1000flux 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 KE""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                    input pruned ntp                                    output _oxsx fake data ntuple(no .root)   nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 s13 E1orKE

# NOT NEEDED FOR FIT, file produces root file w/ histograms showing results of applying coincidence tagging on a pruned ntuple
def OscEPrompt():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh OscEPrompt.cpp
./OscEPrompt /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root EPromptOut.root 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 MC""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                    input pruned ntp                  output _oxsx fake data ntuple     nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 s13 MCorData looping method

def PruneFlatTree():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
./prune_flat_tree /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root EPromptOut.root 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 MC""").format(ENVIRONMENT_PATH), shell = True)

##############################
#are you using lh2dgaus??????#
#   using right csv file?    #
##############################
lh2d_loc = '/data/snoplus/blakei/antinu/mc/plots/sensitivity/SNOP_1yr_FakeAsimovData_LowStatPdfs_d21_7.4e-5_s12_0.297_0.0215_E2.6_R5500_loglog_tan.root'

textfile_loc = '/data/snoplus/blakei/antinu/temp/fitresulttxtfiles/Results'
textfile_loc += "_"
textfile_name = (textfile_loc.rsplit('/',1)[1]).split('_')[0]
tempfile_loc = '/data/snoplus/blakei/antinu/temp/filltempfiles/temp'
print textfile_name

def LH2Dsubmit():
    Emin = 2.6  #1
    Emax = 8.  #2
    numbins = 16 #3

    infofile = '/home/blakei/code/workshop/oxsx/1000km.ratdb' #4
    constraintcsv = '/home/blakei/antinu_analysis/passes/processed/passes_snop_rat6169_2017_flux1_scintFitter_cleanround1_plots.csv' #1000km 100pass E2.6 R1_R2 5500 #7

    dataconstraints_loc = '/data/snoplus/blakei/antinu/temp/data_constraints.root' #5
    lh2d = sutil.setup_histogram_loglog_tan() #sutil.setup_histogram() #setup_histogram()
    
    #Data:
    #1
    #data_loc = ''  # for use if not using fake data #6
    # OR
    faked21 = 7.4e-5 #4.8e-5 #8
    fakes12 = 0.297 #0.36 #9
    fakes13 = 0.0215 #10
    
    PHWRunoscfile = '/data/snoplus/blakei/antinu/mc/ntuples/rat6169_2017/PHWR_flux1_day360_passcombined100_cleanround1.root' #11
    PWRunoscfile = '/data/snoplus/blakei/antinu/mc/ntuples/rat6169_2017/1000km_no3cad_flux1_day360_passcombined100_cleanround1.root' #12
    
    subprocess.check_call(("""source {0}
cd ~/oxsx/examples/
. ~/oxsx/bin/compile_with_ratIwan.sh lh2dInit2.cpp
./lh2dInit2 {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}
. ~/oxsx/bin/compile_with_ratIwan.sh lh2d2.cpp""").format(ENVIRONMENT_PATH,Emin,Emax,numbins,infofile,dataconstraints_loc.split('.')[0],tempfile_loc+'.root',constraintcsv,faked21,fakes12,fakes13,PHWRunoscfile,PWRunoscfile), shell = True)#"array"), shell = True)
##############################
#are you using lh2dgaus??????#
##############################
    
    bins_x_n = lh2d.GetXaxis().GetNbins()
    bins_y_n = lh2d.GetYaxis().GetNbins()

    submitcommandsarray = []
    
    totalbins = bins_x_n*bins_y_n
    fitspersub = 15
    numsubs = totalbins/fitspersub
    print 'total bins: '+str(totalbins)
    print 'fitspersub: '+str(fitspersub)
    print 'number of subs: '+str(numsubs)
    dividesevenly = True
    if (totalbins % fitspersub != 0):
        dividesevenly = False
        print "total bins and fits per sub doesn't divide evenly!"
        #print "wrong num of fits per batch submit"
        #return 0

    commands = """source {0}
cd ~/oxsx/examples/
""".format(ENVIRONMENT_PATH)
    
    k = 0
    SUB = 1
    for i in range(bins_x_n):
        for j in range(bins_y_n):
            s12 = np.power(np.sin(np.arctan(np.sqrt(lh2d.GetXaxis().GetBinCenter(i+1)))),2)
            d21 = lh2d.GetYaxis().GetBinCenter(j+1)*10**-5
            
##############################
#are you using lh2dgaus??????#
##############################            
            commands += ("""
./lh2d2 {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}
""").format(PHWRunoscfile,dataconstraints_loc.split('.')[0],infofile,d21,s12,fakes13,i+1,j+1,Emin,Emax,numbins,textfile_loc+'{0}.txt'.format(SUB),tempfile_loc+'{0}.root'.format(SUB),PWRunoscfile,faked21,fakes12,fakes13)
            
            k += 1
            #print k
            if (k % fitspersub == 0):
                SUB += 1
                submitcommandsarray.append(commands)
                #print commands
                commands = """source {0}
cd ~/oxsx/examples/
""".format(ENVIRONMENT_PATH)

            if (k == (bins_x_n*bins_y_n) and dividesevenly == False ):
                print 'last k!'
                if (commands != """source {0}
cd ~/oxsx/examples/
""".format(ENVIRONMENT_PATH)):
                    print 'adding last commands!'
                    submitcommandsarray.append(commands)
                
                
                
    
    print SUB, len(submitcommandsarray)
    if (dividesevenly == True):
        if (numsubs != len(submitcommandsarray)):
            print "fits per sub does not match the submit array size!!"
            return 0

    for sub in range(len(submitcommandsarray)):
        shell_script_name = ("/data/snoplus/blakei/antinu/temp/subs/h2sub_{0}of{1}.sh").format(sub+1,len(submitcommandsarray))
        write_shell_script(submitcommandsarray[sub], shell_script_name)
        #subprocess.check_call(submitcommandsarray[sub], shell = True)  # No Submission of jobs
        subprocess.check_call("qsub -l cput=01:59:59 {0}".format(shell_script_name), shell = True)  #Submission of jobs
        
def LH2Dplot():

    lh2d = sutil.setup_histogram_loglog_tan() #sutil.setup_histogram() #setup_histogram()
    
    textfiles = []
    for path, subdirs, files in os.walk(textfile_loc.rsplit('/',1)[0]):
        for filename in files:
            if (filename.split('_')[0] == textfile_name):
                textfiles.append(os.path.join(path,filename))
            
    print textfiles

    for k in range(len(textfiles)):
        textfile = open(textfiles[k], "r")
        textfilearray = textfile.readlines()
        textfile.close()
        for i in range(len(textfilearray)):
            bin_s12 = (textfilearray[i].split('biny')[0]).split('binx_')[1]
            bin_d21 = (textfilearray[i].split('LHval')[0]).split('biny_')[1]
        
            LHval = (textfilearray[i].split('LHval_')[1]).split('\n')[0]

            #print 'bins12: ',int(bin_s12),' bind21: ',int(bin_d21),' LH: ',ROOT.Double(LHval)
            lh2d.SetBinContent(int(bin_s12),int(bin_d21),ROOT.Double(LHval))

        subprocess.check_call(("""
mv {0} {1}""").format(textfiles[k], textfiles[k].rsplit('/fitresulttxtfiles/',1)[0] + "/oldfitresulttxtfiles/" + textfiles[k].rsplit('/',1)[1]), shell = True)

    c = ROOT.TCanvas()
    c.SetLogx(True)
    c.SetLogy(True)
    lh2d.Draw("COLZ")
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.2)
    c.SetBottomMargin(0.2)


    fout = ROOT.TFile(lh2d_loc, 'RECREATE')
    lh2d.Write()
    c.Write()
    fout.Close()

def H2interact():
    """
    fin = ROOT.TFile.Open(lh2d_loc, 'READ')
    
    lh2d = fin.Get("h_2d")
    
    bins_x_n = lh2d.GetXaxis().GetNbins()
    bins_y_n = lh2d.GetYaxis().GetNbins()

    Max = lh2d.GetMaximum()
    lh2d.SetMaximum(Max - 2.)
    lh2d.SetMinimum(Max - 4.5)
    lh2d.SetFillColor(2)
    
    c = ROOT.TCanvas()
    c.SetLogx(True)
    c.SetLogy(True)
    #lh2d.Draw("CONTZ")
    lh2d.Draw("box")#"scat=20")
    
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.2)
    c.SetBottomMargin(0.2)

    #fout = ROOT.TFile(lh2d_loc.rsplit('.root',1)[0] + "Doctored.root" , 'RECREATE')
    fout = ROOT.TFile(lh2d_loc.rsplit('.root',1)[0] + "cont3.root" , 'RECREATE')
    c.Write()
    lh2d.Write()
    
    fout.Close()
    """
    
    fin1 = ROOT.TFile.Open(lh2d_loc, 'READ')
    fin2 = ROOT.TFile.Open(lh2d_loc, 'READ')
    fin3 = ROOT.TFile.Open(lh2d_loc, 'READ')
    
    lh2d1 = fin1.Get("h_2d")
    lh2d2 = fin2.Get("h_2d")
    lh2d3 = fin3.Get("h_2d")
    
    bins_x_n = lh2d1.GetXaxis().GetNbins()
    bins_y_n = lh2d1.GetYaxis().GetNbins()

    c = ROOT.TCanvas()
    c.SetLogx(True)
    c.SetLogy(True)

    Max = lh2d1.GetMaximum()
    lh2d3.SetMaximum(Max - 2.)
    lh2d3.SetMinimum(Max - 4.5)
    lh2d3.SetFillColor(2)
    lh2d3.DrawCopy("same")
    lh2d2.SetMaximum(Max - 0.5)
    lh2d2.SetMinimum(Max - 2.)
    lh2d2.SetFillColor(8)
    lh2d2.DrawCopy("same")
    lh2d1.SetMinimum(Max - 0.5)
    lh2d1.SetFillColor(4)
    lh2d1.DrawCopy("same")
    
    #c.SetLeftMargin(0.2)
    #c.SetRightMargin(0.2)
    #c.SetBottomMargin(0.2)

    #fout = ROOT.TFile(lh2d_loc.rsplit('.root',1)[0] + "Doctored.root" , 'RECREATE')
    fout = ROOT.TFile(lh2d_loc.rsplit('.root',1)[0] + "cont.root" , 'RECREATE')
    c.Write()
    lh2d1.Write()
    lh2d2.Write()
    lh2d3.Write()
    
    fout.Close()
    
    

def H2Contour():
    conf_lev1 = 0.68 #0.95
    conf_lev2 = 0.95 #0.99
    conf_lev3 = 0.9973
    Zscore1 = st.norm.ppf(conf_lev1+((1-conf_lev1)/2))
    Zscore2 = st.norm.ppf(conf_lev2+((1-conf_lev2)/2))
    Zscore3 = st.norm.ppf(conf_lev3+((1-conf_lev3)/2))
    delL1 = ((Zscore1)**2)/2
    delL2 = ((Zscore2)**2)/2
    delL3 = ((Zscore3)**2)/2

    print delL1, delL2, delL3
    
    fin = ROOT.TFile.Open(lh2d_loc, 'READ')
    
    lh2d = fin.Get("h_2d")

    Max = lh2d.GetMaximum()
    bin = lh2d.GetMaximumBin()
    maxxbin = ROOT.Long(0)
    maxybin = ROOT.Long(0)
    binz = ROOT.Long(0)
    #Max = 0
    lh2d.GetBinXYZ(bin,maxxbin,maxybin,binz)
    #contours = np.array([Max - 4.5,Max - 2.,Max - 0.5])
    contours = np.array([Max - delL3,Max - delL2,Max - delL1])
    maxx = lh2d.GetXaxis().GetBinCenter(maxxbin)
    maxy = lh2d.GetYaxis().GetBinCenter(maxybin)
    print maxx, maxy, Max

    
#### from Jeff ########    
    c = ROOT.TCanvas("c1", "c1", int(1000*1.18), 1000)
    c.SetLogx(True)
    c.SetLogy(True)
    lh2d.GetXaxis().SetRangeUser(0.1, 6.05)
    lh2d.GetYaxis().SetRangeUser(1.5, 25)
    lh2d.SetContour(3,contours)
    
    ###### 1 ########
    lh2d.Draw("colz")
    lh2d.SetMinimum(Max - delL3)
    
    ###### 2 #######
    #lh2d.Draw("contz same")
    
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.2)
    c.SetBottomMargin(0.2)
    
    '''
    c = ROOT.TCanvas()
    c.SetLogx(True)
    c.SetLogy(True)
    lh2d.SetContour(3,contours)
    
    #lh2d.SetLineColor(4)
    #### 1 ####
    # alternate way to  contour, comment out if and use Draw method below if only single peak shown in data
    #lh2d.SetMinimum(Max - delL3)
    #lh2d.DrawCopy("colz")
    
    #### 2 ####
    lh2d.Draw("contz same")
    
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.2)
    c.SetBottomMargin(0.2)
    '''
    fout = ROOT.TFile(lh2d_loc.rsplit('.root',1)[0] + "Contoured.root" , 'RECREATE')
    c.Write()
    
    fout.Close()
    

# Currently input data for LH2D is generated within the example, to input data own data, comment out custom produced data part and uncomment out generic datantuple filling part
def LH2DEVindex():
    subprocess.check_call(("""source {0}
cd ~/oxsx/examples/
. ~/oxsx/bin/compile.sh LH2dEVindexPlot.cpp
./LH2dEVindexPlot /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root USUALLYPUTDATA 3CAD.ratdb /data/snoplus/blakei/antinu/Test.root 1yrflux13CADConstrain 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297""").format(ENVIRONMENT_PATH), shell = True)
#arguments:        input pruned unosc ntp                                Data to Fit          reactoinfofile temp file for filling outputLH2dfilename(no .root) nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 (for output naming)


def OscFit():
    subprocess.check_call(("""source {0}
cd ~/oxsx/examples/
. ~/oxsx/bin/compile.sh OscillationFitTest.cpp
./OscillationFitTest /data/snoplus/blakei/antinu/test/Bruce5yr1000fluxKE_oxsx.root /data/snoplus/blakei/antinu/test/3CAD_36000evsKEds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root 3CAD.ratdb Test.root""").format(ENVIRONMENT_PATH), shell = True)
#arguments:        input pruned unosc ntp                 Data to Fit                   reactoinfofile temp file for filling outputLH2dfilename  nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 (for output naming)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Select func to run: (-r) \n prunefortag \n fake data \n fakedatake \n osceprompt \n lh2d \n oscfit")
    parser.add_argument("-r", dest="func", help="need to select which function to run", default=0)
    args = parser.parse_args()

    Func = str(args.func)
    if (Func.find("prunefortag") != -1):
        PruneforCoincTag()
    if (Func.find("fakedata") != -1):
        FakeDataEPrompt()
    if (Func.find("fakedatake") != -1):
        FakeDataKE()
    if (Func.find("osceprompt") != -1):
        OscEPrompt()
    if (Func.find("lh2dsub") != -1):
        LH2Dsubmit()
    if (Func.find("lh2dplot") != -1):
        LH2Dplot()
    if (Func.find("lh2dchange") != -1):
        H2interact()
    if (Func.find("oscfit") != -1):
        OscFit()
    if (Func.find("lh2dcont") != -1):
        H2Contour()
