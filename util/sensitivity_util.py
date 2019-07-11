'''
functions used in the generate and plot scripts
'''
#!/usr/bin/env python
import reactor_list_selector
import os
import sys
from re import sub
import numpy as np # numpy has to be imported before ROOT
import ROOT
ROOT.gROOT.SetBatch(True) # ROOT in batch mode

def setup_histogram():
    '''
    create histogram
    '''

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    bins_x_n = 10 #small = 10, medium = 40, large = 160
    bin_x_min = 0.0
    bin_x_max = 1.0
    bin_x_width = (bin_x_max-bin_x_min)/bins_x_n
    bins_x = [bin_x_min] + [(bin_x_min + (i+1)*bin_x_width) for i in range(bins_x_n)] # len must be +1 of array because... root.
    
    bins_y_n = 15 #small = 15, medium = 60, large = 240
    bin_y_min = 6.0e-5
    bin_y_max = 9.0e-5
    bin_y_width = (bin_y_max-bin_y_min)/bins_y_n
    bins_y = [bin_y_min*y_multiplier] + [(bin_y_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)] # len must be +1 of array because... root.

    h_2d = ROOT.TH2D("h_2d", "h_2d", bins_x_n, np.array(bins_x), bins_y_n, np.array(bins_y))
    
    h_2d.SetTitle("")
    h_2d.GetYaxis().SetTitle("#Deltam_{21}^{2} (x10^{-5} eV^{2})")
    h_2d.GetXaxis().SetTitle("sin^{2} #theta_{12}")
    return h_2d

def setup_histogram_loglin_sin():
    '''
    create histogram
    '''

    ##h_2d = ROOT.TH2D("h_2d", "h_2d", \
                       #nbins_delmsqr21, min_delmsqr21, max_delmsqr21,\
                       #nbins_sinsqrtheta12, min_sinsqrtheta12, max_sinsqrtheta12)

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    bins_x_n = 10 #small = 10, medium = 40, large = 160
    bin_x_min = 0.0
    bin_x_max = 1.0
    bin_x_width = (bin_x_max-bin_x_min)/bins_x_n
    bins_x = [bin_x_min] + [(bin_x_min + (i+1)*bin_x_width) for i in range(bins_x_n)] # len must be +1 of array because... root.

    bins_y_n = 40 #small = 40, medium = 80, large = 160
    bin_y_min = 1e-6 #1.5e-5
    bin_y_max = 2e-3#2e-3+1.5e-5
    bin_logy_min = np.log10(bin_y_min)
    bin_logy_max = np.log10(bin_y_max)
    bin_y_width = (bin_logy_max-bin_logy_min)/bins_y_n
    bins_y = [np.power(10, bin_logy_min)*y_multiplier] + [np.power(10, bin_logy_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)]
    
    h_2d = ROOT.TH2D("h_2d", "h_2d", bins_x_n-2, np.array(bins_x), bins_y_n, np.array(bins_y))  #medium

    # fill = ROOT.TF2("fill","xygaus",0,10,0,10) #test data dumy histogram
    # fill.SetParameters(1,0.5,0.5,6e-5,6e-5)
    # h2.FillRandom("fill")

    h_2d.SetTitle("")
    h_2d.GetYaxis().SetTitle("#Deltam_{21}^{2} (x10^{-5} eV^{2})")
    h_2d.GetXaxis().SetTitle("sin^{2} #theta_{12}")

    #h_2d.GetXaxis().SetNoExponent
    h_2d.GetYaxis().SetNoExponent(True)
    #h_2d.GetXaxis().SetLabelSize(0.025)
    #h_2d.GetYaxis().SetLabelSize(0.025)
    #h_2d.GetXaxis().SetMoreLogLabels()
    h_2d.GetYaxis().SetMoreLogLabels()
    h_2d.GetXaxis().SetTitleOffset(1.2)
    #h_2d.GetYaxis().SetTitleOffset(1.6)
    return h_2d

def setup_histogram_loglog_tan():
    '''
    create histogram
    '''

    ## histogram format:
    ##h_2d = ROOT.TH2D("h_2d", "h_2d", \
                       #nbins_delmsqr21, min_delmsqr21, max_delmsqr21,\
                       #nbins_sinsqrtheta12, min_sinsqrtheta12, max_sinsqrtheta12)

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    bins_x_n = 120 #small = 20, medium = 40, large = 160
    bin_x_min = 0.1
    bin_x_max = 6
    bin_logx_min = np.log10(bin_x_min)
    bin_logx_max = np.log10(bin_x_max)
    bin_x_width = (bin_logx_max-bin_logx_min)/bins_x_n
    bins_x = [np.power(10, bin_logx_min)] + [np.power(10, bin_logx_min + (i+1)*bin_x_width) for i in range(bins_x_n)]

    bins_y_n = 120 #small = 40, medium = 80, large = 160
    bin_y_min = 1.5e-5 #6.5e-5 
    bin_y_max = 2e-4+1.5e-5 #9e-5
    bin_logy_min = np.log10(bin_y_min)
    bin_logy_max = np.log10(bin_y_max)
    bin_y_width = (bin_logy_max-bin_logy_min)/bins_y_n
    bins_y = [np.power(10, bin_logy_min)*y_multiplier] + [np.power(10, bin_logy_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)] #*1e5 #axis scale factor for plotting

    h_2d = ROOT.TH2D("h_2d", "h_2d", bins_x_n, np.array(bins_x), bins_y_n, np.array(bins_y))  #medium

    # fill = ROOT.TF2("fill","xygaus",0,10,0,10) #test data dumy histogram
    # fill.SetParameters(1,0.5,0.5,6e-5,6e-5)
    # h2.FillRandom("fill")

    h_2d.SetTitle("")
    h_2d.GetYaxis().SetTitle("#Deltam_{21}^{2} (x10^{-5} eV^{2})")
    h_2d.GetXaxis().SetTitle("tan^{2} #theta_{12}")

    h_2d.GetXaxis().SetNoExponent(True)
    h_2d.GetYaxis().SetNoExponent(True)
    #h_2d.GetXaxis().SetLabelSize(0.025)
    #h_2d.GetYaxis().SetLabelSize(0.025)
    h_2d.GetXaxis().SetMoreLogLabels()
    h_2d.GetYaxis().SetMoreLogLabels()
    #h_2d.GetXaxis().SetTitleOffset(1.4)
    #h_2d.GetYaxis().SetTitleOffset(1.6)
    return h_2d

