'''
functions used in the generate and plot scripts
'''
#!/usr/bin/env python
import reactor_list_selector
import os
import sys
from re import sub
import numpy as np # numpy has to be imported before ROOT
#import pandas as pd
import ROOT
ROOT.gROOT.SetBatch(True) # ROOT in batch mode

def get_reactor_info(reactor_list_name, position):
    '''
    run the reactor_list_selector utility
    '''

    reactor_list_reactors_file = '/home/lidgard/rat3/data/REACTORS.ratdb'
    reactor_list_status_file = '/home/lidgard/rat3/data/REACTORS_STATUS.ratdb'

    # get the reactors from the specified list
    if reactor_list_name is not "All":
        reactor_list_name = reactor_list_selector.get_reactor_list(reactor_list_name, reactor_list_reactors_file)

    # then get info
    reactor_ratdb_info = reactor_list_selector.get_reactor_ratdb_info(reactor_list_reactors_file)
    reactor_status_ratdb_info = reactor_list_selector.get_reactor_status_ratdb_info(reactor_list_status_file)

    reactor_info = reactor_list_selector.get_reactor_info(reactor_ratdb_info, \
                        reactor_status_ratdb_info, \
                        reactor_list_name, \
                        position)

    reactor_info =  {k.lower(): v for k, v in reactor_info.items()} # lower case
    reactor_info =  {sub(r'\W+','', k): v for k, v in reactor_info.items()} #remove all but alphanumeric chars
    reactor_info =  {k.replace(" ", ""): v for k, v in reactor_info.items()} #remove whitespace

    return reactor_info

def create_reactor_info(output_file, reactor_list_name, position):
    '''
    write the reactor_list_selector output to file
    '''
    reactor_info = get_reactor_info(reactor_list_name, position)
    reactor_list_selector.write_output_file(output_file, reactor_info)
    
def create_single_parameter_info(output_file, n_jobs, del_21, sin_12, sin_13, n_reps):
    '''
    write the reactor_list_selector output to file
    '''

    # get list of all permutations of parameters
    parameter_list = []
    for i in range(n_reps):
        parameter_list.append([del_21, sin_12, sin_13])

    print "Number of submitted jobs: "+str(len(parameter_list))

    # write the csv file for each 'job'
    for job in range(n_jobs):
        output_file_i = output_file+"_"+str(job)+".csv"
        df = pd.DataFrame(parameter_list, columns=['del_21', 'sin_12', 'sin_13'])
        df.to_csv(output_file_i, index=None, header=True)

def create_parameter_info(output_file, h_2d, n_jobs):
    '''
    write the reactor_list_selector output to file
    '''

    y_multiplier = 1e5

    # create lists of parameters # fix sin_sqr_theta_13
    xvalue_12s = [h_2d.GetXaxis().GetBinCenter(i+1) for i in range(h_2d.GetNbinsX())]
    del_21s = [h_2d.GetYaxis().GetBinCenter(i+1)/y_multiplier for i in range(h_2d.GetNbinsY())]
    sin_13 = 0.02303

    sin_12s = [np.power(np.sin(np.arctan(np.sqrt(i))),2) for i in xvalue_12s] #convert from tan^2 nin values to sin^2 parameter values
    #sin_12s = [np.power(np.sin(0.5*np.arcsin(np.sqrt(i))),2) for i in sin_2_12s] #convert from sin^2 nin values to sin2^2 parameter values

    # get list of all permutations of parameters
    parameter_list = []
    for del_21 in del_21s:
        for sin_12 in sin_12s:
            parameter_list.append([del_21, sin_12, sin_13])

    #find how many parameters per job (minimum is 1)
    n_split = int(round(len(parameter_list)/float(n_jobs)))
    if n_split < 1:
        n_split = 1
        n_jobs = len(parameter_list)

    # split the total list of all the permutations of the parameters into 'n_jobs' lists
    # (one sublist of n_split permutations per job)
    parameter_segments = [parameter_list[x:x+n_split] for x in xrange(0, len(parameter_list), n_split)]

    print "Number of submitted jobs: "+str(len(parameter_list))
    print "each submitted containing (number of fits): "+str(n_split)

    # write the csv file for each 'job'
    for job in range(n_jobs):
        output_file_i = output_file+"_"+str(job)+".csv"
        df = pd.DataFrame(np.array(parameter_segments[job]), columns=['del_21', 'sin_12', 'sin_13'])
        df.to_csv(output_file_i, index=None, header=True)

# osc1
#delmsqr21 = 7.4e-05
#sinsqrtheta12 = 0.297
#sinsqrtheta13 = 0.0215

# osc2
#delmsqr21 = 7.58e-05
#sinsqrtheta12 = 0.359
#sinsqrtheta13 = 0.02303

##h_2d = ROOT.TH2D("h_2d", "h_2d", \
                    #nbins_delmsqr21, min_delmsqr21, max_delmsqr21,\
                    #nbins_sinsqrtheta12, min_sinsqrtheta12, max_sinsqrtheta12)
def setup_histogram(name = "h_2d"):
    '''
    create histogram
    '''

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    bins_x_n = 40 #small = 10, medium = 40, large = 160
    bin_x_min = 0.0
    bin_x_max = 1.0
    bin_x_width = (bin_x_max-bin_x_min)/bins_x_n
    bins_x = [bin_x_min] + [(bin_x_min + (i+1)*bin_x_width) for i in range(bins_x_n)] # len must be +1 of array because... root.

    bins_y_n = 60 #small = 15, medium = 60, large = 240
    bin_y_min = 6.0e-5
    bin_y_max = 9.0e-5
    bin_y_width = (bin_y_max-bin_y_min)/bins_y_n
    bins_y = [bin_y_min*y_multiplier] + [(bin_y_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)] # len must be +1 of array because... root.

    h_2d = ROOT.TH2D(name, name, bins_x_n, np.array(bins_x), bins_y_n, np.array(bins_y))

    h_2d.SetTitle("")
    h_2d.GetYaxis().SetTitle("#Deltam_{21}^{2} (x10^{-5} eV^{2})")
    h_2d.GetXaxis().SetTitle("sin^{2} #theta_{12}")
    return h_2d

def setup_histogram_loglin_sin(name = "h_2d"):
    '''
    create histogram
    '''

    ##h_2d = ROOT.TH2D("h_2d", "h_2d", \
                       #nbins_delmsqr21, min_delmsqr21, max_delmsqr21,\
                       #nbins_sinsqrtheta12, min_sinsqrtheta12, max_sinsqrtheta12)

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    bins_x_n = 40 #small = 10, medium = 40, large = 160
    bin_x_min = 0.0
    bin_x_max = 1.0
    bin_x_width = (bin_x_max-bin_x_min)/bins_x_n
    bins_x = [bin_x_min] + [(bin_x_min + (i+1)*bin_x_width) for i in range(bins_x_n)] # len must be +1 of array because... root.

    bins_y_n = 80 #small = 40, medium = 80, large = 160
    bin_y_min = 1.5e-5
    bin_y_max = 2e-3+1.5e-5
    bin_logy_min = np.log10(bin_y_min)
    bin_logy_max = np.log10(bin_y_max)
    bin_y_width = (bin_logy_max-bin_logy_min)/bins_y_n
    bins_y = [np.power(10, bin_logy_min)*y_multiplier] + [np.power(10, bin_logy_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)]

    h_2d = ROOT.TH2D(name, name, bins_x_n-2, np.array(bins_x), bins_y_n, np.array(bins_y))  #medium

    # fill = ROOT.TF2("fill","xygaus",0,10,0,10) #test data dumy histogram
    # fill.SetParameters(1,0.5,0.5,6e-5,6e-5)
    # h2.FillRandom("fill")

    h_2d.SetTitle("")
    h_2d.GetYaxis().SetTitle("#Deltam_{21}^{2} (x10^{-5} eV^{2})")
    h_2d.GetXaxis().SetTitle("sin^{2} #theta_{12}")

    #h_2d.GetXaxis().SetNoExponent(True)
    h_2d.GetYaxis().SetNoExponent(True)
    #h_2d.GetXaxis().SetLabelSize(0.025)
    #h_2d.GetYaxis().SetLabelSize(0.025)
    #h_2d.GetXaxis().SetMoreLogLabels()
    h_2d.GetYaxis().SetMoreLogLabels()
    h_2d.GetXaxis().SetTitleOffset(1.2)
    #h_2d.GetYaxis().SetTitleOffset(1.6)
    return h_2d

def setup_histogram_KL1_loglin_sin(name = "h_2d"):
    '''
    create histogram
    '''

    ##h_2d = ROOT.TH2D("h_2d", "h_2d", \
                       #nbins_delmsqr21, min_delmsqr21, max_delmsqr21,\
                       #nbins_sinsqrtheta12, min_sinsqrtheta12, max_sinsqrtheta12)

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    bins_x_n = 100 #small = 10, medium = 40, large = 160
    bin_x_min = 0.0
    bin_x_max = 1.0
    bin_x_width = (bin_x_max-bin_x_min)/bins_x_n
    bins_x = [bin_x_min] + [(bin_x_min + (i+1)*bin_x_width) for i in range(bins_x_n)] # len must be +1 of array because... root.

    bins_y_n = 180 #small = 40, medium = 80, large = 160
    bin_y_min = 1.e-6
    bin_y_max = 2e-3
    bin_logy_min = np.log10(bin_y_min)
    bin_logy_max = np.log10(bin_y_max)
    bin_y_width = (bin_logy_max-bin_logy_min)/bins_y_n
    bins_y = [np.power(10, bin_logy_min)*y_multiplier] + [np.power(10, bin_logy_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)]

    h_2d = ROOT.TH2D(name, name, bins_x_n, np.array(bins_x), bins_y_n, np.array(bins_y))  #medium

    # fill = ROOT.TF2("fill","xygaus",0,10,0,10) #test data dumy histogram
    # fill.SetParameters(1,0.5,0.5,6e-5,6e-5)
    # h2.FillRandom("fill")

    h_2d.SetTitle("")
    h_2d.GetYaxis().SetTitle("#Deltam_{21}^{2} (x10^{-5} eV^{2})")
    h_2d.GetXaxis().SetTitle("sin^{2} 2*#theta_{12}")#"sin^{2} #theta_{12}")

    #h_2d.GetXaxis().SetNoExponent(True)
    h_2d.GetYaxis().SetNoExponent(True)
    #h_2d.GetXaxis().SetLabelSize(0.025)
    #h_2d.GetYaxis().SetLabelSize(0.025)
    #h_2d.GetXaxis().SetMoreLogLabels()
    h_2d.GetYaxis().SetMoreLogLabels()
    h_2d.GetXaxis().SetTitleOffset(1.2)
    #h_2d.GetYaxis().SetTitleOffset(1.6)

    h_2d.GetXaxis().SetRangeUser(bin_x_min, bin_x_max)
    h_2d.GetYaxis().SetRangeUser(bin_y_min, bin_y_max)
    
    return h_2d

def setup_histogram_loglog_tan(name = "h_2d"):
    '''
    create histogram
    '''

    ## histogram format:
    ##h_2d = ROOT.TH2D("h_2d", "h_2d", \
                       #nbins_delmsqr21, min_delmsqr21, max_delmsqr21,\
                       #nbins_sinsqrtheta12, min_sinsqrtheta12, max_sinsqrtheta12)

    ROOT.gStyle.SetOptStat(0)
    y_multiplier = 1e5

    # delmsqr21 = 7.58e-05
    # sinsqrtheta12 = 0.359
    # sinsqrtheta13 = 0.02303

    bins_x_n = 180 #small = 20, medium = 40, large = 160
    bin_x_min = 0.03 #np.power(np.tan(np.arcsin(np.sqrt(0.355))),2)
    bin_x_max = np.power(np.tan(np.arcsin(np.sqrt(0.365))),2)
    bin_logx_min = np.log10(bin_x_min)
    bin_logx_max = np.log10(bin_x_max)
    bin_x_width = (bin_logx_max-bin_logx_min)/bins_x_n
    bins_x = [np.power(10, bin_logx_min)] + [np.power(10, bin_logx_min + (i+1)*bin_x_width) for i in range(bins_x_n)]

    bins_y_n = 180 #small = 40, medium = 80, large = 160
    bin_y_min = 1.e-6
    bin_y_max = 2e-3
    #bin_y_min = 6e-05
    #bin_y_max = 9e-05
    bin_logy_min = np.log10(bin_y_min)
    bin_logy_max = np.log10(bin_y_max)
    bin_y_width = (bin_logy_max-bin_logy_min)/bins_y_n
    bins_y = [np.power(10, bin_logy_min)*y_multiplier] + [np.power(10, bin_logy_min + (i+1)*bin_y_width)*y_multiplier for i in range(bins_y_n)] #*1e5 #axis scale factor for plotting

    h_2d = ROOT.TH2D(name, name, bins_x_n, np.array(bins_x), bins_y_n, np.array(bins_y))  #medium

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
    
    h_2d.GetXaxis().SetRangeUser(bin_x_min, bin_x_max)
    h_2d.GetYaxis().SetRangeUser(bin_y_min, bin_y_max)
    return h_2d
