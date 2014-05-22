from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np

scale_factor = 0.965

ngc = np.loadtxt('Beutler_et_al_2013_ps_cmass_DR11v1_NGC_111_212_114_Yamamoto_16000000_2000.dat')
sgc = np.loadtxt('Beutler_et_al_2013_ps_cmass_DR11v1_SGC_84_163_86_Yamamoto_5000000_2000.dat')

ngc_theory = np.loadtxt('test_window_1.dat')
sgc_theory = np.loadtxt('test_window_2.dat')

ngc_theory_planck = np.loadtxt('test_planck_window_1.dat')
sgc_theory_planck = np.loadtxt('test_planck_window_2.dat')

for width in [8.8]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    ax = fig.add_subplot(111)
        
    plt.plot(ngc_theory[:,0],scale_factor*ngc_theory[:,1],color='k')
    plt.plot(sgc_theory[:,0],scale_factor*sgc_theory[:,1],color='r')

    plt.plot(ngc_theory[:,0],scale_factor*ngc_theory[:,2],color='k')
    plt.plot(sgc_theory[:,0],scale_factor*sgc_theory[:,2],color='r')

    plt.plot(ngc_theory_planck [:,0],scale_factor*ngc_theory_planck [:,1],color='k',linestyle='--')
    plt.plot(sgc_theory_planck [:,0],scale_factor*sgc_theory_planck [:,1],color='r',linestyle='--')

    plt.plot(ngc_theory_planck [:,0],scale_factor*ngc_theory_planck [:,2],color='k',linestyle='--')
    plt.plot(sgc_theory_planck [:,0],scale_factor*sgc_theory_planck [:,2],color='r',linestyle='--')

    # Data
    plt.errorbar(ngc[:,0],ngc[:,1],ngc[:,4],color='k',fmt='mo',markersize='3',label='NGC')
    plt.errorbar(sgc[:,0],sgc[:,1],sgc[:,4],color='r',fmt='mo',markersize='3',label='SGC')

    plt.errorbar(ngc[:,0],ngc[:,2],ngc[:,5],color='k',fmt='mo',markersize='3',label='')
    plt.errorbar(sgc[:,0],sgc[:,2],sgc[:,5],color='r',fmt='mo',markersize='3',label='')

    # legend
    leg = plt.legend(frameon=True,loc='upper right',numpoints=1)
    # remove box around legend
    leg.get_frame().set_edgecolor("white")
    leg.get_frame().set_alpha(0.8)

    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

    ax.set_yscale('log')
 
    # grid
    plt.grid(True, which="major", axis="both")

    # axes limits
    plt.ylim([500, 100000]); plt.xlim([0, 0.2]);

    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

    plt.xlabel(r'$k [h \, Mpc^{-1}]$')
    plt.ylabel(r'$P_{\ell} (k) [Mpc\, h^{-1}]^3$')

    #plt.savefig("pk_rsd_total_%dmm.pdf" % int(width*10),pad_inches=0.02)    

    plt.savefig("pk_rsd_total_%dmm.pdf" % int(width*10), bbox_inches='tight', pad_inches=0.02)
