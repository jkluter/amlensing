
from astropy.table import Table, MaskedColumn
import astrometry
import numpy as np
import microlensing
import matplotlib.pyplot as plt
from utils import DEG       
import plot_functions  as pf
import random
import os 
os.chdir(os.path.dirname(__file__))
folder = '../Images/Timescale'
os.makedirs(folder+'/', exist_ok=True)
s_circ = np.sin(np.linspace(0,2*np.pi,1000))
c_circ = np.cos(np.linspace(0,2*np.pi,1000))
def plot_timescales(new_source_ra, new_source_dec, max_ra,max_dec,shift, dist2):
    fig = plt.figure()
    index_max = np.where(shift == np.max(shift))[0][0]
    index = np.arange(len(shift))<= index_max
    plt.plot(new_source_ra,new_source_dec, color=pf.blue, )
    plt.plot(new_source_ra[shift>0.1],new_source_dec[shift>0.1],  color=pf.limitcolor, label = r"$t_{\rm aml}$")
    plt.plot(new_source_ra[(dist2<0.1)& index],new_source_dec[(dist2<0.1)& index], color=pf.green, label = r"$t_{0.1{\rm mas}}$")
    plt.plot(new_source_ra[index_max], new_source_dec[index_max], 'kx', label='maximum shift')
    plt.plot(new_source_ra[index_max]+0.1*s_circ, new_source_dec[index_max]+0.1*c_circ, 'k',label = r'$0.1 {\rm mas}$ distance from maximum' )
    plt.plot(0, 0, 'ko', label='unlensed position')
    plt.plot(0.1*s_circ, 0.1*c_circ, 'k--',label = r'$\delta\theta_{+} = 0.1 {\rm mas}$ ' )
    shift_max = np.max(shift)
    #plt.plot((np.max(shift)-0.1)*s_circ, (np.max(shift)-0.1)*c_circ, 'k:' )
    plt.gca().axis('equal')
    plt.ylim(np.min(new_source_dec)-0.15*shift_max, np.max(new_source_dec)+0.15*shift_max)
    plt.xlim(np.min(new_source_ra)-0.15*shift_max, np.max(new_source_ra)+0.15*shift_max)
    plt.legend(loc= 'upper right')
    plt.xlabel(r'$\delta\theta_{+,\,\alpha}$ [mas]')
    plt.ylabel(r'$\delta\theta_{+,\,\delta}$ [mas]')
    return fig



def timescales(tab):

    epoch_predef = np.arange(-150, 150, 1/365)
    pos_sun_predef = astrometry.pos_sun(epoch_predef)
    ts_01 = []
    ts_50pc = []
    ts_05 = []
    for i, row in enumerate(tab,1):
        if i%500==0 : print(i)
        ra = row['ra']
        dec = row['dec']
        pmra = row['pmra']
        pmdec = row['pmdec']
        parallax = row['parallax']
        ThetaE = row['ThetaE']
        ob_ra = row['ob_ra']
        ob_dec = row['ob_dec']
        ob_pmra = row['ob_pmra'] if not np.isnan(row['ob_pmra']) else (row['ob_displacement_ra_doubled'] if not np.isnan(row['ob_displacement_ra_doubled']) else 0)
        ob_pmdec = row['ob_pmdec'] if not np.isnan(row['ob_pmdec']) else (row['ob_displacement_dec_doubled'] if not np.isnan(row['ob_displacement_dec_doubled']) else 0)
        ob_parallax = row['ob_parallax'] if not np.isnan(row['ob_parallax'])  else 4.740/75 *np.sqrt(ob_pmra*ob_pmra +ob_pmdec*ob_pmdec)
        
        pos_lens = astrometry.movePm_parallax(ra, dec, pmra, pmdec, parallax, epoch_predef, 
            earthVec = pos_sun_predef)
        pos_source = astrometry.movePm_parallax(ob_ra, ob_dec, ob_pmra, ob_pmdec, ob_parallax,
            epoch_predef, earthVec = pos_sun_predef)
        pos_lens_max=astrometry.movePm_parallax(ra, dec, pmra, pmdec, parallax, row['TCA']-2016)
        pos_source_max=astrometry.movePm_parallax(ob_ra, ob_dec, ob_pmra, ob_pmdec, ob_parallax,
            row['TCA']-2016)
        lens_ra, lens_dec = astrometry.dirVecToCelCoos(pos_lens)
        source_ra, source_dec = astrometry.dirVecToCelCoos(pos_source)
        lens_max_ra, lens_max_dec = astrometry.dirVecToCelCoos(pos_lens_max)
        source_max_ra, source_max_dec = astrometry.dirVecToCelCoos(pos_source_max)
        dist_max=astrometry.getGCDist_3Vec(pos_lens_max, pos_source_max)*3.6e6

        dist = astrometry.getGCDist_3Vec(pos_lens, pos_source)*3.6e6
        shift_plus = microlensing.calc_shift_plus(dist, ThetaE)

        shift_max = microlensing.calc_shift_plus(dist_max, ThetaE)
        print(ob_ra,ob_dec,ob_pmra,ob_pmdec,ob_parallax)
        print(np.max(shift_plus), np.where(shift_plus == np.max(shift_plus)))
        index_max = np.where(shift_plus == np.max(shift_plus))[0][0]



        new_source_ra = shift_plus * (source_ra-lens_ra)*3.6e6/dist * np.cos(lens_dec*DEG)
        new_source_dec = shift_plus * (source_dec-lens_dec)*3.6e6/dist

        max_source_ra = shift_max * (source_max_ra-lens_max_ra)*3.6e6/dist_max * np.cos(lens_max_dec*DEG)
        max_source_dec = shift_max * (source_max_dec-lens_max_dec)*3.6e6/dist_max


        dist2 = np.sqrt((new_source_ra-max_source_ra)**2+(new_source_dec-max_source_dec)**2)
        # plt.close()
        if np.min(dist2)< 0.1:
            enter_01 = np.where(dist2 < 0.1)[0][0]
            t_enter_01 = epoch_predef[enter_01]
            ts_01.append(row['TCA'] - t_enter_01-2016)
        else: 
            ts_01.append(1/365)
        # plt.plot(epoch_predef,dist2)
        # plt.plot(epoch_predef[dist2<0.1],dist2[dist2<0.1], 'r',label='dist2<0.1')
        # plt.legend
        if ts_01[-1]>100: 
            fig = plot_timescales(new_source_ra, new_source_dec, max_source_ra, max_source_dec,shift_plus, dist2)
            fig.savefig(f'{folder}/{row["ob_source_id"]}.png')
            print(row["shift_plus"], max(shift_plus))
            plt.show()

        plt.close()

        if np.max(dist2)< 0.5: 
            ts_05.append(np.nan) 
        elif np.min(dist2)< 0.5:
            enter_05 = np.where(dist2 < 0.5)[0][0]
            t_enter_05 = epoch_predef[enter_05]
            ts_05.append(row['TCA'] - t_enter_05-2016)
        else:
            ts_05.append(1/365)
        if min(dist2) < 0.5*max(shift_max):
            enter_50pc = np.where(dist2 < 0.5*shift_max)[0][0]
            t_enter_50pc = epoch_predef[enter_50pc]
            ts_50pc.append(row['TCA'] - t_enter_50pc-2016)
        else:
            ts_50pc.append(1/365)
    print(len(tab.colnames))
    tab = add_collums(ts_01,ts_05,ts_50pc,tab)  
    print(len(tab.colnames))
    plot_hist(ts_01,ts_05,ts_50pc,tab)
    return tab
def add_collums(ts_01,ts_05,ts_50pc,tab):
    ts_01 = np.minimum(ts_01, 100)
    ts_01 = np.maximum(ts_01, 1/365)
    ts_50pc = np.maximum(ts_50pc, 1/365)
    ts_50pc = np.minimum(ts_50pc, 100)

    T01 = MaskedColumn(ts_01, 
        dtype = 'float64',unit = 'years', 
        description = 'Timescale over which the 2D positional shift differs by 0.1mas from its maximum value'
            r'(\lvert$\delta_\boldsymbol{theta}_{+,\,TCA} - $\delta_\boldsymbol{theta}_{+,\,TCA\,-\,t_0.1mas}\rvert = 0.1\mas'
        )
    T50PC = MaskedColumn(ts_50pc, 
        dtype = 'float64',unit = 'years', 
        description = 'Timescale over which the 2D positional shift differs by 50% from its maximum value'
            r'(\lvert$\delta_\boldsymbol{\theta}_{+,\,TCA} - $\delta_\boldsymbol{\theta}_{+,\,TCA\,-\,t_50pc}\rvert = 0.5 \delta_{\theta}'
        )
    tab.add_column(T01, name='t_0.1mas', index=38)
    tab.add_column(T50PC, name='t_50pc', index=39)
    return tab

def plot_hist(ts_01,ts_05,ts_50pc,tab):
    fig = plt.figure()
    plt.hist(ts_01, bins=40, range=(0,40), histtype='step', color='green', label=r'$t_{0.1{ \rm mas}}$')
    #plt.hist(ts_05, bins=40, range=(0,40), histtype='step', color='blue', label=r'$t_{0.5{ \rm mas}}$')
    plt.hist(ts_50pc, bins=40, range=(0,40), histtype='step', color='red', label=r'$t_{50\%}$')
    plt.hist(tab['t_aml'], bins=40, range=(0,40), histtype='step', color='black', label=r'$t_{aml}$')
    plt.gca().set_yscale('log')
    plt.xlim(0,30)
    plt.legend()
    fig.savefig(f'{folder}/timescales.png')
    plt.show()


if __name__ == "__main__":
    infile = '../amlensing_FINAL.fits'
    tab = Table.read(infile) 
    tab = timescales(tab)
    for j,i in enumerate(tab.keys().copy()):
            tab.meta['TCOMM%i'%(j+1)] = tab[i].description
    tab.write(f'{infile[:-5]}.timescales.fits', format = 'fits', overwrite = True)

