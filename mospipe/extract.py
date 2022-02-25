import os
import glob
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import scipy.ndimage as nd

from astropy.table import Table
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

from astropy.modeling import models, fitting
from astropy.modeling.models import Gaussian2D
from scipy.signal import fftconvolve

import mpdaf.obj
from mpdaf.obj import airtovac, vactoair

gau = models.Gaussian1D(mean=0, stddev=1)

from grizli.utils_c import interp
from grizli import utils
utils.set_warnings()

band_lims = {'Y': (9600, 11360), 
             'J': (11440, 13560),
             'H': (14580, 18150), 
             'K': (18880, 24160)}

plt.rcParams['figure.max_open_warning'] = 100
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['image.interpolation'] = 'Nearest'

def optimal_extract(file, show_sn=False, slx=slice(300,-300), prof_model=models.Gaussian1D, prof_keys=['mean','stddev'], fit_profile=False, rescale_errors=True, prof_sigma=3, prof_offset=0, escl=2, zi=0, binf=64, flux_corr=None, to_vacuum=False, limits_2d=[-0.18, 0.18], gwidth=(4,2), clip_edge=5, suffix=''):
    """
    Optimal extraction from 2D fits file
    """
    #print(file)
    im = pyfits.open(file)
    valid = np.isfinite(im[0].data)
    if valid.sum() == 0:
        print('Bad spectrum')
        return False
    
    band = im[0].header['FILTER']
    sh = im[0].data.shape

    if file.endswith('_sp.fits'):
        is_mine = True
        wht = im['WHT'].data
        valid &= wht > 0
        sig = 1/np.sqrt(wht)
        sig[~valid] = 0
        
    else:
        sig = pyfits.open(file.replace('_eps','_sig'))[0].data
        is_mine = False
        
    yp, xp = np.indices(sh)
    oky = valid.sum(axis=0)
    vx = oky > 3
    
    xarr = np.arange(sh[-1])
    if is_mine:
        h = im[0].header
        lam_coeffs = np.array([h[f'LAMCOEF{i}'] for i in range(h['LAMORDER']+1)])
        wave = np.polyval(lam_coeffs, xarr-h['CRPIX1'])
        #print('xx', wave.min(), wave.max(), lam_coeffs)
        wcs = pywcs.WCS(im[0].header)
        wave = wcs.all_pix2world(xarr, xarr*0., 0)[0]*1.e10
        
        #sh = outsci.shape
        #xarr = np.arange(sh[1])
        yp = np.ones(sh[1]) + sh[0]/2
        lam_wcs = pywcs.WCS(im[0].header)
        #wave, _y = lam_wcs.all_pix2world(xarr, yp, 0)
        #wave *= 1.e10 # A

        #pix_lim = np.interp(band_lims[band], wave, xarr, left=0, right=2048)
        
        pix_lim = (np.array(band_lims[band])-lam_coeffs[-1])/lam_coeffs[-2] + h['CRPIX1']
        
    else:
        wcs = pywcs.WCS(im[0].header)
        wave = wcs.all_pix2world(xarr, xarr*0., 0)[0]*1.e10
        pix_lim = wcs.all_world2pix(np.array(band_lims[band])/1.e10, [0,0], 0)[0] #wave[iw], xarr[iw])

    if flux_corr is not None:
        try:
            to_flam = flux_corr[0](flux_corr[1], wave)
        except:
            to_flam = flux_corr(wave)
        
        vx &= np.isfinite(to_flam)
        if band in ['K']:
            vx &= to_flam/np.nanmin(to_flam) < 6
        else:
            vx &= to_flam/np.nanmin(to_flam) < 6
            
    else:
        to_flam = 1.0
    
    vx = np.where(vx)[0]
    xvalid = vx.min(), vx.max()
    valid &= (xp > xvalid[0]) & (xp < xvalid[1])
    
    if slx is None:
        slx = slice(*xvalid)
        
    if is_mine:
        targname = os.path.basename(file.split('_sp.fits')[0]) + suffix
    else:
        targname = os.path.basename(file.split('_eps.fits')[0]) + suffix
    
    if os.path.exists(file.replace('_eps','_itime')):
        itime = pyfits.open(file.replace('_eps','_itime'))
        exptime = np.nanpercentile(itime[0].data[itime[0].data > 0], [50,90])
    else:
        if is_mine:
            exptime = im[0].header['EXPTIME'], im[0].header['EXPTIME']
        else:
            exptime = [0,0]
    
    y0 = im[0].header['CRPIX2']-im[0].header['CRVAL2']
    
    if gwidth is not None:
        ivar = 1/sig**2
        ivar[~np.isfinite(ivar) | ~valid] = 0

        gau = Gaussian2D(x_mean=0, x_stddev=gwidth[0], y_mean=0, y_stddev=gwidth[1])
        xgarr = np.arange(-4*gwidth[0], 4.1*gwidth[0], 1)
        ygarr = np.arange(-4*gwidth[1], 4.1*gwidth[1], 1)

        xp, yp = np.meshgrid(xgarr, ygarr)
        gm = gau(xp, yp)


        sci = im[0].data*1
        sci[~valid] = 0

        num = fftconvolve(sci*ivar, gm, mode='same')
        den = fftconvolve(ivar, gm**2, mode='same')

        smoothed = num/den*valid
        
        if show_sn:
            smoothed *= np.sqrt(den)
        
        yarr = np.arange(smoothed.shape[0])
        ysl = np.abs(yarr-y0) < 10
        
        perc = np.nanpercentile(smoothed[ysl,:][valid[ysl,:]], [16,50,84])
        limits_2d = perc[1] - 3*np.diff(perc)[0], perc[1] + 3*np.diff(perc)[1]
        
        if show_sn:
            lmax = np.clip(limits_2d[1], 5, 40)
            #print('xxx', lmax)
            limits_2d = [-lmax, lmax]
            
    else:
        smoothed = im[0].data
    
    smoothed[~valid] = 0
    
    figs = []
    
    fig, axes = plt.subplots(1,2, figsize=(12,3), gridspec_kw={'width_ratios':[3,1]}, sharey=True)
    ax = axes[0]
    figs.append(fig)
    
    ax.imshow(smoothed, origin='lower', vmin=limits_2d[0], vmax=limits_2d[1], cmap='gray')
    ax.set_aspect('auto')
    
    if slx.stop < 0:
        ax.vlines([slx.start, sh[1]+slx.stop], 0, sh[0]-1, color='w', linewidth=3, alpha=0.35)
        ax.vlines([slx.start, sh[1]+slx.stop], 0, sh[0]-1, color='r', linewidth=1, alpha=0.35)
    else:
        ax.vlines([slx.start, slx.stop], 0, sh[0]-1, color='w', linewidth=3, alpha=0.35)
        ax.vlines([slx.start, slx.stop], 0, sh[0]-1, color='r', linewidth=1, alpha=0.35)
    
    ax.hlines(y0+np.array([-10,10]), 0, sh[1]-1, 
               color='w', linewidth=3, alpha=0.35)
    ax.hlines(y0+np.array([-10,10]), 0, sh[1]-1, 
               color='r', linewidth=1, alpha=0.35)
    
    ivar = 1/sig**2
    sci = im[0].data*1
    imask = (sig == 0) | ~np.isfinite(ivar) | ~np.isfinite(sci) | (~valid)
    ivar[imask] = 0
    sci[imask] = 0
    
    iw = np.where(np.isfinite(wave))[0]
    #print('x limits: ', sh, pix_lim)
    ax.set_xlim(*pix_lim)
    ax.xaxis.set_major_locator(MultipleLocator(200))
    #ax.set_xticklabels([])
    xt = ax.get_xticks()
    ax.set_xticks(xt[2:-2])
    ax.set_xticklabels(np.cast[int](xt[2:-2]))
    
    yt = ax.get_yticks()
    for j in [-3, -2]:
        ax.text(0.01*(pix_lim[1]-pix_lim[0])+pix_lim[0], yt[j], f'{int(yt[j])}', ha='left', va='center', 
            fontsize=7, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
    
    ax.set_yticklabels([])
    
    new_sci, new_ivar = sci, ivar
    ax.text(0.98, 0.98, targname, ha='right', va='top', 
            transform=ax.transAxes, fontsize=8, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
    
    #fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax = axes[1]
    #figs.append(fig)

    prof = (new_sci*new_ivar)[:,slx].sum(axis=1)/(new_ivar[:,slx].sum(axis=1))
    yarr = np.arange(len(prof))*1.
    ax.plot(prof, yarr)
    
    y0 = im[0].header['CRPIX2']-im[0].header['CRVAL2']+prof_offset
    keys = {prof_keys[0]:y0, prof_keys[1]:prof_sigma}
    
    prof_mask = np.isfinite(prof) & (np.abs(yarr-y0) < 10)
    #print('Prof: ', yarr.shape, prof_mask.sum())
    
    gau = prof_model(amplitude=prof[prof_mask].max(), **keys)
    gau.bounds['mean'] = (y0-8,y0+8)
    gau.bounds['stddev'] = (1, 4)
    
    fit_status = True
    
    if fit_profile:
        fitter = fitting.LevMarLSQFitter()
        try:
            gau = fitter(gau, yarr[prof_mask & (prof > 0)], prof[prof_mask & (prof > 0)])
        except:
            fit_status = False
            
    if fit_status:
        prof_offset = gau.parameters[gau.param_names.index(prof_keys[0])] + im[0].header['CRVAL2']
    
    ax.plot(gau(yarr), yarr)
    
    ymax = 1.5*gau.amplitude.value
    ax.set_xlim(-ymax, ymax)
    ax.set_xticklabels([])
    ax.text(gau.amplitude.value, 0.05*sh[0], f'{gau.amplitude.value:.4f}', 
            ha='center', va='center', bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
    
    #ax.plot(yarr[prof_mask], gau(yarr)[prof_mask])

    ax.hlines(y0+np.array([-10,10]),*plt.xlim(),  
               color='r', linewidth=1, alpha=0.35)
    
    ax.grid()
        
    fig.tight_layout(pad=0.5)
    #fig.savefig(file.replace('_eps.fits','_extract.png'))
    
    if not fit_status:
        return {'fig_extract':fig}
    
    yfull = np.linspace(yarr[0], yarr[-1], 1024)
    gnorm = np.trapz(gau(yfull), yfull)
    gprof = gau(yarr)/gnorm

    num = (new_sci*new_ivar/escl**2*gprof[:,None]).sum(axis=0)
    den = (gprof[:,None]**2*new_ivar/escl**2).sum(axis=0)
    opt_flux = num/den
    opt_err = 1/np.sqrt(den)
    
    if rescale_errors:
        ok = np.isfinite(opt_flux + opt_err) & (den > 0)
        sn = (opt_flux/opt_err)[ok]
        if np.median(sn) < 10:
            df = np.diff(opt_flux[ok])
            de = np.sqrt(opt_err[ok][1:]**2+opt_err[ok][:-1]**2)
            scl = utils.nmad(df/de)
            print(f'Rescale uncertainties: {scl:.3f}')
            opt_err *= scl
        else:
            print(f'Rescale uncertainties (med SN={np.median(sn):.2f})')
            
    #fig, ax = plt.subplots(1,1, figsize=(9, 3))
    fig, axes = plt.subplots(1,2, figsize=(12,3), gridspec_kw={'width_ratios':[3,1]}, sharey=True)
    ax = axes[0]
    axes[1].axis('off')
    
    figs.append(fig)

    #ax.plot(wave/(1+zi), sp[0].data)
    #ax.plot(wave/(1+zi), sp[1].data)

    opt_ivar = 1/opt_err**2
    
    ax.set_ylim(-0.05, 0.1)
    #ax.set_xlim(0.98e4/(1+zi), 1.04e4/(1+zi))
    
    # Lines
    xline = opt_ivar < 0.7*np.median(opt_ivar)
    opt_ivar[xline] *= 0. #0.05
        
    opt_flux *= to_flam
    opt_err *= to_flam
    opt_ivar /= to_flam**2

    ok_idx = np.where(np.isfinite(opt_ivar + opt_flux + to_flam))[0]
    if len(ok_idx) > 2*clip_edge:
        # print('Clip edge', len(ok_idx), len(opt_ivar))
        
        ok_idx = ok_idx[clip_edge:-clip_edge]
    
    opt_mask = np.ones(len(opt_ivar), dtype=bool)
    opt_mask[ok_idx] = False
    
    opt_ivar[opt_mask] = 0
    opt_err[opt_mask] = 1e8
    opt_flux[opt_mask] = 0

    ax.plot(wave/(1+zi), opt_flux, alpha=0.4, color='0.5')
    ax.plot(wave/(1+zi), opt_err, alpha=0.5, color='pink')
    
    bkern = np.ones(binf)
    bnum = nd.convolve1d(opt_flux*opt_ivar, bkern)[binf//2::binf]
    bwnum = nd.convolve1d(wave*opt_ivar, bkern)[binf//2::binf]
    bden = nd.convolve1d(opt_ivar, bkern)[binf//2::binf]
    bflux = bnum/bden
    berr = 1/np.sqrt(bden)
    bwave = bwnum/bden

    ymax = np.percentile(bflux[np.isfinite(bflux)], 90)*5
    ax.set_ylim(-0.5*ymax, ymax)
    #ax.set_ylim(-0.05, 0.11)
    ax.set_xlim(*(bwave[np.isfinite(bwave)][np.array([0,-1])]/(1+zi)))
    xl = ax.get_xlim()
    
    if (zi > 0):
        ax.set_xlim(6500, 6800)

        ax.vlines([3727., 4102.9, 4341.7, 4862., 4960., 5008., 6302, 6563., 6548, 6584, 
                   6679., 6717, 6731, 7137.77, 7321.94, 7332.17], 
                   ax.get_ylim()[0], 0, color='r', linestyle=':')
        ax.set_xlabel(f'rest wave, z={zi:.4f}')
        
    ax.set_xlim(*xl)
    ax.text(0.98, 0.98, targname, ha='right', va='top', 
            transform=ax.transAxes, fontsize=8, bbox={'edgecolor':'None', 'facecolor':'w'})
    
    if (zi > 6):
        ax.vlines([1216.], *ax.get_ylim(), color='r', linestyle=':')

    ax.errorbar(bwave/(1+zi), bflux, berr, color='k', alpha=0.4, linestyle='None', marker='.')
    ax.plot(wave/(1+zi), wave*0., color='k', linestyle=':')
    ok = np.isfinite(bflux+bwave)
    utils.fill_between_steps(bwave[ok]/(1+zi), bflux[ok], bflux[ok]*0., 
                             ax=ax, color='orange', alpha=0.5, zorder=-1)
    
    ax.set_xlim(*band_lims[band])
    
    yt = ax.get_yticks()
    for j in [0, yt[-3]]:
        if j == 0:
            labl = '0'
        else:
            if j < 0.1:
                labl = f'{j:.2f}'
            elif j < 1:
                labl = f'{j:.1f}'
            else:
                labl = f'{int(j)}'
                
        ax.text(0.01*(band_lims[band][1]-band_lims[band][0])+band_lims[band][0], 
                j, labl, ha='left', va='center', 
                fontsize=7, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
    
    ax.set_yticklabels([])
    fig.tight_layout(pad=0.5)

    #fig.savefig(file.replace('_eps.fits','_spec.png'))

    if to_vacuum:
        try:
            wave = airtovac(wave)
            bwave = airtovac(bwave)
        except:
            pass
    
    spec = {'wave':wave, 'opt_flux':opt_flux, 'opt_err':opt_err, 
            'wave_bin':bwave, 'bin_flux':bflux, 'bin_err':berr}
        
    spec['yarr'] = yarr
    spec['prof_model'] = gau
    spec['gprof'] = gau(yarr)
    spec['prof'] = prof
    spec['prof_offset'] = prof_offset
    spec['fig_extract'] = figs[0]
    spec['fig_1d'] = figs[1]
    
    spec['to_flam'] = to_flam
    spec['targname'] = targname
    spec['im'] = im
    spec['file'] = file
    spec['filter'] = band
    spec['xarr'] = xarr
    spec['shape'] = sh
    
    tab = utils.GTable()
    tab['wave'] = wave.astype(np.float32)
    tab['flux'] = opt_flux.astype(np.float32)
    tab['err'] = opt_err.astype(np.float32)
    tab.meta['ny'], tab.meta['nx'] = sh
    tab['ny'] = oky
    
    tab.meta['slx0'] = slx.start, '2D slice start'
    tab.meta['slx1'] = slx.stop, '2D slice stop'
    
    tab['to_flam'] = np.cast[np.float32](to_flam)
    
    tab.meta['itime50'] = exptime[0], 'Median exposure time in 2D'
    tab.meta['itime90'] = exptime[1], '90th percentile 2D exposure time'
    
    tab.meta['wmin'] = tab['wave'][opt_ivar > 0].min(), 'Min valid wavelength'
    tab.meta['wmax'] = tab['wave'][opt_ivar > 0].max(), 'Max valid wavelength'
    snperc = np.nanpercentile((tab['flux']/tab['err'])[opt_ivar > 0], [16, 50, 84, 99])
    tab.meta['sn16'] = snperc[0], 'SN 16th percentile'
    tab.meta['sn50'] = snperc[1], 'SN median'
    tab.meta['sn84'] = snperc[2], 'SN 84th percentile'
    tab.meta['sn99'] = snperc[3], 'SN 99th percentile'
    
    tab.meta['prof_amp'] = spec['prof_model'].amplitude.value, 'Profile model amplitude'
    tab.meta['prof_sig'] = spec['prof_model'].stddev.value, 'Profile model sigma'
    tab.meta['prof_mu'] = spec['prof_model'].mean.value, 'Profile model mean'
    
    ima = np.nanargmax(prof)
    tab.meta['prof_yma'] = yarr[ima], 'Location of profile max'
    tab.meta['prof_ma'] = prof[ima], 'Profile max'
    imi = np.nanargmin(prof)
    tab.meta['prof_ymi'] = yarr[imi], 'Location of profile min'
    tab.meta['prof_mi'] = prof[imi], 'Profile min'
    
    for k in ['prof_offset','file','filter','targname']:
        tab.meta[k] = spec[k]
    
    stats = {}
    cols = ['SKYPA3','AIRMASS','GUIDFWHM']
    tr = {'SKYPA3':'pa','AIRMASS':'airm','GUIDFWHM':'fwhm'}
    
    for k in cols:
        stats[k] = []

    for ki in spec['im'][0].header:
        if '_img' not in ki:
            continue

        ks = ki.split('_img')[0]
        if ks in cols:
            stats[ks].append(spec['im'][0].header[ki])

            
    for k in cols:
        if len(stats[k]) == 0:
            stats[k].append(0)
            
    for k in stats:
        #print(k, len(stats[k]), np.median(stats[k]))
        tab.meta['{0}_min'.format(tr[k])] = np.nanmin(stats[k]), f'Min {k}'
        tab.meta['{0}'.format(tr[k])] = np.nanmedian(stats[k]), f'Median {k}'
        tab.meta['{0}_max'.format(tr[k])] = np.nanmax(stats[k]), f'Max {k}'
    
    full_path = os.path.join(os.getcwd(), file)
    tab.meta['file'] = os.path.basename(full_path), 'Extraction filename'
    tab.meta['path'] = os.path.dirname(full_path), 'File path'
    tab.meta['datemask'] = im[0].header['DATEMASK'], 'Unique mask identifier'
    
    spec['opt_spec'] = tab
    
    return spec


##################
## Find peak 
def find_max(file, gwidth=(5,2), pad=10, erode=10, suffix=''):
    """
    Find peak S/N in 2D spectrum file
    """
    import scipy.ndimage as nd
    
    im = pyfits.open(file)
    valid = np.isfinite(im[0].data)
    if erode:
        valid = nd.binary_erosion(valid, iterations=erode)
        
    if valid.sum() == 0:
        return (-1, (0,0), 0)
    
    if file.endswith('_sp.fits'):
        targname = os.path.basename(file.split('_sp.fits')[0]) + suffix
        is_mine = True
        wht = im['WHT'].data
        valid &= wht > 0
        sig = 1/np.sqrt(wht)
        sig[~valid] = 0
        
    else:
        sig = pyfits.open(file.replace('_eps','_sig'))[0].data
        targname = os.path.basename(file.split('_eps.fits')[0]) + suffix
        is_mine = False

    #sig = pyfits.open(file.replace('_eps','_sig'))

    if os.path.exists(file.replace('_eps','_itime')):
        itime = pyfits.open(file.replace('_eps','_itime'))
        exptime = np.nanpercentile(itime[0].data[itime[0].data > 0], [50,90])
    else:
        exptime = [0,0]

    if gwidth is not None:
        ivar = 1/sig**2
        ivar[~np.isfinite(ivar) | ~valid] = 0

        gau = Gaussian2D(x_mean=0, x_stddev=gwidth[0], y_mean=0, y_stddev=gwidth[1])
        xarr = np.arange(-4*gwidth[0], 4.1*gwidth[0], 1)
        yarr = np.arange(-4*gwidth[1], 4.1*gwidth[1], 1)

        xp, yp = np.meshgrid(xarr, yarr)
        gm = gau(xp, yp)

        sci = im[0].data*1
        sci[~valid] = 0

        num = fftconvolve(sci*ivar, gm, mode='same')
        den = fftconvolve(ivar, gm**2, mode='same')

        smoothed = num/den*valid
        smoothed_sn = smoothed *np.sqrt(den)

        perc = np.percentile(smoothed[valid], [16,50,84])
        limits_2d = perc[1] - 3*np.diff(perc)[0], perc[1] + 3*np.diff(perc)[1]
    else:
        smoothed = im[0].data

    yp, xp = np.indices(ivar.shape)
    msk = (np.abs(yp-(im[0].header['CRPIX2']-im[0].header['CRVAL2'])) <= pad) & (valid)

    #fig, ax = plt.subplots(1,1,figsize=(16,4))
    #ax.imshow(smoothed, vmin=-0.05, vmax=0.05, origin='lower')
    #ax.imshow(smoothed_sn*msk, vmin=-10, vmax=10, origin='lower')
    #print(smoothed_sn[msk].max())

    ix = np.nanargmax(smoothed_sn*msk)
    ij = np.unravel_index(ix, smoothed_sn.shape)
    #ax.set_aspect('auto')
    #xp.shape
    off = ij[0]-(im[0].header['CRPIX2']-im[0].header['CRVAL2'])
    
    return smoothed_sn[ij], ij, off

##############
# Show emission line
def show_line(spec, wave=1.735e4, dv=800, dy=20, vm=(-0.2, 1), vstep=100, yscale=6, aspect=0.8, ax=None):
    """
    Show cutout of line in 2D
    """
    x0 = np.interp(wave, spec['wave'], spec['xarr'])
    y0 = spec['prof_model'].mean.value
    
    xlim = np.interp(wave + np.array([-1,1])*dv/3.e5*wave, spec['wave'], spec['xarr'])
    #print(x0, xlim, y0)
    ylim = y0-dy, y0+dy
    
    sh = spec['im'][0].data.shape
    
    slx = slice(*np.clip(np.cast[int](xlim), 0, sh[1]))
    sly = slice(*np.clip(np.cast[int](ylim), 0, sh[0]))
    
    box = spec['im'][0].data[sly,slx]
    okb = (box != 0) & (np.isfinite(box))
    
    perc = np.percentile(box[okb], [16,50,84])
    limits_2d = perc[1] - yscale/3*np.diff(perc)[0], perc[1] + yscale*np.diff(perc)[1]

    if ax is None:
        newaxis = True
        fig, ax = plt.subplots(1,1, figsize=(4,4))
    else:
        newaxis = False
        
    ax.imshow(spec['im'][0].data, vmin=limits_2d[0], vmax=limits_2d[1])
    ax.set_xlim(*xlim)
    
    #vstep = 100
    vmax = np.floor(dv/vstep)*vstep
    vticks = np.arange(-vmax, vmax+0.01, vstep)
    
    xticks = np.interp(vticks, (spec['wave']-wave)/wave*3.e5, spec['xarr'])
    ax.set_xticks(xticks)
    if newaxis:
        ax.set_xticklabels([int(vt) for vt in vticks])
    else:
        ax.set_xticklabels([])
        
    if newaxis:
        ax.set_xlabel(r'$\Delta v$, km/s' + '\n' + '$\lambda$=' + f'{wave:.1f}')
        ax.text(0.98, 0.98, spec['targname'], ha='right', va='top', 
            transform=ax.transAxes, fontsize=8, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
    else:
        ax.text(0.98, 0.98, r'$\lambda$=' + f'{wave:.1f}', ha='right', va='top', 
            transform=ax.transAxes, fontsize=8, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
        
    ax.grid(linestyle=':')
    
    ax.set_ylim(*ylim)
    
    amax = np.floor(dy*0.18)
    aticks = np.arange(-amax, amax+0.01, 1)
    yarr = np.arange(spec['im'][0].data.shape[0])
    ya = (yarr-y0)*0.18
    
    yticks = np.interp(aticks, ya, yarr)
    ax.set_yticks(yticks)
    if newaxis:
        ax.set_yticklabels(aticks)
        ax.set_ylabel(r'$\Delta y$, arcsec')
    else:
        ax.set_yticklabels([])
        yt = ax.get_yticks()
        xt = ax.get_xticks()
        xli = ax.get_xlim()
        yli = ax.get_ylim()
        for j in [0, -1]:
            ax.text(0.01*(xli[1]-xli[0])+xli[0], yt[j], f'{int(aticks[j])}', ha='left', va='center', 
                fontsize=7, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
            ax.text(xt[j], 0.01*(yli[1]-yli[0])+yli[0], f'{int(vticks[j])}', ha='center', va='bottom', 
                fontsize=7, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
        
    ax.set_aspect(aspect) #'auto')
    
    if newaxis:
        fig.tight_layout(pad=0.5)
        spec['fig_line2d'] = fig
    else:
        fig = plt.gcf()
    
    return fig


#######
# Find Gaussian emission lines
def find_gaussian_lines(sp, binw=50, linew=120, zi=None, left_pad=10, peak_threshold=3.5, ch_order=7, 
                        show_cutout=False, cutout_args=dict(dv=400, dy=20, vm=(-0.2, 1)), 
                        zlines=[6564.6, 5008.24, 3728.4], sens_err=0.5):
    """
    Find emission lines with a gaussian filter
    """
    from grizli import utils
    from mpdaf.obj.spectrum import airtovac, vactoair
    
    import peakutils

    #binw = 50 # km/s
    #logw = utils.log_zgrid([sp['wave'][0], sp['wave'][-1]], binw/3.e5)
    logw = utils.log_zgrid([3000, 2.5e4], binw/3.e5)
    clip = (logw > sp['wave'][0]) & (logw < sp['wave'][-1])
    logw = logw[clip]
    
    band = sp['filter']
    
    # Rebin
    valid = np.isfinite(sp['opt_flux']+sp['opt_err']+sp['to_flam'])
    sp['opt_err'] *= (sp['to_flam']/np.nanmin(sp['to_flam'][valid]))**sens_err
    sp['opt_flux'][~valid] = 0
    sp['opt_err'][~valid] = 1e8

    rebin = interp.rebin_weighted_c(logw, sp['wave'], sp['opt_flux'], sp['opt_err'])
    lw, lf, le = rebin
    
    ok = (lw > 0) & (le < 8) & (le > 0)
    ok &= lf/le > -2
    
    ivar = 1/le**2
    ivar[~ok] = 0
    lf[~ok] = 0

    line_step = linew / binw

    nx = np.ceil(5*line_step)

    xarr = np.arange(-nx,nx+1,1)
    yg = 1/np.sqrt(2*np.pi*line_step**2)*np.exp(-xarr**2/2/line_step**2)
    
    #lfg = lf - nd.gaussian_filter(lf, 25)
    if ch_order >= 0:
        nn = len(ok)
        xch = np.linspace(0,1,nn)
        perc = np.percentile((lf/le)[ok], [10, 90])
        ok2 = ok & (lf/le > perc[0]) & (lf/le < perc[1])
        
        ch = np.polynomial.chebyshev.chebfit(xch[ok2], lf[ok2], ch_order, w=1/le[ok2])
        lf_ch = np.polynomial.chebyshev.chebval(xch, ch)
        lfg = lf - lf_ch
    else:
        lfg = lf - np.median(lf[ok])
        lf_ch = None
        
    num = nd.convolve1d(lfg*ivar, yg)
    den = nd.convolve1d(ivar, yg**2)

    line_err = np.sqrt((1/den))
    line_flux = num/den
    line_sn = (num/den/line_err)
    
    #print('xx', logw.shape, lw.shape)
    
    log_spec = utils.GTable()
    log_spec['wave'] = (logw[:-1]+np.diff(logw)/2.).astype(np.float32)
    log_spec['awave'] = lw.astype(np.float32)
    log_spec['flux'] = lf.astype(np.float32)
    log_spec['err'] = le.astype(np.float32)
    log_spec['line_flux'] = line_flux.astype(np.float32)
    log_spec['line_err'] = line_err.astype(np.float32)
    log_spec.meta['binw'] = binw
    log_spec.meta['linew'] = linew
    log_spec.meta['z'] = zi
    
    for k in sp['opt_spec'].meta:
        log_spec.meta[k] = sp['opt_spec'].meta[k]
        
    fig, axes = plt.subplots(1,2,figsize=(12,3), gridspec_kw={'width_ratios':[3,1]})
    ax = axes[0]
    
    _ = ax.scatter(lw[ok], lf[ok]/le[ok], marker='.', color='k', alpha=0.2)

    ax.plot(lw[ok], line_sn[ok] - np.median(line_sn[ok]))
    ax.plot(lw[ok], line_err[ok], color='pink', alpha=0.5)
    #plt.plot(lw[ok], err[ok])
    
    if lf_ch is not None:
        ax.plot(lw[ok], lf_ch[ok]/np.median(le[ok]), color='orange')
        
    ax.set_ylim(-2,8)
    #ax = plt.gca()
    ax.grid(linestyle=':')

    ax.text(0.98, 0.98, sp['targname'], ha='right', va='top', 
            transform=ax.transAxes, fontsize=8)
    
    ax.set_xlim(*band_lims[sp['filter']])
    xl = ax.get_xlim()
    
    lwv = lw #airtovac(lw) # already vacuu
    
    #line_sn[~ok] = 0
    peak_sn = line_sn*1.
    peak_sn[~ok] = 0
    peak_sn[np.where(ok)[0][:left_pad]] = 0
    peak_sn[np.where(ok)[0][-left_pad:]] = 0
    
    lcolors = ['r','purple','orange','g','b','magenta']
    
    all_lines = np.array([1215.6701, 1240.81, 1549.48, 1640.42, 1661.241, 1666.15, 1906.68, 1908.73, 
                          3727.1, 3869.86, 3968.59, 4102.892, 4341.692, 4862.71, 4960.295, 5008.24, 
                          #5193.27, 5201.705, 5519.242, 5539.411, # ArIII, NI, ClIII, ClIII
                          6302.046, 6365.535, 6564.6, 6549.86, 6585.27, 6718.294, 6732.673, 7137.77,
                          9071.1, 9533.2, 10833., 12821.57])
    
    perc = np.nanpercentile(peak_sn[peak_sn != 0], [50, 84])
    thresh = perc[0] + peak_threshold*(perc[1]-perc[0])
    #print('xx', peak_threshold)
    
    ax.hlines(thresh, *ax.get_xlim(), color='r', linestyle=':')
    
    indexes = peakutils.indexes(np.maximum(peak_sn, -1), thres=thresh, min_dist=4, thres_abs=True)
    
    log_spec.meta['pthresh'] = peak_threshold
    log_spec.meta['lthresh'] = thresh
    sp['line_tab'] = None
    
    if (len(indexes) > 0) & (len(indexes) < 10):
        textkw = dict(fontsize=7, ha='center', va='center', bbox={'facecolor':'w', 'alpha':0.8, 'edgecolor':'None'})
        
        lso = indexes[np.argsort(line_sn[indexes])[::-1]]
        log_spec.meta['nline'] = len(indexes)
        
        lwx = peakutils.interpolate(lw, line_flux, 
                                    ind=indexes, width=10)

        line_tab = utils.GTable()
        line_tab['wave'] = lw[indexes]/1.e4
        sp['line_tab'] = line_tab
        
        for i, ii in enumerate(lso):
            
            log_spec.meta[f'linew{i:02d}'] = lwv[ii]
            log_spec.meta[f'linesn{i:02d}'] = line_sn[ii]
            log_spec.meta[f'linef{i:02d}'] = line_flux[ii]
            log_spec.meta[f'linee{i:02d}'] = line_err[ii]
            
            ax.text(lw[ii], (0.85 - 0.07*(i % 2 == 1))*ax.get_ylim()[1], f'{lwv[ii]:.1f}', **textkw)
            
            if i == 0:
                print(i, lwv[ii], lwv[ii]/np.array(zlines)-1)
                
                for j, ll in enumerate(zlines):
                    zl = lwv[ii]/ll-1
                    ax.text(lw[ii], (0.85 - 0.07*(j+2))*ax.get_ylim()[1], 
                        f'{zl:.4f}', color=lcolors[j], **textkw)
                    
                    ax.vlines((all_lines)*(1+zl), 
                                        (0.85 - 0.07*(j+1.5))*ax.get_ylim()[1], 
                                        (0.85 - 0.07*(j+2.5))*ax.get_ylim()[1], color=lcolors[j], linewidth=0.5)
                
                if show_cutout:
                    #print(line_flux[ii], yg.max())
                    lfig = show_line(sp, wave=lw[ii], yscale=np.maximum(line_sn[ii]*yg.max()*2, 5), 
                                     ax=axes[1], **cutout_args)    
                    figs = (fig, lfig)
                    lfig.tight_layout(pad=0.5)

                else:
                    figs = fig
                    axes[1].axis('off')
                    
                #ax.text(lw[ii], (0.95 - 0.07*2)*ax.get_ylim()[1], 
                #    f'{lwv[ii]/5008.24-1:.4f}', **textkw)
                

    else:
        log_spec.meta['nline'] = len(indexes)
        axes[1].axis('off')
        
        figs = (fig)
    
    sp['log_spec'] = log_spec
    sp['fig_findline'] = fig
    
    if zi is not None:
        ax.vlines((all_lines)*(1+zi), 
                   ax.get_ylim()[0], 0, color='r', linestyle=':')

        ax.set_xlim(*xl)

    yt = ax.get_yticks()
    for j in [0, yt[-2]]:
        ax.text(0.01*(band_lims[band][1]-band_lims[band][0])+band_lims[band][0], 
                j, f'{int(j)}', ha='left', va='center', 
                fontsize=7, bbox=dict(edgecolor='None', facecolor='w', alpha=0.9))
    
    ax.set_yticklabels([])

    fig.tight_layout(pad=0.5)

    return log_spec, figs


######
# Save outputs
def save_data(sp, path=None, suffix=''):
    if 'opt_spec' not in sp:
        return []
    
    files = []
    
    if path is None:
        outdir = sp['opt_spec'].meta['path']
        if len(outdir) == 2:
            outdir = outdir[0]
    else:
        outdir = path
    
    outroot = sp['opt_spec'].meta['targname'] + suffix
    #print(outdir, outroot)
    
    if 'opt_spec' in sp:
        files.append(f'{outroot}_1d.fits')
        sp['opt_spec'].write(f'{outdir}/{outroot}_1d.fits', overwrite=True)
    
    if 'log_spec' in sp:
        files.append(f'{outroot}_log1d.fits')
        sp['log_spec'].write(f'{outdir}/{outroot}_log1d.fits', overwrite=True)

    if 'fig_extract' in sp:
        files.append(f'{outroot}_2d.png')
        sp['fig_extract'].savefig(f'{outdir}/{outroot}_2d.png')
    
    if 'fig_1d' in sp:
        files.append(f'{outroot}_1d.png')
        sp['fig_1d'].savefig(f'{outdir}/{outroot}_1d.png')

    if 'fig_findline' in sp:
        files.append(f'{outroot}_line1d.png')
        sp['fig_findline'].savefig(f'{outdir}/{outroot}_log1d.png')
    
    if 'fig_line2d' in sp:
        files.append(f'{outroot}_line2d.png')
        sp['fig_line2d'].savefig(f'{outdir}/{outroot}_line2d.png')
        
    return files


def read_telluric_correction(basename='MOSFIRE_telluric_correction_long2pos_20170509_{b}_HIP105437.fits'):
    """
    Read telluric files and initialize interpolation objects
    """
    from scipy.interpolate import interp1d
    
    telluric_path = os.path.join(os.path.dirname(__file__), 'data')
    
    flux_corr_interp = {}
    for b in 'YJHK':
        telluric_file = basename.format(b=b)
        hip_corr = Table.read(f'{telluric_path}/{telluric_file}')
        flux_corr_interp[b] = interp1d(hip_corr['wave'], 
                                       hip_corr['flux_corr'], 
                                       bounds_error=False, fill_value=1e10)
    
    return flux_corr_interp
    
flux_corr_interp = read_telluric_correction()

def runit(file, zi=0, save=True):
    """
    Run the full extraction
    """
    full_file = os.path.join(os.getcwd(), file)
    #band = full_file.split('/')[-2]
    im = pyfits.open(file)
    band = im[0].header['FILTER']
    im.close()
    
    #fmax = find_max(file)
    try:
        fmax = find_max(file, gwidth=(5,2), erode=10)
    except:
        fmax = (-1, (0,0), 0)
    
    try:
        redo = False
        
        if fmax[0] > 5:
            try:
                slx = slice(np.maximum(fmax[1][1]-20, 0), fmax[1][1]+20)
                test = optimal_extract(file, 
                                       slx=slx,
                                       fit_profile=1, escl=1, zi=zi, binf=2,
                                       prof_sigma=1.8, prof_offset=fmax[-1],
                                       show_sn=True,
                                       flux_corr=flux_corr_interp[band], 
                                       gwidth=(5, 1), clip_edge=2)
                redo = False
            except ValueError:
                plt.close()
                redo = True
        
        if (fmax[0] < 5) | redo:
            test = optimal_extract(file, slx=None, fit_profile=1, 
                         escl=1, zi=zi, binf=2, prof_sigma=1.8, prof_offset=0, 
                         flux_corr=flux_corr_interp[band], gwidth=(5, 1), clip_edge=2)
   
    except ValueError:
        print('Failed')
        return test
        
    if test:
        cutout_args = dict(vstep=400, dv=1198, dy=20,
                           vm=(-0.2, 1), aspect='auto')
                           
        if 'wave' in test:
            find_gaussian_lines(test, peak_threshold=3.5, 
                                zlines=[6564.6, 5008.24, 3728.4],
                                left_pad=1, 
                                show_cutout=True, cutout_args=cutout_args, 
                                sens_err=0.5)
        else:
            return test
            
        if 'opt_spec' not in test:
            return test
            
        meta = test['opt_spec'].meta
        key = '{0}/{1}'.format(meta['datemask'][0], meta['targname'])
        #spectra[key] = test
        test['key'] = key
        
        if save:
            xfiles = save_data(test, suffix='_sp')
            #plt.close('all')
    
    return test