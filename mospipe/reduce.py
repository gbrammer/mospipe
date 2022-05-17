"""
Extract slits from a MOSFIRE slitmask
"""
import glob
import os
import traceback

import astropy.io.fits as pyfits
import numpy as np
from grizli import utils
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import drizzlepac

import scipy.ndimage as nd
import peakutils
from skimage.feature import match_template
from skimage.registration import phase_cross_correlation
from tqdm import tqdm

utils.LOGFILE = 'mospipe.log'

utils.set_warnings()

def grating_dlambda(band):
    """
    returns the dlambda/dpixel in angstrom for a band
    
    (From MosfireDRP)
    """
    orders = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = orders[band]
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length
    dlambda = scale * d / order * 10000

    return dlambda

# grating_summary = {'Y': {'edge':[9612, 11350]}, 
#                    'J': {'edge':[11550, 13623]}, 
#                    'H': {'edge':[14590, 18142]},
#                    'K': {'edge':[19118, 24071]}}
grating_summary = {'Y': {'edge':[9612, 11350]}, 
                  'J': {'edge':[11450, 13550]}, 
                  'H': {'edge':[14590, 18142]},
                  'K': {'edge':[18900, 24150]}}
                  
for k in grating_summary:
    edge = grating_summary[k]['edge']
    grating_summary[k]['dlam'] = dlam = grating_dlambda(k)
    grating_summary[k]['N'] = int(np.ceil(edge[1]-edge[0])/dlam)
    grating_summary[k]['ref_wave'] = (edge[1]+edge[0])/2


def get_grating_loglam(filter):
    """
    Get polynomial and WCS coefficients that approximate logarithmic 
    wavelength spacing
    """
        
    gr = grating_summary[filter]
    edge, dlam, N = gr['edge'], gr['dlam'], gr['N']

    loglam = np.logspace(np.log10(edge[0]), np.log10(edge[1]), N)

    xarr = np.arange(N)
    xref = N//2-1

    # Polynomial representation of logarithmic wavelength
    c = np.polyfit(xarr-xref, loglam, 3)

    # WAVE-LOG WCS
    w = np.log(loglam/loglam[xref])*loglam[xref]
    #plt.plot(mask.xarr-1024, w)

    #plt.plot(mask.xarr, w)

    cl = np.polyfit(xarr-xref, w, 1)
    
    #plt.plot(loglam, (np.polyval(c, xarr-xref) - loglam)/loglam)

    #print(N, xref, c, cl)
    return N, loglam[xref], c, cl


def fit_wavelength_from_sky(sky_row, band, order=3, make_figure=True, nsplit=5, plot_axis=None, debug=False, use_oliva=True, **kwargs):
    
    lines = {}
    line_intens = {}
    
    xarr = np.arange(2048)
    
    # OH Sky lines from MosfireDRP
    lines['Y'] = np.array([9793.6294, 9874.84889, 9897.54143, 9917.43821,
                           10015.6207, 10028.0978, # 10046.7027, 10085.1622
                           10106.4478, 10126.8684, 10174.623, 10192.4683,
                           10213.6107, 10289.3707, 10298.7496, 10312.3406,
                           10350.3153, 10375.6394, 10399.0957, 10421.1394,
                           10453.2888, 10471.829, 10512.1022, 10527.7948,
                           10575.5123, 10588.6942, 10731.6768, 10753.9758,
                           10774.9474, 10834.1592, 10844.6328, 10859.5264,
                           10898.7224, 10926.3765, 10951.2749, 10975.3784,
                           11029.8517, 11072.4773, 11090.083, 11140.9467,
                           11156.0366 ])

    lines['J'] = np.array([11538.7582, 11591.7013, 11627.8446, 11650.7735, 
                           11696.3379, 11716.2294, 11788.0779, 11866.4924,
                           11988.5382, 12007.0419, 12030.7863, 12122.4957,
                           12135.8356, 12154.9582, 12196.3557, 12229.2777,
                           12257.7632, 12286.964, 12325.9549, 12351.5321,
                           12400.8893, 12423.349, 12482.8503, 12502.43,
                           12905.5773, 12921.1364, 12943.1311, 12985.5595,
                           13021.6447, 13052.818, 13085.2604, 13127.8037,
                           13156.9911, 13210.6977, 13236.5414, 13301.9624,
                           13324.3509, 13421.579])

    lines['H'] = np.array([14605.0225, 14664.9975, 14698.7767, 14740.3346, 
                           14783.7537, 14833.029, 14864.3219, 14887.5334,
                           14931.8767, 15055.3754, 15088.2599, 15187.1554,
                           15240.922, 15287.7652, 15332.3843, 15395.3014,
                           15432.1242, 15570.0593, 15597.6252, 15631.4697,
                           15655.3049, 15702.5101, 15833.0432, 15848.0556,
                           15869.3672, 15972.6151, 16030.8077, 16079.6529,
                           16128.6053, 16194.6497, 16235.3623, 16317.0572,
                           16351.2684, 16388.4977, 16442.2868, 16477.849,
                           16502.395, 16553.6288, 16610.807, 16692.2366,
                           16708.8296, 16732.6568, 16840.538, 16903.7002,
                           16955.0726, 17008.6989, 17078.3519, 17123.5694,
                           17210.579, 17248.5646, 17282.8514, 17330.8089,
                           17386.0403, 17427.0418, 17449.9205, 17505.7497,
                           17653.0464, 17671.843, 17698.7879, 17811.3826,
                           17880.341, 17993.9600, 18067.9500 ])

    lines['K'] = np.array([19518.4784, 19593.2626, 19618.5719, 19642.4493, 
                           19678.046, 19701.6455, 19771.9063, 19839.7764,
                           20008.0235, 20193.1799, 20275.9409, 20339.697,
                           20412.7192, 20499.237, 20563.6072, 20729.032,
                           20860.2122, 20909.5976, 21176.5323, 21249.5368,
                           21279.1406, 21507.1875, 21537.4185, 21580.5093,
                           21711.1235, 21802.2757, 21873.507, 21955.6857,
                           22125.4484, 22312.8204, 22460.4183, 22517.9267,
                           22690.1765, 22742.1907, 22985.9156, 23914.55,
                           24041.62])

    for b in lines:
        line_intens[b] = np.ones_like(lines[b])

    line_intens['K'] = np.array([0.05, 0.1, 0.1, 0.25, 0.1,
                                 0.35, 0.4, 0.15, 0.7, 0.1, 0.7, 0.35, 1,
                                 0.25, 0.65, 0.45, 0.2, 0.25, 0.25, 0.3, 0.05,
                                 0.75, 0.3, 0.1, 0.2, 0.6, 0.25, 0.75, 0.5,
                                 0.3, 0.05, 0.1, 0.05, 0.05, 0.05, 0.15,
                                 0.05])
    
    # Oliva et al. sky lines
    # https://ui.adsabs.harvard.edu/abs/2015A%26A...581A..47O
    if use_oliva:
        data_dir = os.path.dirname(__file__)
        oh_file = os.path.join(data_dir, 'data', 'oliva_oh_lines.vot')
        oh = utils.read_catalog(oh_file)
        for b in lines:
            lines[b] = (oh['lambda1'] + oh['lambda2'])/2
            line_intens[b] = oh['Flux']/1000.
        
    dlam = grating_dlambda(band)
    msg = f'Band={band}, dlambda={dlam:.3f} A/pix'
    utils.log_comment(utils.LOGFILE, msg, verbose=True)
    
    band_lims = {'Y': (9701, 11280), 
                 'J': (11501, 13650),
                 'H': (14510, 17700), 
                 'K': (19100, 23990)}
    
    for k in band_lims:
        ok = (lines[k] > band_lims[k][0]) & (lines[k] < band_lims[k][1])
        lines[k] = lines[k][ok]
        line_intens[k] = line_intens[k][ok]
        
    if band == 'K':
        x0 = 19577.35927
    elif band == 'J':
        x0 = 11536.0
    elif band == 'H':
        x0 = 14600.
    elif band == 'Y':
        x0 = 9700.
    
    ###############
    # First guess linear from simple cross correlation
    for it in range(3):
        wave = xarr*dlam+x0
        yline = wave*0.
        for il, l in enumerate(lines[band]):
            yline += line_intens[band][il]*np.exp(-(wave-l)**2/2/dlam**2)

        shift, error, diffphase = phase_cross_correlation(sky_row[None,:], yline[None,:], upsample_factor=100, reference_mask=~np.isnan(sky_row[None,:]), moving_mask=~np.isnan(yline[None,:]), return_error=True)
        msg = f'  phase shift {it} {shift[1]:.2f}  x0={x0-shift[1]*dlam:.3f}'
        utils.log_comment(utils.LOGFILE, msg, verbose=True)
        
        x0 -= shift[1]*dlam

    wave = xarr*dlam+x0
    
    if debug:
        # Template match cross correlation
        xres = np.squeeze(match_template(sky_row[None,:], yline[None,:],
                                         pad_input=True))
        fig, ax = plt.subplots(1,1)
        ax.plot(xarr, xres)

    indexes = peakutils.indexes(yline, thres=0.08, min_dist=12)
    base = peakutils.baseline(sky_row/sky_row.max(), 4)
    
    # Peak centers
    peaks_x = peakutils.interpolate(wave, sky_row/sky_row.max()-base, 
                                    ind=indexes, width=10)
    peaks_pix = peakutils.interpolate(xarr, sky_row/sky_row.max()-base, 
                                      ind=indexes, width=10)
                                      
    peaks_model = peakutils.interpolate(wave, yline, ind=indexes, width=10)

    peaks_y = np.interp(peaks_x, wave, sky_row/sky_row.max())
    
    if debug:
        fig, axes = plt.subplots(nsplit,1,figsize=(12,12*nsplit/5))
        for i, ax in enumerate(axes):
            #ax.scatter(wave[indexes], yline[indexes], marker='x', color='r')
            ax.scatter(peaks_x, peaks_y, marker='x', color='r')

            ax.plot(wave, sky_row/sky_row.max())
            ax.plot(wave, yline*0.5)
            ax.set_xlim(wave[i*2048//nsplit], wave[(i+1)*2048//nsplit-1])
    
    #######
    # Fit polynomial dispersion
    lam = peaks_x*1

    for it in range(3):
        ok = (np.abs(peaks_x-peaks_model) < 10) 
        ok &= (peaks_model > 9000) & (peaks_x > 9000)
        ok &= (peaks_model < 2.3e4)
        
        if it > 0:
            ok &= np.abs(lam-peaks_model) < 0.5

        lam_coeffs = np.polyfit((peaks_pix[ok]-1023), peaks_model[ok], 
                                order, w=yline[indexes][ok])

        lam = np.polyval(lam_coeffs, peaks_pix-1023)
    
    wave_fit = np.polyval(lam_coeffs, xarr-1023)
    
    if plot_axis is not None:
        ax = plot_axis
        make_figure = False
        
        ax.scatter(peaks_model[ok]/1.e4, (peaks_x-peaks_model)[ok], 
                   label='Linear')
        ax.plot(peaks_model[ok]/1.e4, (lam-peaks_model)[ok],
                label=f'poly deg={order}')
        ax.grid()
        ax.set_ylim(-3*dlam, 3*dlam)
        ax.set_ylabel(r'$\Delta\lambda$ [$\mathrm{\AA}$]')
        ax.legend(ncol=2, fontsize=8)

    if make_figure | debug:
        
        fig1, ax = plt.subplots(1,1,figsize=(6,6))
        ax.scatter(peaks_model[ok]/1.e4, (peaks_x-peaks_model)[ok],
                   label='Linear')
        ax.plot(peaks_model[ok]/1.e4, (lam-peaks_model)[ok],
                label=f'poly deg={order}')
        ax.grid()
        ax.set_ylim(-3, 3)
        ax.set_ylabel(r'$\Delta\lambda$ [A]')
        ax.legend(ncol=2, fontsize=8)
        fig1.tight_layout(pad=0.1)
        
        ynew = wave_fit*0.
        for il, l in enumerate(lines[band]):
            ynew += line_intens[band][il]*np.exp(-(wave_fit-l)**2/2/dlam**2)

        fig, axes = plt.subplots(nsplit,1,figsize=(12,9*nsplit/5))
        for i, ax in enumerate(axes):
            #ax.scatter(wave[indexes], yline[indexes], marker='x', color='r')
            #ax.scatter(peaks_x, peaks_y, marker='x', color='r')

            ax.plot(wave_fit, sky_row/sky_row.max())
            ax.plot(wave, sky_row/sky_row.max(), alpha=0.4)
            ax.plot(wave_fit, ynew/2)
            ax.set_xlim(wave_fit[i*2048//nsplit],
                        wave_fit[(i+1)*2048//nsplit-1])
            ax.set_yticklabels([])
            ax.grid()
        
        fig.tight_layout(pad=0.3)
        
        figs = (fig1, fig)
    else:
        figs = None
        
    return lam_coeffs, wave_fit, figs


# Header keys to pull out for a single exposure explaining the mask
OBJ_KEYS = ['OBJECT', 'FRAMDESC', 'MASKNAME','OBSERVER','PATTERN',
            'DTHCOORD', 'SKYPA3','PROGID','PROGPI','PROGTL1','SEMESTER', 
            'WAVERED','WAVEBLUE']

# Heaader keys for the exposure sequence
SEQ_KEYS = ['AIRMASS','GUIDFWHM','MJD-OBS']

utils.set_warnings()

def show_ls_targets(path):
    """
    Show Long2pos targets in a specified path
    """
    offset_files = glob.glob(os.path.join(path, 'Offset*txt'))
    ls_targets = []
    for file in offset_files:
        if '_Pos' in file:
            ls_targets.append(file.split('_')[-2])
    
    if len(ls_targets) == 0:
        print(f'No LS targets found in {path}')
        return ['']
    
    ls_targets = np.unique(ls_targets).tolist()
    for i, t in enumerate(ls_targets):
        print(f'{i}: {t}')
    
    return ls_targets

class MosfireMask(object):
    """
    A group os MOSFIRE exposures for a single mask / night / filter
    """
    def __init__(self, path='mask12_13_new_20200225/Reduced/mask12_13_new/2020feb26/H', min_nexp=3, ls_target='', use_pixel_flat=True):
        
        if path.endswith('/'):
            self.path = path[:-1]
        else:
            self.path = path
            
        self.filter = self.path.strip('/').split('/')[-1]
        
        self.datemask = os.path.join(os.getcwd(), self.path).split('/')[-5]

        logfile = f'{self.datemask}-{self.filter}.log'
        
        utils.LOGFILE = os.path.join(os.getcwd(), logfile)
        self.logfile = utils.LOGFILE

        utils.log_comment(self.logfile, f'Start mask {self.path}', 
                          verbose=True, show_date=True)
        
        flat_file = glob.glob(os.path.join(self.path, 'combflat_2d*fits'))[0]
        self.flat_file = flat_file
        self.flat = pyfits.open(flat_file)
        
        pixel_flat_file = flat_file.replace('combflat','pixelflat')
        
        if use_pixel_flat & os.path.exists(pixel_flat_file):
            self.pixel_flat = pyfits.open(pixel_flat_file)[0].data
        else:
            self.pixel_flat = 1
            
        self.offset_files = glob.glob(os.path.join(self.path, f'Offset*{ls_target}*txt'))

        self.groups = {}
        self.read_exposures(min_nexp=min_nexp)
        
        #self.target_names = []
        #self.target_keys = []
        self.read_ssl_table()
    
        self.nslits = 0 #len(self.ssl)
        
        # Info
        self.info(log=True)
        
        self.slit_hdus = {}
        self.slit_info = {}
        self.plan_pairs = []
        
    @property 
    def keys(self):
        """
        Keys of the exposure groups, like 'A','B','Ap','Bp', etc.
        """
        return list(self.groups.keys())
    
    @property 
    def namestr(self):
        """
        Descriptive name
        """
        return f'{self.datemask} {self.filter}'
        
    def __repr__(self):
        return self.namestr
    
    @property
    def meta(self):
        """
        Metadata dictionary
        """
        meta = {}
        for i, gr in enumerate(self.groups):
            grm = self.groups[gr].meta
            if i == 0:
                for k in OBJ_KEYS:
                    meta[k] = grm[k]
                for k in SEQ_KEYS:
                    meta[k] = grm[k]
            else:
                for k in SEQ_KEYS:
                    meta[k].extend(grm[k])
        
        for k in SEQ_KEYS:
            try:
                meta[f'{k}_MIN'] = np.nanmin(meta[k])
                meta[f'{k}_MED'] = np.nanmedian(meta[k])
                meta[f'{k}_MAX'] = np.nanmax(meta[k])
            except:
                meta[f'{k}_MIN'] = 0.
                meta[f'{k}_MED'] = 0.
                meta[f'{k}_MAX'] = 0.
                
        return meta
    
    
    @property
    def plans(self):
        """
        Dither groups ['A','B'], ['Ap','Bp'], etc.
        """
        keys = self.keys
        plans = []
        if ('A' in keys) & ('B' in keys):
            plans.append(['A','B'])
        
        for i in range(5):
            pp = 'p'*(i+1)
            if (f'A{pp}' in keys) & (f'B{pp}' in keys):
                plans.append([f'A{pp}',f'B{pp}'])
        
        return plans
    
    
    def get_plan_pairs(self, tpad=90, show=False):
        """
        Paired exposures in a given offset "plan"
        """
        plan_pairs = []
        if len(self.plans) == 0:
            msg = f'{self.namestr} # No plans found!'
            utils.log_comment(self.logfile, msg, verbose=True)
            return [], None

        for ip, plan in enumerate(self.plans):

            pa, pb = plan

            ta = np.array(self.groups[pa].meta['MJD-OBS'])
            if ip == 0:
                t0 = ta[0]

            ta = (ta - t0)*86400
            tb = np.array(self.groups[pb].meta['MJD-OBS'])
            tb = (tb - t0)*86400

            ia = []
            ib = []
            npl = len(self.plans)
            for i, a in enumerate(ta):
                dt = np.abs(tb - a)
                tstep = (self.groups[pa].truitime[i]+tpad)*npl
                ok = np.where(dt < tstep)[0]
                if len(ok) > 0:
                    for j in ok:
                        if (j not in ib) & (i not in ia):
                            ia.append(i)
                            ib.append(j)

            pd = {'plan':plan, 'ia':ia, 'ib':ib,
                  'ta':ta, 'tb':tb, 't0':t0,
                  'n':len(ia), 
                  'fwhm':np.ones(len(ia)), 'shift':np.zeros(len(ia)), 
                  'scale':np.ones(len(ia))}

            plan_pairs.append(pd)

        self.plan_pairs = plan_pairs

        if show & (len(plan_pairs) > 0):
            fig, ax = plt.subplots(1,1,figsize=(12,3))
            #ax.plot(ta[ia], tb[ib]-ta[ia], marker='o', label='B - A')

            for i, pd in enumerate(plan_pairs):
                pa, pb = pd['plan']
                ta, tb, ia, ib = pd['ta'], pd['tb'], pd['ia'], pd['ib']

                p = 6
                y0 = 2*p*i
                ax.scatter(ta[ia], self.groups[pa].truitime[ia] + y0, 
                           color='r', zorder=100)

                for j in range(len(ia)):
                    ax.plot([ta[ia[j]], tb[ib[j]]],
                            np.ones(2)*self.groups[pa].truitime[ia[j]] + y0, 
                            marker='o', color='0.5')

                ax.vlines(ta, self.groups[pa].truitime-p+y0,
                          self.groups[pa].truitime-2+y0, color='r', alpha=0.3,
                          linestyle='--')
                ax.vlines(ta[ia], self.groups[pa].truitime[ia]-p+y0,
                          self.groups[pa].truitime[ia]-2+y0, color='r',
                          alpha=0.9)
                ax.vlines(tb, self.groups[pb].truitime+p+y0,
                          self.groups[pb].truitime+2+y0, color='0.5',
                          alpha=0.3, linestyle='--')
                ax.vlines(tb[ib], self.groups[pb].truitime[ib]+p+y0,
                          self.groups[pb].truitime[ib]+2+y0, color='0.5',
                          alpha=0.9)

            ax.set_xlabel(r'$\Delta t$, $\mathrm{MJD}_0 = $'+ '{0:.2f}   {1}'.format(pd['t0'], self.namestr))

            xl = ax.get_xlim()
            xlab = 0.05*(xl[1]-xl[0])
            for i, pd in enumerate(plan_pairs):
                pa, pb = pd['plan']
                ta, tb, ia, ib = pd['ta'], pd['tb'], pd['ia'], pd['ib']
                y0 = 2*p*i
                yi = np.interp(self.groups[pb].truitime[ib][0] + y0, ax.get_ylim(), [0,1])
                ax.text(0.01, yi, 
                        f"{pa} - {pb}  {pd['n']}", rotation=90, 
                        ha='left', va='center', transform=ax.transAxes)

            ax.set_yticklabels([])
            ax.set_yticks(ax.get_ylim())

            fig.tight_layout(pad=0.1)

        else:
            fig = None


        return plan_pairs, fig


    def plan_pairs_info(self):
        """
        Print a summary of the plain_pairs data (shifts, fwhm, etc)
        """
        for pd in self.plan_pairs:
            #print(pd)
            row = '# fileA             fileB                    dt    fwhm    shift    scale\n'
            row += '# {0}\n'.format(self.namestr)
            row += '# plan: {0} {1}\n'.format(*pd['plan'])

            pa, pb = pd['plan']
            gra = self.groups[pa]
            grb = self.groups[pb]

            for i, (ia, ib) in enumerate(zip(pd['ia'], pd['ib'])):
                row += '{0}   {1}   {2:>7.1f}  {3:6.2f}  {4:6.2f}   {5:6.2f}\n'.format(gra.files[ia], grb.files[ib], 
                                                                              pd['ta'][i], pd['fwhm'][i], 
                                                                              pd['shift'][i], pd['scale'][i])
            utils.log_comment(self.logfile, row, verbose=True)


    @property
    def exptime(self):
        """
        Total exposure time across groups
        """
        return np.sum([self.groups[k].truitime.sum() for k in self.groups])
    
    
    @property
    def nexp(self):
        """
        Total number of exposure across groups
        """
        return np.sum([self.groups[k].N for k in self.groups])

    
    def read_ssl_table(self):
        """
        Read the attached Science_Slit_List table
        """
        from astropy.coordinates import SkyCoord

        img = self.groups[self.keys[0]].img[0]
        #img.info()
        ssl = utils.GTable.read(img['Science_Slit_List'])
        #msl = utils.GTable.read(img['Mechanical_Slit_List'])
        
        valid_slits = np.where([t.strip() != '' for t in ssl['Target_Name']])[0]
        self.ssl = ssl = ssl[valid_slits][::-1]
        
        # Get mag from Target_List
        tsl = utils.GTable.read(img['Target_List'])
        
        tsl_names = [n.strip() for n in tsl['Target_Name']]
        ssl_names = [n.strip() for n in ssl['Target_Name']]
        
        self.ssl['Magnitude'] = -1.
        self.ssl['target_ra'] = -1.
        self.ssl['target_dec'] = -1.
        
        for i, n in enumerate(ssl_names):
            if n in tsl_names:
                ti = tsl_names.index(n)
                ras = '{0}:{1}:{2}'.format(tsl['RA_Hours'][ti],
                                            tsl['RA_Minutes'][ti],
                                            tsl['RA_Seconds'][ti])
                des = '{0}:{1}:{2}'.format(tsl['Dec_Degrees'][ti],
                                            tsl['Dec_Minutes'][ti],
                                            tsl['Dec_Seconds'][ti])
                
                target_rd = SkyCoord(ras, des, unit=('hour','deg'))
                self.ssl['target_ra'][i] = target_rd.ra.value
                self.ssl['target_dec'][i] = target_rd.dec.value
                if 'Magnitude' in tsl.colnames:
                    try:
                        self.ssl['Magnitude'][i] = float(tsl['Magnitude'][ti])
                    except ValueError:
                        self.ssl['Magnitude'][i] = 99.
                    
        # Coords
        ras = []
        des = []
        for i in range(len(ssl)):
            ras.append('{0}:{1}:{2}'.format(ssl['Slit_RA_Hours'][i],
                                            ssl['Slit_RA_Minutes'][i],
                                            ssl['Slit_RA_Seconds'][i]))
            des.append('{0}:{1}:{2}'.format(ssl['Slit_Dec_Degrees'][i],
                                            ssl['Slit_Dec_Minutes'][i],
                                            ssl['Slit_Dec_Seconds'][i]))

        slit_rd = SkyCoord(ras, des, unit=('hour','deg'))
        self.ssl['slit_ra'] = slit_rd.ra.value
        self.ssl['slit_dec'] = slit_rd.dec.value
        

    @property 
    def ssl_stop(self):
        sl = np.cast[float](self.ssl['Slit_length'])
        ssl_stop = np.cumsum(sl/0.1799+5.35)-9
        return np.minimum(ssl_stop, 2045)


    @property
    def ssl_start(self):
        sl = np.cast[float](self.ssl['Slit_length'])
        return np.maximum(self.ssl_stop - sl/0.1799, 4)


    @property 
    def target_names(self):
        target_names = [t.strip().replace('/','-').replace(' ','_') 
                        for t in self.ssl['Target_Name']]
        return target_names


    @property 
    def target_slit_numbers(self):
        slit_numbers = [int(n.strip()) for n in self.ssl['Slit_Number']]
        return slit_numbers


    @property 
    def target_keys(self):
        target_keys = [f'{self.datemask}-{self.filter}-slit_{n:02d}-{t}'
                       for n, t in zip(self.target_slit_numbers, 
                                       self.target_names)]
        return target_keys


    def info(self, log=True):
        """
        Print summary info of the mask
        """
        
        msg = '\n============================\n'
        msg += f'{self.namestr}     path={self.path}\n'
        msg += '============================\n'
        
        meta = self.meta
        msg += f"{self.namestr} {meta['SEMESTER']} {meta['PROGID']} "
        msg += f" {meta['PROGPI']} | {meta['PROGTL1']}\n"
        
        for k in self.keys:
            msg += f'{self.namestr} {self.groups[k].namestr}\n'

        for i in range(len(self.ssl)):
            msg += f"{self.namestr} {i:>4} {self.target_names[i]} "
            ri, di = self.ssl['target_ra'][i], self.ssl['target_dec'][i]
            msg += f"{ri:.6f} {di:.6f} {self.ssl['Magnitude'][i]:.1f}\n"
            
        if log:
            utils.log_comment(self.logfile, msg, verbose=True, 
                              show_date=True)
        else:
            print(msg)
            

    def make_region_file(self, regfile=None, make_figure=True, region_defaults='global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman"'):
        """
        Make a region file showing the slits
        """
        
        if regfile is None:
            regfile = f'{self.datemask}-{self.filter}.reg'
            
        sq = np.array([[-1, 1, 1, -1], [-1,-1,1,1]]).T

        #print('PA', mask.meta['SKYPA3'])
        
        pa = self.meta['SKYPA3']
        theta = pa/180*np.pi

        _mat = np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta), np.cos(theta)]])

        if make_figure:
            fig, ax = plt.subplots(1,1,figsize=(10,10))

        with open(regfile,'w') as fp:
            rows = [f'# {self.namestr}\n {region_defaults}\nfk5\n']
            fp.write(rows[0])
            for i in range(len(self.ssl)):
                ra, dec = self.ssl['slit_ra'][i], self.ssl['slit_dec'][i]
                cosd = np.cos(dec/180*np.pi)

                dim = np.cast[float]([self.ssl['Slit_width'][i], 
                                      self.ssl['Slit_length'][i]])/2.
                dim2 = np.cast[float]([self.ssl['Slit_width'][i], 
                                       self.ssl['Slit_width'][i]])/2.
                #dim[0] *= 500
                dy = float(self.ssl['Target_to_center_of_slit_distance'][i])

                yoff = np.array([0, dy])

                scl = np.array([1./cosd, 1])/3600.

                rot = (sq*dim).dot(_mat)*scl + np.array([ra, dec])
                rot2 = (sq*dim2+yoff).dot(_mat)*scl + np.array([ra, dec])

                if make_figure:
                    pl = ax.plot(*rot.T)
                    pl = ax.plot(*rot2.T)
                    ax.scatter(ra, dec, color=pl[0].get_color(), marker='x')

                row = 'polygon('
                row += ','.join([f'{v:.6f}' for v in rot.flatten()]) + ')\n'
                row += 'polygon('
                row += ','.join([f'{v:.6f}' for v in rot2.flatten()]) + ')'
                sw = float(self.ssl['Slit_width'][i])/2
                reg_label = ' # text=<<<{0} {1}>>>\n'
                row += reg_label.format(self.target_slit_numbers[i], 
                                        self.target_names[i])
                row = row.replace('<<<','{').replace('>>>','}')
                fp.write(row)
                rows.append(row)

        if make_figure:
            ax.set_aspect(1/cosd)
            ax.set_xlim(ax.get_xlim()[::-1])
        else:
            fig = None
        
        return rows, fig

    def search_targets(self, query, ignore_case=True):
        """
        String search on the target list
        """
        if ignore_case:
            has_q = [query.lower() in t.lower() for t in self.target_names]
        else:
            has_q = [query in t for t in self.target_names]
        
        if np.sum(has_q) > 0:
            qidx = np.where(has_q)[0]
            for q in qidx:
                print(f'{q} : {self.target_names[q]}')
            
            return qidx
        else:
            return None
        
        
    def read_exposures(self, min_nexp=4):
        """
        Read exposure groups from "Offset_{shift}.txt" files produced by 
        mospy
        """
        self.groups = {}
        for offset_file in self.offset_files:
            grp = ExposureGroup(offset_file=offset_file, min_nexp=min_nexp, 
                                flat=self.pixel_flat)
            if grp.frameid is None:
                continue
                
            key = grp.frameid + ''
            while key in self.groups:
                key += 'p'
            
            grp.frameid = key
            self.groups[key] = grp
    
    
    def flag_cosmic_rays(self, **kwargs):
        """
        Flag CRs in individual groups
        """
        for k in self.groups:
            self.groups[k].flag_cosmic_rays(**kwargs)
    
    
    def interpolate_slits(self, order=1, debug=False):
        """
        Fit slit coeffs and insert for neighboring slits found in the 
        Science_Slit_List table if not the same number as the slits found 
        from the flat
        """
        if hasattr(self, 'tr_fit'):
            tr_fit = self.tr_fit
        else:
            tr_coeffs = np.array([self.trace_stop[j]
                                  for j in range(self.nslits-1)])
            tr_fit = []
            for i in range(tr_coeffs.shape[1]):
                tc = np.polyfit(self.slit_stop[:-1]-1024, tr_coeffs[:,i], 
                                order)
                tr_fit.append(tc)

            tr_fit = np.array(tr_fit)
            self.tr_fit = tr_fit
            self.trace_stop_orig = self.trace_stop
            self.trace_start_orig = self.trace_start
            
        if len(self.ssl_stop) > self.nslits:
            utils.log_comment(self.logfile, 
                              f'{self.namestr}: Found extra slit!', 
                              verbose=True)
                            
            si = np.interp(self.ssl_stop, self.slit_stop,
                           np.arange(self.nslits))
            new = np.where(np.abs(si-np.round(si)) > 0.15)[0][::-1]
            
            if debug:
                print('ssl interp', si)
                
            for j in new:
                xi = int(si[j])+1
                utils.log_comment(self.logfile,
                                  f'{self.namestr}: Insert slit in {xi}',
                                  verbose=True)
                
                self.slit_stop = np.insert(self.slit_stop, xi,
                                           self.ssl_stop[j])
                self.slit_start = np.insert(self.slit_start, xi+1,
                                            self.ssl_stop[j]+5.35)
        
        elif ((self.nslits - len(self.ssl)) > 0) & ((self.nslits - len(self.ssl)) <= 2):
            pop = []
            keep = []
            for j in range(len(self.slit_stop)):
                ds = self.slit_stop[j] - self.ssl_stop
                if np.abs(ds).min() < 10:
                    keep.append(j)
                else:
                    msg = f'Pop extra slit {j} at {self.slit_stop[j]}'
                    msg += f' (ds={ds[np.nanargmin(np.abs(ds))]:.2f})'
                    utils.log_comment(self.logfile, msg, verbose=True)
                    
                    pop.append(j)
                    
            self.slit_stop = self.slit_stop[keep]
            self.slit_start = self.slit_start[keep]
            self.trace_stop = [self.trace_stop[k] for k in keep]
            self.trace_start = [self.trace_start[k] for k in keep]
            #self.trace_stop_orig = [self.trace_stop_orig[k] for k in keep]
            #self.trace_start_orig = [self.trace_start_orig[k] for k in keep]
            
        self.nslits = len(self.slit_stop)
        
        if len(self.ssl_stop) > self.nslits:
            if self.ssl_stop[0] < self.slit_start[0]+5:
                utils.log_comment(self.logfile,
                                  f'{self.namestr}: Remove first empty slit',
                                  verbose=True)

                self.ssl = self.ssl[1:]
        
        self.trace_stop = np.array([np.polyval(tr_fit[j,:], self.slit_stop-1024)
                                    for j in range(tr_fit.shape[0])]).T
        
        self.trace_mid = np.array([np.polyval(tr_fit[j,:], (self.slit_stop + self.slit_stop)/2.-1024)
                                    for j in range(tr_fit.shape[0])]).T

        self.trace_start = np.array([np.polyval(tr_fit[j,:], self.slit_start-1024)
                                     for j in range(tr_fit.shape[0])]).T

        
    def find_slits(self, x0=1024, initial_thresh=800, thresh_step=1.3, make_figure=True, interpolate_coeffs=True, verbose=True, max_iter=5, use_ssl=False):
        """
        Find slits based on gradient of the combflat image
        """
        import peakutils

        grad = (np.gradient(self.flat[0].data, axis=0))
        self.grad = grad

        xarr = np.arange(2048)
        self.xarr = xarr

        if initial_thresh is None:
            gx = np.abs(grad[:,x0])
            thres = 0.2*np.nanmax(gx)
        else:
            thres = initial_thresh*1
        
        it = 0
        pindexes, nindexes = [0], [0,0]
        pn, nn = len(pindexes), len(nindexes)
        self.use_ssl_slit = use_ssl
        if use_ssl:
            msg = f'{self.namestr} # Use SSL table for slit definition'
            utils.log_comment(self.logfile, msg, verbose=True)
            
            self.slit_stop = np.cast[int](np.round(self.ssl_stop))
            self.slit_start = np.cast[int](np.round(self.ssl_start))
            
            gx = np.abs(grad[:,x0])
            thres = np.nanmax(gx)*0.2
            pindexes = peakutils.indexes(grad[:,x0], thres=thres,
                                         thres_abs=True, min_dist=12)
            nindexes = peakutils.indexes(-grad[:,x0], thres=thres,
                                         thres_abs=True, min_dist=12)
            
            self.nslits = len(self.slit_stop)
            for j in range(self.nslits):
                ds = self.slit_stop[j] - nindexes
                if np.nanmin(np.abs(ds)) < 8:
                    isl = np.nanargmin(np.abs(ds))
                    self.slit_stop[j] = nindexes[isl]
                #
                ds = self.slit_start[j] - pindexes
                if np.nanmin(np.abs(ds)) < 8:
                    isl = np.nanargmin(np.abs(ds))
                    self.slit_start[j] = pindexes[isl]
                
            self.tr_fit = np.array([[-2.76466436e-09,  2.70016208e-08],
                                    [ 2.91466864e-07,  1.29269336e-03],
                                    [ 1.00002258e+00,  1.02418277e+03]])
            self.interpolate_slits()
            
        else:
            while ((pn != nn) | (np.maximum(nn, pn) > len(self.ssl))) & (it < max_iter):
        
                pindexes = peakutils.indexes(grad[:,x0], thres=thres,
                                             thres_abs=True, min_dist=12)
                nindexes = peakutils.indexes(-grad[:,x0], thres=thres,
                                             thres_abs=True, min_dist=12)
                pn, nn = len(pindexes), len(nindexes)
        
                msg = f'{self.namestr} # gradient iteration {it},'
                msg += f' thresh={thres:.0f} nn={nn} np={pn} nssl={len(self.ssl)}'
                utils.log_comment(self.logfile, msg, verbose=verbose)
        
                thres *= thresh_step
                it += 1

            self.slit_stop = nindexes
            self.slit_start = pindexes
            self.nslits = len(self.slit_stop)
            if self.nslits != (len(self.ssl)):
                raise ValueError
                
            ############
            # Fit curved slits
            slitp = []
            slitn = []
            for pi in pindexes:
                yi = grad[pi-1:pi+2,x0]
                xi = xarr[pi-1:pi+2]
                c = np.polyfit(xi, yi, 2)
                ymax = -c[1]/(2*c[0])
                slitp.append([[x0,ymax]])

            for pi in nindexes:
                yi = -grad[pi-1:pi+2,x0]
                xi = xarr[pi-1:pi+2]
                c = np.polyfit(xi, yi, 2)
                ymax = -c[1]/(2*c[0])
                slitn.append([[x0,ymax]])

            for x in range(16, 2048-16, 16):
                pi = peakutils.indexes(grad[:,x], thres=thres/2, 
                                       thres_abs=True, min_dist=8)
                ni = peakutils.indexes(-grad[:,x], thres=thres/2,
                                       thres_abs=True, min_dist=8)
                for j in range(self.nslits):
                    dp = pindexes[j]-pi
                    dn = nindexes[j]-ni

                    for k in np.where(np.abs(dp) < 5)[0]:
                        yi = grad[pi[k]-1:pi[k]+2,x]
                        xi = xarr[pi[k]-1:pi[k]+2]
                        c = np.polyfit(xi, yi, 2)
                        ymax = -c[1]/(2*c[0])
                        slitp[j].append([x,ymax])

                    for k in np.where(np.abs(dn) < 5)[0]:
                        yi = -grad[ni[k]-1:ni[k]+2,x]
                        xi = xarr[ni[k]-1:ni[k]+2]
                        c = np.polyfit(xi, yi, 2)
                        ymax = -c[1]/(2*c[0])
                        slitn[j].append([x,ymax])

            ###########
            # Fit them
            trp = []
            trn = []
            for i in range(self.nslits):
                # msg = f'{self.datemask}  Fit slit {i:>2}: {self.slit_start[i]:>4}'
                #msg += f' - {self.slit_stop[i]:>4} | {self.target_names[i]}'
                #utils.log_comment(self.logfile, msg, verbose=verbose)

                ap = np.array(slitp[i])
                cp = np.polyfit(ap[:,0]-1024, ap[:,1], 2)
                vp = np.polyval(cp, ap[:,0]-1024)
                ok = np.abs(vp-ap[:,1]) < 1
                cp = np.polyfit(ap[ok,0]-1024, ap[ok,1], 2)

                an = np.array(slitn[i])
                cn = np.polyfit(an[:,0]-1024, an[:,1], 2)
                vn = np.polyval(cn, an[:,0]-1024)
                ok = np.abs(vn-an[:,1]) < 1
                cn = np.polyfit(an[ok,0]-1024, an[ok,1], 2)

                trp.append(cp)
                trn.append(cn)

            self.trace_start = trp
            self.trace_stop = trn
        
            if interpolate_coeffs:
                self.interpolate_slits()
        
        ##########
        # Make the figure
        if make_figure:
            fig, ax = plt.subplots(1,1,figsize=(12, 5))

            ax.plot(grad[:, x0])
            #ax.plot(grad[:, 1200])

            ax.scatter(xarr[self.slit_start], grad[self.slit_start,x0],
                       color='r')

            ax.scatter(xarr[self.slit_stop], grad[self.slit_stop,x0],
                       color='r')

            for pi, ni in zip(self.slit_start, self.slit_stop):
                ax.plot(xarr[[pi, ni]], grad[[pi, ni], x0], color='r',
                        alpha=0.5)

            yl = ax.get_ylim()
            dlab = 0.1*(yl[1]-yl[0])

            for j, pi in enumerate(self.slit_start):
                ax.text(xarr[pi], grad[pi, x0]+dlab, f'{j}', ha='center',
                        va='bottom', fontsize=8)

            for j, pi in enumerate(self.slit_stop):
                ax.text(xarr[pi], grad[pi, x0]-dlab, f'{j}', ha='center',
                        va='top', fontsize=8)

            ax.set_ylim(yl[0]-3*dlab, yl[1]+3*dlab)
            
            # Target table expected stop 
            ax.vlines(self.ssl_stop, yl[0], -200, color='orange', alpha=0.3)
            
            ax.grid()
            ax.text(0.05, 0.95, self.namestr, ha='left', va='top',
                    transform=ax.transAxes)
            ax.set_xlabel('y pixel')
            ax.set_ylabel('gradient')
            fig.tight_layout(pad=0.3)

        else:
            fig = None

        return fig


    def find_longpos_trace(self, use_plan=None, thresh=100, make_figure=True):
        import peakutils

        x0 = 1024
        if use_plan is None:
            use_plan = self.plans[0]

        pa, pb = use_plan
        diff = self.groups[pb].sci[0,:,:] - self.groups[pa].sci[0,:,:]

        # Mask out extra exposures
        for p in [pa, pb]:
            if self.groups[p].nexp > 1:
                self.groups[p].var[1:,:,:] = 0
                self.groups[p].flag_cosmic_rays(minmax_nrej=0, sigma=1000)

        prof = diff[:,1024]

        pindexes = []
        thres = thresh*1
        it = -1

        while (len(pindexes) != 1) & (it < 10):
            it += 1
            thres *= 1.2
            pindexes = peakutils.indexes(diff[:,x0], thres=thres,
                                 thres_abs=True, min_dist=12)

            print(thres, len(pindexes))

        peak_flux = diff[pindexes[0], x0]

        nindexes = peakutils.indexes(-diff[:,x0], thres=thres,
                                             thres_abs=True, min_dist=12)

        pn, nn = len(pindexes), len(nindexes)
        #plt.plot(prof)

        self.xarr = np.arange(2048)
        yarr = np.arange(2048)
        slitp = []
        slitn = []

        for x in range(16, 2048-16, 16):
            px = peakutils.indexes(diff[:,x], thres=thres/2., 
                               thres_abs=True, min_dist=8)
            if len(px) == 1:
                peaky = peakutils.interpolate(yarr+0.5, diff[:,x], 
                                        ind=px, width=10)
                slitp.append([x, peaky[0]])

            nx = peakutils.indexes(-diff[:,x], thres=thres/2., 
                               thres_abs=True, min_dist=8)
            if len(nx) == 1:
                peaky = peakutils.interpolate(yarr+0.5, -diff[:,x], 
                                        ind=nx, width=10)
                slitn.append([x, peaky[0]])

        # Fit the curved traces
        
        i = 0
        ap = np.array(slitp)
        an = np.array(slitn)

        #print('xx', an.shape, ap.shape)

        ok = (np.abs(ap[:,1]-np.nanmedian(ap[:,1])) < 20)
        cp = np.polyfit(ap[ok,0]-1024, ap[ok,1], 2)
        vp = np.polyval(cp, ap[:,0]-1024)
        ok = (np.abs(vp-ap[:,1]) < 1) & (np.abs(ap[:,1]-np.median(ap[:,1])) < 10)
        cp = np.polyfit(ap[ok,0]-1024, ap[ok,1], 2)
        ap = ap[ok,:]

        ok = (np.abs(an[:,1]-np.nanmedian(an[:,1])) < 20)
        cn = np.polyfit(an[ok,0]-1024, an[ok,1], 2)
        vn = np.polyval(cn, an[:,0]-1024)
        ok = (np.abs(vn-an[:,1]) < 1) & (np.abs(an[:,1]-np.median(an[:,1])) < 10)
        cn = np.polyfit(an[ok,0]-1024, an[ok,1], 2)
        an = an[ok,:]
        
        #print('xx', nindexes, pindexes)
        
        start = np.minimum(nindexes, pindexes)[0]-20
        stop = np.maximum(nindexes, pindexes)[0]+20
        ## print('xxx', start, stop)
        
        #stop = start+30
        #start = start+10
        #start = stop-30
        #stop = stop-10

        self.slit_start = np.array([start])
        self.slit_stop = np.array([stop])
        
        dither_offset = 0.1799*(stop-start-40)/2.
        for p in [pa, pb]:
            off = self.groups[p].yoffset*1.
            for j, o in enumerate(off):
                off_i = dither_offset*(1-2*(o < 0))
                print(f'Set {p} yoffset {off_i:.2f}')
                self.groups[p].yoffset[j] = off_i
            
            #print('xx off', self.groups[p].yoffset)
        
        cpp = cp*1
        cnn = cp*1
        
        if nindexes[0] > pindexes[0]:
            cpp[-1] = cp[-1] - 20
            cnn[-1] = cn[-1] + 20
            
            self.trace_start = [cpp]
            self.trace_stop = [cnn]
        else:
            cpp[-1] = cp[-1] + 20
            cnn[-1] = cn[-1] - 20

            self.trace_start = [cnn]
            self.trace_stop = [cpp]

        self.ssl['Target_to_center_of_slit_distance'] = 0.
        self.nslits = 1

        targ = self.groups['A'].img[0][0].header['TARGNAME'].strip().replace(' ','') + f'-{pa}'
        self.ssl.remove_column('Target_Name')
        self.ssl['Target_Name'] = targ
        self.ssl = self.ssl[:1]

        #self.target_names = [t.strip() for t in self.ssl['Target_Name']]
        #self.target_keys = [f'{self.datemask}-{self.filter}-{t}' for t in self.target_names]

        if make_figure:
            sli = slice(start, stop)

            fig, ax = plt.subplots(1,1,figsize=(14,5))
            arr = diff
            #arr = self.groups[pa].sci[0,:,:]
            ax.imshow(arr[sli,:], extent=(0, 2048, start, stop), origin='lower', vmin=-peak_flux, vmax=peak_flux)
            ax.set_aspect('auto')
            ax.plot(*an.T, color='r')
            ax.plot(*ap.T, color='r')
            ax.set_ylim(start, stop)

            ax.text(0.5, 0.5, self.target_keys[0], ha='center', va='center', color='w', transform=ax.transAxes)
            ax.grid()

            fig.tight_layout(pad=0.5)

        else:
            fig = None

        return fig

    
    def get_slit_params(self, slit, img_data=None, pad=16, skip=4, plan=None, xy_order=3, wave_order=3, verbose=True, show=True, wave_kwargs={}):
        """
        Traace parameters of a single slit
        """
        start_trace = np.polyval(self.trace_start[slit], self.xarr-1024)
        stop_trace = np.polyval(self.trace_stop[slit], self.xarr-1024)

        i0 = np.maximum(int(start_trace.min()) - pad, 0)
        i1 = np.minimum(int(stop_trace.max()) + pad, 2048)

        istop = np.polyval(self.trace_stop[slit], 0)
        istart = np.polyval(self.trace_start[slit], 0)

        imed = int((i0+i1)/2)

        msg = f'{self.datemask}  Limits for slit {slit}: {i0} - {i1} ({imed})'
        utils.log_comment(self.logfile, msg, verbose=True)

        if plan is None:
            plan = self.plans[0]

        pa, pb = plan

        if img_data is None: 
            img_data = self.groups[pa].sci[0,:,:]

        ####### median along center of the slit with sky lines
        sky_row = np.nanmedian(img_data[imed-5:imed+5,:], axis=0)

        ############
        # Cross-correlate xshifts
        y0 = int(start_trace.max())
        y1 = int(stop_trace.min())

        ysh = np.arange(y0+skip//2, y1-skip//2, skip)
        xsh = ysh*0.
        for i, yi in tqdm(enumerate(ysh)):
            row = np.nanmedian(img_data[yi-skip//2:yi+skip//2+1,:],
                               axis=0)[None, :]
            cc = phase_cross_correlation(row, sky_row[None,:],
                                         upsample_factor=100, 
                                reference_mask=~np.isnan(row), 
                                moving_mask=~np.isnan(sky_row[None,:]), 
                                return_error=True)
            (_, xsh[i]), error, diffphase = cc
        
        if (i1-i0 > 200):
            dxmax = 300
        else:
            dxmax = 50
            
        xok = np.isfinite(xsh)
        if len(xsh) > 1:
            xok &= np.abs(np.gradient(xsh)) < 3
            xok &= np.abs(xsh) < dxmax
        
        for it in range(3):
            xy_coeffs = np.polyfit(ysh[xok]-imed, xsh[xok], xy_order)
            xfit = np.polyval(xy_coeffs, ysh-imed)
            xok = np.isfinite(xsh)
            if len(xsh) > 1:
                xok &= (np.abs(np.gradient(xsh)) < 3) & (np.abs(xfit-xsh) < 3)
        
        #targname =  self.ssl['Target_Name'][slit].strip()
        #slit_num = int(self.ssl['Slit_Number'][slit])
        targname = self.target_names[slit]
        slit_num = self.target_slit_numbers[slit]
         
        slit_info = {}
        slit_info['slit'] = slit
        slit_info['i0'] = i0
        slit_info['i1'] = i1
        slit_info['sky_row'] = sky_row
        slit_info['xy_coeffs'] = xy_coeffs
        slit_info['filter'] = self.filter
        slit_info['wave_order'] = wave_order
        slit_info['slice'] = slice(i0, i1)
        slit_info['trace_coeffs'] = self.trace_stop[slit]
        slit_info['pad'] = pad
        slit_info['istart'] = istart
        slit_info['istop'] = istop
        
        slit_info['xsh'] = xsh
        slit_info['ysh'] = ysh
        slit_info['xok'] = xok
                
        slit_info['target_name'] = targname
        
        slit_info['width'] = float(self.ssl['Slit_width'][slit])
        slit_info['length'] = float(self.ssl['Slit_length'][slit])
        slit_info['target_offset'] = float(self.ssl['Target_to_center_of_slit_distance'][slit])
        
        # Center of target in slit cutout
        yoff = float(slit_info['target_offset'])/0.1799
        half = slit_info['istop'] - slit_info['istart']
        ty = slit_info['istart'] - slit_info['i0'] + half/2. + yoff
        slit_info['target_y'] = ty
        
        slit_info['slit_ra'] = self.ssl['slit_ra'][slit]
        slit_info['slit_dec'] = self.ssl['slit_dec'][slit]
        slit_info['target_ra'] = self.ssl['target_ra'][slit]
        slit_info['target_dec'] = self.ssl['target_dec'][slit]
        slit_info['target_mag'] = self.ssl['Magnitude'][slit]
        slit_info['target_orig_slit'] = slit_num
        slit_info['datemask'] = self.datemask
        
        ############
        # Fit wavelength
        if show:
            
            # Show difference in slit figure
            #pa, pb = self.plans[0]
            diff = self.groups[pb].sci[0,:,:] - self.groups[pa].sci[0,:,:]
            diff /= np.sqrt(self.groups[pa].var[0,:,:])
            
            fig, axes = plt.subplots(3,2,figsize=(12,7), 
                                     gridspec_kw={'height_ratios':[2,2,1],
                                                  'width_ratios':[3.5,1]})
            
            # Difference
            for ia, _data in enumerate([diff, img_data]):
                ax = axes[ia][0]
                perc = np.nanpercentile(_data[i0:i1,:], [5, 90])
                ax.imshow(_data, vmin=perc[0], vmax=perc[1])

                if ia == 0:
                    ax.text(0.02, 0.96, 
                      f'{self.datemask} {self.filter} {slit_num}: {targname}', 
                            ha='left', va='top', transform=ax.transAxes, 
                            bbox=dict(facecolor='w', edgecolor='None'))
            
                ax.set_aspect('auto')
                ax.plot(self.xarr, start_trace, color='pink')
                ax.plot(self.xarr, stop_trace, color='pink')
                ax.set_ylim(i0, i1)
                ax.hlines(imed, 0, 200, color='r')
                ax.set_xticklabels([])

            # x(y) shift
            ax = axes[1][1]
            ax.scatter(xsh[xok], ysh[xok])
            ax.scatter(xsh, ysh, alpha=0.3)
            ax.set_yticklabels([])
            ax.grid()
            ax.set_ylim(*axes[0][0].get_ylim())
            yy = np.linspace(y0, y1, 256)
            xx = np.polyval(xy_coeffs, yy-imed)
            ax.plot(xx, yy, alpha=0.4)
            #ax.set_xlim(xsh[xok].min()-1, xsh[xok].max()+1)
            ax.set_xlim(xsh.min()-1, xsh.max()+1)
            
            ax.set_xticklabels([])
               
            #ax.set_xlabel(r'$\Delta x$')

            _lfit = fit_wavelength_from_sky(sky_row, self.filter, 
                                            order=wave_order,
                                            make_figure=False, nsplit=5,
                                            plot_axis=axes[2][0],
                                            **wave_kwargs)
            
            lam_coeffs, wave_fit, figs = _lfit
            
            axes[2][0].set_xlim(wave_fit.min()/1.e4, wave_fit.max()/1.e4)
            axes[2][0].legend(loc='upper left')
            axes[0][1].axis('off')
            axes[2][1].axis('off')

            fig.tight_layout(pad=0.5)

        else:
            fig = None
            _lfit = fit_wavelength_from_sky(sky_row, self.filter, 
                                            order=wave_order,
                                            make_figure=False,
                                            nsplit=5, **wave_kwargs)
            lam_coeffs, wave_fit, figs = _lfit
            
        slit_info['lam_coeffs'] = lam_coeffs
        slit_info['wave'] = wave_fit

        msg = ('Slit {0}: {1} {2}x{3} {4:.5f} {5:.5f}'.format(slit_num,
                    slit_info['target_name'], slit_info['width'],
                    slit_info['length'], slit_info['target_ra'],
                    slit_info['target_dec']))
                    
        utils.log_comment(self.logfile, msg, verbose=verbose)
        
        return slit_info, fig    


    def drizzle_slit_plan_single(self, slit_info, plan_i=0, linearize_wave=False, log_wave=False, kernel='point', pixfrac=1., sig_clip=(3,3), mask_trace=True, mask_offset=True, mask_overlap=False, mask_single=False, **kwargs):
        """
        Drizzle a rectified slit
        """
        from drizzlepac import adrizzle
        import astropy.wcs as pywcs

        plan = self.plans[plan_i]
        pd = self.plan_pairs[plan_i]

        pa, pb = plan
        gra = self.groups[pa]
        grb = self.groups[pb]

        ia = pd['ia'][0]
        ib = pd['ib'][0]
        
        ysl = slit_info['slice']
        cutout = (grb.sci[ib,:,:] - gra.sci[ia,:,:])[ysl,:]
        cutout_var = (grb.var[ib,:,:] + gra.var[ia,:,:])[ysl,:]
        
        exptime = grb.truitime[pd['ib']].sum() + gra.truitime[pd['ia']].sum()
        nexp = pd['n']*2

        cutout_wht = 1/cutout_var
        cutout_wht[cutout_var == 0] = 0
        
        # Mask out-of-slit pixels
        slit = slit_info['slit']
        start_trace = np.polyval(self.trace_start[slit], self.xarr-1024)
        stop_trace = np.polyval(self.trace_stop[slit], self.xarr-1024)
        yp, xp = np.indices(cutout.shape)
        cutout_mask = (yp + ysl.start >= start_trace) & (yp + ysl.start <= stop_trace)
        
        rel_offset = np.abs(self.groups[pa].yoffset[0] -
                            self.groups[pb].yoffset[0])
        
        #########
        # Two headers for slit distortions
        h0 = pyfits.Header()
        hdist = pyfits.Header()
        
        trx = slit_info['xy_coeffs']*1
        trn = slit_info['trace_coeffs']*1
        trw = slit_info['lam_coeffs']*1
        
        for h in [h0, hdist]:
            h['NAXIS'] = 2
            h['NAXIS1'] = cutout.shape[1]
            h['NAXIS2'] = cutout.shape[0]
            h['CRPIX1'] = h['NAXIS1']/2.
            h['CRPIX2'] = h['NAXIS2']/2.
            h['CRVAL1'] = trw[-1]#/1.e10
            h['CRVAL2'] = 0.
            h['CD1_1'] = trw[-2]#/1.e10
            h['CD2_2'] = 0.1799 #/3600

            #h['CTYPE1'] = 'RA---TAN-SIP'
            #h['CTYPE2'] = 'DEC--TAN-SIP'
            h['CTYPE1'] = '-SIP'
            h['CTYPE2'] = '-SIP'
            h['A_ORDER'] = 3
            h['B_ORDER'] = 3

        #hdist['CTYPE1'] = 'RA---TAN-SIP'
        #hdist['CTYPE2'] = 'DEC--TAN-SIP'


        ##########
        # Slit distortion as SIP coefficients
        # ToDo: full y dependence of the trace curvature

        # Tilted x(y)
        ncoeff = len(trx)
        for ai in range(ncoeff-1):
            hdist[f'A_0_{ai+1}'] = -trx[ncoeff-2-ai]
        
        # Curved trace y(x)
        ncoeff = len(trn)
        if hasattr(self, 'tr_fit'):
            # Full distorted trace
            print('distorted trace')
            dy = (slit_info['i0'] + slit_info['i1'])/2. - 1024
            for ai in range(ncoeff-1):
                bcoeff = self.tr_fit[ncoeff-2-ai,:]*1.
                bcoeff[1] += bcoeff[0]*dy
                hdist[f'B_{ai+1}_0'] = -bcoeff[1]
                hdist[f'B_{ai+1}_1'] = -bcoeff[0]
        else:
            for ai in range(ncoeff-1):
                hdist[f'B_{ai+1}_0'] = -trn[ncoeff-2-ai]

        # Wavelength
        ncoeff = len(trw)
        for ai in range(ncoeff-2):
            h0[f'A_{ai+2}_0'] = trw[ncoeff-3-ai]/trw[-2]
            hdist[f'A_{ai+2}_0'] = trw[ncoeff-3-ai]/trw[-2]
        
        in_wcs = pywcs.WCS(hdist)
        in_wcs.pscale = 0.1799

        if linearize_wave:
            print('Linearize wavelength array')
            
            gr = grating_summary[self.filter]
            
            # Recenter and set log wavelength spacing
            sh = [cutout.shape[0], gr['N']]
            
            # Set poly keywords
            h0['NAXIS1'] = gr['N']
            h0['CRPIX1'] = gr['N']//2
            h0['CRVAL1'] = gr['ref_wave']#/1.e10
            h0['CD1_1'] = gr['dlam']#/1.e10
            
            if log_wave: 
                loglam_pars = get_grating_loglam('J')
                #h0['NAXIS1'] = loglam_pars[0]
                #h0['CRPIX1'] = loglam_pars[0]/2
                #h0['CRVAL1'] = loglam_pars[1]
                #h0['CD1_1'] = loglam_pars[3][0]
                h0['CTYPE1'] = '-SIP'
                coeffs = loglam_pars[2][::-1]
                h0['CRVAL1'] = coeffs[0]
                h0['CD1_1'] = coeffs[1]
                h0['A_2_0'] = coeffs[2]
                h0['A_3_0'] = coeffs[3]
                h0['A_ORDER'] = 3
                h0['B_ORDER'] = 3

            else:
                # Strip SIP keywords
                h0.remove('CTYPE1')
                h0.remove('CTYPE2')
                #h0['CTYPE1'] = 'WAVE'
                
                for k in list(h0.keys()):
                    test = k.startswith('A_') | k.startswith('B_')
                    test |= k in ['A_ORDER','B_ORDER']
                    if test:
                        h0.remove(k)
                    
            lam_coeffs = np.array([gr['dlam'], gr['ref_wave']])
            
        else:
            sh = cutout.shape
            logwcs = {}
            lam_coeffs = slit_info['lam_coeffs']*1
            
        outsci = np.zeros((pd['n']*2, *sh), dtype=np.float32)
        outwht = np.zeros((pd['n']*2, *sh), dtype=np.float32)
        outctx = np.zeros((pd['n']*2, *sh), dtype=np.int16)

        npl = pd['n']
        slit = slit_info['slit']
        targname = slit_info['target_name']

        msg = f'{self.namestr} # Drizzle N={npl} {pa}-{pb} exposure pairs for slit {slit}: {targname}'
        utils.log_comment(self.logfile, msg, verbose=True)
        
        self.distorted_header = hdist
        
        # Do the drizzling
        for i, (ia, ib) in tqdm(enumerate(zip(pd['ia'], pd['ib']))):
            
            if False: #pd['n'] > 1:
                if i == 0:
                    nb = pd['ib'][i+1]
                    na = pd['ia'][i+1]
                else:
                    nb = pd['ib'][i-1]
                    na = pd['ia'][i-1]
                
                skyb = (grb.sci[ib,:,:] + grb.sci[nb,:,:])/2.
                svarb = (grb.var[ib,:,:] + grb.var[nb,:,:])/2.
                skya = (gra.sci[ia,:,:] + gra.sci[na,:,:])/2.
                svara = (gra.var[ia,:,:] + gra.var[na,:,:])/2.
                
            else:
                skyb = grb.sci[ib,:,:]
                svarb = grb.var[ib,:,:]
                skya = gra.sci[ia,:,:]
                svara = gra.var[ia,:,:]
            
            diff_b = (grb.sci[ib,:,:] - skya)[ysl,:].astype(np.float32, copy=False)
            diff_a = (gra.sci[ia,:,:] - skyb)[ysl,:].astype(np.float32, copy=False)
            
            wht_b = 1./(grb.var[ib,:,:] + svara)[ysl,:].astype(np.float32, copy=False)
            wht_b[~np.isfinite(wht_b)] = 0

            wht_a = 1./(gra.var[ia,:,:] + svarb)[ysl,:].astype(np.float32, copy=False)
            wht_a[~np.isfinite(wht_a)] = 0
                        
            if mask_trace:
                #cutout_wht *= cutout_mask
                wht_b *= cutout_mask
                wht_a *= cutout_mask
                
            #sci = cutout.astype(np.float32, copy=False)
            #wht = cutout_wht.astype(np.float32, copy=False)

            # Trace drift
            hdist['CRPIX2'] = hdist['NAXIS2']/2. + pd['shift'][i]
            in_wcs = pywcs.WCS(hdist)
            in_wcs.pscale = 0.1799

            ### B position, positive
            h0['CRPIX2'] = h['NAXIS2']/2. + self.groups[pb].yoffset[0]/0.1799
            out_wcs = pywcs.WCS(h0)
            out_wcs.pscale = 0.1799
            
            adrizzle.do_driz(diff_b, in_wcs, wht_b, out_wcs, 
                             outsci[i*2,:,:], outwht[i*2,:,:], 
                             outctx[i*2,:,:],
                             1., 'cps', 1,
                             wcslin_pscale=0.1799, uniqid=1,
                             pixfrac=pixfrac, kernel=kernel, fillval='0',
                             wcsmap=None)

            ### A position, negative
            h0['CRPIX2'] = h['NAXIS2']/2. + self.groups[pa].yoffset[0]/0.1799
            out_wcs = pywcs.WCS(h0)
            out_wcs.pscale = 0.1799
            
            adrizzle.do_driz(diff_a, in_wcs, wht_a, out_wcs, 
                             outsci[i*2+1,:,:], outwht[i*2+1,:,:], 
                             outctx[i*2+1,:,:],
                             1., 'cps', 1,
                             wcslin_pscale=0.1799, uniqid=1,
                             pixfrac=pixfrac, kernel=kernel, fillval='0',
                             wcsmap=None)

        scale_fwhm = pd['fwhm'] / np.nanmin(pd['fwhm'])
        scale_flux = np.nanmax(pd['scale']) / pd['scale']
        
        # weight by inverse FWHM
        if 0:
            msg = f'{self.namestr} # Weight by inverse FWHM'
            utils.log_comment(self.logfile, msg, verbose=True)
            scale_weight = 1./scale_fwhm
        else:
            # weight by flux rather than fwhm as per MOSDEF
            if np.allclose(scale_flux, 1):
                msg = f'{self.namestr} # No scales found, weight by FWHM'
                utils.log_comment(self.logfile, msg, verbose=True)
                scale_weight = 1./scale_fwhm
            else:
                msg = f'{self.namestr} # Weight by sum x fwhm'
                utils.log_comment(self.logfile, msg, verbose=True)
                scale_weight = 1./(scale_flux*scale_fwhm)
        
        pd['scale_weight'] = scale_weight
        
        for i in range(pd['n']):
            outsci[i*2,:,:] *= scale_flux[i]
            outwht[i*2,:,:] *= 1./scale_flux[i]**2 * scale_weight[i]
            outsci[i*2+1,:,:] *= scale_flux[i]
            outwht[i*2+1,:,:] *= 1./scale_flux[i]**2 * scale_weight[i]
        
        if mask_single:
            sing = (outwht[0::2,:,:] > 0).sum(axis=0) > 0
            sing &= (outwht[1::2,:,:] > 0).sum(axis=0) > 0
            outwht *= sing
        
        avg = 0    
        if sig_clip is not None:
            clip = np.isfinite(outwht) & (outwht > 0)
            c0 = clip.sum()
            for it in range(sig_clip[0]):    
                if it > 0:
                    resid = (outsci - avg)*np.sqrt(outwht)
                    clip = (np.abs(resid) < sig_clip[1]) & (outwht > 0)

                msg = f'{self.namestr} # Drizzle {slit} sigma clip {it} {(1-clip.sum()/c0)*100:.2f} %'
                utils.log_comment(self.logfile, msg, verbose=True)
                num = (outsci*outwht*clip).sum(axis=0)
                den = (outwht*clip).sum(axis=0)
                avg = num/den

            outwht[~clip] = 0

        msk = (outwht <= 0) | ~np.isfinite(outsci+outwht)
        outsci[msk] = 0
        outwht[msk] = 0

        h0['SIGCLIPN'] = sig_clip[0], 'Sigma clipping iterations'
        h0['SIGCLIPV'] = sig_clip[1], 'Sigma clipping level'

        h0['CRPIX2'] = h['NAXIS2']/2.

        h0['EXPTIME'] = exptime, 'Integration time, seconds'
        h0['NEXP'] = nexp, 'Number of raw exposures'
        h0['KERNEL'] = kernel, 'Drizzle kernel'
        h0['PIXFRAC'] = pixfrac, 'Drizzle pixfrac'

        h0['PLAN'] = f'{pa}{pb}', 'Dither plan'
        h0['OFFSETA'] = self.groups[pa].yoffset[0]/0.1799, f'Offset {pa}, pixels'
        h0['OFFSETB'] = self.groups[pb].yoffset[0]/0.1799, f'Offset {pb}, pixels'

        h0['TARGNAME'] = slit_info['target_name'], 'Target name from slit table'
        h0['RA_SLIT'] = slit_info['slit_ra'], 'Target RA from slit table'
        h0['DEC_SLIT'] = slit_info['slit_dec'], 'Target declination from slit table'
        h0['RA_TARG'] = slit_info['target_ra'], 'Target RA from slit table'
        h0['DEC_TARG'] = slit_info['target_dec'], 'Target declination from slit table'
        h0['MAG_TARG'] = slit_info['target_mag'], 'Magnitude from slit table'
        h0['FILTER'] = self.filter, 'MOSFIRE filter'
        h0['DATEMASK'] = slit_info['datemask'], 'Unique mask identifier'

        h0['SLITIDX'] = slit_info['slit'], 'Slit number, counting from y=0'
        h0['SLITNUM'] = (slit_info['target_orig_slit'], 
                          'Slit number in mask table')
        h0['Y0'] = slit_info['i0'], 'Bottom of the slit cutout'
        h0['Y1'] = slit_info['i1'], 'Top of the slit cutout'
        h0['YSTART'] = slit_info['istart'], 'Bottom of the slit'
        h0['YSTOP'] = slit_info['istop'], 'Top of the slit'
        h0['YPAD'] = slit_info['pad'], 'Cutout padding'
        
        h0['CUNIT2'] = 'arcsec'
        
        if linearize_wave:
            h0['CTYPE1'] = 'WAVE'
            h0['CUNIT1'] = 'Angstrom'
            
            if log_wave:
                h0['CTYPE1'] = 'WAVE-LOG'
                for k in list(h0.keys()):
                    test = k.startswith('A_') | k.startswith('B_')
                    test |= k in ['A_ORDER','B_ORDER']
                    if test:
                        h0.remove(k)
                
        h0['MASKOFF'] = mask_offset, 'Mask pixels outside of offsets'

        if mask_offset:
            yp = np.arange(sh[0])
            ymsk = (yp < np.abs(h0['OFFSETA'])+slit_info['pad']*mask_overlap) 
            ymsk |= (yp > sh[0]-np.abs(h0['OFFSETB'])-slit_info['pad']*mask_overlap)
            outsci[:,ymsk,:] = 0
            outwht[:,ymsk,:] = 0

        # Target position
        h0['TARGOFF'] = slit_info['target_offset'], 'Target offset to slit center, arcsec'
        h0['TARGYPIX'] = slit_info['target_y'], 'Expected central y pixel of target'

        h0['CRPIX2'] = h0['TARGYPIX']

        h0['TRAORDER'] = len(trn), 'Order of curved trace fit'
        for i, c in enumerate(trn):
            h0[f'TRACOEF{i}'] = c, 'Trace coefficient'
            
        h0['LAMORDER'] = len(lam_coeffs)-1, 'Order of wavelength solution'
        for i, c in enumerate(lam_coeffs):
            h0[f'LAMCOEF{i}'] = c, 'Wavelength solution coefficient'

        meta = self.meta
        for k in OBJ_KEYS:
            h0[k] = meta[k]

        for k in SEQ_KEYS:
            for ext in ['_MIN','_MED','_MAX']:
                key = f'{k}{ext}'
                h0[key] = meta[key]
        
        h0['MJD-OBS'] = h0['MJD-OBS_MED']
        return h0, outsci, outwht


    def drizzle_all_plans(self, slit_info, skip_plans=[], **kwargs):
        """
        Drizzle and combine all available plans
        
        Todo: 
            1) Find max S/N within the trace window
            2) Separate extensions for combinations by offset position 
            3) combination across multiple plans
        """
        #drz = {}
        num = None
        fi = 0
        
        self.drizzled_plans = []
        
        for i, plan in enumerate(self.plans):
            if plan in skip_plans:
                continue
                
            key = ''.join(plan)
            drz_i = self.drizzle_slit_plan_single(slit_info, plan_i=i, 
                                                     **kwargs)
            
            self.drizzled_plans.append(drz_i)
            
            if num is None:
                num = (drz_i[1]*drz_i[2]).sum(axis=0)
                wht = drz_i[2].sum(axis=0)
                head = drz_i[0]

            else:
                num += (drz_i[1]*drz_i[2]).sum(axis=0)
                wht += drz_i[2].sum(axis=0)
                
                head['EXPTIME'] += drz_i[0]['EXPTIME']
                head['NEXP'] += drz_i[0]['NEXP']
                head[f'PLAN{i+1}'] = key
            
            pd = self.plan_pairs[i]
            pa, pb = plan
            gra = self.groups[pa]
            grb = self.groups[pb]
            for j in range(pd['n']):
                fi += 1
                head[f'FILEA{fi}'] = gra.files[pd['ia'][j]], f'File from plan {pa}'
                head[f'FILEB{fi}'] = grb.files[pd['ib'][j]], f'File from plan {pb}'
                head[f'SHIFT{fi}'] = np.float32(pd['shift'][j]),  'Shift [pix]'
                head[f'FWHM{fi}'] = np.float32(pd['fwhm'][j]), 'FWHM [pix]'
                head[f'SUM{fi}'] = np.float32(pd['scale'][j]), 'Profile sum'
                head[f'WSCALE{fi}'] = np.float32(pd['scale_weight'][j]), 'Weight scaling'
                
                
        outsci = num/wht
        outwht = wht
        msk = (wht == 0) | (~np.isfinite(outsci+outwht))
        outsci[msk] = 0
        outwht[msk] = 0

        hdu = pyfits.HDUList()
        hdu.append(pyfits.ImageHDU(data=outsci, header=head, name='SCI'))
        hdu.append(pyfits.ImageHDU(data=outwht, header=head, name='WHT'))

        return hdu
    
    
    def get_plan_drift(self, slit_info, plan_i=0, ax=None, fwhm_bounds=[2,10],
                       driz_kwargs=dict(sig_clip=(5, 3), mask_single=True, mask_offset=False, mask_trace=True), profile_model=None, use_peak=True):
        """
        Get drift from bright target
        """
        from astropy.modeling.fitting import LevMarLSQFitter
        from astropy.modeling.models import Lorentz1D, Gaussian1D, Moffat1D
        if profile_model is None:
            profile_model = Moffat1D()

        pd = self.plan_pairs[plan_i]
        pd['fwhm'] = np.ones_like(pd['fwhm'])
        pd['shift'] = np.zeros_like(pd['shift'])
        pd['scale'] = np.ones_like(pd['fwhm'])
        pd['scale_weight'] = np.ones_like(pd['fwhm'])

        h0, outsci, outwht = self.drizzle_slit_plan_single(slit_info, 
                                                plan_i=plan_i, **driz_kwargs)

        mu = slit_info['target_y']
        if 'x_0' in profile_model.param_names:
            profile_model.x_0 = mu
        else:
            profile_model.mean = mu
        
        # bounds on fwhm
        for k in ['stddev', 'fwhm','gamma']:
            if k in profile_model.param_names:
                pi = profile_model.param_names.index(k)
                fscl = profile_model.parameters[pi] / profile_model.fwhm
                profile_model.bounds[k] = [fwhm_bounds[0]*fscl, fwhm_bounds[1]*fscl]
                
        stacksci = outsci
        stackwht = stacksci*0.

        profs = []

        (pa, pb) = plan = self.plans[plan_i]

        fwhm = self.groups[pa].meta['GUIDFWHM']

        npair = pd['n']
        fit_fwhm = np.zeros(npair)
        fit_x = np.zeros(npair)
        fit_sum = np.zeros(npair)

        # Combine A-B
        num = outsci[0::2,:,:]*outwht[0::2,:,:] + outsci[1::2,:,:]*outwht[1::2,:,:]
        den = outwht[0::2,:,:] + outwht[1::2,:,:]
        ab_avg = num/den
        ab_avg[den <= 0] = 0

        for i in range(npair):
            kws = dict(alpha=0.5, color=plt.cm.jet((i+1)/npair))
            yprof = drizzled_profile((h0, ab_avg[i,:,:], den[i,:,:]),
                                     ax=ax, plot_kws=kws)

            xprof = np.arange(len(yprof))
            #ok = (yprof > 0) & (np.abs(xprof-slit_info['target_y']) < 10)
            if use_peak:
                xmax = xprof[np.nanargmax(yprof)]
                #print('xxx', xmax, slit_info['target_y'])
            else:
                xcut = np.abs(xprof-slit_info['target_y']) < 10
                xmax = xprof[np.nanargmax(yprof*xcut)]
                
            ok = (yprof > 0) & (np.abs(xprof-xmax) < 10)
            if 'x_0' in profile_model.param_names:
                profile_model.x_0 = xmax
            else:
                profile_model.mean = xmax
            
            profile_model.amplitude = yprof[ok].max()

            #print(xprof.shape, yprof.shape, ok.sum())

            mfit = LevMarLSQFitter()(profile_model, xprof[ok], yprof[ok])
            fit_fwhm[i] = mfit.fwhm
            fit_x[i] = mfit.x_0.value
            fit_sum[i] = np.trapz(mfit(xprof), xprof)
            #fit_sum[i] = mfit(xprof).max()
            
            #if ax is not None:
            #    ax.plot(xprof, mfit(xprof), color='r')
        
        ok = np.abs(fit_x - np.median(fit_x)) < 2.5
        ok &= (fit_fwhm < 17) & (fit_fwhm > 0.8)
        ok &= np.isfinite(fit_x + fit_fwhm + fit_sum)
        
        bad = np.where(~ok)[0]
        if len(bad) > 0.5*pd['n']:
            msg = f'{self.namestr} # Too many bad exposures found in drift'
            msg += f" ({len(bad)} / {pd['n']})\n"
            utils.log_comment(self.logfile, msg, verbose=True)
            pd['shift'] = np.zeros(pd['n'])
            pd['fwhm'] = np.array(self.groups[pd['plan'][0]].meta['GUIDFWHM'])[pd['ia']]/0.1799
            pd['scale'] = np.ones(pd['n'])
            pd['scale_weight'] = np.ones(pd['n'])

            return pd['fwhm'], pd['shift'], pd['scale']

        if len(bad) > 0:
            for i in bad[::-1]:
                msg = f'{self.namestr} # Remove bad exposure from list: '
                msg += f'{i} fwhm={fit_fwhm[i]:.2f}'
                msg += f' shift={fit_x[i]:.2f} scale={fit_sum[i]:.2f}'
                utils.log_comment(self.logfile, msg, verbose=True)
                
                #print('xxx', i, len(pd['ia']), len(pd['ib']))
                pd['ia'].pop(i)
                pd['ib'].pop(i)
                
            #pd['ia'] = pd['ia'][ok]
            #pd['ib'] = pd['ib'][ok]
            #pd['ta'] = pd['ta'][ok]
            pd['n'] =  len(pd['ia'])
        
        #print('xxx set fwhm')       
        pd['fwhm'] = np.maximum(fit_fwhm[ok], 1.1)
        pd['shift'] = fit_x[ok] - slit_info['target_y'] #fit_x[ok][0]
        pd['scale'] = fit_sum[ok]
        
        #self.plan_pairs[plan_i] = pd
        
        return fit_fwhm, fit_x, fit_sum


    def get_target_drift(self, slit, use_peak=True, profile_model=None):
        """
        Get drifts of all offset plans and make a figure
        """
        fig, axes = plt.subplots(1,2,figsize=(14,5))

        slit_info, fig = self.get_slit_params(slit=slit, xy_order=3, pad=16, 
                                              show=False)
        
        for i in range(len(self.plans)):
            self.get_plan_drift(slit_info, plan_i=i, ax=axes[0], 
                                use_peak=use_peak, 
                                profile_model=profile_model)

            td = self.plan_pairs[i]
            gra = self.groups[self.plans[i][0]]
            airm = np.array(gra.meta['AIRMASS'])[td['ia']]
            guid = np.array(gra.meta['GUIDFWHM'])[td['ia']]
            
            ia = td['ia']
            if i == 0:
                axes[1].plot(td['ta'][ia], td['fwhm']*0.1799, label='FWHM', color='b')
                axes[1].plot(td['ta'][ia], guid, label='Guider FWHM', color='lightblue')
                axes[1].plot(td['ta'][ia], td['shift'], label='Shift', color='g')
                axes[1].plot(td['ta'][ia], td['scale']/td['scale'].max(), label='Scale', color='orange')
                axes[1].plot(td['ta'][ia], airm, label='AIRMASS', color='r')
            else:
                axes[1].plot(td['ta'][ia], td['fwhm']*0.1799, color='b')
                axes[1].plot(td['ta'][ia], guid, color='lightblue')
                axes[1].plot(td['ta'][ia], td['shift'], color='g')
                axes[1].plot(td['ta'][ia], td['scale']/td['scale'].max(), color='orange')
                axes[1].plot(td['ta'][ia], airm, color='r')

        for ax in axes:
            ax.grid()
        
        txt = '{0} slit#{1} {2}'.format(self.namestr, 
                                        slit_info['target_orig_slit'],
                                        slit_info['target_name'])
                                        
        axes[0].text(0.05, 0.95, txt, ha='left', va='top', 
                     transform=axes[0].transAxes)
        
        xc = slit_info['target_y']
        axes[0].vlines([xc-10, xc+10], *axes[0].get_ylim(), color='k',
                       alpha=0.5, linestyle='--')

        axes[1].legend()
        
        return fig


    def align_bright_star(self, mag_lim=20.5, use_peak=False):
        mag = self.ssl['Magnitude']*1
        mag[mag < 0] = 99
        slit = np.argmin(mag)
    
        mag_ok = (self.ssl['Magnitude'][slit] > 10) 
        mag_ok &= (self.ssl['Magnitude'][slit] < mag_lim)
        
        if mag_ok:
        
            msg = f'{self.namestr} # Align on bright object on slit {slit} ("{self.target_names[slit]}") m={mag[slit]:.1f}'
            utils.log_comment(self.logfile, msg, verbose=True)
        
            fig = self.get_target_drift(slit, use_peak=use_peak)
            plt.gcf().savefig(os.path.join(self.path, f'{self.datemask}-{self.filter}_drift.png'))
            
            key = self.target_keys[slit]
            
            try:
                show = 1
                slit_info, fig = self.get_slit_params(slit=slit, xy_order=3, 
                                                      pad=16, show=show)
                if show:
                    fig.savefig(os.path.join(self.path, f'{key}_slit.png'))

                print('Drizzle')

                hdu = self.drizzle_all_plans(slit_info, kernel='point', 
                                             sig_clip=(3, 3), 
                                             linearize_wave=True,
                                             mask_single=False, 
                                             mask_overlap=False, 
                                             mask_offset=False)

                fig = show_drizzled_hdu(hdu, vm=None, vmp=[0.5,98], 
                                  vm1=1.1, smooth=1, ny=3,
                                  cmap=plt.cm.inferno, 
                                  zoom=20, wpower=0.15, use_sn=False, 
                                  clip_wave=0)

            except TypeError:
                pass
            except IndexError:
                pass


        else:
            msg = f'{self.namestr} # No bright stars found, use GUIDFWHM header keywords'
            utils.log_comment(self.logfile, msg, verbose=True)

            for pd in self.plan_pairs:
                pd['shift'] = np.zeros(pd['n'])
                pd['fwhm'] = np.array(self.groups[pd['plan'][0]].meta['GUIDFWHM'])[pd['ia']]/0.1799
                pd['scale'] = np.ones(pd['n'])
        
        self.plan_pairs_info()


    def reset_plan_shifts(self):
        """
        Reset shifts, FWHM, scale of plan exposures
        """
        msg = f'{self.namestr} # Reset plan shifts and scaling'
        utils.log_comment(self.logfile, msg, verbose=True)

        for pd in self.plan_pairs:
            pd['shift'] = np.zeros(len(pd['shift']), dtype=np.float32)
            pd['scale'] = np.ones(len(pd['shift']), dtype=np.float32)
            pd['fwhm'] = np.array(self.groups[pd['plan'][0]].meta['GUIDFWHM'])[pd['ia']]/0.1799
    
    
    def pipeline(self, initial_thresh=500, max_iter=5, align_peak=False, tpad=60, use_ssl_slits=False, **kwargs):
        """
        Run the full pipeline to extract spectra
        """
        self.pipeline_status = False

        # Find slits
        try:
            fig = self.find_slits(x0=1024, initial_thresh=initial_thresh,
                                  use_ssl=use_ssl_slits,
                                  max_iter=max_iter, interpolate_coeffs=True)
        except IndexError:
            try:
                fig = self.find_slits(x0=1024, initial_thresh=1800, 
                                  interpolate_coeffs=True)
            except:
                fig = self.find_slits(x0=1024, initial_thresh=1800, 
                                  use_ssl=True, 
                                  interpolate_coeffs=True)
        except ValueError:
            fig = self.find_slits(x0=1024, initial_thresh=1800, 
                              use_ssl=True, 
                              interpolate_coeffs=True)
                    
        fig.savefig(os.path.join(self.path,
                                 f'{self.datemask}-{self.filter}_slits.png'))
        
        # Region file
        try:
            regs, fig = self.make_region_file(make_figure=True)
        except:
            utils.log_exception(self.logfile, traceback)
            return False
            
        # Find plans
        try:
            plan_pairs, fig = self.get_plan_pairs(tpad=tpad, show=True)
        except:
            utils.log_exception(self.logfile, traceback)
            return False
            
        fig.savefig(os.path.join(self.path, 
                                 f'{self.datemask}-{self.filter}_plan.png'))
        
        # Drifts
        try:
            self.align_bright_star(use_peak=align_peak)
        except:
            utils.log_exception(self.logfile, traceback)
            plan_pairs, fig = self.get_plan_pairs(tpad=tpad, show=True)
            self.reset_plan_shifts()
    
        self.pipeline_status = True


    def extract_all_slits(self, slit_numbers=None, orig_slit_numbers=None, show_trace=False, show_drizzled=False, pad=8, sig_clip=(3, 3), mask_offset=False, mask_trace=False, mask_overlap=False, mask_single=False, kernel='point', pixfrac=1, linearize_wave=True, log_wave=False, zoom=35, save_full_drizzled=False, **kwargs):    
        """
        Extract 2D spectra for all slits
        """
        if slit_numbers is None:
            slit_numbers = range(self.nslits)
        
        if orig_slit_numbers is not None:
            slit_numbers = [self.target_slit_numbers.index(s)
                            for s in orig_slit_numbers]
                            
        for slit in slit_numbers:
            key = self.target_keys[slit]
            try:                
                slit_info, fig = self.get_slit_params(slit=slit, xy_order=3, 
                                                      pad=pad, 
                                                      show=show_trace)
                if show_trace:
                    fig.savefig(os.path.join(self.path, f'{key}_slit.png'))
                    
            except TypeError:
                utils.log_exception(self.logfile, traceback)
                continue
            except IndexError:
                utils.log_exception(self.logfile, traceback)
                continue
            except SystemError:
                utils.log_exception(self.logfile, traceback)
                continue

            try:
                hdu = self.drizzle_all_plans(slit_info, 
                                         sig_clip=sig_clip,
                                         linearize_wave=True,
                                         mask_trace=mask_trace, 
                                         mask_single=mask_single,
                                         mask_overlap=mask_overlap, 
                                         mask_offset=mask_offset, 
                                         kernel=kernel, 
                                         pixfrac=pixfrac, 
                                         log_wave=log_wave)
            except:
                utils.log_exception(self.logfile, traceback)
                continue
                
            self.slit_hdus[slit] = hdu
            self.slit_info[slit] = slit_info

            file = os.path.join(self.path, f'{key}_sp.fits')
            msg = f'{self.namestr} # slitidx {slit} / {self.nslits} {file}' 
            utils.log_comment(self.logfile, msg, verbose=True)
            self.slit_hdus[slit].writeto(file, overwrite=True)
            
            if save_full_drizzled:
                xfile = os.path.join(self.path, f'{key}_xsp.fits')

                xhdu = pyfits.HDUList()
                xhead, xsci, xwht = self.drizzled_plans[0]
                xhead['EXTNAME'] = 'SCI'
                xhdu.append(pyfits.PrimaryHDU(data=xsci, header=xhead))
                xhead['EXTNAME'] = 'WHT'
                xhdu.append(pyfits.ImageHDU(data=xwht, header=xhead))

                xmsg = f'{self.namestr} # Full drizzled file {xfile}' 
                utils.log_comment(self.logfile, xmsg, verbose=True)
                xhdu.writeto(xfile, overwrite=True)
                del(xhdu)
                
            if show_drizzled:
                try:
                    fig = show_drizzled_hdu(hdu, vm=None,
                                  vmp=[0.5,99], vm1=1.5, smooth=1, ny=1,
                                  cmap=plt.cm.Spectral_r,
                                  zoom=zoom, wpower=0.15, use_sn=True, 
                                  clip_wave=0, spix=(6, 2))
                    fig.savefig(file.replace('_sp.fits', '_sp2d.png'))
                    
                    if zoom > 0:
                        fig = show_drizzled_hdu(hdu, vm=None,
                                      vmp=[0.5,99], vm1=1.5, smooth=1, ny=1,
                                      cmap=plt.cm.Spectral_r, 
                                      ys=4,
                                      zoom=-1, wpower=0.15, use_sn=True, 
                                      clip_wave=0, spix=(6, 2))
                        fig.savefig(file.replace('_sp.fits', 
                                                 '_sp2d_full.png'))
                        
                except IndexError:
                    continue
                    
                
        return True


class ExposureGroup(object):
    """
    Group of exposures at a given dither position
    """
    def __init__(self, offset_file='HighzT7_20130117/Reduced/HighzT7/2013jan17/Y/Offset_1.25.txt', min_nexp=4, verbose=True, flat=1):
        
        self.logfile = utils.LOGFILE
        
        msg = f'==========\nRead exposures: {offset_file}\n==========='
        utils.log_comment(self.logfile, msg, verbose=True)
        
        key = os.path.basename(offset_file)
        self.offset_file = offset_file
        with open(offset_file) as fp:
            lines = fp.readlines()

        self.rawdir = lines[1].split()[0]
        self.img = []
        self.minmax_nrej = 0
        self.truitime = []
        self.frameid = None
        self.flat = flat
        
        self.files = []
        for i in range(2,len(lines)):
            exp_file = lines[i].split()[0]
            self.files.append(exp_file)
        
        if self.nexp > min_nexp:
            self.read_fits(verbose=verbose)
        else:
            utils.log_comment(self.logfile,
                              f'Skip: nexp {self.nexp} < {min_nexp}', 
                              verbose=True)


    @property 
    def nexp(self):
        """
        Number of exposures
        """
        return len(self.files)
    
    
    def read_fits(self, verbose=True):
        """
        Read exposure FITS files
        """
        for i, exp_file in enumerate(self.files): 
            img = pyfits.open(self.rawdir + '/' + exp_file)
            #offset_groups[key]['img'].append(img)
            self.img.append(img)
            
            h = img[0].header
            msg = f'{i:2} {exp_file} '
            
            h['AIRM'] = '{0:.2f}'.format(h['AIRMASS'])
            
            for k in ['FILTER', 'XOFFSET', 'YOFFSET', 'FRAMEID', 'FRAME',
                      'REPEAT', 'REPEATS', 'AIRM', 'GUIDFWHM']:
                if k in h:
                    msg += f' {h[k]}'
                else:
                    msg += ' *'
                                                 
            utils.log_comment(self.logfile, msg, verbose=verbose)
                
        self.N = len(self.img)

        self.sci = np.zeros((self.N, 2048, 2048))
        self.var = self.sci*0.
        self.RN = 21

        for i in range(self.N):
            img = self.img[i]
            itime = img[0].header['TRUITIME']
            self.sci[i,:,:] = img[0].data * img[0].header['SYSGAIN'] / itime / self.flat
            RN_i = self.RN**2 / img[0].header['NUMREADS']
            self.var[i,:,:] = np.maximum(self.sci[i,:,:], 0) * itime + RN_i
            self.var[i,:,:] /= itime**2
        
        self.var[:,:,:16] = 0
        self.var[:,:,-16:] = 0
        
        self.frameid = img[0].header['FRAMEID'].replace("'", 'p')
        
        self.truitime = np.array([img[0].header['TRUITIME']*1.
                                  for img in self.img], dtype=np.float32)
        self.xoffset = np.array([img[0].header['XOFFSET']*1.
                                 for img in self.img], dtype=np.float32)
        self.yoffset = np.array([img[0].header['YOFFSET']*1.
                                 for img in self.img], dtype=np.float32)
        
        self.fit_fwhm = np.ones(self.N)
        self.fit_yoffset = np.zeros(self.N)


    @property
    def meta(self):
        """
        metadata
        """        
        meta = {}
        for i, img in enumerate(self.img):
            if i == 0:
                for k in OBJ_KEYS:
                    meta[k] = img[0].header[k]
                
                for k in SEQ_KEYS:
                    if k in img[0].header:
                        meta[k] = [img[0].header[k]]
                    else:
                        meta[k] = [1.0]
            else:
                for k in SEQ_KEYS:
                    if k in img[0].header:
                        meta[k].append(img[0].header[k])
                    else:
                        meta[k].append(1.0)
                        
        
        return meta


    @property 
    def exptime(self):
        """
        Total exposure time
        """
        return self.truitime.sum()


    def flag_cosmic_rays(self, sigma=5, use_bpix=True, verbose=True, grow_bpm=1, minmax_nrej=-1, **kwargs):
        """
        Sigma clipping + minmax rejection to clean up exposure arrays
        """
        
        clip = ~np.isfinite(self.sci)
        
        minmax_clip = clip & True
        
        if minmax_nrej < 0:
            if self.nexp < 4:
                minmax_nrej = 0
            elif self.nexp < 10:
                minmax_nrej = 1
            elif self.nexp < 20:
                minmax_nrej = 2
            else:
                minmax_nrej = 4

        msg = f'{self.namestr}: flag crs, sigma={sigma}, minmax_nrej={minmax_nrej}, use_bpix={use_bpix}'
        utils.log_comment(self.logfile, msg, verbose=verbose)
        
        if minmax_nrej > 0:                
            minmax_mask = self.sci*1.
            for it in range(minmax_nrej):
                mi = np.nanmin(minmax_mask, axis=0)
                ma = np.nanmax(minmax_mask, axis=0)
                minmax_mask[(minmax_mask <= mi) | (minmax_mask >= ma)] = np.nan
            
            minmax_clip |= ~np.isfinite(minmax_mask)
        
        self.minmax_nrej = minmax_nrej
        
        if use_bpix:
            bpix = pyfits.open('badpix_10sep2012.fits')
            bpm = bpix[0].data > 0
            if grow_bpm:
                bpm = nd.binary_dilation(bpm, iterations=grow_bpm*1)

            #bpix.info()
            clip |= bpm
        
        # With minmax clipping
        self.ivar = 1./self.var
        self.ivar[clip | minmax_clip | (self.var <= 0)] = 0
        num = np.nansum(self.sci*self.ivar, axis=0)
        den = np.nansum(self.ivar, axis=0)
        
        avg = num/den
        avg[den == 0] = 0
        
        # Do the sigma clipping
        clip |= np.abs(self.sci - avg)*np.sqrt(self.ivar) > sigma
        
        # Recompute average
        self.ivar[clip] = 0
        num = np.nansum(self.sci*self.ivar, axis=0)
        den = np.nansum(self.ivar, axis=0)
        
        avg = num/den
        avg[den <= 0] = 0
        
        avg_msk = (den <= 0) | ~np.isfinite(num+den)
        
        self.avg_sci = avg
        self.avg_wht = den
        self.avg_var = 1/self.avg_wht
        self.grow_bpm = grow_bpm

        self.avg_sci[avg_msk] = 0
        self.avg_wht[avg_msk] = 0
        self.avg_var[avg_msk] = 0

        self.avg_npix = (self.avg_wht > 0).sum(axis=0)

    @property 
    def namestr(self):
        """
        Descriptive identifier
        """
        return (f'{self.frameid:>5} {self.N:>3} '
                f'{os.path.basename(self.offset_file)}')
    
    
    def __repr__(self):
        return self.namestr


def show_drizzled_hdu(hdu, header=None, xs=14, ys=2, vm=None, vm1=1, vmp=[2,98.9], smooth=False, spix=(1,1), ny=5, zoom=25, use_sn=True, wpower=0.15, clip_wave=True, use_title=False, cmap=plt.cm.viridis):
    """
    Make figure of a drizzled spectrum
    """
    from astropy.modeling.models import Gaussian2D
    import astropy.wcs as pywcs
    
    from scipy.signal import fftconvolve
    
    if hasattr(hdu, 'info'):
        outsci = hdu['SCI'].data*1
        outwht = hdu['WHT'].data*1
        if header is None:
            header = hdu['SCI'].header
    else:
        outsci, outwht = hdu
        
    diff = outsci*1. #(outsci - first_pos)/2.
    dwht = outwht*1.
    if wpower is not None:
        wp = np.nanpercentile(dwht[dwht > 0], 98)
        wscl = (np.minimum(dwht, wp)/wp)**wpower
        dwht *= wscl
            
    if smooth:
        g2 = Gaussian2D(x_mean=0, y_mean=0, x_stddev=spix[0], y_stddev=spix[1])

        sx = np.ceil(5*spix[0])
        sy = np.ceil(5*spix[1])
        xar = np.arange(-sx, sx+0.1)
        yar =np.arange(-sy, sy+0.1)
        gx, gy = np.meshgrid(xar, yar)
        kern = g2(gx, gy)
        kern /= kern.sum()
        
        #dnum = nd.gaussian_filter(diff*dwht, 1) 
        #dden = nd.gaussian_filter(dwht, 1) 
        dnum = fftconvolve(diff*dwht, kern, mode='same')
        dden = fftconvolve(dwht, kern**2, mode='same')
        
        if use_sn:
            diff = dnum/np.sqrt(dden)
    else:
        if use_sn:
            diff = diff*np.sqrt(dwht)
    
    diff[dwht == 0] = 0
    
    #wave = slit_info['wave']/1.e4
    
    sh = outsci.shape
    xarr = np.arange(sh[1])
    yp = np.ones(sh[1]) + sh[0]/2
    lam_wcs = pywcs.WCS(header)
    wave, _y = lam_wcs.all_pix2world(xarr, yp, 0)
    #wave *= 1.e6 # microns
    if wave.max() < 1.e-3:
        wave *= 1.e6 # microns
    elif wave.max() > 1000:
        wave *= 1.e-4
    
    yi = np.clip(np.interp([header['WAVEBLUE'], header['WAVERED']], wave, xarr), 0, sh[1])
    #print('interp: ', yi, header['WAVEBLUE'], header['WAVERED'], wave.min(), wave.max())
    
    if clip_wave:
        sl = slice(*np.cast[int](yi))
    else:
        sl = slice(0, sh[1])
            
    wave = wave[sl]
    diff = diff[:,sl]
    xarr = xarr[sl]
    
    #yoff = header['CRPIX2']
    
    #yoff = slit_info['target_offset']/0.1799
    sh = diff.shape

    #half = slit_info['istop'] - slit_info['istart']
    #center = slit_info['istart'] - slit_info['i0'] + half/2. + yoff
    #half = header['YSTOP'] - header['YSTART']
    #center = header['YSTART'] - header['Y0'] + half/2. + yoff
    center = header['CRPIX2']
    
    yarr = np.arange(sh[0])
    sly = np.abs(yarr-center) < 10
    
    if vm is None:
        vm = np.nanpercentile(diff[sly,:], vmp)
    
    if use_sn:
        # print('xx', vm[1])
        vm[1] = np.clip(vm[1], 5, 20)
        
    vm[1] *= vm1
    vm[0] = -0.8*vm[1]
    
    print(header['SLITNUM'], header['TARGNAME'], vm)
    
    fig, axes = plt.subplots(ny, 1, figsize=(xs, ys*ny))
    
    if ny == 1:
        axes = [axes]
        
    wsplit = np.linspace(wave.min(), wave.max(), ny+1)
    for ia, ax in enumerate(axes):

        ax.imshow(diff, #extent=(wave.min(), wave.max(), 0, diff.shape[0]), 
                  cmap=cmap, vmin=vm[0], vmax=vm[1], origin='lower', 
                  interpolation='Nearest')
        
        ax.set_aspect('auto')
        
        xl = ax.get_xlim()
        yl = ax.get_ylim()
        if ia == 0:
            label = f"{header['DATEMASK']}  {header['FILTER']}  {header['SLITNUM']}  {header['TARGNAME']}"
            if use_title:
                ax.set_title(label)
            else:
                left = f"{header['DATEMASK']}  {header['FILTER']}"
                right = f"{header['SLITNUM']}  {header['TARGNAME']}"
                kws = dict(transform=ax.transAxes,
                           bbox={'facecolor':'w', 'edgecolor':'None', 
                                 'alpha':1.0}, 
                           fontsize=8)
                
                ax.text(0.005, 0.03, left, ha='left', va='bottom', **kws)
                ax.text(0.995, 0.03, right, ha='right', va='bottom', **kws)
                        
        ax.hlines(center + np.array([-11.1,11.1]), *ax.get_xlim(), 
                  color='k', alpha=0.5, linestyle=':')

        if zoom > 0:
            ax.set_ylim(center + np.array([-zoom,zoom]))
        
        ax.set_yticks([])
        ax.set_yticklabels([])
        if 0:
            ax.text(0.006, center/sh[0], f'{center:.1f}', ha='left', va='center', 
                    fontsize=8, transform=ax.transAxes,
                    bbox=dict(facecolor='w', edgecolor='None', alpha=0.9))
                    
        if ny == 1:
            if header['FILTER'] == 'Y':
                # ax.xaxis.set_major_locator(MultipleLocator(0.02))
                # ax.xaxis.set_minor_locator(MultipleLocator(0.005))
                xv = np.arange(0.965, 1.136, 0.005)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi, minor=True)
                xv = np.arange(0.98, 1.121, 0.02)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi)
                #print(xi, wave.min(), wave.max())
                ax.set_xticklabels([f'{f:.2f}' for f in xv])
            
            elif header['FILTER'] == 'J':
                xv = np.arange(1.15, 1.355, 0.01)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi, minor=True)
                xv = np.arange(1.16, 1.34, 0.02)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi)
                #print(xi, wave.min(), wave.max())
                ax.set_xticklabels([f'{f:.2f}' for f in xv])
                
                #ax.xaxis.set_major_locator(MultipleLocator(0.02))
                #ax.xaxis.set_minor_locator(MultipleLocator(0.01))
                
            elif header['FILTER'] == 'H':
                # ax.xaxis.set_major_locator(MultipleLocator(0.05))
                # ax.xaxis.set_minor_locator(MultipleLocator(0.01))
                xv = np.arange(1.46, 1.811, 0.01)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi, minor=True)
                xv = np.arange(1.5, 1.81, 0.05)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi)
                #print(xi, wave.min(), wave.max())
                ax.set_xticklabels([f'{f:.2f}' for f in xv])

            elif header['FILTER'] == 'K':
                # ax.xaxis.set_major_locator(MultipleLocator(0.05))
                # ax.xaxis.set_minor_locator(MultipleLocator(0.01))
                xv = np.arange(1.90, 2.41, 0.01)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi, minor=True)
                xv = np.arange(1.90, 2.41, 0.05)
                xi = np.interp(xv, wave, xarr)                
                ax.set_xticks(xi)
                #print(xi, wave.min(), wave.max())
                ax.set_xticklabels([f'{f:.2f}' for f in xv])
        
        xpi = np.interp((wsplit[ia], wsplit[ia+1]), wave, xarr)
        ax.set_xlim(*xpi)
        ax.grid(axis='x')
        
        ##print('xx', ax.get_xlim(), wave.min(), wave.max())
        
    fig.tight_layout(pad=0.1)
    
    return fig


def drizzled_profile(hdu, wpower=0.15, ax=None, plot_kws={}, clip_wave=True, wlim=None):
    """
    Cross-dispersion profile of a drizzled spectrum
    """
    import astropy.wcs as pywcs
    
    if hasattr(hdu, 'info'):
        outsci = hdu['SCI'].data*1
        outwht = hdu['WHT'].data*1
        head = hdu['SCI'].header
    else:
        head, outsci, outwht = hdu
        
    diff = outsci*1. #(outsci - first_pos)/2.
    dwht = outwht*1.
    if wpower is not None:
        wp = np.nanpercentile(dwht[dwht > 0], 98)
        wscl = (np.minimum(dwht, wp)/wp)**wpower
        dwht *= wscl

    sh = outsci.shape
    xarr = np.arange(sh[1])

    yp = np.ones(sh[1]) + sh[0]/2
    lam_wcs = pywcs.WCS(head)
    wave, _y = lam_wcs.all_pix2world(xarr, yp, 0)

    if wave.max() < 1.e-3:
        wave *= 1.e6 # microns
    elif wave.max() > 1000:
        wave *= 1.e-4
        
    wlim0 = [head['WAVEBLUE'], head['WAVERED']]
    yi = np.clip(np.interp(wlim0, wave, xarr), 0, 2048)
    
    if wlim is not None:
        yi = np.clip(np.interp(wlim, wave, xarr), 0, 2048)
        
    if clip_wave:
        sl = slice(*np.cast[int](yi))
    else:
        sl = slice(0, sh[1])
    
    den = dwht[:,sl].sum(axis=1)
    yprof = (diff*dwht)[:,sl].sum(axis=1)/den
    yprof[den == 0] = 0
    
    
    if ax is None:
        ax = plt
    
    #slit_info['drizzled_wave'] = wave
    ax.plot(yprof, **plot_kws)
    #ax.plot(1/np.sqrt(den), **plot_kws)

    return yprof


filt_params = {
    'Y': {'filter':'Y','df':11, 'fwhm':0.0025, 'dv':0}, 
    'J': {'filter':'J','df':13, 'fwhm':0.0028, 'dv':0}, 
    'H': {'filter':'H','df':11, 'fwhm':0.0038, 'dv':10}, 
    'K': {'filter':'K','df':11, 'fwhm':0.008, 'dv':0},     
}


def fit_vega_spectrum(hip_wave, hip_flux, hip_err, hip_mag=8.503, filter='K', df=11, fwhm=0.0025, dv=0, ax=None, make_figure=True, mask=None):
    
    import astropy.units as u
    from synphot import SourceSpectrum, SpectralElement, units
    from scipy.modeling import models
    
    vega = SourceSpectrum.from_vega()
    
    flam = 1*u.erg/u.second/u.erg/u.cm**2/u.Angstrom
    
    hip_vega = vega(hip_wave*u.micron).to(flam, equivalencies=u.spectral_density(hip_wave*u.micron)) 
    
    xmed = np.median(hip_wave)
    
    hip_scale = np.nanmax(hip_flux)/np.interp(xmed, hip_wave, hip_vega)

    xclip = np.isfinite(hip_wave + hip_err) & (hip_flux > 0.05*np.nanmax(hip_flux))
    
    ix = np.where(xclip)[0]
    xclip[ix[:4]] = False
    xclip[ix[-4:]] = False
    
    if mask.filter == 'K':
        xclip &= np.abs(hip_wave-1.96) > 0.01
        xclip &= np.abs(hip_wave-2.01) > 0.02
        xclip &= np.abs(hip_wave-2.063) > 0.02
    elif mask.filter == 'Y':
        xclip &= hip_wave < 1.125
        
    print('NPIX', xclip.sum())
    bspl = utils.bspline_templates(wave=hip_wave, degree=3, df=df, get_matrix=True)

    #linew = np.array([getattr(line, f'x_0_{i}').value for i in range(nline)])

    #Rh = 1.0973731568160e7
    Rh = 1.09678e7

    # Lines for H recombination
    if mask.filter == 'Y':
        n = np.arange(6,9)
        nr = 3
    elif mask.filter == 'J':
        n = np.arange(5,6)
        nr = 3
    elif mask.filter == 'H':
        n = np.arange(9, 20)
        nr = 4
    else:
        n = np.arange(7, 9)
        nr = 4

    linew = 1./(Rh * (1./nr**2 - 1/n**2))/1.e-6
    
    #linew = np.append(linew, [1.0675, 1.0828])
    
    lines = []
    vlines = []

    for li in linew:
        if (li > hip_wave[xclip].min()) & (li < hip_wave[xclip].max()):
            print(f'Add line: {li}')
            lm = models.Lorentz1D(x_0=li*(1+dv/3.e5), fwhm=fwhm, amplitude=-0.03)
            lines.append(lm(hip_wave))
            lm = models.Lorentz1D(x_0=li*(1+dv/3.e5), fwhm=0.004, amplitude=-0.03)
            vlines.append(lm(hip_wave))

    templ = np.vstack([bspl.T, np.array(lines)])
    vtempl = np.vstack([bspl.T, np.array(vlines)])
    p = 0

    fcorr = 1. #resamp #* mfi

    _x = np.linalg.lstsq(templ[:,xclip].T, (hip_flux/(fcorr)**p)[xclip])
    _m = templ.T.dot(_x[0])

    nb = bspl.shape[1]
    _mc = templ.T[:,:nb].dot(_x[0][:nb])
    _ml = templ.T[:,nb:].dot(_x[0][nb:])

    _xv = np.linalg.lstsq(vtempl[:,xclip].T, (hip_vega)[xclip])
    _mv = vtempl[:,xclip].T.dot(_xv[0])
    _mcv = vtempl.T[:,:nb].dot(_xv[0][:nb])

    hip_model = _mcv*(1+_ml/_mc)

    if make_figure:
        if ax is None:
            fig, ax = plt.subplots(1,1,figsize=(14,8))

        ax.plot(hip_wave, (hip_flux/(fcorr)**p), color='lightblue', alpha=0.5)

        ax.plot(hip_wave[xclip], (hip_flux/(fcorr)**p)[xclip], color='b', alpha=0.5)
        ax.plot(hip_wave[xclip], ((hip_flux/(fcorr)**p)[xclip] - _m[xclip] ), color='k', alpha=0.3)

        ax.plot(hip_wave[xclip], _m[xclip], color='pink', alpha=0.5)
        ax.plot(hip_wave[xclip], _mc[xclip], color='brown', alpha=0.5)
        ax.plot(hip_wave[xclip], _mcv[xclip]*hip_scale, color='brown', alpha=0.5)
        ax.plot(hip_wave, hip_vega*hip_scale, color='pink', alpha=0.5)
    
    mag_scl = 10**(-0.4*hip_mag)
    
    to_flam = hip_model.to(1.e-17*u.erg/u.second/u.cm**2/u.Angstrom)*mag_scl/hip_flux
    to_flam[~xclip] = 0
    
    return to_flam


def run_mask(flat, skip=True, **kwargs):
    import time
    
    path = os.path.dirname(flat)
    
    datemask = path.split('/')[-5]
    filt = path.split('/')[-1]
    
    if os.path.exists(f'{datemask}-{filt}.log'):
        print(f'Extractions found for {datemask}-{filt}.log')
        if skip:
            return True
        
    mask = MosfireMask(path=path, min_nexp=0)
    
    mask.pipeline(align_peak=True, **kwargs)
    
    if len(mask.plans) == 0:
        msg = f'{mask.namestr} # No plans found'
        utils.log_comment(mask.logfile, msg, verbose=True)
        return False
        
    mask.extract_all_slits(show_drizzled=True, zoom=35, **kwargs)
    
    plt.close('all')
    
    return mask
    
    
def dump_all():
        
    with open('source_table/counts.txt','w') as fp:
        fp.write('# i n ncols file\n')
    
    istop = -1
    with open('new_flats.txt') as fp:
        lines = fp.readlines()
    
    N = 20
    istop = -1
    os.system('ls -ltr *_20*log | sed "s/-[YJHK].log//" | awk \'{print $9}\' | tail -NN > logfiles'.replace('NN',f'{N}'))        
    with open ('logfiles') as fp:
        lines = fp.readlines()
    
    flatfiles = []
    for line in lines:
        flatfiles.extend( glob.glob(f'{line.strip()}/Reduced/*/20*/[YJHK]/combflat_2d_[YJHK].fits'))
        
    #flatfiles = [l.strip() for l in lines]
    for istop in range(istop+1,len(flatfiles)):
        file=flatfiles[istop]
        if 'lamps_o' in file:
            continue
            
        #print(istop, file)
        try:
            tab = dumpheaders(file)
        except:
            print(f'Failed to read {file}')
            continue
            
        if hasattr(tab, 'colnames'):
            msg = f'{istop:4} {len(tab):3} {len(tab.colnames):3} {file}'
        else:
            msg = f'{istop:4} {0:3} {60:3} {file}'
        
        #print(msg)
        
        utils.log_comment('source_table/counts.txt', msg, verbose=True)
    
    
    bash = """
    wc mosfire_extractions.ecsv
    cat source_table/*ecsv > all_raw.ecsv
    grep -e "#" -e "file" source_table/1138_final2_20150126-J.ecsv > mosfire_extractions.ecsv
    grep -v -e "long2pos" -e "#" -e "file" all_raw.ecsv >> mosfire_extractions.ecsv
    wc mosfire_extractions.ecsv
    """
    os.system(bash)
    
def read_all():
    """
    cd /home/idies/workspace/Temporary/gbrammer/scratch/MOSFIRE
    
    cat source_table/*ecsv > all_raw.ecsv
    grep -e "#" -e "file" source_table/1138_final2_20150126-J.ecsv > mosfire_extractions.ecsv
    grep -v -e "long2pos" -e "#" -e "file" all_raw.ecsv >> mosfire_extractions.ecsv
    
    for file in `ls mosfire_[me][ax][st]*`; do 
        echo $file
        aws s3 cp ${file} s3://grizli-v1/MosfireSpectra/ --acl public-read
    done
    
    """
    istart = 0
        
    files = glob.glob('source_table/*ecsv')
    
    for istart in range(istart, len(files))[:5]:
        tab = utils.read_catalog(files[istart], format='ascii.ecsv')
        msg = f'{istart:4} {len(tab):3} {len(tab.colnames):3} {files[istart]}'
        print(msg)
        #utils.log_comment('source_table/counts.txt', msg, verbose=True)
        
def dumpheaders(flatfile):
    import os
    import glob
    from grizli import utils
    import astropy.io.fits as pyfits
    import numpy as np
    
    path=os.path.dirname(flatfile)
    files = glob.glob(os.path.join(path, '*sp.fits'))
    if len(files) == 0:
        with open('empty.log','a') as fp:
            fp.write(flatfile+'\n')
        
        return False
    else:
        with open('notempty.log','a') as fp:
            fp.write(f'{flatfile} {len(files)}\n')
        
    rows = []
    for i, file in enumerate(files):
        im = pyfits.open(file)
        if i == 0:
            keys = [k.lower() for k in im[0].header]
        row = [im[0].header[k] for k in im[0].header]
        #print(i, len(keys), len(row))
        rows.append(row)
    
    # ik = keys.index('observer')
    # for r in rows:
    #     r[ik] = r[ik].replace(',','')
        
    tab = utils.GTable(rows=rows, names=keys)
    cols = ['file', 'naxis1', 'naxis2', 'crpix1', 'crpix2', 'crval1', 'crval2', 'cd1_1', 'cd2_2', 'ctype1', 'ctype2', 'sigclipn', 'sigclipv', 'exptime', 'nexp', 'kernel', 'pixfrac', 'offseta', 'offsetb', 'targname', 'ra_slit', 'dec_slit', 'ra_targ', 'dec_targ', 'mag_targ', 'filter', 'datemask', 'slitidx', 'slitnum', 'y0', 'y1', 'ystart', 'ystop', 'ypad', 'maskoff', 'targoff', 'targypix', 'lamorder', 'lamcoef0', 'lamcoef1', 'object', 'maskname', 'observer', 'pattern', 'dthcoord', 'skypa3', 'progid', 'progpi', 'progtl1', 'semester', 'wavered', 'waveblue', 'airmass_min', 'airmass_med', 'airmass_max', 'guidfwhm_min', 'guidfwhm_med', 'guidfwhm_max', 'mjd_obs_min', 'mjd_obs_med', 'mjd_obs_max']
    
    cols = ['file', 'naxis1', 'naxis2', 'crpix1', 'crpix2', 'crval1', 'crval2', 'cd1_1', 'cd2_2', 'sigclipn', 'sigclipv', 'exptime', 'nexp', 'kernel', 'pixfrac', 'plan', 'offseta', 'offsetb', 'targname', 'ra_slit', 'dec_slit', 'ra_targ', 'dec_targ', 'mag_targ', 'filter', 'datemask', 'slitidx', 'slitnum', 'y0', 'y1', 'ystart', 'ystop', 'ypad', 'cunit2', 'ctype1', 'cunit1', 'maskoff', 'targoff', 'targypix', 'traorder', 'tracoef0', 'tracoef1', 'tracoef2', 'lamorder', 'lamcoef0', 'lamcoef1', 'object', 'framdesc', 'maskname', 'observer', 'pattern', 'dthcoord', 'skypa3', 'progid', 'progpi', 'progtl1', 'semester', 'wavered','waveblue', 'airmass_min', 'airmass_med', 'airmass_max', 'guidfwhm_min', 'guidfwhm_med', 'guidfwhm_max', 'mjd_obs_min', 'mjd_obs_max', 'mjd_obs']
    
    tab['file'] = files
    
    so = np.argsort(tab['slitnum'])
    for k in tab.colnames:
        if tab[k].dtype == np.float64:
            tab[k] = tab[k].astype(np.float32)
    
    for k in tab.colnames:
        if '-' in k:
            tab.rename_column(k, k.replace('-','_'))
            
    dm = tab['datemask'][0]
    fi = tab['filter'][i]
    tab['observer'] = [o.replace(',','') for o in tab['observer']]
    tfile = f'source_table/{dm}-{fi}.ecsv'
    tab[cols][so].write(tfile, overwrite=True)
    print(f'    >> {len(tab)}   {file}')
    
    return tab[cols][so]
    
        
if __name__ ==  '__main__':
    
    bash = """
    cd /home/idies/workspace/Temporary/gbrammer/scratch/MOSFIRE
    
    files=`shuf all_flats.txt | grep -v -e long2pos -e align`
    files=`shuf empty.log | grep -v -e long2pos -e align`
    
    # New masks
    rm new_flats.txt
    dirs=`grep "Oct 2[34]" auto.log | awk '{print $3}'`
    for dir in $dirs; do 
        ls ${dir}/Reduced/*/*/?/combfl*fits |grep -v _lamps_o >> new_flats.txt
    done
    
    ls *_20*/Reduced/*/20*/[YJHK]/combflat_2d*fits > all_flats.txt
    
    files=`cat all_flats.txt |grep -v -e _lamps_o -e long2pos`
    
    for file in $files; do 
        python mfpipe.py ${file}
    done
        
    """
    import sys
    flat_file = sys.argv[1]
    print(flat_file)
    
    SKIP = True
    
    run_mask(flat_file, skip=SKIP, initial_thresh=None, max_iter=50, use_ssl_slits=False)
    
        
