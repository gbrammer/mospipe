import os
import sys
import glob
import time

import numpy as np
from grizli import utils

def install_dfits():
    """
    Install dfits/fitsort tools and put them in the python/bin directory
    """
    import os
    import sys
    binpath = os.path.dirname(sys.executable)
    
    os.system('dfits > dfits.log')
    with open('dfits.log') as fp:
        lines = fp.readlines()
    
    os.remove('dfits.log')
    
    repo = 'https://github.com/granttremblay/eso_fits_tools.git'
    
    if len(lines) == 0:
        print(f'Install {repo} to {binpath}')
        os.system(f'git clone {repo}')
        os.chdir('eso_fits_tools')
        os.system('make > make.log.txt')
        os.system(f'cp dfits fitsort {binpath}/')
        os.chdir('../')
        os.system('rm -rf eso_fits_tools')


def run_pipeline(extra_query="AND progpi like '%%obash%%' AND progid='U190' and maskname='gs'", csv_file='mosfire.{hash}.csv', pwd='/GrizliImaging/', skip=True, min_nexp=10, sync=True, **kwargs):
    """
    Run the pipeline to download files and extract 2D spectra
    """
    import secrets
    from astropy.table import Table
    import matplotlib.pyplot as plt
    
    import mospipe.reduce
    
    binpath = os.path.dirname(sys.executable)
    
    install_dfits()
    
    if csv_file.startswith('s3://'):
        os.system(f'aws s3 cp {csv_file} . ')
        csv_file = os.path.basename(csv_file)
        
    if 'KOA_USERNAME' in os.environ:
        from pykoa.koa import Koa
        cookiepath = 'koa.kookie'
        Koa.login(cookiepath, userid=os.environ['KOA_USERNAME'],
                  password=os.environ['KOA_PASSWORD'], debugfile='koa.debug')
    else:
        from pykoa.koa import Archive
        cookiepath = ''
        Koa = Archive(debugfile='koa.debug')
    
    hash_key = secrets.token_urlsafe(16)[:6]
    hash_file = csv_file.format(hash=hash_key)
    
    if not os.path.exists(hash_file):
        cols = ['koaid', 'ofname', 'instrume as instrument', 'targname',
                'koaimtyp', 'frame', 'frameid', 'frameno', 'pattern', 
                'ra', 'dec', 
                "to_char(date_obs,'YYYY-MM-DD') as date_obs", 'ut', 
                'elaptime', 'waveblue', 'wavered', 'gratmode', 'pscale',
                'filter', 'mgtname', 'maskname', 'sampmode', 'numreads',
                'coadds', 'truitime', 'semid', 'proginst', 'progid', 'progpi',
                'progtitl', 'filehand', 'airmass', 'guidfwhm', 'mjd_obs', 
                "CONCAT(CONCAT(TRIM(maskname), '_'), " + 
                        "SUBSTR(koaid, 4, 8)) as datemask"]
    
        query =  f"""select {', '.join(cols)}
                    from koa_mosfire
                    WHERE gratmode='spectroscopy' AND koaimtyp='object'
                    AND maskname NOT LIKE '%%(align)%%'
                    {extra_query}
                    order by utdatetime"""
    
        print(f'======= Query ======= \n{query}\n ===============')
    
        Koa.query_adql(query, hash_file, overwrite=True, format='csv', 
                       cookiepath=cookiepath)
                       
    mfx = utils.read_catalog(hash_file)
    print(f'{len(mfx)} exposures found in {hash_file}')
    
    if len(mfx) < min_nexp:
        return False
        
    # mfx['datemask'] = [f"{mask}_{file.split('.')[1]}"
    #                    for mask, file in zip(mfx['maskname'], mfx['koaid'])]
    # if 'fileurl' not in mfx.colnames:
    #     mfx['fileurl'] = [f'filehand={f}' for f in mfx['filehand']]
    #     
    # mfx.write(hash_file, overwrite=True)
    
    ####### Run flats
    ONLY_FLAT = True
    
    un = utils.Unique(mfx['datemask'])
    
    masks = un.values
    
    pop = []
    
    for mi, mask in enumerate(masks):

        outdir = os.path.join(pwd, mask)

        if not os.path.exists(outdir):
            spl = outdir.split('/')
            #print(spl)
            for i in range(2,len(spl)+1):
                d = '/'+os.path.join(*spl[:i])
                if not os.path.exists(d):
                    print('mkdir '+d)
                    os.mkdir(d)
        else:
            if skip:
                print(f'{mask} exists, skip')
                pop.append(mi)
                continue

        with open(f'{pwd}/auto.log','a') as fp:
            fp.write(f'auto  - {mask} - {time.ctime()}\n')

        sel = (mfx['datemask'] == mask)
        tmp = mfx[sel]
        print(f'{mask}  N={len(tmp)}')
        tmp.write(os.path.join(pwd, f'{mask}_exposures.csv'), overwrite=True)
        
        if sel.sum() < min_nexp:
            pop.append(mi)
            print(f'{mask}: too few exposures found ({sel.sum()} < {min_nexp}), skipping')
            continue

        tmp['instrume'] = tmp['instrument']
        #tmp['filehand'] = [f.split('filehand=')[1] for f in tmp['fileurl']]

        mask_table = os.path.join(outdir, f'{mask}.tbl')
        tmp['koaid','instrume','filehand'].write(mask_table, 
                                                 format='ascii.ipac', 
                                                 overwrite=True)

        for d in ['Raw','Reduced']:
            dd = os.path.join(outdir, d)
            if not os.path.exists(dd):
                print(dd)
                os.mkdir(dd)

        ##### Download files
        outdir = os.path.join(pwd, mask)

        os.chdir(outdir)

        print(f'\n{mask}: Download\n')
        rawdir = os.path.join(pwd, mask, 'Raw')
        
        # Manual download files if missing
        fitsfiles = glob.glob(os.path.join(rawdir, '*fits'))
        wget = 'wget https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand={0} -O {1}'
        if len(fitsfiles) > 0:
            Koa.download(mask_table, 'ipac', rawdir, calibfile=0)   
            # Still missing?
            for f in tmp['filehand']:
                if not os.path.exists(os.path.join(rawdir, 
                                      os.path.basename(f))):
                    print('File still missing, try one more time')
                    Koa.download(mask_table, 'ipac', rawdir, calibfile=0)   
        else:
            if cookiepath:
                Koa.login(cookiepath, userid=os.environ['KOA_USERNAME'],
                          password=os.environ['KOA_PASSWORD'], 
                          debugfile='koa.debug')

            Koa.download(mask_table, 'ipac', rawdir, 
                         calibfile=1, cookiepath=cookiepath)   

        ##### Check for aborted exposures, which break the pipeline
        rawdir = os.path.join(pwd, mask, 'Raw')
        os.chdir(rawdir)

        files = glob.glob('MF*fits')
        if len(files) == 0:
            print(f'No downloaded files found for mask {mask}')
            continue

        os.system('dfits MF*fits | fitsort ABORTED > aborted.txt')

        info = Table.read('aborted.txt', format='ascii')
        bad = info['ABORTED'] == 'T'

        if bad.sum() > 0:
            for file in info['FILE'][bad]:
                print(f'{mask}: remove aborted file {file}')
                os.remove(file)
        else:
            print(f'{mask}: no aborted files')

        ###### Run the whole thing
        redpath = os.path.join(pwd, mask, 'Reduced')
        rawpath = os.path.join(pwd, mask, 'Raw')

        # Move files around
        os.chdir(redpath)
        os.system('rm ../Raw/translate.csh')

        print(f'\n {redpath}: translator\n')
        os.system(f'koa_translator_mod {rawpath}')

        # "handle" masks/filters
        print(f'\n {redpath}: handle\n')

        os.system(f'{sys.executable} {binpath}/mospy_handle.py '+ 
                  f'{redpath}/MOSFIRE/*/*/*fits > handle.log')

        # Put calibs and exposures in same directory if they were split up
        all_dirs = []
        for grat in 'YJHK':
            dirs = glob.glob(f'*/20*/{grat}')
            dirs.sort()
            if len(dirs) > 1:
                # Put all "txt" files in last directory
                for d in dirs[:-1]:
                    txt_files = glob.glob(f'{d}/*txt')
                    if len(txt_files) > 0:
                        for f in txt_files:
                            cmd = f'mv {f} {dirs[-1]}'
                            print(cmd)
                            os.system(cmd)

                    cmd = f'rm -rf {d}'
                    print(cmd)
                    os.system(cmd)

                all_dirs.append(dirs[-1])

            elif len(dirs) == 1:
                all_dirs.append(dirs[0])
            else:
                continue

        # Run it on all filters for a given mask
        for dir in all_dirs:
            if not os.path.exists(dir):
                continue

            os.chdir(dir)
            print(f'===========\nProcess mask {dir}\n============')
            os.system(f'{sys.executable} {binpath}/AutoDriver.py')

            longfiles = glob.glob('Long*py')

            # Don't run extractions, which are slow
            if os.path.exists('Driver.py'):

                if ONLY_FLAT:
                    print(f'\n#####\n Only flats! log={os.getcwd()}/mospy.log \n ######\n')

                    flat_files = glob.glob('combflat*')
                    if len(flat_files) > 0:
                        print(f'Flat files found: {flat_files}, skip')
                        continue

                    # Only up to flats
                    print('Run only flat')
                    with open('Driver.py') as fp:
                        lines = fp.readlines()

                    with open('RunFlat.py','w') as fp:
                        for line in lines:
                            fp.write(line)

                            if 'Flats.handle_flats' in line:
                                break

                    os.system(f'{sys.executable} RunFlat.py > mospy.log')

                else:
                    print(f'\n#####\n Running Driver.py, log={os.getcwd()}/mospy.log \n ######\n')
                    os.system('perl -pi -e "s/Extract.extract_spectra/# Extract.extract_spectra/" Driver.py')
                    os.system(f'{sys.executable} Driver.py > mospy.log')

            elif len(longfiles) > 0:
                ### Stop at flats for LongSlit reductions
                pyfile = longfiles[0]

                # Only up to flats
                print('Run only flat')
                with open(pyfile) as fp:
                    lines = fp.readlines()

                with open('RunFlat.py','w') as fp:
                    for line in lines:
                        fp.write(line)

                        if 'Flats.handle_flats' in line:
                            break

                os.system(f'{sys.executable} RunFlat.py > mospy.log')

            os.chdir(redpath)
                    
        # Extractions
        os.chdir(pwd)
        flat_files = glob.glob(f'{mask}/*/*/*/*/*combflat_2d*fits')
        for flat_file in flat_files:
            os.chdir(pwd)
            msk = mospipe.reduce.run_mask(flat_file, skip=skip, 
                                    initial_thresh=None, max_iter=50, 
                                    use_ssl_slits=False)
            
            plt.close('all')
            
        slit_summary(mask, outfile='slit_objects.csv')
        
        if sync:
            sync_results(mask)


def sync_results(mask, bucket='mosfire-pipeline', prefix='Spectra'):
    """
    Send files to S3 and update database
    """
    import pandas as pd
    from grizli.aws import db
    engine = db.get_db_engine()
    
    owd = os.getcwd()
    
    obj_file = f'{mask}_slit_objects.csv'
    
    if not os.path.exists(obj_file):
        return False
    
    obj = utils.read_catalog(obj_file)
    df = obj.to_pandas()
    
    db.execute_helper('DELETE FROM mosfire_extractions WHERE '
                       f"datemask='{mask}'", engine)

    df.to_sql('mosfire_extractions', engine, index=False, 
              if_exists='append', method='multi')
    
    print(f'{mask}_slit_objects > `mosfire_extractions`')
    
    # Exposures / Masks
    exp = utils.read_catalog(f'{mask}_exposures.csv')
    
    mask_cols = ['instrument', 'targname', 'koaimtyp', 'pattern', 'date_obs', 
                 'mgtname', 'maskname', 'semid', 'proginst', 'progid',
                 'progpi', 'progtitl', 'datemask']
    
    exp_cols = ['datemask', 'koaid', 'ofname', 'frame', 'frameid', 'frameno', 
                'ra', 'dec', 'ut', 'filehand', 'airmass', 'guidfwhm',
                'mjd_obs', 'elaptime', 'filter', 'waveblue', 'wavered',
                'gratmode', 'pscale', 'sampmode', 'numreads', 'coadds',
                'truitime']
    
    df_mask = exp[mask_cols][0:1].to_pandas()
    df_exp = exp[exp_cols].to_pandas()
    
    # Mask
    db.execute_helper('DELETE FROM mosfire_datemask WHERE '
                       f"datemask='{mask}'", engine)
    
    df_mask.to_sql('mosfire_datemask', engine, index=False, 
              if_exists='append', method='multi')
    
    # Exposures        
    db.execute_helper('DELETE FROM mosfire_exposures WHERE '
                       f"datemask='{mask}'", engine)

    df_exp.to_sql('mosfire_exposures', engine, index=False, 
              if_exists='append', method='multi')
    
    print(f'{mask}_exposures > `mosfire_exposures`, `mosfire_datemask`')

    #os.chdir(mask)
    os.system(f'aws s3 rm s3://{bucket}/{prefix}/{mask}/ --recursive')
    os.system(f'cd {owd}/{mask}; '+ 
              f'aws s3 sync ./ s3://{bucket}/{prefix}/{mask}/ ' + 
              '--exclude "*" --include "Reduced/*/*/[YJHK]/*"')
    
    files = glob.glob(f'{mask}*.*g')
    files.sort()
    for file in files:
        os.system(f'aws s3 cp {file} s3://{bucket}/Log/')
            
    os.chdir(owd)
    return True


def slit_summary(mask, outfile='slit_objects.csv'):
    """
    Summary of *extracted* slit spectra
    """
    import astropy.io.fits as pyfits
    
    files = glob.glob(f'{mask}/*/*/*/*/*-slit_*sp.fits')
    files.sort()
    
    if len(files) == 0:
        return None
        
    rows = []
    keys = ['SLITNUM','DATEMASK','TARGNAME','FILTER', 
            'NEXP', 'EXPTIME',
            'RA_SLIT','DEC_SLIT','RA_TARG','DEC_TARG','SKYPA3',
            'TARGOFF', 'TARGYPIX', 
            'TRAORDER', 'TRACOEF0', 'TRACOEF1', 'TRACOEF2', 
            'LAMORDER', 'LAMCOEF0', 'LAMCOEF1', 
            'Y0', 'Y1', 'YSTART', 'YSTOP', 'YPAD', 
            'MJD-OBS']
    
    colnames = ['file']
    colnames += [k.lower() for k in keys]
    colnames += ['slit_width', 'slit_length']

    for file in files:
        sp = pyfits.open(file)
        row = [file]
        for k in keys:
            row.append(sp[0].header[k])

        slit_length = (sp[0].header['YSTOP'] - sp[0].header['YSTART'])*0.1799

        row.extend([0.7, slit_length])

        rows.append(row)

    tab = utils.GTable(rows=rows, names=colnames)
    for k in ['RA_SLIT','DEC_SLIT','RA_TARG','DEC_TARG']:
        tab[k.lower()].format = '.6f'

    for k in ['EXPTIME','slit_width','slit_length']:
        tab[k.lower()].format = '.1f'

    tab['skypa3'].format = '6.1f'

    tab['slitnum'].format = '2d'

    for k in ['datemask','targname']:
        tab[k].format = '24'
    
    so = np.argsort(tab['slitnum'])
    tab = tab[so]
    
    if outfile:
        tab.write(f'{mask}_{outfile}', overwrite=True)
        print(f'Slit summary to {mask}_{outfile}')
    
    return tab
    
    
if __name__ == '__main__':
    argv = sys.argv
    kws = {}
    for arg in sys.argv:
        if arg.startswith('--') & ('=' in arg):
            kw = arg.split('=')[0][2:]
            val = '='.join(arg.split('=')[1:])
                        
            if val.lower() == 'false':
                val = False
            elif val.lower() == 'true':
                val = True

            if isinstance(val, str):
                try:
                    val = int(val)
                except ValueError:
                    pass
        
            kws[kw] = val
            
    print(f'Run pipeline with kwargs: {kws}')
    
    if ('extra_query' in kws) | ('csv_file' in kws):
        run_pipeline(**kws)
    else:
        print(f"Abort: kwargs must include 'extra_query' or 'csv_file'")

