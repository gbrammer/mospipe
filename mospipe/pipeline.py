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


def run_pipeline(extra_query="AND progpi like '%%obash%%' AND progid='U190' and maskname='gs'", csv_file='mosfire.{hash}.csv', pwd='/GrizliImaging/', skip=True, min_nexp=10):
    """
    Run the pipeline to download files and extract 2D spectra
    """
    import secrets
    from astropy.table import Table
    import matplotlib.pyplot as plt
    
    import mospipe.reduce
    
    binpath = os.path.dirname(sys.executable)
    
    install_dfits()
    
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
        
        if 'lont2pos' in mask:
            continue
            
        # Extractions
        os.chdir(pwd)
        flat_files = glob.glob(f'{mask}/*/*/*/*/*combflat*fits')
        for flat_file in flat_files:
            os.chdir(pwd)
            mospipe.reduce.run_mask(flat_file, skip=False, 
                                    initial_thresh=None, max_iter=50, 
                                    use_ssl_slits=False)
            plt.close('all')