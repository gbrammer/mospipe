import os
import sys
import glob
import time

import numpy as np
import astropy.time

from grizli import utils
utils.set_warnings()

MJDREF = astropy.time.Time('2022-01-01').mjd

def install_dfits():
    """
    Install dfits/fitsort tools and put them in the python/bin directory
    """
    import os
    import sys
    import secrets
    
    binpath = os.path.dirname(sys.executable)
    
    hash_key = secrets.token_urlsafe(16)[:6]
    hash_file = f'dfits.{hash_key}.log'
    
    os.system(f'dfits > {hash_file}')
    with open(hash_file) as fp:
        lines = fp.readlines()
    
    os.remove(hash_file)
    
    repo = 'https://github.com/granttremblay/eso_fits_tools.git'
    
    if len(lines) == 0:
        print(f'Install {repo} to {binpath}')
        os.system(f'git clone {repo}')
        os.chdir('eso_fits_tools')
        os.system('make > make.log.txt')
        os.system(f'cp dfits fitsort {binpath}/')
        os.chdir('../')
        os.system('rm -rf eso_fits_tools')


def master_query():
    """
    """
    import os
    import pandas as pd
    import astropy.time
    
    from grizli.aws import db
    engine = db.get_db_engine()
    
    if 'KOA_USERNAME' in os.environ:
        from pykoa.koa import Koa
        cookiepath = 'koa.kookie'
        Koa.login(cookiepath, userid=os.environ['KOA_USERNAME'],
                  password=os.environ['KOA_PASSWORD'], debugfile='koa.debug')
    else:
        from pykoa.koa import Archive
        cookiepath = ''
        Koa = Archive(debugfile='koa.debug')
    
    csv_file = 'master.csv'
    
    cols = ['maskname', 'filter', 'SUBSTR(koaid, 4, 8)']

    query =  f"""select maskname, filter, SUBSTR(koaid, 4, 8) night, COUNT(maskname) count
                from koa_mosfire
                WHERE gratmode='spectroscopy' AND koaimtyp='object'
                AND maskname NOT LIKE '%%(align)%%'
                AND maskname NOT LIKE '%%long2pos%%'                
                AND (progid = 'U190' OR progid = 'N097')
                GROUP BY maskname, filter, SUBSTR(koaid, 4, 8)
                """
    
    cos = "(contains(point('J2000',ra ,dec), circle('J2000', 150.1 2.0, 1.6))=1)"
    uds = "(contains(point('J2000',ra ,dec), circle('J2000', 34.409 -5.163, 2.0))=1)"
    gds = "(contains(point('J2000',ra ,dec), circle('J2000', 53.110 -27.830, 1.0))=1)"
    egs = "(contains(point('J2000',ra ,dec), circle('J2000', 214.8288 52.8067, 1.8))=1)"
    gdn = "(contains(point('J2000',ra ,dec), circle('J2000', 189.236	62.257, 1.8))=1)"
    
    query =  f"""select maskname, filter, progpi, SUBSTR(koaid, 4, 8) night, COUNT(maskname) count
                from koa_mosfire
                WHERE gratmode='spectroscopy' AND koaimtyp='object'
                AND maskname NOT LIKE '%%(align)%%'
                AND maskname NOT LIKE '%%long2pos%%'                
                AND maskname NOT LIKE '%%LONGSLIT%%'                
                AND ({egs} OR {cos} OR {uds} OR {gds} OR {gdn})
                AND (filter = 'Y' OR filter = 'J' OR filter = 'H' OR filter = 'K')
                GROUP BY maskname, filter, progpi, SUBSTR(koaid, 4, 8)
                """
    #
    query =  f"""select maskname, gratmode, koaimtyp, filter, progpi, progid, SUBSTR(koaid, 4, 8) night, ra, dec, object
    from koa_mosfire
    WHERE maskname = 'ud_van2'
    AND (filter = 'Y' OR filter = 'J' OR filter = 'H' OR filter = 'K')
                """
    #
    Koa.query_adql(query, csv_file, overwrite=True, format='csv', 
                   cookiepath=cookiepath)
    
    x = utils.read_catalog(csv_file)
    
    # Has flat
    query =  f"""select maskname, filter, progpi, progid, SUBSTR(koaid, 4, 8) night, AVG(dra) dra, AVG(ddec) ddec, AVG(ra) ra, AVG(dec) dec, COUNT(maskname) count
                from koa_mosfire
                WHERE maskname NOT LIKE '%%(align)%%'
                AND UPPER(maskname) NOT LIKE '%%LONG%%'                
                AND UPPER(maskname) NOT LIKE 'NGC%%'       
                AND UPPER(maskname) NOT LIKE 'HATP%%'       
                AND UPPER(maskname) NOT LIKE 'KEPL%%'       
                AND UPPER(maskname) NOT LIKE 'WASP%%'       
                AND UPPER(maskname) NOT LIKE 'GJ%%'       
                AND maskname = 'cr7_mosfire'
                AND (object LIKE 'Flat:%%')
                AND (filter = 'Y' OR filter = 'J' OR filter = 'H' OR filter = 'K')
                GROUP BY maskname, filter, progpi, progid, SUBSTR(koaid, 4, 8)
                """
    
    print(f'======= Query ======= \n{query}\n ===============')

    Koa.query_adql(query, csv_file, overwrite=True, format='csv', 
                   cookiepath=cookiepath)
    
    rflat = utils.read_catalog(csv_file)
    rflat['datemask'] = [f'{m}_{d}' for m, d in zip(rflat['maskname'], rflat['night'])]
    
    query =  f"""select maskname, filter, progpi, progid, SUBSTR(koaid, 4, 8) night, AVG(dra) dra, AVG(ddec) ddec, AVG(ra) ra, AVG(dec) dec, COUNT(maskname) count
                from koa_mosfire
                WHERE gratmode='spectroscopy' AND koaimtyp='object'
                AND maskname NOT LIKE '%%(align)%%'
                AND UPPER(maskname) NOT LIKE '%%LONG%%'                
                AND UPPER(maskname) NOT LIKE 'NGC%%'       
                AND UPPER(maskname) NOT LIKE 'HATP%%'       
                AND UPPER(maskname) NOT LIKE 'KEPL%%'       
                AND UPPER(maskname) NOT LIKE 'WASP%%'       
                AND UPPER(maskname) NOT LIKE 'GJ%%'       
                AND progpi NOT LIKE '%%ngineer%%'         
                AND progpi NOT LIKE 'mosfireeng'
                AND ((object NOT LIKE 'Flat:%%' AND object NOT LIKE 'Arc:%%') OR object IS NULL)
                AND (filter = 'Y' OR filter = 'J' OR filter = 'H' OR filter = 'K')
                GROUP BY maskname, filter, progpi, progid, SUBSTR(koaid, 4, 8)
                """
    
    print(f'======= Query ======= \n{query}\n ===============')

    Koa.query_adql(query, csv_file, overwrite=True, format='csv', 
                   cookiepath=cookiepath)
    
    res = utils.read_catalog(csv_file)
    
    res = res[res['count'] >= 8]
    res = res[res['count'] < 200]
    
    ############### Full query
    query =  f"""SELECT instrume, targname, koaimtyp, pattern, date_obs, mgtname, maskname, semid, proginst, progid, progpi, progtitl, SUBSTR(koaid, 4, 8) night
                from koa_mosfire
                WHERE gratmode='spectroscopy' AND koaimtyp='object'
                AND maskname NOT LIKE '%%(align)%%'
                AND UPPER(maskname) NOT LIKE '%%LONG%%'                
                AND UPPER(maskname) NOT LIKE 'NGC%%'       
                AND UPPER(maskname) NOT LIKE 'HATP%%'       
                AND UPPER(maskname) NOT LIKE 'KEPL%%'       
                AND UPPER(maskname) NOT LIKE 'WASP%%'       
                AND UPPER(maskname) NOT LIKE 'GJ%%'       
                AND progpi NOT LIKE '%%ngineer%%'         
                AND progpi NOT LIKE 'mosfireeng'
                AND ((object NOT LIKE 'Flat:%%' AND object NOT LIKE 'Arc:%%') OR object IS NULL)
                AND (filter = 'Y' OR filter = 'J' OR filter = 'H' OR filter = 'K')
                """
    
    #query += "AND date_obs > '2020-08-01'"
    
    Koa.query_adql(query, csv_file, overwrite=True, format='csv', 
                   cookiepath=cookiepath)
    
    res = utils.read_catalog(csv_file)
    
    res['datemask'] = [f'{m}_{d}' for m, d in zip(res['maskname'], res['night'])]
    
    from astropy.coordinates import SkyCoord
    coo = SkyCoord(res['ra'], res['dec'], unit=('deg','deg'))
    res['gal_lat'] = coo.galactic.b
    res['gal_lon'] = coo.galactic.l
    
    so = np.argsort(res['night'])[::-1]
    res = res[so]
    
    done = db.from_sql('select datemask, progid from mosfire_datemask', engine)
    #skip = False
    #for d in done['datemask']:
    #    skip |= res['datemask'] == d
    skip = np.in1d(res['datemask'], done['datemask'])
        
    res['status'] = 0
    
    hyp = np.array(['Hyperi' in d for d in res['datemask']])
    new = np.array([d > 20220000 for d in res['night']])
    highz = np.array(['z8' in d.lower() for d in res['datemask']])
    highz |= np.array(['z7' in d.lower() for d in res['datemask']])
    highz |= np.array(['z5' in d.lower() for d in res['datemask']])
    highz |= np.array(['z6' in d.lower() for d in res['datemask']])
    highz |= np.array(['z9' in d.lower() for d in res['datemask']])
    highz |= np.array(['z10' in d.lower() for d in res['datemask']])
    highz |= np.array(['red' in d.lower() for d in res['datemask']])
    highz |= np.array(['nug' in d.lower() for d in res['datemask']])
    
    pi = np.array(['illing' in d.lower() for d in res['progpi']])
    pi |= np.array(['ellis' in d.lower() for d in res['progpi']])
    pi |= np.array(['ellis' in d.lower() for d in res['progpi']])
    pi |= np.array(['finkels' in d.lower() for d in res['progpi']])
    pi |= np.array(['oesc' in d.lower() for d in res['progpi']])
    pi |= np.array(['ilson' in d.lower() for d in res['progpi']])
    pi |= np.array(['tanaka' in d.lower() for d in res['progpi']])
    pi |= np.array(['glazebrook' in d.lower() for d in res['progpi']])
    pi |= np.array(['apovich' in d.lower() for d in res['progpi']])
    
    keep = (~skip) #& (hyp | new | highz | pi)
    
    un = utils.Unique(res['datemask'][keep], verbose=False)
    df = pd.DataFrame()
    df['datemask'] = un.values
    df['status'] = 0
    df['updtime'] = astropy.time.Time.now().mjd - MJDREF
    
    df.to_sql('mosfire_datemask', engine, index=False, 
              if_exists='append', method='multi')
    
    # Full extractions
    run_one(datemask='cosmos_ID3_20220223', delete_from_s3=False, orig_slit_numbers=[17], save_full_drizzled=True, clean=False, skip=False)

    run_one(datemask='C2020_maskID2B_20220123', delete_from_s3=False, orig_slit_numbers=[13], save_full_drizzled=True, clean=False, skip=False)

    run_one(datemask='C2020_maskID2_20220114', delete_from_s3=False, orig_slit_numbers=[12], save_full_drizzled=True, clean=False, skip=False)

    run_one(datemask='C2020_maskID1B_20220123', delete_from_s3=False, orig_slit_numbers=[28], save_full_drizzled=True, clean=False, skip=False)

    run_one(datemask='C2020_maskID1_20220114', delete_from_s3=False, orig_slit_numbers=[3], save_full_drizzled=True, clean=False, skip=False)
    
    #
    ############### Full query
    query =  f"""SELECT instrume, targname, koaimtyp, pattern, date_obs, mgtname, maskname, semid, proginst, progid, progpi, progtitl, SUBSTR(koaid, 4, 8) night
                from koa_mosfire
                WHERE gratmode='spectroscopy' AND koaimtyp='object'
                AND maskname NOT LIKE '%%(align)%%'
                AND UPPER(maskname) NOT LIKE '%%LONG%%'                
                AND UPPER(maskname) NOT LIKE 'NGC%%'       
                AND UPPER(maskname) NOT LIKE 'HATP%%'       
                AND UPPER(maskname) NOT LIKE 'KEPL%%'       
                AND UPPER(maskname) NOT LIKE 'WASP%%'       
                AND UPPER(maskname) NOT LIKE 'GJ%%'       
                AND progpi NOT LIKE '%%ngineer%%'         
                AND progpi NOT LIKE 'mosfireeng'
                AND ((object NOT LIKE 'Flat:%%' AND object NOT LIKE 'Arc:%%') OR object IS NULL)
                AND (filter = 'Y' OR filter = 'J' OR filter = 'H' OR filter = 'K')
                """
    #
    Koa.query_adql(query, csv_file, overwrite=True, format='csv', 
                   cookiepath=cookiepath)
    
    res = utils.read_catalog(csv_file)
    
    res['datemask'] = [f'{m}_{d}' for m, d in zip(res['maskname'], res['night'])]
    
    un = utils.Unique(res['datemask'], verbose=False)
    ix = []
    res['count'] = 0
    
    for i, v in enumerate(un.values):
        res['count'][un[v]] = un.counts[i]
        if un.counts[i] < 200:
            ix.append(np.where(un[v])[0][0])
    
    res = res[ix]
    res.rename_column('instrume','instrument')
    res.remove_column('night')
    
    so = np.argsort(res['date_obs'])
    res = res[so]
    
    done = db.SQL('select datemask, progid, status from mosfire_datemask')
    skip = np.in1d(res['datemask'], done['datemask'])
    
    abell = np.array(['2744' in t for t in res['targname']])
    
    aold = np.array(['2744' in t for t in done['datemask']])
    
    
def run_pipeline(extra_query="AND progpi like '%%obash%%' AND progid='U190' and maskname='gs'", csv_file='mosfire.{hash}.csv', pwd='/GrizliImaging/', skip=True, min_nexp=10, sync=True, query_only=False, download_only=False, skip_long2pos=True, **kwargs):
    """
    Run the pipeline to download files and extract 2D spectra
    """
    import secrets
    from astropy.table import Table
    import matplotlib.pyplot as plt
    from MOSFIRE import IO
    
    import mospipe.reduce
    
    binpath = os.path.dirname(sys.executable)
    
    install_dfits()
    
    os.chdir(pwd)
    
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
                    WHERE gratmode='spectroscopy'
                    AND (koaimtyp='object' OR koaimtyp='flatlamp' OR koaimtyp='flatlampoff')
                    AND maskname NOT LIKE '%%(align)%%'
                    {extra_query}
                    order by utdatetime"""
    
        print(f'======= Query ======= \n{query}\n ===============')
    
        Koa.query_adql(query, hash_file, overwrite=True, format='csv', 
                       cookiepath=cookiepath)
        
    mfx = utils.read_catalog(hash_file)
    if query_only:
        return mfx
        
    print(f'{len(mfx)} exposures found in {hash_file}')
    
    if len(mfx) < min_nexp:
        return False
    
    ####### Run flats
    ONLY_FLAT = True
    
    un = utils.Unique(mfx['datemask'])
    
    datemasks = un.values
    
    pop = []
    
    for mi, datemask in enumerate(datemasks):
        
        LOGFILE = f'/GrizliImaging/{datemask}.pipeline.log'
        
        if ('long2pos' in datemask) & (skip_long2pos):
            print(f'Skip {datemask}')
            
        outdir = os.path.join(pwd, datemask)

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
                msg = f'{datemask} exists, skip'
                utils.log_comment(LOGFILE, msg, verbose=True, 
                                  show_date=True, mode='a')
                    
                pop.append(mi)
                continue

        with open(f'{pwd}/auto.log','a') as fp:
            fp.write(f'auto  - {datemask} - {time.ctime()}\n')

        sel = (mfx['datemask'] == datemask)
        tmp = mfx[sel]
        msg = f'{datemask}  N={len(tmp)}'
        utils.log_comment(LOGFILE, msg, verbose=True, 
                          show_date=True, mode='a')
            
        tmp.write(os.path.join(pwd, f'{datemask}_exposures.csv'), 
                  overwrite=True)
        
        if sel.sum() < min_nexp:
            pop.append(mi)
            msg = (f'{datemask}: too few exposures found '
                  f'({sel.sum()} < {min_nexp}), skipping')

            utils.log_comment(LOGFILE, msg, verbose=True, 
                              show_date=True, mode='a')
                
            continue

        tmp['instrume'] = tmp['instrument']
        #tmp['filehand'] = [f.split('filehand=')[1] for f in tmp['fileurl']]

        mask_table = os.path.join(outdir, f'{datemask}.tbl')
        tmp['koaid','instrume','filehand'].write(mask_table, 
                                                 format='ascii.ipac', 
                                                 overwrite=True)

        for d in ['Raw','Reduced']:
            dd = os.path.join(outdir, d)
            if not os.path.exists(dd):
                print(dd)
                os.mkdir(dd)

        ##### Download files
        outdir = os.path.join(pwd, datemask)

        os.chdir(outdir)

        msg = f'\n{datemask}: Download\n'
        utils.log_comment(LOGFILE, msg, verbose=True, 
                          show_date=True, mode='a')
            
        rawdir = os.path.join(pwd, datemask, 'Raw')
        
        # Try to sync from S3 staging
        s3stage = f's3://mosfire-pipeline/RawFiles/{datemask}/'
        os.system(f'aws s3 sync {s3stage} {rawdir}/')

        os.system(f'aws s3 sync s3://mosfire-pipeline/Spectra/{datemask}/ {pwd}/{datemask}/ --exclude "*" --include "*_sp.fits"')
        
        # Manual download files if missing
        fitsfiles = glob.glob(os.path.join(rawdir, '*fits'))
        wget = 'wget https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand={0} -O {1}'
        if len(fitsfiles) > 0:
            if cookiepath:
                Koa.login(cookiepath, userid=os.environ['KOA_USERNAME'],
                          password=os.environ['KOA_PASSWORD'], 
                          debugfile='koa.debug')

            Koa.download(mask_table, 'ipac', rawdir, calibfile=0, 
                         cookiepath=cookiepath)  
                          
            # Still missing?
            for f in tmp['filehand']:
                if not os.path.exists(os.path.join(rawdir, 
                                      os.path.basename(f))):
                    print('File still missing, try one more time')
                    Koa.download(mask_table, 'ipac', rawdir, calibfile=0, 
                                 cookiepath=cookiepath)   
        else:
            if cookiepath:
                Koa.login(cookiepath, userid=os.environ['KOA_USERNAME'],
                          password=os.environ['KOA_PASSWORD'], 
                          debugfile='koa.debug')

            Koa.download(mask_table, 'ipac', rawdir, 
                         calibfile=1, cookiepath=cookiepath)   

        ##### Check for aborted exposures, which break the pipeline
        rawdir = os.path.join(pwd, datemask, 'Raw')
        os.chdir(rawdir)

        # sync to s3 staging
        os.system(f'aws s3 sync {rawdir}/ {s3stage}')
        
        # Download complete
        update_mask_db_status(datemask, 2, verbose=True)
        
        if download_only:
            update_mask_db_status(datemask, -2, verbose=True)
            continue
            
        files = glob.glob('MF*fits')
        if len(files) == 0:
            msg = f'No downloaded files found for mask {datemask}'
            utils.log_comment(LOGFILE, msg, verbose=True, 
                              show_date=True, mode='a')
                
            continue

        os.system('dfits MF*fits | fitsort ABORTED > aborted.txt')

        info = Table.read('aborted.txt', format='ascii')
        bad = info['ABORTED'] == 'T'

        if bad.sum() > 0:
            for file in info['FILE'][bad]:
                msg = f'{datemask}: remove aborted file {file}'
                utils.log_comment(LOGFILE, msg, verbose=True, 
                                  show_date=True, mode='a')
                os.remove(file)
        else:
            msg = f'{datemask}: no aborted files'
            utils.log_comment(LOGFILE, msg, verbose=True, 
                              show_date=True, mode='a')
                
        ###### Run the whole thing
        redpath = os.path.join(pwd, datemask, 'Reduced')
        rawpath = os.path.join(pwd, datemask, 'Raw')

        # Move files around
        os.chdir(redpath)
        os.system('rm ../Raw/translate.csh')

        msg = f'\n {redpath}: translator\n'
        utils.log_comment(LOGFILE, msg, verbose=True, 
                          show_date=True, mode='a')
            
        os.system(f'koa_translator_mod {rawpath}')

        # "handle" masks/filters
        msg = f'\n {redpath}: handle\n'
        utils.log_comment(LOGFILE, msg, verbose=True, 
                          show_date=True, mode='a')
        
        # Remove "Improper" files
        
        raw_files = glob.glob(f'{redpath}/MOSFIRE/*/*/*fits')
        raw_files.sort()
        for file in raw_files:
            try:
                _ = IO.readmosfits(file, {})
            except:
                msg = f'remove invalid file {file}'
                utils.log_comment(LOGFILE, msg, verbose=True, 
                                  show_date=True, mode='a')
                
                os.remove(file)
        
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
        
        utils.log_comment(LOGFILE, f'all_dirs: {all_dirs}', verbose=True, 
                          show_date=True, mode='a')
        
        # Run it on all filters for a given mask
        for dir in all_dirs:
            if not os.path.exists(dir):
                continue

            os.chdir(dir)
            msg = f'===========\nProcess mask {dir}\n============'
            utils.log_comment(LOGFILE, msg, verbose=True, 
                              show_date=True, mode='a')
            
            # ToDo: check that offset pairs exist                  
            offset_files = glob.glob('Offset*txt')
            if len(offset_files) < 2:
                msg = f'{dir} only found {len(offset_files)} Offset*txt files'
                utils.log_comment(LOGFILE, msg, verbose=True, 
                                  show_date=True, mode='a')
                os.chdir(redpath)
                continue
                
            os.system(f'{sys.executable} {binpath}/AutoDriver.py')

            longfiles = glob.glob('Long*py')
                
            # Don't run extractions, which are slow
            if os.path.exists('Driver.py'):

                if ONLY_FLAT:
                    msg = f'Only flats! log={os.getcwd()}/mospy.log'
                    utils.log_comment(LOGFILE, msg, verbose=True, 
                                      show_date=True, mode='a')

                    flat_files = glob.glob('combflat*')
                    if len(flat_files) > 0:
                        msg = f'Flat files found: {flat_files}, skip'
                        utils.log_comment(LOGFILE,
                            msg, verbose=True, show_date=True, mode='a')
                            
                        continue
                    
                    # Only up to flats
                    msg = 'Run only flat'
                    utils.log_comment(LOGFILE, msg, verbose=True, 
                                      show_date=True, mode='a')
                                      
                    with open('Driver.py') as fp:
                        lines = fp.readlines()

                    with open('RunFlat.py','w') as fp:
                        for line in lines:
                            fp.write(line)

                            if 'Flats.handle_flats' in line:
                                break

                    os.system(f'{sys.executable} RunFlat.py > mospy.log')

                else:
                    msg = f'Running Driver.py, log={os.getcwd()}/mospy.log'
                    utils.log_comment(LOGFILE, msg, verbose=True, 
                                      show_date=True, mode='a')
                                      
                    os.system('perl -pi -e "s/Extract.extract_spectra/# Extract.extract_spectra/" Driver.py')
                    os.system(f'{sys.executable} Driver.py > mospy.log')

            elif len(longfiles) > 0:
                ### Stop at flats for LongSlit reductions
                pyfile = longfiles[0]

                # Only up to flats
                msg = 'Run only flat'
                utils.log_comment(LOGFILE, msg, verbose=True, 
                                  show_date=True, mode='a')
                                  
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
        flat_files = glob.glob(f'{datemask}/*/*/*/*/*combflat_2d*fits')
        
        if len(flat_files) > 0:
            update_mask_db_status(datemask, 3, verbose=True)
        else:
            update_mask_db_status(datemask, 9, verbose=True)
            
        
        for flat_file in flat_files:
            os.chdir(pwd)
            msk = mospipe.reduce.run_mask(flat_file, skip=skip, 
                                          initial_thresh=None, max_iter=50, 
                                          use_ssl_slits=False, **kwargs)
            
            plt.close('all')
            
        slit_summary(datemask, outfile='slit_objects.csv')
        
        if sync:
            sync_results(datemask, **kwargs)


def sync_results(datemask, bucket='mosfire-pipeline', prefix='Spectra', delete_from_s3=False, make_public=True, **kwargs):
    """
    Send files to S3 and update database
    """
    import pandas as pd
    from astropy.time import Time
    from grizli.aws import db
    engine = db.get_db_engine()
    
    owd = os.getcwd()
    
    obj_file = f'{datemask}_slit_objects.csv'
    
    if not os.path.exists(obj_file):
        return False
    
    obj = utils.read_catalog(obj_file)
    df_obj = obj.to_pandas()
        
    # Exposures / Masks
    exp = utils.read_catalog(f'{datemask}_exposures.csv')
    #exp['status'] = 4
    #exp['updtime'] = Time.now().mjd
    
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

    # Send to tables
    files = glob.glob(f'{datemask}/*/*/*/*/*log1d.fits')
    files.sort()    
    filters = np.unique([f.split('/')[-2].lower() for f in files])
    for f in filters:
        db.execute_helper(f'DELETE FROM mosfire_spectra_{f} WHERE '
                       f"datemask='{datemask}'", engine)
                       
    db.execute_helper('DELETE FROM mosfire_extractions WHERE '
                       f"datemask='{datemask}'", engine)

    db.execute_helper('DELETE FROM mosfire_exposures WHERE '
                       f"datemask='{datemask}'", engine)
    
    db.execute_helper('DELETE FROM mosfire_datemask WHERE '
                       f"datemask='{datemask}'", engine)
        
    df_mask.to_sql('mosfire_datemask', engine, index=False, 
              if_exists='append', method='multi')
    
    update_mask_db_status(datemask, 4, verbose=True)
    
    df_exp.to_sql('mosfire_exposures', engine, index=False, 
              if_exists='append', method='multi')
    
    print(f'{datemask}_exposures > `mosfire_exposures`, `mosfire_datemask`')

    df_obj.to_sql('mosfire_extractions', engine, index=False, 
              if_exists='append', method='multi')
    
    # Done, remove exposure files
    # os.system(f'aws s3 rm --recursive s3://mosfire-pipeline/RawFiles/{datemask}/')
    
    print(f'{datemask}_slit_objects > `mosfire_extractions`')
    
    # 1D spectra
    for file in files:
        spec = utils.read_catalog(file)
        
        df = pd.DataFrame()
        df['datemask'] = [spec.meta['DATEMASK']]
        df['filter'] = [spec.meta['FILTER']]
        df['slitnum'] = [spec.meta['SLITNUM']]
        df['flux'] = [spec["flux"].astype(np.float32).tolist()]
        df['err'] = [spec["err"].astype(np.float32).tolist()]
        df['lineflux'] = [spec["line_flux"].astype(np.float32).tolist()]
        df['lineerr'] = [spec["line_err"].astype(np.float32).tolist()]
        
        oned_table = f"mosfire_spectra_{spec.meta['FILTER']}".lower()
        msg = f'{os.path.basename(file)} > `{oned_table}`'
        utils.log_comment(f'/GrizliImaging/{datemask}.pipeline.log',
            msg, verbose=True, show_date=True, mode='a')
        
        df.to_sql(oned_table, 
                  engine, index=False, 
                  if_exists='append', method='multi')
    
    if delete_from_s3:
        os.system(f'aws s3 rm s3://{bucket}/{prefix}/{datemask}/ --recursive')
    
    pub_str = ' --acl public-read ' if make_public else ''
    
    os.system(f'cd {owd}/{datemask}; '+ 
              f'aws s3 sync ./ s3://{bucket}/{prefix}/{datemask}/ ' + 
              f'--exclude "*" --include "Reduced/*/*/[YJHK]/*" {pub_str}')
    
    files = glob.glob(f'{datemask}*.*g')
    files.sort()
    for file in files:
        os.system(f'aws s3 cp {file} s3://{bucket}/Log/ {pub_str}')
            
    os.chdir(owd)
    return True


def slit_summary(datemask, outfile='slit_objects.csv'):
    """
    Summary of *extracted* slit spectra
    """
    import astropy.io.fits as pyfits
    from astropy.time import Time
    from . import extract
    
    files = glob.glob(f'{datemask}/*/*/*/*/*-slit_*_sp.fits')
    files.sort()
    
    if len(files) == 0:
        return None
        
    keys = ['SLITNUM', 'DATEMASK','TARGNAME','FILTER', 
            'NEXP', 'EXPTIME',
            'RA_SLIT','DEC_SLIT','RA_TARG','DEC_TARG','SKYPA3',
            'TARGOFF', 'TARGYPIX', 
            'TRAORDER', 'TRACOEF0', 'TRACOEF1', 'TRACOEF2', 
            'LAMORDER', 'LAMCOEF0', 'LAMCOEF1', 
            'Y0', 'Y1', 'YSTART', 'YSTOP', 'YPAD', 
            'MJD-OBS']
    
    colnames = ['file', 'modtime']
    colnames += [k.lower() for k in keys]
    colnames += ['slit_width', 'slit_length']
    
    oned_cols = ['binw','linew','wmin','wmax','sn16','sn50','sn84',
                 'sn99', 'prof_amp', 'prof_sig', 'prof_mu',
                 'prof_yma', 'prof_ma', 'prof_ymi', 'prof_mi', 
                 'prof_offset', 'pthresh', 'lthresh', 'nline']
    
    oned_ints = ['binw','linew','nline', 'prof_mi', 'prof_ma']
    
    for j in range(4):
        oned_cols += [f'linew{j:02d}', f'linesn{j:02d}', f'linef{j:02d}', f'linee{j:02d}']
    
    if False:
        ## Add columns to database
        from grizli.aws import db
        engine = db.get_db_engine()
        
        for k in oned_cols:
            if k in oned_ints:
                dtype = 'INT'
            else:
                dtype = 'REAL'
             
            cmd = f'ALTER TABLE mosfire_extractions ADD COLUMN {k} {dtype}'  
            print(cmd)
            engine.execute(cmd)
                
    colnames += oned_cols
    
    rows = []
    
    for file in files:
        
        sp = pyfits.open(file)
        modtime = Time(os.path.getmtime(file), format='unix').mjd - MJDREF
        row = [file, modtime]
        for k in keys:
            row.append(sp[0].header[k])

        slit_length = (sp[0].header['YSTOP'] - sp[0].header['YSTART'])*0.1799

        row.extend([0.7, slit_length])
        
        # Oned extraction
        try:
            #if 1:
            spec = extract.runit(file)
            meta = spec["log_spec"].meta
            
            oned_row = []
            for k in oned_cols:
                if k in meta:
                    if hasattr(meta[k], '__len__'):
                        val = meta[k][0]
                    else:
                        val = meta[k]
                    
                    if k in oned_ints:
                        try:
                            val = int(val)
                        except ValueError:
                            val = -999
                else:
                    if k in oned_ints:
                        val = -999
                    else:
                        val = -999.
                        
                oned_row.append(val)
            
            msg = f'{datemask}: 1D extraction for {file}'
            utils.log_comment(f'/GrizliImaging/{datemask}.pipeline.log',
                msg, verbose=True, show_date=True, mode='a')
            
        except:
            msg = f'{datemask}: 1D extraction *failed* for {file}'
            utils.log_comment(f'/GrizliImaging/{datemask}.pipeline.log',
                msg, verbose=True, show_date=True, mode='a')
            
            oned_row = []
            for k in oned_cols:
                if k in oned_ints:
                    val = -9999
                else:
                    val = -9999.
                        
                oned_row.append(val)
               
        rows.append(row + oned_row)

    tab = utils.GTable(rows=rows, names=colnames)
    for k in ['RA_SLIT','DEC_SLIT','RA_TARG','DEC_TARG']:
        tab[k.lower()].format = '.6f'

    for k in ['EXPTIME','slit_width','slit_length']:
        tab[k.lower()].format = '.1f'

    tab['skypa3'].format = '6.1f'

    tab['slitnum'].format = '2d'
    tab['modtime'].format = '.3f'
    #tab['slitidx'].format = '2d'
    
    for k in oned_cols:
        if k not in oned_ints:
            tab[k].format = '.2f'
            
    for k in ['datemask','targname']:
        tab[k].format = '24'
    
    so = np.argsort(tab['slitnum'])
    tab = tab[so]
    
    tab.rename_column('targname', 'target_name')
    
    if outfile:
        tab.write(f'{datemask}_{outfile}', overwrite=True)
        msg = f'Slit summary to {datemask}_{outfile}'
        utils.log_comment(f'/GrizliImaging/{datemask}.pipeline.log',
            msg, verbose=True, show_date=True, mode='a')
            
    return tab
    
    
def setup_db_tables():
    """
    Set up table indices
    """
    from grizli.aws import db
    engine = db.get_db_engine()
    
    engine.execute('ALTER TABLE mosfire_datemask ADD PRIMARY KEY (datemask)')
    
    SQL = """
    ALTER TABLE mosfire_exposures 
    ADD CONSTRAINT fk_exp_datemask 
    FOREIGN KEY (datemask) 
    REFERENCES mosfire_datemask (datemask);
    """
    engine.execute(SQL)

    SQL = """
    ALTER TABLE mosfire_extractions 
    ADD CONSTRAINT fk_ext_datemask 
    FOREIGN KEY (datemask) 
    REFERENCES mosfire_datemask (datemask);
    """
    engine.execute(SQL)
    
    SQL = """
    ALTER TABLE mosfire_extractions ADD PRIMARY KEY (datemask, filter, slitnum);"""
    engine.execute(SQL)
    
    for filt, n in zip('yjhk', [996, 1009, 1305, 1469]):
        
        SQL = f"""
    CREATE TABLE IF NOT EXISTS mosfire_spectra_{filt} (
        datemask VARCHAR, 
        filter VARCHAR, 
        slitnum INT, 
        flux REAL[{n}], 
        err REAL[{n}], 
        lineflux REAL[{n}], 
        lineerr REAL[{n}], 
        FOREIGN KEY (datemask, filter, slitnum)
            REFERENCES mosfire_extractions (datemask, filter, slitnum)
        );
        """
    
        engine.execute(SQL)
    
    SQL = """
    CREATE TABLE IF NOT EXISTS mosfire_spectra_k (
        datemask VARCHAR, 
        filter VARCHAR, 
        slitnum INT, 
        flux REAL[1469], 
        err REAL[1469], 
        lineflux REAL[1469], 
        lineerr REAL[1469], 
        FOREIGN KEY (datemask, filter, slitnum)
            REFERENCES mosfire_extractions (datemask, filter, slitnum)
        );
        """
    
    engine.execute(SQL)
    
    engine.execute('GRANT ALL PRIVILEGES ON ALL TABLEs IN SCHEMA public TO db_iam_user')
    engine.execute('GRANT SELECT ON ALL TABLEs IN SCHEMA public TO redshift_fit_public')
    
    # Test


def get_oned_wavelengths(binw=50, filter='K'):
    """
    Get log-spaced wavelengths of 1d spectra in the database tables
    """
    
    wlim = {'Y':( 9614.2, 11348.8),
            'J':(11453.2, 13548.7),
            'H':(14594.2, 18137.4),
            'K':(18905.9, 24146.9)}
    
    logw = utils.log_zgrid([3000, 2.5e4], binw/3.e5)

    logw = (logw[:-1]+np.diff(logw)/2.).astype(np.float32)
    clip = (logw > wlim[filter][0]) & (logw < wlim[filter][1])
    logw = logw[clip]
    
    return logw


def update_mask_db_status(datemask, status, verbose=True):
    """
    Set status flag of a mask in the `mosfire_datemask` table
    """
    import pandas as pd
    from astropy.time import Time
    
    from grizli.aws import db
    engine = db.get_db_engine()
    
    NOW = Time.now().mjd - MJDREF
    
    table = 'mosfire_datemask'
    
    sqlstr = f"""UPDATE {table}
        SET status = {status}, updtime = '{NOW}'
        WHERE (datemask = '{datemask}');"""

    try:
        db.execute_helper(sqlstr, engine)
        if verbose:
            msg = (f'Update status = {status} for {datemask} '
                   f'on `{table}` ({NOW})')
            print(msg)
    except:
        msg = (f'FAILED Update status = {status} for {datemask} '
               f'on `{table}` ({NOW})')
        print(msg)
        

def get_random_mask(extra='', status='status = 0', **kwargs):
    """
    Find a mask that needs processing
    """
    from grizli.aws import db
    engine = db.get_db_engine()
        
    all_masks = db.from_sql('SELECT DISTINCT(datemask) FROM mosfire_datemask' 
                             f' WHERE {status} ' + extra, engine)
    
    if len(all_masks) == 0:
        return None
    
    random_mask = all_masks[np.random.randint(0, len(all_masks))][0]
    return random_mask


def run_one(datemask=None, clean=True, **kwargs):
    """
    Process a with status=0
    """
    import os
    import time
    from grizli.aws import db

    engine = db.get_db_engine()

    if datemask is None:
        datemask = get_random_mask(**kwargs)
    
    if datemask is None:
        with open('/GrizliImaging/mosfire_finished.txt','w') as fp:
            fp.write(time.ctime() + '\n')
    else:
        print(f'==============  Run datemask  ===============')
        print(f'{time.ctime()}   {datemask}')
        
        maskname = '_'.join(datemask.split('_')[:-1])
        koaid = 'MF.{0}'.format(datemask.split('_')[-1])
        query = f"AND koaid LIKE '{koaid}%%' AND maskname='{maskname}'"
        kws = dict(extra_query=query, 
                   csv_file=f'mosfire.{datemask}.csv', 
                   **kwargs)
        
        update_mask_db_status(datemask, 1, verbose=True)
        
        with open('/GrizliImaging/mosfire_processed.txt','a') as fp:
            fp.write(f'{datemask} {time.ctime()}\n')
            
        run_pipeline(**kws)
        
        os.chdir('/GrizliImaging/')
        
        if clean & (len(datemask) > 0):
            os.system(f'rm -rf /GrizliImaging/{datemask}')


def parse_sys_argv(argv):
    """
    Parse keywords and boolean switches from sys.argv
    
    Parameters
    ----------
    argv : list
    
        ```
        ./test.py xxx --word=test --boolean=True --integer=3 --float=4.345 \
                      --quoted="[(*this is a string = 3.4)" \
                      --nan_val=nan --none_val=None --inf_val=inf \
                      --neg_int=-3 --neg-float=-3.1415 \
                      -false_switch +true_switch
        ```
    
    Returns
    -------
    kws : dict
        
        ```
        {'word': 'test', 'boolean': True, 'integer': 3, 'float': 4.345, 
         'quoted': '[(*this is a string = 3.4)',
         'nan_val': nan, 'none_val': None, 'inf_val': inf,
         'neg_int': -3, 'neg-float': -3.1415,
         'false_switch': False, 'true_switch': True}
        ```
    
    """
    import numpy as np
    
    kws = {}
    for arg in sys.argv:
        if arg.startswith('--') & ('=' in arg):
            kw = arg.split('=')[0][2:]
            val = '='.join(arg.split('=')[1:])
            
            if val.lower() == 'false':
                val = False
            elif val.lower() == 'true':
                val = True
            elif val.lower() == 'none':
                val = None
            elif val.lower() == 'inf':
                val = np.inf
            elif val.lower() == 'nan':
                val = np.nan
                
            if isinstance(val, str):
                try:
                    val = int(val)
                except ValueError:
                    try:
                        val = float(val)
                    except ValueError:
                        pass
        
            kws[kw] = val
            
        elif arg.startswith('-') & (arg[1] != '-'):
            # False switch
            kw = arg[1:]
            kws[kw] = False
        
        elif arg.startswith('+'):
            # True switch
            kw = arg[1:]
            kws[kw] = True
            
    return kws


if __name__ == '__main__':
    argv = sys.argv
    kws = parse_sys_argv(argv)
    print(f'Run pipeline with kwargs: {kws}')
    
    if ('extra_query' in kws) | ('csv_file' in kws):
        run_pipeline(**kws)
    else:
        print(f"Abort: kwargs must include 'extra_query' or 'csv_file'")
    
    point_query = """
    where (contains(point('J2000',ra ,dec), circle('J2000', 215.050566 53.007441, 0.2))=1)
     """
