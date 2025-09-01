from warnings import warn
import re

import pandas as pd

from phot import make_expmap
from run import heasoftpy_run, cmd_run, update_regions_table, update_pix_regions_table
from utils import  _generic_xselect, gen_pix_sel_str, generate_pixreg_file
from constants import GRADES, INST_NAMES, INST_MODES
from classes import EventList
from conf import ANALYSIS_DIR
from timing import make_gti

def make_spectra(evt: EventList, min_energy: int = None, max_energy: int = None, 
                 include_pix: list[int] = None, exclude_pix: list[int] = None, 
                 pix_name: str = None, regmode: str = 'DET', regionfile: str = None, 
                 regionname: str = None, grades: list[str] = None) -> list[str]:
    """
    Make spectra from XRISM event lists. Note about the region filtering: for resolve, set the 
    variables: include_pix OR exclude_pix and pix_name, for xtend set the variables regmode and 
    regionfile and regionname. You always need to set either pix_name or region_name if you want to 
    do any region filtering. Note that this code doesn't check that regionfile contents are in 
    actually in the coords specified by regmode.

    Parameters
    ---------
    evt : EventList
        The EventList object to be read into xselect
    min_energy : str
        in pha, the min energy cut
    max_energy : str
        in pha, the max energy cut
    include_pix : list[int]
        list of pixel indicies to include in products.
    exclude_pix : list[int]
        list of pixel indicies to exclude from products.
    pix_name : str
        this is like a nickname for the selected pixel regions, this will be appended to file and 
        dir names so it acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    regmode : str
        coordinate system that the regionfile is in
    regionfile : str
        absolute path to regionfile used to filter the product
    regionname : str
        this is like a nickname for the region, this will be appended to file and dir names so it 
        acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    grades : list[str]
        list of grades to include when generating products. ie, should be given as ['Hp', 'Mp'] not 
        [0, 1]. 

    Returns
    -------
    listdir : list[str]
        list of filenames that were created by the command run
    """
    
    return _generic_xselect(evt, 'spectrum', min_energy=min_energy,
                     max_energy=max_energy, include_pix=include_pix, exclude_pix=exclude_pix, 
                     pix_name=pix_name, regmode=regmode, regionfile=regionfile, 
                     regionname=regionname, grade=grades)

@heasoftpy_run
def make_arf(evt: EventList, ra: float, dec: float, expmap: str, rmf: str, regmode: str = None, 
             regionfile: str = None, regionname: str = None, numphoton: int = 60000, 
             include_pix: list[int] = None, exclude_pix: list[int] = None, pix_name: str = None):
    """
    Makes a XRISM arf.
    
    Note that for Resolve, only regions in DET coordinates should be used because for spectra 
    extracted in SKY coordinates, photons cannot be traced back to the individual pixels that they 
    landed in due to the large size of the pixels compared to the size of the extraction region.

    The regmode and regionfile options should match the source regions used for the spectrum 
    extraction.

    Paramters
    ---------
    evt : EventList
        The EventList object to generate an arf for. In this function this is only used to see which
        instrument the arf needs to be generated for.
    ra : float
        source ra in deg
    dec : float 
        source dec in deg
    expmap:
        absolute path to the exposure map needed for arf generation
    rmf:
        absolute path to the rmf needed for arf generation
    regmode : str
        coordinate system that the regionfile is in
    regionfile : str
        absolute path to regionfile used to filter the product
    regionname : str
        this is like a nickname for the region, this will be appended to file and dir names so it 
        acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.  
    numphoton : int
        the number of photons to used in the raytrace simulation
    include_pix : list[int]
        list of pixel indicies to include in products.
    exclude_pix : list[int]
        list of pixel indicies to exclude from products.
    pix_name : str
        this is like a nickname for the selected pixel regions, this will be appended to file and 
        dir names so it acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    """
    params = {
        'source_ra':ra,
        'source_dec': dec,
        'telescope': 'XRISM',
        'instrume': evt.inst,
        'emapfile': expmap,
        'rmffile': rmf,
        'numphoton': numphoton,
        'sourcetype': 'POINT',
        'cleanup': 'no',
        'history': 'no'
    }

    region_info = ''
    heasoft_cmd = 'xaarfgen'

    prefix = f'{evt.inst}'

    if evt.inst.upper() == 'RESOLVE':
        if pix_name is None:
            pix_name = 'ALLPIX'
        
        use_regionfile = generate_pixreg_file('to_move_' + pix_name + '.reg', include_pix, exclude_pix)

        params['regmode'] = 'DET'
        params['regionfile'] = use_regionfile
        prefix += f"_pixreg-{pix_name}"
        
        if pix_name != 'ALLPIX':
            region_info = {}
            region_info['include_pix'] = include_pix
            region_info['exclude_pix'] = exclude_pix
            region_info['pix_name'] = pix_name
            region_info['pix_reg'] = True

    else:
        if regionname is not None:
            params['regmode'] = regmode
            params['regionfile'] = regionfile
            prefix += f"_reg-{regionname}"
            region_info = {}
            # If region is provided, read its content
            with open(regionfile, 'r') as f:
                region_content = f.read()
            
            region_info['name'] = regionname
            region_info['file'] = regionfile
            region_info['content'] = region_content
            region_info['pix_reg'] = False

    raytracefile = 'rayt_' + prefix + '.fits'
    outfile = prefix + '.arf'

    params['xrtevtfile'] = raytracefile
    params['outfile'] = outfile

    return params, heasoft_cmd, prefix, region_info, 'arf'
        

@heasoftpy_run
def make_rmf(evt: EventList, whichrmf: str, splitrmf: str = 'no', splitcomb: str = 'no', 
             regmode: str = None, grades: list[str] = None, regionfile: str = None, 
             regionname: str = None, include_pix: list[int] = None, exclude_pix: list[int] = None, 
             pix_name: str = None):
    """ 
    Makes a XRISM RMF. Note input event file should contain events of all X-ray grades.
    splitcomb = yes --> the two response matrices (corresponding to the core and ELC components) 
    will be placed into a single file.

    Paramters
    ---------
    evt : EventList
        The EventList object to generate an rmf for.
    whichrmf : str
        either 'S', 'M', 'L', 'X'
    splitrmf : str
        If 'splitrmf=yes', split the RMF into core and ELC responses.
    splitcomb : str
        If 'splitrmf=yes' and 'splitcomb=yes', the two response matrices (corresponding to the core 
        and ELC components) will be placed into a single file.
    grades : list[str]
        grades of events to include in rmf. by default all will be included. should be input using
        names of grades, eg. ['Hp', 'Mp']
    regmode : str
        coordinate system that the regionfile is in
    regionfile : str
        absolute path to regionfile used to filter the product
    regionname : str
        this is like a nickname for the region, this will be appended to file and dir names so it 
        acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names. 
    include_pix : list[int]
        list of pixel indicies to include in products.
    exclude_pix : list[int]
        list of pixel indicies to exclude from products.
    pix_name : str
        this is like a nickname for the selected pixel regions, this will be appended to file and 
        dir names so it acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    """
    region_info = ''

    if evt.inst == 'resolve':
        params = {
            'infile': evt.path,
            'splitrmf': splitrmf,
            'splitcomb': splitcomb,
            'whichrmf': whichrmf,
            'cleanup': 'no',
            'history': 'no'
        }

        short_inst = INST_NAMES[evt.inst]
        obs_md = INST_MODES[short_inst][evt.obs_mode]

        prefix = f'{evt.obs_id}{short_inst}_{obs_md}'
        # That means use all pixels
        if grades is None:
            params['resolist'] = 'ALL'
        else:
            params['resolist'] = ",".join([str(GRADES[grd]) for grd in grades])
            prefix += f"_grd-{''.join([str(GRADES[grd]) for grd in grades])}"

        if include_pix is None and exclude_pix is None:
            params['regionfile'] = 'ALLPIX'
            params['regmode'] = 'DET'
        else:
            params['pixlist'] = gen_pix_sel_str(include_pix, exclude_pix, separator='-')
            params['regionfile'] = 'NONE'
            prefix += f"_pixreg-{pix_name}"
            
            region_info = {}

            region_info['include_pix'] = include_pix
            region_info['exclude_pix'] = exclude_pix
            region_info['pix_name'] = pix_name
            region_info['pix_reg'] = True

        if regionname is not None:
            params['regionfile'] = regionfile
            params['regmode'] = regmode
            warn(f"Double check that the coord system given by regmode:{regmode} match the one "
                f"used by the regionfile provided: {regionfile}")
            prefix += f"_reg-{regionname}"

            region_info = {}
            # If region is provided, read its content
            with open(regionfile, 'r') as f:
                region_content = f.read()
            
            region_info['name'] = regionname
            region_info['file'] = regionfile
            region_info['content'] = region_content
            region_info['pix_reg'] = False
        
        params['outfileroot'] = prefix
        heasoft_cmd = 'rslmkrmf'

    else:
        prefix = f'{evt.obs_id}{short_inst}_{obs_md}'
        outfile = prefix + '.rmf'
        params = {'infile': evt.path,
                'outfile': outfile,
                'cleanup': 'no',
                'history': 'no'}
        
        heasoft_cmd = 'xtdrmf'

    return params, heasoft_cmd, prefix, region_info, 'rmf' 


def make_spectra_inc_rsp_and_nxb(evt: EventList, ra: float, dec: float, whichrmf: str, 
                                 splitrmf: str = 'no', splitcomb: str = 'no', 
                                 min_energy: int = None, max_energy: int = None, 
                                 include_pix: list[int] = None, exclude_pix: list[int] = None, 
                                 pix_name: str = None, regmode: str = 'DET', 
                                 regionfile: str = None, regionname: str = None, 
                                 grades: list[str] = None, numphoton: int = 60000, 
                                 nxbfile: str = None, nxbehkfile: str = None, cortime: int = None,
                                 apply_rsl_rise_time: bool = True, status: list[int] = None,
                                 timefirst: int = 150, timelast: int = 150) -> dict:
    """
    Makes a spectra including all the calibration products (arf and rmf). Can also generate NXB if
    argument nxbfile is set.

    Paramters
    ---------
    evt : EventList
        The EventList object to generate spectrum and calibration products from
    ra : float
        source ra in deg
    dec : float 
        source dec in deg
    whichrmf : str
        either 'S', 'M', 'L', 'X'
    splitrmf : str
        If 'splitrmf=yes', split the RMF into core and ELC responses.
    splitcomb : str
        If 'splitrmf=yes' and 'splitcomb=yes', the two response matrices (corresponding to the core 
        and ELC components) will be placed into a single file.
    min_energy : str
        in pha, the min energy cut
    max_energy : str
        in pha, the max energy cut
    include_pix : list[int]
        list of pixel indicies to include in products.
    exclude_pix : list[int]
        list of pixel indicies to exclude from products.
    pix_name : str
        this is like a nickname for the selected pixel regions, this will be appended to file and 
        dir names so it acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    regmode : str
        coordinate system that the regionfile is in
    regionfile : str
        absolute path to regionfile used to filter the product
    regionname : str
        this is like a nickname for the region, this will be appended to file and dir names so it 
        acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    grades : list[str]
        list of grades to include when generating products. ie, should be given as ['Hp', 'Mp'] not 
        [0, 1]. 
    numphoton : int
        the number of photons to used in the raytrace simulation
    nxbfile : str
        absolute path to the nxb eventlist
    nxbehkfile : str
        absolute path to the nxb ehk file
    cortime : int
        minmum cortime allowed in ehk files when creating gtis to filter event lists by.
    apply_rsl_rise_time: bool
        apply recommended rise time screening to resolve event lists
    status: list[int]
        event statuses to exclude from eventlists that spectra and calibration products are made
        from
    timefirst : int
        days before the observation to include from the nxb eventlist in the final spectrum
    timelast : int
        days after the observation to include from the nxb eventlist in the final spectrum
    
    Return
    ------
    files_made : dict
        dictionary containing the absolute paths to relevant products made. They keys are 
        'spectrum', 'expmap', 'arf', 'rmf' and if the nxbfile argument was set there will also be 
        and key for 'nxb'
    """
    files_made = {}

    if nxbfile is not None:
        evt = make_nxb_spec_inc_filtering(evt, nxbfile, nxbehkfile, cortime=cortime, 
                                 apply_rsl_rise_time=apply_rsl_rise_time, min_energy=min_energy, 
                                 max_energy=max_energy, status=status, 
                                 grades=grades, timefirst=timefirst, timelast=timelast,
                                 include_pix=include_pix, exclude_pix=exclude_pix, 
                                 pix_name=pix_name)
        # Constructing the name of the nxb spectrum made
        short_inst = INST_NAMES[evt.inst]
        obs_md = INST_MODES[short_inst][evt.obs_mode]
        pix_name_str = pix_name if pix_name is not None else 'ALLPIX'
        cortime_name = f'CORT-{cortime}' if cortime is not None else 'noCORFILT'
        prefix = f'{evt.obs_id}{short_inst}_{obs_md}_pix-{pix_name_str}_{cortime_name}'
        nxb_abs_path = ANALYSIS_DIR + '/nxb/' + prefix + '/' + prefix + '_nxb.pi' 
        files_made['nxb'] = nxb_abs_path

    spec_files = make_spectra(evt, min_energy=min_energy, max_energy=max_energy, 
                              include_pix=include_pix, exclude_pix=exclude_pix, 
                 pix_name=pix_name, regmode=regmode, regionfile=regionfile, regionname=regionname, 
                 grades=grades)
    spectrum_file = [f for f in spec_files if f.endswith('spectrum.fits')][0]
    spectrum_abs_path = ANALYSIS_DIR + '/spectrum/' + spectrum_file.removesuffix('_spectrum.fits') \
                        + '/' + spectrum_file
    files_made['spectrum'] = spectrum_abs_path
    
    expmap_files = make_expmap(evt, outmaptype='EXPOSURE', stopsys='SKY')
    expmap_file = [f for f in expmap_files if f.endswith('.fits')][0]
    expmap_abs_path = ANALYSIS_DIR + '/expmap/' +  expmap_file.removesuffix('.fits') + '/' \
                      + expmap_file 
    files_made['expmap'] = expmap_abs_path

    rmf_files = make_rmf(evt, whichrmf, splitrmf=splitrmf, splitcomb=splitcomb, regmode=regmode, 
                        grades=grades, regionfile=regionfile, regionname=regionname, 
                        include_pix=include_pix, exclude_pix=exclude_pix, pix_name=pix_name)
    rmf_file = [f for f in rmf_files if f.endswith('.rmf')][0]
    rmf_abs_path = ANALYSIS_DIR + '/rmf/' + rmf_file.removesuffix('.rmf') + '/' + rmf_file
    files_made['rmf'] = rmf_abs_path
    
    arf_files = make_arf(evt, ra, dec, expmap_abs_path, rmf_abs_path, regmode=regmode, 
                         regionfile=regionfile, regionname=regionname, numphoton=numphoton, 
                         include_pix=include_pix, exclude_pix=exclude_pix, pix_name=pix_name)
    arf_file = [f for f in arf_files if f.endswith('.arf')][0]
    arf_abs_path = ANALYSIS_DIR + '/arf/' + arf_file.removesuffix('.arf') + '/' + arf_file
    files_made['arf'] = arf_abs_path

    return files_made

#@heasoftpy_run
@cmd_run
def make_nxb_spec(evt: EventList, nxbfile: str, nxbehkfile: str, timefirst: int = 150, 
                  timelast: int = 150, include_pix: list[int] = None, exclude_pix: list[int] = None, 
                  pix_name: str = None, regmode: str = 'DET', regionfile: str = None, 
                  regionname: str = None):
    """
    Generates a NXB spectrum. This assumes that screening and gtifiltering has already taken place.
    Use function make_nxb_spec_inc_filtering if you want to do that all in one go.

    Parameters
    ----------
    evt : EventList
        the eventlist to generate a nxb spectrum for
    innxbfile : str
        absolute path to the nxb eventlist
    innxbehkfile : str
        absolute path to the nxb ehk file
    timefirst : int
        days before the observation to include from the nxb eventlist in the final spectrum
    timelast : int
        days after the observation to include from the nxb eventlist in the final spectrum
    include_pix : list[int]
        list of pixel indicies to include in products.
    exclude_pix : list[int]
        list of pixel indicies to exclude from products.
    pix_name : str
        this is like a nickname for the selected pixel regions, this will be appended to file and 
        dir names so it acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    regmode : str
        coordinate system that the regionfile is in
    regionfile : str
        absolute path to regionfile used to filter the product
    regionname : str
        this is like a nickname for the region, this will be appended to file and dir names so it 
        acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    """
    # XTEND
    # Define two regions one for actual source region nxb, and another for a larger region, as the
    # nxb is spatially varying, but also need sufficient counts for good statistics
    # can use the whole area minus calibration sources as the region that is the region to get counts from
    # nxb region needs to be in DET coords
    # want two regions to be in DET coords as well
    
    
    region_info = ''

    # Will be overwritten if there was a cortime filter applied
    cortime_str = '0,4,5,6,7,8,9,10,11,12,13,99'
    cortime_name = 'noCORFILT'
    # Getting the cortime info
    if evt.gti_applied:
        df = pd.read_csv(ANALYSIS_DIR + '/gti/gti_lookup.csv')

        matched = df[
            (df['gtifile'].apply(lambda x: x.endswith(evt.gtifile))) &
            (df['ehkfile'] == evt.ehkfile)
        ]

        expr = matched['expr'].values[0]

        match = re.search(r'CORTIME>(-?\d+\.?\d*)', expr)
        if match:
            number = int(match.group(1)) # This will always be a lower limit
            cortime_str = f'{number}, 99'
            cortime_name = f'CORT-{number}'

    # Then we just need to format the command for the event list
    short_inst = INST_NAMES[evt.inst]
    obs_md = INST_MODES[short_inst][evt.obs_mode]

    prefix = f'{evt.obs_id}{short_inst}_{obs_md}_{cortime_name}'
    outpifile = prefix + '_nxb.pi'
    outehkfile = prefix + '.ehk'


    if evt.inst.upper() == 'RESOLVE':
        heasoft_cmd = 'rslnxbgen'
        if pix_name is not None:
            pix_str = gen_pix_sel_str(include_pix=include_pix, exclude_pix=exclude_pix)
            pix_name_str = pix_name

            region_info = {}
            region_info['include_pix'] = include_pix
            region_info['exclude_pix'] = exclude_pix
            region_info['pix_name'] = pix_name
            region_info['pix_reg'] = True

        else:
            pix_str = '-'
            pix_name_str = 'ALLPIX'  

        prefix = prefix + f'_pix-{pix_name_str}'

        # there is a bug in heasoftpy which means this cant be used, but im keeping it incase the 
        # bug is fixed
        params = {
            'infile': evt.path,
            'ehkfile': evt.ehkfile,
            'regfile': 'NONE',
            'pixels': pix_str,
            'innxbfile': nxbfile,
            'innxbehk': nxbehkfile,
            'timefirst': f'-{timefirst}',
            'timelast': f'+{timelast}',
            'SORTCOL': 'CORTIME',
            'sortbin': cortime_str, # needs to match cor filtering on science data
            'outpifile': outpifile, 
            'outnxbehk': outehkfile,
            'database': 'LOCAL',
            'db_location': './'

        }

        cmd = (f'rslnxbgen infile={evt.path} ehkfile={evt.ehkfile} regfile=NONE pixels="{pix_str}" '
               f'innxbfile={nxbfile} innxbehk={nxbehkfile} timefirst=-{timefirst} '
               f'timelast=+{timelast} SORTCOL=CORTIME sortbin="{cortime_str}" outpifile={outpifile} '
               f'outnxbehk={outehkfile} database=LOCAL db_location=./ clobber=yes')

        cmd = 'source ~/headas.sh ;' + cmd
    
    else:
        heasoft_cmd = 'xtdnxbgen'

        if regionname is not None:
            use_regfile = regionfile
            warn(f"Double check that the coord system given by regmode:{regmode} match the one "
                f"used by the regionfile provided: {regionfile}")
            prefix += f"_reg-{regionname}"

            region_info = {}
            # If region is provided, read its content
            with open(regionfile, 'r') as f:
                region_content = f.read()
            
            region_info['name'] = regionname
            region_info['file'] = regionfile
            region_info['content'] = region_content
            region_info['pix_reg'] = False
        
        else:
            use_regfile = 'NONE'

        cmd = ('source ~/headas.sh; '
               f'xtdnxbgen infile={evt.path} ehkfile={evt.ehkfile} regmode={regmode} '
               f'regfile={use_regfile} innxbfile={nxbfile} innxbehk={nxbehkfile} apply_xtdtools=no '
               f'database=LOCAL db_location=./ timefirst=-{timefirst} timelast=+{timelast} '
               f'SORTCOL=CORTIME sortbin="{cortime_str}" '
               f'outpifile={outpifile} outnxbehk={outehkfile}')
    
        cmd = 'source ~/headas.sh ;' + cmd

    return cmd, prefix, region_info, 'nxb'

def make_nxb_spec_inc_filtering(evt: EventList, nxbfile: str, nxbehkfile: str, cortime: int = None, 
                                 apply_rsl_rise_time: bool = True, min_energy: int = None, 
                                 max_energy: int = None, status: list[int] = None, 
                                 grades: list[str] = None, remove_old_CI_rows: bool = False, 
                                 timefirst: int = 150, timelast: int = 150, 
                                 include_pix: list[int] = None, exclude_pix: list[int] = None, 
                                 pix_name: str = None, regmode: str = 'DET', regionfile: str = None, 
                                 regionname: str = None):
    """
    Applies earth night and optionally CORTIME fitering to make GTIs. GTIs are applied to 
    eventlists. The applies screening conditions to eventlists which are then used to generate a 
    NXB spectrum. 

    Parameters
    ----------
    evt : EventList
        the eventlist to generate a nxb spectrum for
    nxbfile : str
        absolute path to the nxb eventlist
    nxbehkfile : str
        absolute path to the nxb ehk file
    cortime : int
        minmum cortime allowed in ehk files when creating gtis to filter event lists by.
    apply_rsl_rise_time: bool
        (RESOLVE only) apply recommended rise time screening to resolve event lists
    min_energy : str
        in pha, the min energy cut for the spectrum (nxb and source)
    max_energy : str
        in pha, the max energy cut for the spectrum (nxb and source)
    status: list[int]
        event statuses to exclude from eventlists that spectra and calibration products are made
        from
    grades : list[str]
        (RESOLVE only) Accepts the strings 'Hp', 'Mp', 'Ms', 'Lp', 'Ls'.
    remove_old_CI_rows : bool
        (XTEND only), for NXB generation purposes, current NXBDB contains NTE data only before 
        the change of CI rows, if this argument is set to True, this removes all CI rows and 
        preceding/trailing rows to make source spectrum compatible with the NXB spectrum
    timefirst : int
        days before the observation to include from the nxb eventlist in the final spectrum
    timelast : int
        days after the observation to include from the nxb eventlist in the final spectrum
    include_pix : list[int]
        list of pixel indicies to include in products.
    exclude_pix : list[int]
        list of pixel indicies to exclude from products.
    pix_name : str
        this is like a nickname for the selected pixel regions, this will be appended to file and 
        dir names so it acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    regmode : str
        coordinate system that the regionfile is in
    regionfile : str
        absolute path to regionfile used to filter the product
    regionname : str
        this is like a nickname for the region, this will be appended to file and dir names so it 
        acts like a quick lookup so you know the region filtering that has happened
        just from reading the file names.
    """
    
    # Need to apply additional screening to NXB data, so that it doesnt include earth day 
    # also need to apply same screening that you do to your source data
    # Produce GTI that gives you night earth also could then optionally screen for CORTIME this 
    # should be applied to NXB and source
    nxbgtifile = make_gti(nxbehkfile, inst=evt.inst, source=False)
    gtifile = make_gti(evt.ehkfile, inst=evt.inst, source=True, cortime=cortime)

    # Apply GTI and screening to source
    gti_filt_evt = evt.apply_gti(gtifile)
    screened_evt = gti_filt_evt.screen(apply_rsl_rise_time=apply_rsl_rise_time, min_pi=min_energy, 
                                    max_pi=max_energy, status=status, grades=grades,
                                    remove_old_CI_rows=remove_old_CI_rows)

    # Apply GTI and screening to NXB
    nxb_evt = EventList(nxbfile, obs_id=evt.obs_id, inst=evt.inst, obs_mode=evt.obs_mode,
                        root_dir=evt.root_dir)

    nxb_evt.ehkfile = nxbehkfile
    gti_filt_nxb_evt = nxb_evt.apply_gti(nxbgtifile)
    screened_nxb_evt = gti_filt_nxb_evt.screen(apply_rsl_rise_time=apply_rsl_rise_time, 
                                                min_pi=min_energy, max_pi=max_energy, status=status, 
                                                grades=grades)

    make_nxb_spec(screened_evt, screened_nxb_evt.path, nxbehkfile, timefirst=timefirst, 
                  timelast=timelast, include_pix=include_pix, exclude_pix=exclude_pix, 
                  pix_name=pix_name, regmode=regmode, regionfile=regionfile, regionname=regionname)  
    
    return screened_evt