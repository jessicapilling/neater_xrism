import os

from constants import INST_NAMES, INST_MODES
from run import heasoftpy_run
from utils import _generic_xselect
from classes import EventList

def make_image(evt: EventList, min_energy: int = None, max_energy: int = None, regmode: str = 'DET', 
               regionfile: str = None, regionname: str = None, include_pix: list[int] = None, 
               exclude_pix: list[int] = None, pix_name: str = None) -> list[str]:
    """
    Make images from XRISM event lists. Note about the region filtering: for resolve, set the 
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
    
    Returns
    -------
    listdir : list[str]
        list of filenames that were created by the command run
    """  
    return _generic_xselect(evt, 'image', min_energy=min_energy,
                     max_energy=max_energy, regmode=regmode, regionfile=regionfile, 
                     regionname=regionname, include_pix=include_pix, exclude_pix=exclude_pix, 
                     pix_name=pix_name)

@heasoftpy_run
def make_expmap(evt: EventList, outmaptype: str = 'EXPOSURE', stopsys: str = 'SKY') -> list[str]:
    """
    Makes exposure maps for XRISM eventlists.

    Parameters
    ----------
    evt : EventList
        The EventList object to generate an exposure map for.
    outmaptype : str
        the type of exposure map to be generated. The option 'outmaptype=EXPOSURE' is needed for 
        the xaarfgen task, because this option creates additional extensions that contain pixel 
        partial exposure information for each attitude bin. The 'outmaptype=EFFICIENCY' option is 
        more appropriate for correcting images for detector and exposure time effects.
    stopsys : str
        Output coordinate system for the exposure map image. Options are: DET, FOC, SKY

    Returns
    -------
    listdir : list[str]
        list of filenames that were created by the command run
    """
    region_info = ''

    # Then we just need to format the command for the event list
    short_inst = INST_NAMES[evt.inst]
    obs_md = INST_MODES[short_inst][evt.obs_mode]

    pixgtifile = f'xa{evt.obs_id}{short_inst}_{obs_md}_exp.gti.gz'
    pixgtipath = os.path.join(evt.root_dir, evt.obs_id, evt.inst, 'event_uf', pixgtifile)
    ehkfile = f'xa{evt.obs_id}.ehk.gz'
    ehkpath = os.path.join(evt.root_dir, evt.obs_id, 'auxil', ehkfile)
    prefix = f'{evt.obs_id}{short_inst}_{obs_md}_regmode-{stopsys}_type-{outmaptype}'
    expmapfile = prefix + '.fits'
    logfile = prefix + '.log'
    
    params = {
        'ehkfile': ehkpath,
        'gtifile': evt.path,
        'instrume': evt.inst,
        'badimgfile': 'NONE',
        'pixgtifile': pixgtipath,
        'outfile': expmapfile,
        'stopsys': stopsys,
        'logfile': logfile,
        'cleanup': 'no',
        'history': 'no'}
    
    heasoft_cmd = 'xaexpmap'
    cmd_dir = prefix
    
    return params, heasoft_cmd, cmd_dir, region_info, 'expmap'

        



