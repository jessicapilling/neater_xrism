from run import heasoftpy_run
from classes import EventList
from constants import INST_MODES, INST_NAMES

@heasoftpy_run
def run_rslbratios(evt: EventList, filetype: str = 'cl', lcurve: bool = False):
    """
    Runs rslbratios.
    """

    region_info = ''
    heasoft_cmd = 'rslbratios'

    short_inst = INST_NAMES[evt.inst]
    obs_md = INST_MODES[short_inst][evt.obs_mode]

    get_lcurve = 'yes' if lcurve else 'no'

    prefix = f'{evt.obs_id}{short_inst}_{obs_md}_{filetype}_lc-{int(lcurve)}'
    params = {
        'infile': evt.path,
        'filetype': filetype,
        'outroot': prefix,
        'lcurve': get_lcurve
    }

    return params, heasoft_cmd, prefix, region_info, 'rslbratio'
