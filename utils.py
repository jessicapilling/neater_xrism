import os 

import pandas as pd

from classes import EventList
from constants import PIXEL_REGS, GRADES, INST_MODES, INST_NAMES, HEASOFT_CMDS
from conf import ANALYSIS_DIR
from run import xselect_run, small_heasoftpy

@small_heasoftpy
def extractor(filename: str, outfile: str, timefile: str):
    """
    Used to apply gtifiles to eventlists.

    Parameters
    ----------
    filename : str
        events to have gti applied to
    outfile : str
        new event lists that have had gti applied to them
    timefile : str
        gtis to be applied
    """
    params = {
        'filename': filename,
        'eventsout': outfile,
        'imgfile': 'NONE',
        'phafile': 'NONE',
        'fitsbinlc': 'NONE',
        'regionfile': 'NONE',
        'timefile': timefile,
        'xcolf': 'X',
        'ycolf': 'Y',
        'tcol': 'TIME',
        'ecol': 'PI',
        'xcolh': 'DETX',
        'ycolh': 'DETY',
        'noprompt': True}

    return 'extractor', params

@small_heasoftpy
def maketime(infile: str, outfile: str, expr: str) -> str:
    """
    Used to make gtifiles based upon parameters given by ehk files.

    Parameters
    ----------
    infile : str
        ehkfile to be used to make gti files
    outfile : str
        the name of the gtifile to be made
    expr : str
        the conditions to filter the ehk file by to generate gtis for

    Return
    ------
    outfile : str
        the name of the gtifile made
    """

    params = {'infile': infile,
              'outfile': outfile,
              'expr': expr,
              'compact': 'no',
              'time': 'TIME',
             'noprompt': True}

    return 'maketime', params


def ftcopy(evtpath: str, condition: str, outfile: str) -> str:   
    """
    Used to screen eventlists. Runs the heasoftpy cmd ftcopy to apply a screening to eventlists.

    Parameters
    ----------
    evtpath : str
        absolute path to eventlists that will have screening applied to.
    condition : str
        the expression given to ftcopy that gives the conditions to filter the event lists by
    outfile : str
        the name of the new event file that has had screening applied
    
    Return
    ------
    outfile : str
        the name of the new event file that has had screening applied
    """ 
    og_cwd = os.getcwd()
    os.chdir(ANALYSIS_DIR)

    infile = evtpath + '[EVENTS]' + condition

    params = {'infile': infile,
              'outfile': outfile}
    
    result = HEASOFT_CMDS['ftcopy'](params)

    print(result.stdout)
    if result.stderr is not None:
        os.chdir(og_cwd)
        raise KeyError(result.stderr)

    os.chdir(og_cwd)
    return outfile


def generate_pixreg_file(filename: str, include_indices: list[int] = None, 
                         exclude_indices: str = None) -> str:
    """
    The argument 'regionfile' is input into xaarfgen, and should include all the pixels to be 
    included in the ARF calculation for resolve ARFs. This is a convenience function so you can just 
    give a list of pixel indicies to include or exclude, and the regionfile will be generated. NB: 
    Pixel indices input into this function should be its name (ie 0-35). It's not wise to specify 
    include_indices and exclude, I havent tested what will happen - probably not what you want.

    Parameters
    ----------
    filename : str
        name of the regionfile to be made.
    include_indicies : list[int]
        list of indicies of pixels to be included in analysis
    exclude_indicies : list[int]
        list of indicies of pixels to be excluded from analysis

    Returns
    -------
    rtype : str
        abs path of created region file
    """
    # collecting the region of each pixel in DET coords
    if include_indices is not None:
        selected_regions = [PIXEL_REGS[i] for i in include_indices if i in PIXEL_REGS]
    else:
        selected_regions = list(PIXEL_REGS.values())
    
    if exclude_indices is not None:
        selected_regions = [region for i, region in PIXEL_REGS.items() if i not in exclude_indices]
    
    og_cwd = os.getcwd()
    os.chdir(ANALYSIS_DIR)

    # writing the reg file
    with open(filename, 'w') as f:
        f.write("physical\n")
        for region in selected_regions:
            f.write(region + "\n")
    
    print(f"Region file '{filename}' generated successfully.")

    os.chdir(og_cwd)

    return os.path.join(ANALYSIS_DIR, filename)


def gen_pix_sel_str(include_pix: list[int] = None, exclude_pix: list[int] = None, 
                    separator: str = ':') -> str:
    """
    Lots of product generation functions take an input where you need to provide the pixel selection
    so this function takes the pixels you want to exclude/include and converts them into a string
    that is readable by these product generation functions.

    Parameters
    ----------
    include_pix : list[int]
        pixels to include in analysis
    exclude_pix : list[int]
        pixels to EXclude in analysis >:(
    separator : str
        the character to separate the ranges of pixels with in the return. (some functions use 
        different ones)
    
    Returns
    -------
    rtype : str
        the formatted string, based on given pixel selection
    """
    # Default to full range if include_pix is None
    if include_pix is None:
        include_pix = set(range(36))
    else:
        include_pix = set(include_pix)

    # Default to an empty set if exclude_pix is None
    exclude_pix = set(exclude_pix) if exclude_pix is not None else set()

    # Get valid pixels (0 to 35) that are included but not excluded
    valid_pixels = sorted([p for p in range(36) if p in include_pix and p not in exclude_pix])

    # If no pixels are left, return an empty string
    if not valid_pixels:
        raise ValueError("The include and exclude pixels argument have been set so there are no " 
                         "pixels left!")

    # Find contiguous ranges
    ranges = []
    start = valid_pixels[0]

    for i in range(1, len(valid_pixels)):
        if valid_pixels[i] != valid_pixels[i - 1] + 1:
            # End previous range
            end = valid_pixels[i - 1]
            ranges.append(f"{start}{separator}{end}" if start != end else f"{start}")
            # Start new range
            start = valid_pixels[i]

    # Append final range
    end = valid_pixels[-1]
    ranges.append(f"{start}{separator}{end}" if start != end else f"{start}")

    return ", ".join(ranges)

@xselect_run
def _generic_xselect(evt: EventList, product: str, min_energy: int = None, max_energy: int = None, 
                     exp: str = None, binsize: int = None, include_pix: list[int] = None, 
                     exclude_pix: list[int] = None, pix_name: str = None, regmode: str = None,
                    regionfile: str = None, regionname: str = None, grade: list[str] = None):
    """
    Since scripts to extract images, lightcurves, and spectra using Xselect are very similar, I have
    combined this into one function. Note about the region filtering: for resolve, set the variables:
    include_pix OR exclude_pix and pix_name, for xtend set the variables regmode and regionfile and
    regionname. You always need to set either pix_name or region_name if you want to do any region
    filtering. Note that this code doesn't check that regionfile contents are in actually in the 
    coords specified by regmode.

    Parameters
    ---------
    evt : EventList
        The EventList object to be read into xselect
    product : str
        the type of product that will be generated from this run
    min_energy : str
        in pha, the min energy cut
    max_energy : str
        in pha, the max energy cut
    exp : str
        (lc only) the exposure parameter controls the rejection of bins from the light curve, based 
        on how much of the bin lies within the GTI's. If exposure=1, then a bin is only written to 
        the light curve if it lies wholly within a GTI. If exposure=0, then bins are written to the 
        light curve if any portion of them lies within a GTI. Exposure can have any value between 
        these two extremes. Of course, if a bin does not fall inside any of the GTI's, it will not 
        be written to the light curve, no matter what the value of exposure is.
    binsize : int
        (lc only) sets time binsize
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
    grade : list[str]
        list of grades to include when generating products. ie, should be given as ['Hp', 'Mp'] not 
        [0, 1]. 
    """   

    cmd_prefix = ''' #!/bin/sh
xselect << EOF
{name}
read eve {evt}
{dir}
set image {regmode}
'''
    cmd_suffix='''extr {product} {exp}
save {product} {outfile}
exit
no
EOF
'''
    region_info = ''

    # Then we just need to format the command for the event list
    short_inst = INST_NAMES[evt.inst]
    obs_md = INST_MODES[short_inst][evt.obs_mode]
    
    prefix = f'{evt.obs_id}{short_inst}_{obs_md}'

    if binsize is not None:
        cmd_prefix += f'set binsize {binsize} \n'
        prefix += f'_bin-{binsize}'
    
    if min_energy is not None or max_energy is not None:
        use_min = min_energy if min_energy is not None else 0
        use_max = max_energy if max_energy is not None else 26000
        cmd_prefix += f'filter pha_cutoff {use_min} {use_max} \n'
        prefix += f'_en-{use_min}-{use_max}'
    
    if grade is not None:
        if not all(g in GRADES.keys() for g in grade):
            raise ValueError('Youve entered a dodgy grade.')
        # Convert grade names to numeric values, ignoring invalid keys
        values = sorted([GRADES[g] for g in grade])

        # Find contiguous ranges
        ranges = []
        start = values[0]

        for i in range(1, len(values)):
            if values[i] != values[i - 1] + 1:
                # End previous range
                end = values[i - 1]
                ranges.append(f"{start}:{end}" if start != end else f"{start}")
                # Start new range
                start = values[i]

        # Append final range
        end = values[-1]
        ranges.append(f"{start}:{end}" if start != end else f"{start}")

        # Join into final string
        grade_filt_string =  ", ".join(ranges)
        prefix += f"_grd-{''.join([str(GRADES[grd]) for grd in grade])}"
        cmd_prefix += f'filter GRADE {grade_filt_string} \n'

    if regionname is not None:
        region_info = {}

        if not os.path.exists(regionfile):
            raise ValueError("Regionfile given doesn't exist.")

        cmd_prefix += f'filter region {regionfile} \n'
        prefix += f'_reg-{regionname}'

        # If region is provided, read its content
        with open(regionfile, 'r') as f:
            region_content = f.read()
        
        region_info['name'] = regionname
        region_info['file'] = regionfile
        region_info['content'] = region_content
        region_info['pix_reg'] = False
    
    if pix_name is not None:
        pixel_filt_string = gen_pix_sel_str(include_pix, exclude_pix)
        cmd_prefix +=  f'filter column "pixel={pixel_filt_string}" \n'
        prefix += f'_pixreg-{pix_name}'
        region_info = {}
        region_info['include_pix'] = include_pix
        region_info['exclude_pix'] = exclude_pix
        region_info['pix_name'] = pix_name
        region_info['pix_reg'] = True

    if exp is None:
        exp = ''
    else:
        exp = f'exposure={exp}'

    outfile = prefix + f'_{product}.fits'

    cmd_prefix = cmd_prefix.format(name=prefix,
                                    evt=evt.path.split('/')[-1],
                                    dir = "/".join(evt.path.split('/')[:-1]),
                                    regmode=regmode)
    
    cmd = cmd_prefix + cmd_suffix.format(product=product, exp=exp, outfile=outfile)
    
    return cmd, prefix, region_info, product