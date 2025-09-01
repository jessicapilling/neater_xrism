import os
import shutil

import pandas as pd

from utils import maketime, _generic_xselect
from conf import ANALYSIS_DIR
from classes import EventList

def make_lc(evt: EventList, min_energy: int = None, max_energy: int = None, exp: float = 0.8,
            binsize: int = 128.0, include_pix: list[int] = None, exclude_pix: list[int] = None, 
            pix_name: str = None, regmode: str = None, regionfile: str = None, 
            regionname: str = None, grades: list[str] = None):
    """
    Python wrapper for xselect that makes lightcurves from XRISM event lists.

    Parameters
    ---------
    evt : EventList
        The EventList object to be read into xselect
    min_energy : str
        in pha, the min energy cut
    max_energy : str
        in pha, the max energy cut
    exp : str
        the exposure parameter controls the rejection of bins from the light curve, based 
        on how much of the bin lies within the GTI's. If exposure=1, then a bin is only written to 
        the light curve if it lies wholly within a GTI. If exposure=0, then bins are written to the 
        light curve if any portion of them lies within a GTI. Exposure can have any value between 
        these two extremes. Of course, if a bin does not fall inside any of the GTI's, it will not 
        be written to the light curve, no matter what the value of exposure is.
    binsize : int
        sets time binsize
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

    """
    return _generic_xselect(evt, 'curve', min_energy=min_energy,
                    max_energy=max_energy, exp=exp, binsize=binsize, include_pix=include_pix,          
                    exclude_pix=exclude_pix, pix_name=pix_name, regmode=regmode, 
                    regionfile=regionfile, regionname=regionname, grade=grades)

def make_gti(ehkfile: str, inst: str, source: bool = True, cortime: int = None):
    """
    Making a .gti file based on a selection of parameters from an ehkfile.
    This is meant to be used in the context of making a NXB.

    Parameters 
    ----------
    ehkfile: str
        absolute path to ehk file
    source : bool
        True if generating a gti from source ehk file, False if generating a gti from a nxb ehk 
        file.
    cortime: int
        cortime above which will be included in gtis. 
    """
    if inst.upper() == 'RESOLVE':
        # if we are extracting a source gti or an nxb gti we need slightly different expressions
        expr =  "T_SAA_SXS>0 && DYE_ELV>5" if source else "T_SAA_SXS>0 && ELV<-5 && DYE_ELV>5"
    else:
        # For Xtend there is a different instrument
        expr = "SAA==0 && T_SAA_SXI>300" if source else "SAA==0 && ELV<-5 && DYE_ELV>100 && T_SAA_SXI>300"

    if cortime is not None:
        expr = expr + f" && CORTIME>{cortime}"

    # if this dir already exists, it means that gti making has happened before, so need to check
    # that the name we use for gti hasnt been used before
    if os.path.exists(ANALYSIS_DIR + '/gti'):
        gtis = [e for e in os.listdir(ANALYSIS_DIR + '/gti') if e.endswith('.gti')]
        ident = [int(gti.split('.gti')[0][-1]) for gti in gtis]
        highest_id = sorted(ident)[-1]
        # This is the number that will go on sel{new_id}.gti
        new_id = str(highest_id + 1)

    else:
        os.makedirs(ANALYSIS_DIR + '/gti', exist_ok=True)
        new_id = '1'
    
    outfile = ANALYSIS_DIR + f"/gti/sel{new_id}.gti"

    # if this already exists, it means that gti making has happened before, so need to check
    # that the gti making with current arguments hasnt been run before 
    table_path = ANALYSIS_DIR + "/gti/gti_lookup.csv"
    if os.path.exists(table_path):
        df = pd.read_csv(table_path)

        # Check if row already exists
        existing = df[
        (df["ehkfile"] == ehkfile) &
        (df["expr"] == expr) &
        (df["inst"] == inst)
            ]
        
        # If row does exist, screening has been done before, so no need to perform ftcopy
        if not existing.empty:
            file  = existing["gtifile"].values[0]
            print(f"Gti has already been made: {file}")

            return file

    else:
        df = pd.DataFrame(columns=["ehkfile", "expr", "inst", "gtifile"])
    
    # making the gti
    maketime(ehkfile, outfile, expr)

    # appending this screening to the df
    new_row = {"ehkfile": ehkfile,
            "expr": expr,
            "inst": inst,
            "gtifile": outfile}

    df = df.append(new_row, ignore_index=True)

    # Save updated table
    df.to_csv(table_path, index=False)

    return outfile

