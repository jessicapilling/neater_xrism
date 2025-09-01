import os
import json

import pandas as pd

from constants import INST_MODES, INST_NAMES, GRADES
from conf import ANALYSIS_DIR

class EventList():
    """
    For single event lists (ie. for a particular inst and obs mode)

    Parameters
    ----------
    path : str
        Absoulte path of eventlist.
    obs_id : str
        Ob_id of eventlist.
    inst : str
        either 'resolve' or 'xtend'
    obs_mode : str
        Given in human readable format ie. 'open' instead of 'px1000'. The translation from 
        readable to code is in the INST_MODES dictionary in constants.py
    root_dir : str
        the absolute path to the root_dir, ie the dir containing the obs_id dir that contains the 
        auxil, pfiles, log, resolve, xtend dirs
    """
    def __init__(self, path: str, obs_id: str, inst: str, obs_mode: str, root_dir: str):
        self.inst = inst
        self.path = path
        self.obs_id = obs_id
        self.obs_mode = obs_mode
        self.root_dir = root_dir
        self.gti_applied = False
        self.gtifile = None
        self.screening_applied = False
        self.screening_info = None
        self.ehkfile = root_dir + f'/{self.obs_id}/auxil/xa{self.obs_id}.ehk.gz' 
        
    
    def _gen_screen_condition(self, apply_rsl_rise_time: bool = True, min_pi: int = None, 
                              max_pi: int = None, status: list[int] = None, 
                              grades: list[str] = None, remove_old_CI_rows: bool = False):
        """
        Makes a string in the logical filtering format that can be used to filter event lists with 
        ftcopy. ie for _gen_screen_condition(grades=['Hp']) returns [(ITYPE==0)]. Link to what
        different status numbers mean:
        https://heasarc.gsfc.nasa.gov/docs/xrism/analysis/abc_guide/XRISM_Data_Specifics.html#SECTION00710000000000000000

        Parameters
        ----------
        apply_rsl_rise_time : bool
            There is a reccomended filtering on rise time for resolve event lists which is applied
            if this is set to True.
        min_pi : int
            minimum PI for events
        max_pi : int 
            maximum PI for events
        status : list[int]
            should be input as a list, even if there is just one status to be filtered. (see the 
            link in the docstring for proper info about what status number means what)
        grades : list[str]
            Accepts the strings 'Hp', 'Mp', 'Ms', 'Lp', 'Ls'.
        remove_old_CI_rows : bool
            (XTEND only), for NXB generation purposes, current NXBDB contains NTE data only before 
            the change of CI rows, if this argument is set to True, this removes all CI rows and 
            preceding/trailing rows to make source spectrum compatible with the NXB spectrum
        """
        # Conditions to be appended here 
        conditions = []

        # Handle PI conditions
        if min_pi is not None:
            conditions.append(f"(PI >= {min_pi})")
        if max_pi is not None:
            conditions.append(f"(PI <= {max_pi})")

        # Handle rise time filter
        if apply_rsl_rise_time:
            rise_time_cond = "(((((RISE_TIME + 0.00075*DERIV_MAX) > 46) && ((RISE_TIME + 0.00075*DERIV_MAX) < 58))&&ITYPE<4)||(ITYPE==4))"
            conditions.append(rise_time_cond)

        # Handle ITYPE conditions from grades
        if grades:
            # Translating the strings into ITYPE format
            grade_conditions = [f"(ITYPE == {GRADES[g]})" for g in grades]
            if len(grade_conditions) > 1:
                itype_cond =  " || ".join(grade_conditions)
            else:
                itype_cond = grade_conditions[0]
            
            itype_cond = "(" + itype_cond + ")"
            conditions.append(itype_cond)

        # Add STATUS condition if provided
        if status is not None:
            status_conditions = [f"(STATUS[{s}] == b0)" for s in status]
            if len(status_conditions) > 1:
                status_cond =  " && ".join(status_conditions)
            else:
                status_cond = status_conditions[0]
            
            conditions.append(status_cond)
        
        if remove_old_CI_rows:
            ci_cond = "(abs(rawy-0)>=3 && abs(rawy-80)>=3 && abs(rawy-160)>=3 && abs(rawy-240)>=3 && abs(rawy-320)>=3 && abs(rawy-400)>=3 && abs(rawy-480)>=3 && abs(rawy-560)>=3) && (abs(rawy-0)>=3 && abs(rawy-35)>=3 && abs(rawy-115)>=3 && abs(rawy-195)>=3 && abs(rawy-275)>=3 && abs(rawy-355)>=3 && abs(rawy-435)>=3 && abs(rawy-515)>=3 && abs(rawy-595)>=3) && (grade==0 || grade==2 || grade==3 || grade==4 || grade==6) && DETY>360 && DETY<1440"
            conditions.append(ci_cond)

        # Join everything with &&
        final_condition = " && ".join(conditions)
        final_condition = f"[{final_condition}]"

        return final_condition

    def _get_outfile_name(self):
        """
        Internal method to read the events dir and see how to name a new screened event. Basically
        every screened or filtered eventlist will have the suffix 'cl{i}' where i is an int that
        increases by one everytime a new filtering is made.
        """
        if os.path.exists(ANALYSIS_DIR + '/event'):
            events = [e for e in os.listdir(ANALYSIS_DIR + '/event') if e.endswith('.evt.gz')]
            ident = [int(ev.split('.evt')[0][-1]) for ev in events]
            highest_id = sorted(ident)[-1]
            # This is the number that will go on cl{new_id}.evt.gz
            new_id = str(highest_id + 1)
        
        else:
            os.makedirs(ANALYSIS_DIR + '/event')
            new_id = '1'

        # formatting name of filtered event
        short_inst = INST_NAMES[self.inst]
        obs_md = INST_MODES[short_inst][self.obs_mode]
        outfile = ANALYSIS_DIR + f'/event/{self.obs_id}{short_inst}_{obs_md}_cl{new_id}.evt.gz'
        
        return outfile

    def _gen_filtered_class(self, file, screening_applied, screening_info, gti_applied, gtifile):
        """
        Internal method to instantiate a new Eventlist object with the filtering information
        included.

        Parameters
        ----------
        file : str
            relative path of new eventlist to be instantiated
        screening_applied : bool
            True if screening has been applied
        screening_info : dict
            dictionary containing the configuration of the screening that has taken place
        gti_applied : bool
            True if gti filtering has been applied
        gtifile: str
            relative path to gtifile that has been applied to events
        
        Returns
        -------
        filtered_evt : EventList
            new eventlist object with filtering information included in attributes.
        """
        filtered_evt = EventList(file, self.obs_id, self.inst, self.obs_mode, self.root_dir)
        filtered_evt.screening_applied = screening_applied
        filtered_evt.screening_info = screening_info
        filtered_evt.gti_applied = gti_applied
        filtered_evt.gtifile = gtifile
    
        return filtered_evt


    def screen(self, apply_rsl_rise_time: bool = True, min_pi: int = None, max_pi: int = None, 
               status: list[int] = None, grades: list[str] = None, remove_old_CI_rows: bool = False):
        """
        Screens/ cleans event lists based on conditions given in the arguments.

        Parameters
        ----------
        apply_rsl_rise_time : bool
            There is a reccomended filtering on rise time for resolve event lists which is applied
            if this is set to True.
        min_pi : int
            minimum PI for events
        max_pi : int 
            maximum PI for events
        status : list[int]
            should be input as a list, even if there is just one status to be filtered. (see the 
            link in the docstring for proper info about what status number means what)
        grades : list[str]
            Accepts the strings 'Hp', 'Mp', 'Ms', 'Lp', 'Ls'.
        
        Returns
        -------
        rtype : EventList
            The screened eventlist loaded into an Eventlist object.
        """
        # Getting the string that interprets the arguments to this function and parses them into 
        # a logical statement that can be used with ftcopy
        screen_condition = self._gen_screen_condition(apply_rsl_rise_time=apply_rsl_rise_time, 
                                                      min_pi=min_pi, max_pi=max_pi, status=status,
                                                      grades=grades, 
                                                      remove_old_CI_rows=remove_old_CI_rows)
        # Sorting these so that when we check if these conditions have been parsed before 
        # we dont need to worry about the order
        if status is not None:
            status = sorted(status)
        
        # Sorting these so that when we check if these conditions have been parsed before 
        # we dont need to worry about the order
        if grades is not None:
            grades = sorted(grades, key=lambda g: GRADES[g])
        
        outfile = self._get_outfile_name()
    
        # appending this screening to the df
        screening_info = {"evtfile": outfile,
                   "base_evtfile": self.path,
                "rise_time_applied": int(apply_rsl_rise_time),
                "min_pi": min_pi,
                "max_pi": max_pi,
                "status": json.dumps(status) if status is not None else None,
                "grades": json.dumps(grades) if grades is not None else None,
                "remove_old_CI_rows": int(remove_old_CI_rows),
                "gtifile": self.gtifile}

        # if this already exists, it means that screening has happened before, so need to check
        # that the screening with current arguments hasnt been run before 
        table_path = ANALYSIS_DIR + "/event/screen_lookup.csv"
        if os.path.exists(table_path):
            df = pd.read_csv(table_path) 
            # annoying stuff happens if you try to compare with a Nan, so we have to do all this
            min_pi_condition = (
                df['min_pi'].isnull() if min_pi is None else
                df['min_pi'] == min_pi
            )

            max_pi_condition = (
                df['max_pi'].isnull() if max_pi is None else
                df['max_pi'] == max_pi
            )
            status_condition = (
                df["status"].isnull() if status is None else
                df["status"].apply(lambda x: tuple(json.loads(x)) if pd.notnull(x) else None) == tuple(status)
            )

            grades_condition = (
                df["grades"].isnull() if grades is None else
                df["grades"].apply(lambda x: tuple(json.loads(x)) if pd.notnull(x) else None) == tuple(grades)
            )

            # Check if row already exists
            existing = df[
            (df["rise_time_applied"] == int(apply_rsl_rise_time)) &
            min_pi_condition &
            max_pi_condition &
            status_condition &
            grades_condition &
            (df["remove_old_CI_rows"] == int(remove_old_CI_rows)) &
            (df["gtifile"] == self.gtifile) &
            (df['base_evtfile'] == self.path)
             ]

            # If row does exist, screening has been done before, so no need to perform ftcopy
            if not existing.empty:
                file  = existing["evtfile"].values[0]
                print(f"Screening has already been applied, made: {file}")

                return self._gen_filtered_class(file, screening_applied=True, screening_info=screening_info,
                                        gti_applied=self.gti_applied, gtifile=self.gtifile)
    
        else:
            df = pd.DataFrame(columns=["evtfile", "base_evtfile", "rise_time_applied", "min_pi", 
                                       "max_pi", "status", "grades", "remove_old_CI_rows", 
                                       "gtifile"])
        from utils import ftcopy
        # actually performing the screening
        ftcopy(self.path, screen_condition, outfile)

        df = df.append(screening_info, ignore_index=True)
        df.to_csv(table_path, index=False)

        return self._gen_filtered_class(outfile, screening_applied=True, screening_info=screening_info,
                                        gti_applied=self.gti_applied, gtifile=self.gtifile)

    def apply_gti(self, gtifile):
        """
        Applies GTI to events.

        Parameters
        ----------
        gtifile: str
            path to gtifile that will be used to filter events.
        
        Returns
        -------
        filtered_events : Eventlist
            new eventlist object of events with gti applied.
        """
        outfile = self._get_outfile_name()

        # Check if gtifile has been applied before
        # if this already exists, it means that screening has happened before, so need to check
        # that the screening with current arguments hasnt been run before 
        table_path = ANALYSIS_DIR + "/event/screen_lookup.csv"
        if os.path.exists(table_path):
            df = pd.read_csv(table_path)
        
            existing = df[
                (df['base_evtfile'] == self.path) &
                (df['gtifile'] == gtifile)
            ]

            if not existing.empty:
                file  = existing["evtfile"].values[0]
                print(f"GTI {gtifile} has already been applied to this event list, making {file}.")

                return self._gen_filtered_class(file, screening_applied=self.screening_applied,
                                                screening_info=self.screening_info, gti_applied=True,
                                                gtifile=gtifile)
    
        else:
            df = pd.DataFrame(columns=["evtfile", "base_evtfile", "rise_time_applied", "min_pi", 
                                       "max_pi", "status", "grades", "remove_old_CI_rows", "gtifile"])
        
        from utils import extractor
        # actually performing the screening
        extractor(self.path, outfile, gtifile)

        if self.screening_applied:
            # appending this screening to the df
            new_row = {"evtfile": outfile,
                    "base_evtfile": self.path,
                    "rise_time_applied": self.screening_info['rise_time_applied'],
                    "min_pi": self.screening_info['min_pi'],
                    "max_pi": self.screening_info['max_pi'],
                    "status": self.screening_info['status'],
                    "grades": self.screening_info['grades'],
                    "remove_old_CI_rows": self.screening_info['grades'],
                    "gtifile": gtifile}
        else:
            # appending this screening to the df
            new_row = {"evtfile": outfile,
                    "base_evtfile": self.path,
                    "rise_time_applied": 0,
                    "min_pi": None,
                    "max_pi": None,
                    "status": None,
                    "grades": None,
                    "remove_old_CI_rows": 0,
                    "gtifile": gtifile}

        df = df.append(new_row, ignore_index=True)
        df.to_csv(table_path, index=False)

        return self._gen_filtered_class(outfile, screening_applied=self.screening_applied,
                                                screening_info=self.screening_info, gti_applied=True,
                                                gtifile=gtifile)
            
class XRISMSource():
    """
    Class that can load in all eventlists given a root_dir where the data is given.
    """
    def __init__(self, obs_id, root_dir, obs_modes):
        self.obs_id = obs_id
        self.obs_modes = obs_modes
        self.root_dir = root_dir

        self.events = {}
        for long_inst, inst in INST_NAMES.items():
            self.events[long_inst] = {}
            for obs_mode in obs_modes[long_inst]:
                evt_list = f"xa{obs_id}{inst}_p0{INST_MODES[inst][obs_mode]}_cl.evt.gz"
                evt_path = os.path.join(root_dir, obs_id, long_inst, 'event_cl', evt_list)
                self.events[long_inst][obs_mode] = EventList(evt_path, self.obs_id, long_inst, 
                                                             obs_mode, root_dir)
