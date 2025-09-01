# Im just dumping all the general utility functions in here
from functools import wraps
import subprocess
from subprocess import Popen, PIPE
import os
import shutil
import pandas as pd
import yaml

from constants import *
from conf import ANALYSIS_DIR

def xselect_run(xselect_func):
    """
    Decorator for functions that need to use xselect. Xselect functions output a string for a script,
    the dir name for products to be stored in, any information about region filtering used, and the
    type of the product. This decorator will write the script file and run it, then organise the 
    products into an expected dir structure.
    """
    @wraps(xselect_func)
    def wrapper(*args, **kwargs):
        cmd, cmd_dir, region_info, product_type = xselect_func(*args, **kwargs)

        # i will change dirs later so need to know where to go back to
        og_cwd = os.getcwd()

        # if data product has already been made we dont need to run the rest of the func
        if os.path.exists(os.path.join(ANALYSIS_DIR, product_type, cmd_dir)):
            print(f'{product_type} with config {cmd_dir} already made.')
            return os.listdir(os.path.join(ANALYSIS_DIR, product_type, cmd_dir))

        os.makedirs(cmd_dir, exist_ok=True)
        os.chdir(cmd_dir)

        script_name = cmd_dir + '.sh'

        # Writing the script
        with open(script_name, 'w') as file:
            file.write(cmd)

        # Make the script executable
        os.chmod(script_name, 0o755)

        # running the script
        print(f"Starting make {product_type}.")
        result = subprocess.run(f"source ~/headas.sh; source {script_name}", shell=True, 
                                capture_output=True)
        write_res = result.stdout
        if result.stderr is not None:
            write_res += result.stderr
        # i wanted to save everything to a file
        with open(cmd_dir + '_output.txt', 'w') as file:
            file.write(str(write_res))
        print(f"{product_type} made.")
        print(f"Files made: {os.listdir(os.getcwd())}")

        # depending on what product we are making ie, image/expmap/spectra etc, we will move this
        # all into a dir named after the product type, just to organise the files in a clearer way
        product_type_dir = os.path.join(ANALYSIS_DIR, product_type)
        if not os.path.exists(product_type_dir):
            os.makedirs(product_type_dir)
        
        # returning to the og cwd and moving the dir containing the data products and generation
        # files into the product dir
        os.chdir(og_cwd)
        shutil.move(cmd_dir, product_type_dir)

        # if products were made using a region, this information is stored
        if region_info != '':
            if region_info['pix_reg']:
                # pixel regions behave differently to regular regions
                update_pix_regions_table(region_info)
            else:
                update_regions_table(region_info)
        
        return os.listdir(os.path.join(ANALYSIS_DIR, product_type, cmd_dir))
    return wrapper

def heasoftpy_run(heasoftpy_func):
    """
    Decorator for functions that need to use heasoftpy. Heasoftpy functions output a dictionary of
    parameters to be parsed into heasoftypy tasks, the relevant heasoft command to be run, the dir 
    name for products to be stored in, any information about region filtering used, and the
    type of the product. This decorator will run heasoftpy with the right task with the parameters 
    parsed to the decorator, then organise any output products into an expected dir structure.
    """
    @wraps(heasoftpy_func)
    def wrapper(*args, **kwargs):
        params, heasoft_cmd, cmd_dir, region_info, product_type = heasoftpy_func(*args, **kwargs)
        
        og_cwd = os.getcwd()
        # This accepts any defaults for arguments we didnt include in the params dict
        params['noprompt'] = True

        # if data product has already been made we dont need to run the rest of the func
        if os.path.exists(os.path.join(ANALYSIS_DIR, product_type, cmd_dir)):
            print(f'{product_type} with config {cmd_dir} already made.')
            return os.listdir(os.path.join(ANALYSIS_DIR, product_type, cmd_dir))

        os.makedirs(cmd_dir, exist_ok=True)
        os.chdir(cmd_dir)

        # since params are in a dict, easy to use the yml format to save this neatly
        with open(cmd_dir + '_input.yml', 'w') as file:
            yaml.dump(params, file, default_flow_style=False)
        
        # Actually running the heasoft command
        print(f"Starting {product_type}.")
        result = HEASOFT_CMDS[heasoft_cmd](params)
        write_res = result.stdout
        # saving everything to a file
        if result.stderr is not None:
            write_res += result.stderr
        with open(cmd_dir + '_output.txt', 'w') as file:
            file.write(write_res)
        print(f"{product_type} ended.")
        print(f"Files made: {os.listdir(os.getcwd())}")

        # depending on what product we are making ie, image/expmap/spectra etc, we will move this
        # all into a dir named after the product type, just to organise the files in a clearer way
        product_type_dir = os.path.join(ANALYSIS_DIR, product_type)
        if not os.path.exists(product_type_dir):
            os.makedirs(product_type_dir)

        # returning to the og cwd and moving the dir containing the data products and generation
        # files into the product dir
        os.chdir(og_cwd)
        shutil.move(cmd_dir, product_type_dir)

        # For ARFs I write pixel region files in the heasoftpy_func, which is in the og_cwd, so 
        # need to move this to the dir where the rest of the command was run
        if product_type == 'arf':
            # I give it the 'to_move' prefix so it is easy to find            
            file_to_move = [f for f in os.listdir(ANALYSIS_DIR) if 'to_move' in f]
            print(file_to_move)

            if len(file_to_move) != 0:
                file_to_move = file_to_move[0]
                shutil.move(os.path.join(ANALYSIS_DIR, file_to_move), os.path.join(ANALYSIS_DIR, product_type, cmd_dir))
                # removing the prefix
                os.rename(os.path.join(ANALYSIS_DIR, product_type, cmd_dir, file_to_move), 
                        os.path.join(ANALYSIS_DIR, product_type, cmd_dir, file_to_move.removeprefix('to_move_')))
        
        # if products were made using a region, this information is stored
        if region_info != '':
            if region_info['pix_reg']:
                # pixel regions behave differently to regular regions
                update_pix_regions_table(region_info)
            else:
                update_regions_table(region_info)

        return os.listdir(os.path.join(ANALYSIS_DIR, product_type, cmd_dir))
    return wrapper

def cmd_run(func):
    """
    Decorator for functions that need to use heasoft via the command line. Functions output the cmd
    to be run, the dir name for products to be stored in, any information about region filtering 
    used, and the type of the product. This decorator will run the cmd, then organise any output 
    products into an expected dir structure.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        cmd, cmd_dir, region_info, product_type = func(*args, **kwargs)
        
        og_cwd = os.getcwd()
        # This accepts any defaults for arguments we didnt include in the params dict

        # if data product has already been made we dont need to run the rest of the func
        if os.path.exists(os.path.join(ANALYSIS_DIR, product_type, cmd_dir)):
            print(f'{product_type} with config {cmd_dir} already made.')
            return os.listdir(os.path.join(ANALYSIS_DIR, product_type, cmd_dir))

        os.makedirs(cmd_dir, exist_ok=True)
        os.chdir(cmd_dir)

        # since params are in a dict, easy to use the yml format to save this neatly
        with open(cmd_dir + '_input.txt', 'w') as file:
            file.write(cmd)
        
        # Actually running the heasoft command
        print(f"Starting {product_type}.")
        out, err = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()
        out = out.decode("UTF-8", errors='ignore')
        err = err.decode("UTF-8", errors='ignore')

        # saving everything to a file
        write_res = out + err if err is not None else out
        with open(cmd_dir + '_output.txt', 'w') as file:
            file.write(write_res)
        print(f"{product_type} ended.")
        print(f"Files made: {os.listdir(os.getcwd())}")

        # depending on what product we are making ie, image/expmap/spectra etc, we will move this
        # all into a dir named after the product type, just to organise the files in a clearer way
        product_type_dir = os.path.join(ANALYSIS_DIR, product_type)
        if not os.path.exists(product_type_dir):
            os.makedirs(product_type_dir)

        # returning to the og cwd and moving the dir containing the data products and generation
        # files into the product dir
        os.chdir(og_cwd)
        shutil.move(cmd_dir, product_type_dir)

        # For ARFs I write pixel region files in the heasoftpy_func, which is in the og_cwd, so 
        # need to move this to the dir where the rest of the command was run
        if product_type == 'arf':
            # I give it the 'to_move' prefix so it is easy to find            
            file_to_move = [f for f in os.listdir(ANALYSIS_DIR) if 'to_move' in f]
            print(file_to_move)

            if len(file_to_move) != 0:
                file_to_move = file_to_move[0]
                shutil.move(os.path.join(ANALYSIS_DIR, file_to_move), os.path.join(ANALYSIS_DIR, product_type, cmd_dir))
                # removing the prefix
                os.rename(os.path.join(ANALYSIS_DIR, product_type, cmd_dir, file_to_move), 
                        os.path.join(ANALYSIS_DIR, product_type, cmd_dir, file_to_move.removeprefix('to_move_')))
        
        # if products were made using a region, this information is stored
        if region_info != '':
            if region_info['pix_reg']:
                # pixel regions behave differently to regular regions
                update_pix_regions_table(region_info)
            else:
                update_regions_table(region_info)

        return os.listdir(os.path.join(ANALYSIS_DIR, product_type, cmd_dir))
    return wrapper

def small_heasoftpy(func):
    """
    Ok so this is just a simplified version of heasoftpy_run but for uses where I dont need to do 
    any directory restructuring or sorting of region files. So that's why it's called 
    small_heasoftpy - it is just change dir, run heasoftpy, change back to og dir.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        og_cwd = os.getcwd()
        os.chdir(ANALYSIS_DIR)

        heasoft_cmd, params = func(*args, **kwargs)

        result = HEASOFT_CMDS[heasoft_cmd](params)

        print(result.stdout)
        if result.stderr is not None:
            os.chdir(og_cwd)
            raise KeyError(result.stderr)
        
        os.chdir(og_cwd)
    return wrapper

def update_pix_regions_table(region_info: dict):
    """
    Updates or creates a table storing pixel region information. Each pixel gets its own column 
    with True if included. his table can then be used to check what files have had which region 
    filtering applied. (products will have the naming scheme prefix_pixreg-{pixel_name}_suffix, so 
    you can look in this table for the pixel_name, and see which pixels it corresponds to.)

    Parameters
    ----------
    region_info : dict
      pix_name : A nickname for the region, so you can look at a file name and understand what
        the region is. ie. inner_pix, if it is a region that only selects central pixels, then files 
        and naming structures will contain 'inner_pix' in. These should be unique for each pixel
        selection.
      include_pix: Set or list of pixels to include.
      exclude pix: list of pixels to exclude.
    """
    og_cwd = os.getcwd()
    os.chdir(ANALYSIS_DIR)

    table_path = "pixel_region_lookup.csv"
    include_pix = region_info['include_pix']
    exclude_pix = region_info['exclude_pix']
    pix_name = region_info['pix_name']

    # Convert to sets for fast lookup
    include_pix = set(include_pix) if include_pix is not None else set(range(36))
    exclude_pix = set(exclude_pix) if exclude_pix is not None else set()

    # Create row data with pixel columns named pixel_0 to pixel_35
    row_data = {f'pixel_{i}': (i in include_pix and i not in exclude_pix) for i in range(36)}
    row_data['pix_name'] = pix_name

    # Convert to DataFrame
    new_row_df = pd.DataFrame([row_data])

    # If file exists, check if identical row exists
    if os.path.exists(table_path):
        existing_df = pd.read_csv(table_path)

        # Reindex to ensure all columns are present
        all_columns = [f'pixel_{i}' for i in range(36)] + ['pix_name']
        existing_df = existing_df.reindex(columns=all_columns, fill_value=False)
        new_row_df = new_row_df.reindex(columns=all_columns, fill_value=False)

        # Check if an identical row already exists
        duplicate = (existing_df == new_row_df.iloc[0]).all(axis=1).any()
        if duplicate:
            print("Identical pixel region row already exists. No update made.")
            os.chdir(og_cwd)
            return

        # Append new data
        updated_df = pd.concat([existing_df, new_row_df], ignore_index=True)
    else:
        updated_df = new_row_df

    # Write to file
    updated_df.to_csv(table_path, index=False)

    os.chdir(og_cwd)


def update_regions_table(region_info: dict):
    """
    Updates or creates a table storing region information. This table can then be used to check 
    what files have had which region filtering applied. (products will have the naming scheme
    prefix_reg-{regionname}_suffix, so you can look in this table for the regionname, and see
    which regionfile it corresponds to.)

    Parameters
    ----------
    region_info : dict
        contains:
        name (str): A nickname for the region, so you can look at a file name and understand what
        the region is. ie. no_cal_srcs, if it is a region that removes the cal sources from xtend.
        then files and naming structures will contain 'no_cal_srcs' in. These should be unique for
        each region used.
        file (str): Path to the region file.
        content (str): Contents of the region file.
    """
    og_cwd = os.getcwd()
    os.chdir(ANALYSIS_DIR)

    region_name = region_info['name']
    region_path = region_info['file']
    region_content = region_info['content']

    table_path = 'region_lookup.csv'
    # Load existing table if it exists
    if os.path.exists(table_path):
        df = pd.read_csv(table_path)
    else:
        df = pd.DataFrame(columns=["Region Name", "Region Path", "Region Content"])
    
    # Check if row already exists
    existing = df[
        (df["Region Name"] == region_name) &
        (df["Region Path"] == region_path) &
        (df["Region Content"] == region_content)
    ]

    if not existing.empty:
        print("Matching row already exists. No update made.")
        os.chdir(og_cwd)
        return
    
    # Append new row
    new_row = {"Region Name": region_name,
            "Region Path": region_path,
            "Region Content": region_content}
    df = df.append(new_row, ignore_index=True)

    # Save updated table
    df.to_csv(table_path, index=False)

    os.chdir(og_cwd)
