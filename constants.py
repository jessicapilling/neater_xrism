import sys
sys.path.insert(0, "/Users/jp735/Desktop/software/heasoft-6.35/x86_64-apple-darwin22.6.0/lib/python")
import heasoftpy as hsp

# Translating the long names to the short ones
INST_NAMES = {'resolve': 'rsl', 'xtend': 'xtd'}

# Translating the obsmodes into something more readable
INST_MODES = {'rsl': {'undefined': 'px000',
                     'open': 'px1000',
                     'al/polymide': 'px2000',
                     'neutral_density': 'px3000',
                     'be': 'px4000',
                     'fe55': 'px5000'},
             'xtd': {'ccd1234_full_wdw': '30000010',
                     'ccd12_1/8_wdw': '31100010',
                     'ccd12_full_wdw_brst_md': '31200010',
                     'ccd12_1/8_wdw_brst_md': '31300010',
                     'ccd34_full_wdw': '32000010'}}

# Grades to ITYPE translation
GRADES = {'Hp': 0,
          'Mp': 1,
          'Ms': 2,
          'Lp': 3,
          'Ls': 4}

# For use by wrapper @heasoftpy_run
HEASOFT_CMDS = {'xaexpmap': hsp.xaexpmap,
                'rslmkrmf': hsp.rslmkrmf,
                'xtdrmf': hsp.xtdrmf,
                'xaarfgen': hsp.xaarfgen,
                'ftcopy': hsp.ftcopy,
                'maketime': hsp.maketime,
                'extractor': hsp.extractor,
                'rslnxbgen': hsp.rslnxbgen,
                'rslbratios': hsp.rslbratios}

# These are regions in DET coords for each resolve pixel
PIXEL_REGS = {0: 'box(4,3,1,1,0)',  # pixel 0
              1: 'box(6,3,1,1,0)', # pixel 1
              2: 'box(5,3,1,1,0)', # pixel 2
              3: 'box(6,2,1,1,0)', # pixel 3
              4: 'box(5,2,1,1,0)', # pixel 4
              5: 'box(6,1,1,1,0)', # pixel 5
              6: 'box(5,1,1,1,0)', # pixel 6
              7: 'box(4,2,1,1,0)', # pixel 7
              8: 'box(4,1,1,1,0)', # pixel 8
              9: 'box(1,3,1,1,0)', # pixel 9
              10: 'box(2,3,1,1,0)', # pixel 10
              11: 'box(1,2,1,1,0)', # pixel 11
              13: 'box(2,2,1,1,0)', # pixel 13
              14: 'box(2,1,1,1,0)', # pixel 14
              15: 'box(3,2,1,1,0)', # pixel 15
              16: 'box(3,1,1,1,0)', # pixel 16
              17: 'box(3,3,1,1,0)', # pixel 17
              18: 'box(3,4,1,1,0)', # pixel 18
              19: 'box(1,4,1,1,0)', # pixel 19
              20: 'box(2,4,1,1,0)', # pixel 20
              21: 'box(1,5,1,1,0)', # pixel 21
              22: 'box(2,5,1,1,0)', # pixel 22
              23: 'box(1,6,1,1,0)', # pixel 23
              24: 'box(2,6,1,1,0)', # pixel 24
              25: 'box(3,5,1,1,0)', # pixel 25
              26: 'box(3,6,1,1,0)', # pixel 26
              27: 'box(6,4,1,1,0)', # pixel 27
              28: 'box(5,4,1,1,0)', # pixel 28
              29: 'box(6,5,1,1,0)', # pixel 29
              30: 'box(6,6,1,1,0)', # pixel 30
              31: 'box(5,5,1,1,0)', # pixel 31
              32: 'box(5,6,1,1,0)', # pixel 32
              33: 'box(4,5,1,1,0)', # pixel 33
              34: 'box(4,6,1,1,0)', # pixel 34
              35: 'box(4,4,1,1,0)' } # pixel 35