import numpy as np # version >= 1.19.5
import os
import sys
import time
from astropy.table import vstack, Table # version >= 4.2



sys.path.append(__file__[:-11]+'src')
import setup as SH

argv = sys.argv[1:]
keys = {}
i = 0
order = 1
cpy = 0
blacklist = 1

download_data = False
if len(argv)>0:
    if argv[0] == "download_data": 
        download_data = True
        argv = argv[1:]
# read command imput
while i < len(argv):

    if argv[i][0] == "-":
        order = 0
    if argv[i] == "--raw_cands_table" or argv[i] == "-r" or order*(i+1) == 1 :
        SH.raw_cands_table = argv[i+1-order]
    elif argv[i] == "--good_hpms_file" or argv[i] == "-H" or order*(i+1) == 2:
        SH.good_hpms_file = argv[i+1-order]
        if SH.good_hpms_file == "py":
            SH.good_hpms_file = "good_HPMS.py"
            if argv[i+2-order][0] != '-':
                SH.hpms_eDR3_file = argv[i+1-order]
                i+=1
                cpy +1
    elif argv[i] == "--good_bgs_file" or argv[i] == "-B" or order*(i+1) == 4:
        SH.good_bgs_file = argv[i+1-order]
        if SH.good_bgs_file == "py":
            SH.good_bgs_file = "good_BGS.py"
            if argv[i+2][0] != '-':
                SH.bgs_eDR3_file = argv[i+2-order]
                i+=1
                cpy +1
            if argv[i+2][0] != '-':
                SH.DR2_good_bgs_file = argv[i+2-order]
                i+=1
                cpy +1
    elif argv[i] == "--n_core" or argv[i] == "-n" or order*(i+1) == 7:
        SH.n_core = int(argv[i+1-order])

    elif argv[i] == "--make_plots" or argv[i] == "-m" or order*(i+1) == 8:
        if order: SH.make_plots = bool(argv[i]) 
        else:
            try:
                if argv[i+1][0] == "-":
                    SH.make_plots = True
                else: SH.make_plots = bool(argv[i+1])
            except:
                SH.make_plots = True
    elif argv[i] == "--save_table" or argv[i] == "-s" or order*(i+1) == 9:
        if order: SH.save_table_process = bool(argv[i]) 
        else:
            try:
                if argv[i+1][0] == "-":
                    SH.save_table_process = True
                else: SH.save_table_process = bool(argv[i+1])
            except:
                SH.save_table_process = True
    elif argv[i] == "--filter" or argv[i] == "-f" or order*(i+1) == 9:
        if order: SH.do_filter = int(argv[i])
        else:
            try:
                if argv[i+1][0] == "-":
                    SH.do_filter = 1
                else: SH.do_filter = int(argv[i+1])
            except:
                SH.do_filter = 1
    elif argv[i] == "--prefix" or argv[i] == "-p" or order*(i+1) == 10:
        SH.prefix = argv[i+1-order]
        if SH.prefix[0] == "'" or SH.prefix[0] == '"':
            SH.prefix = SH.prefix[1:-1]
        SH.prefix = "."+SH.prefix 
    elif argv[i] == "--blacklist" or argv[i] == "-b" or order*(i+1) == 11:
        SH.Blacklist_file = argv[i+1-order]
    elif argv[i][:1] == "-n":pass
    elif argv[i][:1] == "-n":pass

    elif argv[i][0] == "-":
        print("unknown keyword: " + argv[i])
    i+=1

if download_data:
    import query_data
    query_data.main(query_data.input())
    sys.exit()
import main
main.broadcast()

try:
  result = main.main(n = keys.get('n_core', SH.n_core))
except FileNotFoundError as ve:
  print(ve)
  print('If you intend to use the default files, you might want to run "python %s download_data" first.'%sys.argv[0])