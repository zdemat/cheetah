#!/usr/bin/env python
import configobj
import os,sys
import time,datetime
import socket
import numpy as np
import run,terminal,google
import getpass

# Are we on the right host machine?
hostname = socket.gethostname()
if not ("psexport" in hostname or "pslogin" in hostname):
    sys.exit("ERROR: Cannot submit jobs from this host. You are logged in to %s but you need to be logged in to psexport of pslogin to submit jobs to the queue." % hostname)

# Read configuration
if len(sys.argv) < 2:
    sys.exit("Usage: busb_cheetah.py bsub_cheetah.cfg")

debug = False
#print sys.argv
#sys.exit()
if len(sys.argv) >= 3:
    if sys.argv[2] == "DEBUG":
        debug = True

def main():
    configfilename = sys.argv[1]
    C = configobj.ConfigObj(configfilename)
    Cl = C["locations"]
    E = {}
    for k,l in Cl.items():
        E[k] = os.path.expandvars(l)
    Cs = C["spreadsheet"]
    email = Cs["email"]
    if "password" not in Cs.keys():
        print "Gmail password [%s]:" % email
        password = getpass.getpass()
    else:
        password = Cs["password"]
    spreadsheet_name = Cs["spreadsheet_name"]
    worksheet_name = Cs["worksheet_name"]
    dry_run = Cs.as_bool("dry_run")

    # Init google table client
    gtab = google.GoogleTable(email,password,spreadsheet_name,worksheet_name)
    ttab = terminal.TerminalTable()
    runs = {}
    
    counter = 0
    
    while True:
        gtab.update_table_dict()
        valid_runs = set(gtab.get_valid_runs())
        runs_to_update = set(gtab.get_runs_to_update())
        runs_to_add = valid_runs - set(runs.keys())

        ns = list(runs_to_add)
        ns.sort()
        for n in ns:
            ttab.note("%s - adding" % n)
            runs[n] = run.Run(n,gtab,E)

        ns = list(runs_to_update)
        ns.sort()
        for n in ns:
            r = runs[n]
            ttab.note("%s - updating" % r.run_name)
            r.update()
            if r.status == "new":
                ttab.note("%s - initialization" % r.run_name)
                r.init_process()
                if r.prepared:
                    r.start()
                    ttab.note("%s - started" % r.run_name)
                    if Cs["mode"] == "swmr":
                        r.start_swmr()                                                     
                        ttab.note("%s - started (swmr)" % r.run_name)
                else:
                    ttab.note("%s - processing postponed." % r.run_name)
    
        ttab.set_runs(runs)
    
        if (counter % 10) == 0:
            ttab.note("Writing to google spreadsheet...")
            gtab.write(runs)

    #mill = ["-","\\","|","/"]
    #sys.stdout.write("\r"+mill[i%len(mill)] + " Waiting for new XTC files to appear...")
    
        time.sleep(0.1)
    #except:
    #    print "ERROR OCCURED - RESTARTING IN A FEW SECONDS..."
    #    time.sleep(5)

if debug:
    main()
else:
    while True:
        try:
            main()
        except:
            print "WARNING: run() crashed, restarting..."