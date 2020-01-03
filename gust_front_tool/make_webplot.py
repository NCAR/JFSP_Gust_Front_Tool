#!/usr/bin/env python

import sys, time, os
from webplot import webPlot, readGrid

def log(msg): print time.ctime(time.time()),':', msg

log('Begin Script')
stime = time.time()

newPlot = webPlot()
        
log('Reading Data')
newPlot.readEnsemble()

file_not_created, num_attempts = True, 0
while file_not_created and num_attempts <= 3:
        newPlot.createFilename()
        fname = newPlot.outfile
 
        log('Creating Map')
        newPlot.loadMap()

        log('Plotting Data')
        newPlot.plotFields()
        newPlot.plotTitleTimes()

       #newPlot.plotReports() # CSS added

        log('Writing Image')
        newPlot.saveFigure()

        if os.path.exists(fname):
            file_not_created = False 
            log('Created %s, %.1f KB'%(fname,os.stat(fname).st_size/1000.0))
    
        num_attempts += 1

etime = time.time()
log('End Plotting (took %.2f sec)'%(etime-stime))

