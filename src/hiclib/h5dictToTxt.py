"""
Saves mirnylab.h5dict.h5dict file to a folder with different datsets. 

Usage
-----

python h5dictToTxt.py in_hdf5_file out_folder

out_folder may not exist, but may not be file. 
Existing files will be overwritten. 

Each key of the h5dict dataset is converted to a text file. 
If key is an array, a matlab-compatible array is returned. 
If key is not an array, python's repr command is used 
to get an exact representation of an object.  
"""

from mirnylib.h5dict import h5dict
from mirnylib.numutils import generalizedDtype
import os,sys
import numpy

if len(sys.argv) != 3:
    print "Converts h5dict file to a bunch of txt files"
    print 'Each key of the array is converted to a separate file'
    print "Numpy.arrays are converted to a matlab-loadable txt files"
    print "Other keys are converted using python's repr command"
    print "Usage: python h5dictToTxt h5dictFile folderName"
    print "Folder will be created if not exists"
    exit() 
    
filename = sys.argv[1]
folder = sys.argv[2]

if not os.path.exists(filename):
    raise IOError("File not found: %s" % filename)
if os.path.isfile(folder):
    raise IOError("Supplied folder is a file! ")
if not os.path.exists(folder):
    os.mkdir(folder)
     
mydict = h5dict(filename,'r')
for i in mydict.keys():
    data = mydict[i] 
    savefile = os.path.join(folder,i)
    if issubclass(type(data),numpy.ndarray):
        print "saving numpy array", i, "to", savefile  
        if len(data.shape) > 0:
            numpy.savetxt(savefile,data)
            continue
        
    if type(data) == str:
        datarepr = data
    else:
        datarepr = repr(data) 
    print "saving data", i, "to", savefile
    with open(savefile,'w') as f:
        f.write(datarepr)
        
    
            
        
              