import systemutils
systemutils.setExceptionHook() 
from copy import  deepcopy
import matplotlib
import matplotlib.pyplot as plt 

import numpy 



def cmap_map(function,cmap,mapRange = [0,1]):
    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
    
    Also trims the "range[0]:range[1]" fragment from the colormap - use this to cut the part of the "jet" colormap! 
    """
    cdict = cmap._segmentdata
    
    for key in cdict.keys():
        print cdict[key]
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
    
    step_list = sum(step_dict.values(), [])
    array = numpy.array
    step_list = array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : array(cmap(step)[0:3])
    old_LUT = array(map( reduced_cmap, mapRange[0] + step_list * (mapRange[1] - mapRange[0]) ))
    new_LUT = array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def showPolymerRasmol(x,y=None,z=None):
    """
    Shows the polymer using rasmol.
    Can't properly treat continuous chains (they will have linkers of 5 balls between chains) 
    Accepts data as x,y,z or as an array of Nx3 or 3xN shape"
    Shows system  by drawing spheres
    draws 4 spheres in between any two points (5 * N spheres total)
    
    Parameters
    ----------
    x : array
        Nx3 or 3xN array, or N-long array
    y,z : array or None
        if array, it corresponds to y coordinate. If none, x is assumed to have 3 coordinates
    
    """
    
    import os, tempfile      
    #if you want to change positions of the spheres along each segment, change these numbers
    #e.g. [0,.1, .2 ...  .9] will draw 10 spheres, and this will look better
    shifts = [0.,0.2,0.4,0.6,0.8]
    
    if y == None: data = numpy.array(x)
    else: data = numpy.array([x,y,z])
    if len(data[0]) != 3: 
        data = numpy.transpose(data)
    if len(data[0]) != 3:
        print "wrong data!"
        return
    
    #determining the 95 percentile distance between particles,  
    meandist = numpy.percentile(numpy.sqrt(numpy.sum(numpy.diff(data,axis = 0)**2,axis = 1)),95)
    #rescaling the data, so that bonds are of the order of 1. This is because rasmol spheres are of the fixed diameter. 
    data /= meandist
    
    #writing the rasmol script. Spacefill controls radius of the sphere. 
    rascript = tempfile.NamedTemporaryFile()
    rascript.write("""wireframe off 
    color temperature
    spacefill 100 
    background white
    """)
    rascript.flush()
    
    
    #creating the array, linearly chanhing from -225 to 225, to serve as an array of colors 
    #(rasmol color space is -250 to 250, but it  still sets blue to the minimum color it found and red to the maximum). 
    colors = numpy.array([int((j*450.)/(len(data)))-225 for j in xrange(len(data))])    
    
    #creating spheres along the trajectory
    #for speedup I just create a Nx4 array, where first three columns are coordinates, and fourth is the color      
    newData = numpy.zeros((len(data) * len(shifts) - (len(shifts) - 1) ,4))  
    for i in xrange(len(shifts)):            
        #filling in the array like 0,5,10,15; then 1,6,11,16; then 2,7,12,17, etc. 
        #this is just very fast
        newData[i:-1:len(shifts),:3] = data[:-1] * shifts[i] + data[1:] * ( 1 - shifts[i])            
        newData[i:-1:len(shifts),3] = colors[:-1]
    newData[-1,:3] = data[-1]
    newData[-1,3] = colors[-1]
                
    towrite = tempfile.NamedTemporaryFile()
    towrite.write("%d\n\n"%(len(newData)))  #number of atoms and a blank line after is a requirement of rasmol
        
    for i in newData:                     
        towrite.write("CA\t%lf\t%lf\t%lf\t%d\n" % tuple(i)) 
    towrite.flush()
    #For windows you might need to change the place where your rasmol file is  
    if os.name == "posix":  #if linux 
        os.system("rasmol -xyz %s -script %s" % (towrite.name, rascript.name))
    else:     #if windows 
        os.system("C:/RasWin/raswin.exe -xyz %s -script %s" % (towrite.name, rascript.name))
        



def scatter3D(x,y,z,color):
    """shows a scatterplot in 3D"""
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if (type(color) == numpy.ndarray) or (type(color) == list):
        color = numpy.array(color,dtype = float)
        color -= color.min()
        color /= float(color.max() - color.min() )
        if len(set(color)) > 20:                        
            for i in xrange(len(x)):
                ax.scatter(x[i], y[i], z[i], c=plt.cm.get_cmap("jet")(color[i]))
        else:
            colors = set(color)
            for mycolor in colors:
                mask = (color == mycolor)
                ax.scatter(x[mask],y[mask],z[mask],c=plt.cm.get_cmap("jet")(mycolor))
    else: ax.scatter(x,y,z,c=color)
    plt.show()
                
                
def removeAxes(mode = "normal",shift = 0, ax = None):
    
    if ax == None: 
        ax = plt.gca()
    for loc, spine in ax.spines.iteritems():
        if mode == "normal":
            if loc in ['left','bottom']:
                if shift != 0: spine.set_position(('outward',shift)) # outward by 10 points
            elif loc in ['right','top']:
                spine.set_color('none') # don't draw spine
            else:
                raise ValueError('unknown spine location: %s'%loc)
        else:
            if loc in ['left','bottom','right','top']:
                spine.set_color('none') # don't draw spine
            else:
                raise ValueError('unknown spine location: %s'%loc)
            
def removeBorder(ax = None):
    removeAxes("all",0,ax = ax)
    if ax == None: ax = plt.gca()    
    for _, line in enumerate(ax.get_xticklines() + ax.get_yticklines()): 
        line.set_visible(False)     
    if ax == None: ax = plt.axes()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    

    

def niceShow(mytype = None):
    if mytype == "log":
        plt.xscale("log")
        plt.yscale("log")

    legend = plt.legend(loc=0,prop={"size":15})
    if legend != None:
        legend.draw_frame(False)
    removeAxes(shift = 0)
    plt.gcf().subplots_adjust(left=0.07, bottom=0.12, top=0.98, right=0.98)
    plt.show()

        
def mat_img(a,cmap="hot_r",trunk = False, **kwargs):
    a = numpy.array(a,float)
    if trunk != False:
        if trunk == True:
            trunk = 0.01
        sa = numpy.sort(a.ravel())
        a[a>sa[(1 - trunk) * len(sa)]] = sa[(1 - trunk) * len(sa)]
        a[a<sa[trunk * len(sa)]] = sa[trunk * len(sa)]
    #plt.ioff() 
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #cax = ax.imshow(a, interpolation = 'nearest', cmap = cmap)
    #cbar = fig.colorbar(cax)
    
    def do_all():
        plt.imshow(a, interpolation = 'nearest', cmap = cmap,**kwargs)
        plt.colorbar()
        plt.show()
    do_all()
    
    #plt.close() 



def plot_occupation(*prop):
    print prop
    try:
        prop[0][0][0]
        labels = prop[1]
        prop = prop[0]
    except:
        labels = [str(i) for i in range(len(prop))]
        pass
    
    global savedplot
    prop = list(prop)
        
    savedplot = deepcopy(prop)    
    #steps = range(1,len(prop[0])+1)

    for j,i in enumerate(prop):
         
       
        plt.plot(range(len(i)),i,label=labels[j])
    plt.legend(loc=0)
    
    plt.show()
    

 
 
    
def pointplot(data,lines=1,labels=None,size=(5,5),fs=15,linewidth=None,markersize=None):
    plt.figure(figsize=size)
    fg = plt.subplot(1,1,1)
    
    if lines >2:        
        plt.xscale("log")
        plt.yscale("log")
    
    #if (numpy.min(numpy.array(data))<0) and loglog == True: 
    #    loglog = False
    #    print "loglog is set but data is below zero"

    if labels == None:
        labels = [str(i) for i in range(len(data))]
    global savedpointplot
    savedpointplot = deepcopy(data)
    if linewidth == None: linewidth = [1 for i in xrange(len(data))]
    if markersize == None: markersize = [3 for i in xrange(len(data))]
    for j,i in enumerate(data):
        
                
        
        if lines%2==0:
            shape = 'o'
        else: shape = '-'
        if len(i) > 2:
            color = i[2]
        else: color = ''
        sstring = "%s%s"%(shape,color)
        fg.plot(i[0],i[1],sstring,label=labels[j],linewidth=linewidth[j],markersize = markersize[j])
        
        
        
    #if loglog == True: 
    #    plt.loglog(*pl)
    #else: plt.plot(*pl)
    plt.xlabel("Rg(N^{2/3})",fontsize=int(1.2*fs))
    plt.ylabel("Log( knot pol)",fontsize=int(1.2*fs))
    ax = plt.axes()
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
            
    
    
    #plt.grid()
    plt.legend(loc=0,prop={"size":int(1.1*fs)})
    plt.show()

