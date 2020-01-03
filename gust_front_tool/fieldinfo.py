import os

def readcm(name):
    '''Read colormap from file formatted as 0-1 RGB CSV'''
    rgb = []
    fh = open(name, 'r')
    for line in fh.read().splitlines(): rgb.append(map(float,line.split()))
    return rgb

def readNCLcm(name):
    '''Read in NCL colormap for use in matplotlib'''
    rgb, appending = [], False
#   fh = open('/glade/apps/opt/ncl/6.2.0/intel/12.1.5/lib/ncarg/colormaps/%s.rgb'%name, 'r')
#   fh = open(os.getenv('NCARG_ROOT','/glade/apps/opt/ncl/6.2.0/intel/12.1.5')+'/lib/ncarg/colormaps/%s.rgb'%name, 'r') # CSS made variable, commented out previous line
    fh = open('./colormaps/%s.rgb'%name, 'r') # CSS made variable, commented out previous line
    for line in fh.read().splitlines():
        if appending: rgb.append(map(float,line.split()))
        if ''.join(line.split()) in ['#rgb',';RGB']: appending = True
    maxrgb = max([ x for y in rgb for x in y ])
    if maxrgb > 1: rgb = [ [ x/255.0 for x in a ] for a in rgb ]
    return rgb

fieldinfo = {
  # gust front stuff
  'binary-gf'    :{ 'levels' : [0.5], 'cmap': [readNCLcm('precip2_17lev')[i] for i in (0,1,2,4,5,6,7,8,10,12,13,14,15)], 'fname':['binary_gust_front'] },
  'cref'         :{ 'levels' : [5,10,15,20,25,30,35,40,45,50,55,60,65,70], 'cmap': readcm('cmap_rad.rgb')[1:14], 'fname': ['REFL_MAX'] },
  'wind10m'      :{ 'fname'  : ['U10', 'V10'], 'skip': 10 },
  'speed10m'     :{ 'levels' : [3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51], 'cmap': readNCLcm('wind_17lev')[1:],'fname'  : ['U10', 'V10']},
  'precip'       :{ 'levels' : [0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3.0], 'cmap': [readNCLcm('precip2_17lev')[i] for i in (0,1,2,4,5,6,7,8,10,12,13,14,15)], 'fname':['tot_precip'] },
}
