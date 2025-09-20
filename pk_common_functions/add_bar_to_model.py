# -*- coding: future_fstrings -*-

from omegaconf import OmegaConf
import sys
import os
import numpy as np
import copy
from dataclasses import dataclass, field
from typing import List,Optional
from omegaconf import MISSING
from astropy.io import fits
from pk_common_functions.functions import get_model_DHI,copy_disk,\
    isiterable,tirific_template,equal_length
from scipy.interpolate import interp1d
#from multiprocessing import cpu_count

from typing import List, Optional
#import support_functions as sf
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Ellipse

@dataclass
class defaults:
    bar_angle: List = field(default_factory=lambda: [[37.]]) #Bar angle in degrees
    bar_length: Optional[float] = None #Bar length in arcsec, default is the estimated ILR
    bar_brightness: List = field(default_factory=lambda: [[None]]) 
    bar_velocities: List = field(default_factory=lambda: [[None]]) 
    bar_thickness: float = 0.1 #Thickness of bar in fraction of the length  
    input_disk: List = field(default_factory=lambda: [1]) #Disk to copy (bar_type ='Free') or add the bar to (bartype ='Harmonics)
    disk_length: Optional[float] = None #length of the optical disk in arcsec, defaults to R_HI
    input_def: Optional[str] = None #The original def file
    output_def: Optional[str] = None #The def file that contains the new model
    bar_type: str = 'Free'  # 'Free'  adds the bar as 2 disk s to the model have complete freedom, great for fitting a known bar.
                      # 'Harmonics' add the bar to the model as 2nd order harmonics as described in http://gigjozsa.github.io/tirific/model_geometry.html#Velocity%20structure%20I
  
def calculate_length(Template,Radii,disk_length,ext):
    velocity = np.array([x for x in Template[f'VROT{ext}'].split()]\
            ,dtype=float)
       
    max_rad = Radii[-1]
    #The pattern speed at a given radius is vrot/radii
    V_Rot = interp1d(Radii, velocity, fill_value="extrapolate")
    # The radius of co-ration can be approximated by the extend of the visible disk ~ Warpstart (Roberts et al. 1975)
    Omega_CR = V_Rot(disk_length)/disk_length
    # From this we can estimate the inner and outer Lindblad resonances (Eq 38 Dobbs & Baba 2012)
    #The epicyclic frequency k^2=R*d/drOmega^2+4*Omega^2
    # f'(x) = (f(x+h)-f(x-h))/2h
    h = disk_length/1000.
    derive = (V_Rot(float(disk_length+h))**2/(disk_length+h)**2
            - V_Rot(float(disk_length-h))**2/(disk_length-h)**2)/(2*h)
    k_CR = (disk_length * derive+4*Omega_CR**2)**0.5
    # So the ILR =
    LLR = Omega_CR-k_CR/2.
    ULR = Omega_CR+k_CR/2.
    Radii[0] = 0.1
    om = interp1d(Radii, velocity/Radii, fill_value="extrapolate")
    Radii[0] = 0.
    r_cur = Radii[1]
    while om(r_cur) > ULR and r_cur < max_rad:
        r_cur += 0.1
    length = 0.75*r_cur #This is the inner lindblad resonance estimate
    return length

def create_bar(Template, disks=[1], length_in=None, thickness = 0.1, \
        disk_length =None, bar_angle= [[37.]],\
        bar_velocities = [[None]], bar_brightness = [[None]] ):
    #disk_length= the size of the visible disk (WarpStart in AGC) in arcsec, 
    # set to last radius if unset

    '''Routine to create a Bar into a model'''
    Template_New =copy.deepcopy(Template)
    Radii = np.array([x for x in Template['RADI'].split()]\
            ,dtype=float)
    if disk_length is None: disk_length = Radii[-1]
    bar_brightness = equal_length(bar_brightness,disks)
    bar_angle = equal_length(bar_angle,disks)
    bar_velocities = equal_length(bar_velocities,disks)
  
    for i,disk in enumerate(disks):
        if disk == 1:
            ext = ''
        else:
            ext = f'_{disk}'
        if length_in is None:  
            length = calculate_length(Template,Radii,disk_length,ext)
        else:
            length=copy.deepcopy(length_in)
        bar_width = thickness*length # the bar thickness is 10% of the n
        # The width has to be 180 when R < width and 180*width/(pi*r)
        width = np.zeros(len(Radii))
        width[Radii <= bar_width] = 180.
        # the angle made up of the radius and width *2.
        width[Radii > bar_width] = 360./np.pi * \
            np.arcsin(bar_width/Radii[Radii > bar_width])
        # We set the full brightness to the maximum of the disk if not given
       
        if bar_brightness[i][0] is None:
            #Has to be float else Omega conf will not accept
            sbr_max = float(np.max([float(x) for x in Template[f'SBR{ext}'].split()]))
            bar_brightness[i] = [sbr_max if rad < length else 0. for rad in Radii]
        
            #bar_brightness[i][list(np.where(Radii < length))] = np.max([float(x) \
            #    for x in Template[f'SBR{ext}'].split()])
        if not isiterable(bar_brightness[i]):
            bar_brightness[i] = [bar_brightness[i] if rad < length else 0. for rad in Radii]
        if not isiterable(bar_angle[i]):
            bar_angle[i] = [bar_angle[i] for rad in Radii]
   
        # Get the number of disks present
        ndisk = int(Template["NDISKS"])
        ndisk += 1
    
        # we also need streaming motions
        if bar_velocities[i][0] is None:
            bar_velocities[i] = [-25. for rad in Radii ] 
        if not isiterable(bar_velocities[i]):
            bar_velocities[i] = [bar_velocities[i] for rad in Radii]
        #copy the input disk
        Template_New = copy_disk(Template_New, disk)
        #Template_New = copy.deepcopy(Template)


        #We offset by 37 deg.

        Template_New.insert(f"VSYS_{ndisk:d}", f"AZ1P_{ndisk:d}",
            f"{' '.join([str(e) for e in bar_angle[i]])}")
        Template_New.insert(f"AZ1P_{ndisk:d}", f"AZ1W_{ndisk:d}"
            , f"{' '.join([str(e) for e in width])}")
       
        Template_New.insert(f"AZ1W_{ndisk:d}", f"AZ2P_{ndisk:d}"
        , f"{' '.join([str(x+180.) for x in bar_angle[i]])}")
        Template_New.insert(f"AZ2P_{ndisk:d}", f"AZ2W_{ndisk:d}"
        ,f"{' '.join([str(e) for e in width])}")
        # And we add streaming motions to the bar km/s
        Template_New.insert(f"VROT_{ndisk:d}", f"VRAD_{ndisk:d}"
        , f"{' '.join([str(e) for e in bar_velocities[i]])}")
        Template_New[f"SBR_{ndisk:d}"] = f"{' '.join(str(e) for e in bar_brightness[i])}"
   

    return Template_New
create_bar.__doc__ = f'''
NAME:
   create_bar

PURPOSE:
    Create a bar

CATEGORY:
   agc

INPUTS:
    velocity,Radii,disk_brightness,Template, disk=1,WarpStart=-1
    velocity = Rotation curve
    Radii= the radii
    disk_brightness = the SBR profile of the disk
    disk=1,WarpStart=-1,Bar="No_Bar"

OPTIONAL INPUTS:
    disk = 1
    Are we creating in the approaching (1) or receding side (2)

    WarpStart=-1
    radius of the start of the warp

    Bar="No_Bar"
    Boolean whether bar is included or not

OUTPUTS:
    phase = the phases of the arms
    brightness = the brightness amplitude of the arms
    width =the used width in deg

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def create_bar_harmonics(Template, disks=[1], length_in=None, \
        disk_length =None, bar_angle= [[37.]],\
        bar_velocities = [[None]], bar_brightness = [[None]]):
    #disk_length= the size of the visible disk (WarpStart in AGC) in arcsec, 
    # set to last radius if unset

    '''Routine to create a Bar into a model'''
    Template_New =copy.deepcopy(Template)
    Radii = np.array([x for x in Template['RADI'].split()]\
            ,dtype=float)
    if disk_length is None: disk_length = Radii[-1]
    bar_brightness = equal_length(bar_brightness,disks)
    bar_angle = equal_length(bar_angle,disks)
    bar_velocities = equal_length(bar_velocities,disks)
  

    for i,disk in enumerate(disks):
        if disk == 1:
            ext = ''
        else:
            ext = f'_{disk:d}'
        if length_in is None:  
            length = calculate_length(Template,Radii,disk_length,ext)
        else:
            length=copy.deepcopy(length_in)
       
        
            #bar_brightness[i][list(np.where(Radii < length))] = np.max([float(x) \
            #    for x in Template[f'SBR{ext}'].split()])
       
        if bar_brightness[i][0] is None:
            #Has to be float else Omega conf will not accept
            sbr_max = float(np.max([float(x) for x in Template[f'SBR{ext}'].split()])) 
            bar_brightness[i] = [sbr_max if rad < length else 0. for rad in Radii]
        if not isiterable(bar_brightness[i]):
            bar_brightness[i] = [bar_brightness[i] if rad < length else 0. for rad in Radii]
        if not isiterable(bar_angle[i]):
            bar_angle[i] = [bar_angle[i] for rad in Radii]   
        if bar_velocities[i][0] is None:
            bar_velocities[i] = [-25. if rad < length else 0. for rad in Radii] 
        if not isiterable(bar_velocities[i]):
            bar_velocities[i] = [bar_velocities[i] if rad < length else 0. for rad in Radii]
        if len(bar_velocities[i]) < len(Radii):
            tmp = []
            for j,rad in enumerate(Radii):
                if j < len(bar_velocities[i]):
                    tmp.append(bar_velocities[i][j])
                else:
                    tmp.append(bar_velocities[i][-1] if rad < length else 0.)     
            bar_velocities[i] = tmp


        velocity_amplitude = [x/2. for x in bar_velocities[i]]
        bar_velocities[i] = [x/2. for x in bar_velocities[i]]
        
        #First we add the velocity Harmonics
        Template_New.insert(f"VROT{ext}", f"RO2P{ext}",
            f"{' '.join([str(e) for e in bar_angle[i]])}") 
        Template_New.insert(f"RO2P{ext}", f"RA2P{ext}",
            f"{' '.join([str(e-45.) for e in bar_angle[i]])}")
        Template_New.insert(f"RO2P{ext}", f"RO2A{ext}",
            f"{' '.join([str(e) for e in velocity_amplitude])}")
        Template_New.insert(f"RA2P{ext}", f"RA2A{ext}",
            f"{' '.join([str(e) for e in velocity_amplitude])}")
        Template_New.insert(f"RO2A{ext}", f"VRAD{ext}"
        , f"{' '.join([str(e) for e in bar_velocities[i]])}")
        
       
        #Then the brightnes harmonics in the same way 
        Template_New.insert(f"SBR{ext}", f"SM2P{ext}",
            f"{' '.join([str(e) for e in bar_angle[i]])}") 
        Template_New.insert(f"SBR{ext}", f"SM2A{ext}",
            f"{' '.join([str(e) for e in bar_brightness[i]])}")
       
        #RA2P,SM2P and RO2P should be tied together
        Template_New['VARY'] = f'{Template_New["VARY"]}, RA2P{ext} 1:{Template_New["NUR"]} RO2P{ext} 1:{Template_New["NUR"]} SM2P{ext} 1:{Template_New["NUR"]}'
    






    return Template_New



create_bar_harmonics.__doc__ = f'''
NAME:
   create_bar

PURPOSE:
    Create a bar by using 2nd order harmonics

CATEGORY:
   agc

INPUTS:
    velocity,Radii,disk_brightness,Template, disk=1,WarpStart=-1
    velocity = Rotation curve
    Radii= the radii
    disk_brightness = the SBR profile of the disk
    disk=1,WarpStart=-1,Bar="No_Bar"

OPTIONAL INPUTS:
    disk = 1
    Are we creating in the approaching (1) or receding side (2)

    WarpStart=-1
    radius of the start of the warp

    Bar="No_Bar"
    Boolean whether bar is included or not

OUTPUTS:
    phase = the phases of the arms
    brightness = the brightness amplitude of the arms
    width =the used width in deg

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


import subprocess
def main():
        argv = sys.argv[1:]
        cfg_defaults = OmegaConf.structured(defaults)
        # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
        inputconf = OmegaConf.from_cli(argv)
        cfg = OmegaConf.merge(cfg_defaults,inputconf)
        if cfg.input_def is None:
            raise InputError(f'''There is no default for the input file please provide one by calling:
add_bar_to_model input_def=<input tirific def file>''') 

        if cfg.bar_type.lower() not in ['free','harmonics']:
            raise InputError(f'''The bar_type has to be either Free or Harmonics''')


        if cfg.output_def is None:
            base = os.path.splitext(cfg.input_def)[0]
            cfg.output_def =f'{base}_bar.def'
        # Read the input_file
        print(f'Reading {cfg.input_def}')
        Template = tirific_template(f'{cfg.input_def}')
        if cfg.disk_length is None:
            cfg.disk_length = get_model_DHI(f'{cfg.input_def}')/2.

        base_files =  os.path.splitext(Template['OUTSET'])[0]
        Template['OUTSET'] = f'{base_files}_bar.fits'
        Template['TIRDEF']= f'{base_files}_bar.def'
        bar_options= ['Free','Harmonics']

        if cfg.bar_type.lower() == 'free':
            Template_New = create_bar(Template,disks=cfg.input_disk,\
                length_in = cfg.bar_length,disk_length=cfg.disk_length,\
                thickness = cfg.bar_thickness,bar_angle=cfg.bar_angle,\
                bar_velocities = cfg.bar_velocities,\
                bar_brightness=cfg.bar_brightness) 
            
        elif  cfg.bar_type.lower() == 'harmonics':
            Template_New = create_bar_harmonics(Template,disks=cfg.input_disk,\
                length_in = cfg.bar_length,disk_length=cfg.disk_length,\
                bar_angle=cfg.bar_angle,\
                bar_velocities = cfg.bar_velocities,\
                bar_brightness=cfg.bar_brightness) 
        else:
            raise InputError(f'''The bar_type {cfg.bar_type} is not valid.
please choose from {','.join(bar_options)}.''')

        
        print(f'Writing {cfg.output_def}')
        with open(f'{cfg.output_def}', 'w') as file:
            for key in Template_New:
                if key[0:5] == 'EMPTY':
                    file.write('\n')
                else:
                    file.write(f"{key}= {Template_New[key]} \n")

      


if __name__ == '__main__':
    #with VizTracer(output_file="FAT_Run_Viztracer.json",min_duration=1000) as tracer:
    main(sys.argv[1:])
