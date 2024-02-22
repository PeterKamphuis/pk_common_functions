

import numpy as np
import copy
import os
import re
import warnings
import traceback
from astropy.wcs import WCS
from collections import OrderedDict #used in Proper_Dictionary

from astropy.io import fits
from scipy.ndimage import rotate, map_coordinates
from scipy.optimize import curve_fit, OptimizeWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Ellipse

class InputError(Exception):
    pass
class FittingError(Exception):
    pass
# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
    def __setitem__(self, key, value):
        if key not in self:
            # If it is a new item we only allow it if it is not Configuration or Original_Cube or if we are in setup_configuration
            try:
                function,variable,empty = traceback.format_stack()[-2].split('\n')
            except ValueError: 
                function,variable = traceback.format_stack()[-2].split('\n')
            function = function.split()[-1].strip()
            variable = variable.split('[')[0].strip()
            if variable == 'Original_Configuration' or variable == 'Configuration':
                if function != 'setup_configuration':
                    raise ProgramError("FAT does not allow additional values to the Configuration outside the setup_configuration in support_functions.")
        OrderedDict.__setitem__(self,key, value)
    #    "what habbens now")
    def insert(self, existing_key, new_key, key_value):
        done = False
        if new_key in self:
            self[new_key] = key_value
            done = True
        else:
            new_orderded_dict = self.__class__()
            for key, value in self.items():
                new_orderded_dict[key] = value
                if key == existing_key:
                    new_orderded_dict[new_key] = key_value
                    done = True
            if not done:
                new_orderded_dict[new_key] = key_value
                done = True
                print(
                    f"----!!!!!!!! YOUR {new_key} was appended at the end as you provided the non-existing {existing_key} to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")

Proper_Dictionary.__doc__=f'''
A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
'''


        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)
def convertRADEC(RAin,DECin,invert=False, colon=False, verbose=False):
    if verbose:
        print(f'''CONVERTRADEC: Starting conversion from the following input.
{'':8s}RA = {RAin}
{'':8s}DEC = {DECin}
''')
    RA = copy.deepcopy(RAin)
    DEC = copy.deepcopy(DECin)
    if not invert:
        try:
            _ = (e for e in RA)
        except TypeError:
            RA= [RA]
            DEC =[DEC]
        for i in range(len(RA)):
            xpos=RA
            ypos=DEC
            xposh=int(np.floor((xpos[i]/360.)*24.))
            xposm=int(np.floor((((xpos[i]/360.)*24.)-xposh)*60.))
            xposs=(((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
            yposh=int(np.floor(np.absolute(ypos[i]*1.)))
            yposm=int(np.floor((((np.absolute(ypos[i]*1.))-yposh)*60.)))
            yposs=(((((np.absolute(ypos[i]*1.))-yposh)*60.)-yposm)*60)
            sign=ypos[i]/np.absolute(ypos[i])
            if colon:
                RA[i]="{}:{}:{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}:{}:{:2.2f}".format(yposh,yposm,yposs)
            else:
                RA[i]="{}h{}m{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}d{}m{:2.2f}".format(yposh,yposm,yposs)
            if sign < 0.: DEC[i]='-'+DEC[i]
        if len(RA) == 1:
            RA = str(RA[0])
            DEC = str(DEC[0])
    else:
        if isinstance(RA,str):
            RA=[RA]
            DEC=[DEC]

        xpos=RA
        ypos=DEC

        for i in range(len(RA)):
            # first we split the numbers out
            tmp = re.split(r"[a-z,:]+",xpos[i])
            RA[i]=(float(tmp[0])+((float(tmp[1])+(float(tmp[2])/60.))/60.))*15.
            tmp = re.split(r"[a-z,:'\"]+",ypos[i])
            if float(tmp[0]) != 0.:
                DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+(float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))
            else:
                DEC[i] = float(np.absolute(float(tmp[0])) + ((float(tmp[1]) + (float(tmp[2]) / 60.)) / 60.))
                if tmp[0][0] == '-':
                    DEC[i] = float(DEC[i])*-1.
        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])
        else:
            RA =np.array(RA,dtype=float)
            DEC = np.array(DEC,dtype=float)
    return RA,DEC

convertRADEC.__doc__ =f'''
 NAME:
    convertRADEC

 PURPOSE:
    convert the RA and DEC in degre to a string with the hour angle

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    RAin = RA to be converted
    DECin = DEC to be converted

 OPTIONAL INPUTS:


    invert=False
    if true input is hour angle string to be converted to degree

    colon=False
    hour angle separotor is : instead of hms

 OUTPUTS:
    converted RA, DEC as string list (hour angles) or numpy float array (degree)

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def columndensity(levels,systemic = 100.,beam=None,\
        channel_width=None,column= False,arcsquare=False,solar_mass_input =False\
        ,solar_mass_output=False, verbose= False, linewidth= None):
    if beam is None:
        if not arcsquare:
            print(f'COLUMNDENSITY: A beam is required to make the proper calculation''')

    if channel_width == None:
        if arcsquare:
            channel_width = 1.
        else:
            print(f'COLUMNDENSITY: A channel width is required to make the proper calculation''')
    if linewidth != None:
        channel_width = linewidth/(np.sqrt(linewidth/channel_width))
    if verbose:
        print(f'''COLUMNDENSITY: Starting conversion from the following input.
{'':8s}Levels = {levels}
{'':8s}Beam = {beam}
{'':8s}channel_width = {channel_width}
''')
    beam=np.array(beam,dtype=float)
    f0 = 1.420405751786E9 #Hz rest freq
    c = 299792.458 # light speed in km / s
    pc = 3.086e+18 #parsec in cm
    solarmass = 1.98855e30 #Solar mass in kg
    mHI = 1.6737236e-27 #neutral hydrogen mass in kg
    if verbose:
        print(f'''COLUMNDENSITY: We have the following input for calculating the columns.
{'':8s}COLUMNDENSITY: level = {levels}, channel_width = {channel_width}, beam = {beam}, systemic = {systemic})
''')

    f = f0 * (1 - (systemic / c)) #Systemic frequency
    if arcsquare:
        #Should we have the (f0/f)**2 factor here????
        HIconv = 605.7383 * 1.823E18 * (2. *np.pi / (np.log(256.)))
        if column:
            # If the input is in solarmass we want to convert back to column densities
            if solar_mass_input:
                levels=np.array([x*solarmass/(mHI*pc**2) for x in levels],dtype=float)
            #levels=levels/(HIconv*channel_width)

            levels = levels/(HIconv*channel_width)
        else:

            levels = HIconv*levels*channel_width
            if solar_mass_output:
                levels=levels*mHI/solarmass*pc*pc
    else:
        if beam.size <2:
            beam= [beam,beam]
        b=beam[0]*beam[1]
        if column:
            if solar_mass_input:
                levels=levels*solarmass/(mHI*pc**2)
            TK = levels/(1.823e18*channel_width)
            levels = TK/(((605.7383)/(b))*(f0/f)**2)
        else:
            TK=((605.7383)/(b))*(f0/f)**2*levels
            levels = TK*(1.823e18*channel_width)
    if not column and solar_mass_input:
        levels = levels*mHI*pc**2/solarmass
    return levels

columndensity.__doc__ =f'''
 NAME:
    columndensity

 PURPOSE:
    Convert the various surface brightnesses to other units

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    levels = the values to convert

 OPTIONAL INPUTS:


    systemic = 100.
    the systemic velocity of the source

    beam  = [-1.,-1.]
    the FWHM of the beam in arcsec, if unset taken from Configuration

    channelwidth = -1. width of a channel in km/s
    channelwidth of the observation if unset taken from Configuration

    column = false
    if True input is columndensities else in mJy

    arcsquare=False
    If true then  input is assumed to be in Jy/arcsec^2.
    If the input is in Jy/arcsec^2*km/s then channelwidth must be 1.
    This is assumed when channelwidth is left unset

    solar_mass_input =False
    If true input is assumed to be in M_solar/pc^2

    solar_mass_output=False
    If true output is provided in M_solar/pc^2

 OUTPUTS:
    The converted values

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Cosmological column densities are taken from Meyer et al 2017
'''
#Convert fluxes based on a set of input and a fits header

def convert_fluxes(array,header,conversion='JB_to_Jykms', IRAM_beam_area = None):
    possible_conversions = ['JB_to_Jykms','JB_to_Jy','JB_to_Tmb','JB_to_Ta',\
                            'JB_to_Tmb_IRAM','Jy_to_Tmb_IRAM','Tmb_to_Jy',\
                            'Tmb_IRAM_to_Jy'  ]
    if conversion not in possible_conversions:
        raise InputError(f'''The conversion {conversion} is not registered
Please use on of  {' ,'.join(possible_conversions)}                         
''')
    conv_array=copy.deepcopy(array)
    if not isiterable(conv_array):
        conv_array = [conv_array]
    #if the input is Janky per beam we first convert to Jansky
    if conversion in ['JB_to_Jykms','JB_to_Jy','JB_to_Tmb','JB_to_Ta']:
        conv_array = JB_to_Jy(conv_array,header)
        if conversion == 'JB_to_Jy':
            return conv_array
        elif conversion == 'JB_to_Jykms':
            raise InputError('Not implemented yet')
    if conversion in ['JB_to_Tmb_IRAM','Jy_to_Tmb_IRAM']: 
        # For converting the Main beam temperature (Yan's) spectra we get (See email 16-02-2024)
        # With the efficiencies as listed on  https://publicwiki.iram.es/Iram30mEfficiencies
        return [x/5.3 for x in conv_array]
    if conversion in ['JB_to_Tmb','Jy_to_Tmb','Tmb_to_Jy']: 
        #Or the theoretical conversion
        #from https://science.nrao.edu/facilities/vla/proposing/TBconv
        # Note that this equation states that I is in mJy/beam but it is not but in mJy
        # See https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php and https://safe.nrao.edu/wiki/pub/Main/UsefulFormulas/planck.pdf
        nu = header['RESTFREQ']/1e9
        const = 1.2222e3
        if IRAM_beam_area == None:
            beam = header['BMAJ']*3600.*header['BMIN']*3600.
        else: 
            #This gives the option to define a beam different from the header
            beam = IRAM_beam_area
        #beam_IRAM = 22.6211340310488**2
        if  conversion in ['JB_to_Tmb','Jy_to_Tmb']: 
            return [x*1000.*const/(nu**2*beam) for x in conv_array] #Temperature in K
        if  conversion in ['Tmb_to_Jy','Tmb_to_JB']:
            conv_array =  [x*nu**2*beam/(1000.*const) for x in conv_array] # convert Temperature to Jy
            if 'Tmb_to_Jy':
                return conv_array
    if conversion in ['Tmb_IRAM_to_Jy']:
        # For converting the Main beam temperature (Yan's) spectra we get (See email 16-02-2024)
        # With the efficiencies as listed on  https://publicwiki.iram.es/Iram30mEfficiencies
        return [x*5.3 for x in conv_array]
 


# function for converting kpc to arcsec and vice versa
def convertskyangle(angle, distance=-1., unit='arcsec', \
        distance_unit='Mpc', physical=False, verbose=False):
    if distance == -1.:
        raise InputError(f'convertskyangle needs a distance')
    if verbose:
        print(f'''CONVERTSKYANGLE: Starting conversion from the following input.
    {'':8s}Angle = {angle}
    {'':8s}Distance = {distance}
''')

    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]

        # if physical is true default unit is kpc
    angle = np.array(angle,dtype=float)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc':
        distance = distance * 10 ** 3
    elif distance_unit.lower() == 'kpc':
        distance = distance
    elif distance_unit.lower() == 'pc':
        distance = distance / (10 ** 3)
    else:
        print(f'''CONVERTSKYANGLE: {distance_unit} is an unknown unit to convertskyangle.
{'':8s}CONVERTSKYANGLE: please use Mpc, kpc or pc.
''')
        raise InputError(f'CONVERTSKYANGLE: {distance_unit} is an unknown unit to convertskyangle.')
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print(f'''CONVERTSKYANGLE: {unit} is an unknown unit to convertskyangle.
{'':8s}CONVERTSKYANGLE: arcsec, arcmin or degree.
''')
            raise InputError(f'CONVERTSKYANGLE: {unit} is an unknown unit to convertskyangle.')


        kpc = 2. * (distance * np.tan(radians / 2.))
    else:
        if unit.lower() == 'kpc':
            kpc = angle
        elif unit.lower() == 'mpc':
            kpc = angle * (10 ** 3)
        elif unit.lower() == 'pc':
            kpc = angle / (10 ** 3)
        else:
            print(f'''CONVERTSKYANGLE: {unit} is an unknown unit to convertskyangle.
{'':8s}CONVERTSKYANGLE: please use Mpc, kpc or pc.
''')
            raise InputError(f'CONVERTSKYANGLE: {unit} is an unknown unit to convertskyangle.')

        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc

convertskyangle.__doc__ =f'''
 NAME:
    convertskyangle

 PURPOSE:
    convert an angle on the sky to a distance in kpc or vice versa

 CATEGORY:
    common_functions

 INPUTS:
    Configuration = Standard FAT configuration
    angle = the angles or lengths to be converted

 OPTIONAL INPUTS:


    distance=1.
    Distance to the galaxy for the conversion

    unit='arcsec'
    Unit of the angle or length options are arcsec (default),arcmin, degree, pc, kpc(default) and Mpc

    distance_unit='Mpc'gauss_parameters = np.array(['NaN','NaN','NaN'],dtype=float)
        gauss_covariance = np.array(['NaN','NaN','NaN'],dtype=float)

    Unit of the distance options are pc, kpc and Mpc

    physical=False
    if true the input is a length converted to an angle

 OUTPUTS:
    converted value or values

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def cutout_cube(filename,sub_cube, outname=None):

    if outname == None:
        outname = f'{os.path.splitext(filename)[0]}_cut.fits'

    Cube = fits.open(filename,uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    hdr = Cube[0].header

    if hdr['NAXIS'] == 3:
        print(sub_cube[0,0],sub_cube[0,1],sub_cube[1,0],sub_cube[1,1],sub_cube[2,0],sub_cube[2,1])
        data = Cube[0].data[sub_cube[0,0]:sub_cube[0,1],sub_cube[1,0]:sub_cube[1,1],sub_cube[2,0]:sub_cube[2,1]]
        hdr['NAXIS1'] = sub_cube[2,1]-sub_cube[2,0]
        hdr['NAXIS2'] = sub_cube[1,1]-sub_cube[1,0]
        hdr['NAXIS3'] = sub_cube[0,1]-sub_cube[0,0]
        hdr['CRPIX1'] = hdr['CRPIX1'] -sub_cube[2,0]
        hdr['CRPIX2'] = hdr['CRPIX2'] -sub_cube[1,0]
        hdr['CRPIX3'] = hdr['CRPIX3'] -sub_cube[0,0]
        #Only update when cutting the cube

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coordinate_frame = WCS(hdr)
            xlow,ylow,zlow = coordinate_frame.wcs_pix2world(1,1,1., 1.)
            xhigh,yhigh,zhigh = coordinate_frame.wcs_pix2world(hdr['NAXIS1'],hdr['NAXIS2'],hdr['NAXIS3'], 1.)
            xlim = np.sort([xlow,xhigh])
            ylim = np.sort([ylow,yhigh])
            zlim =np.sort([zlow,zhigh])/1000.

    elif hdr['NAXIS'] == 2:

        if len(sub_cube) == 3:
            sub_im = sub_cube[1:]
        elif len(sub_cube) == 2:
            sub_im = sub_cube
        else:
            print(f"We don't understand your idea your sub_cube = {len(sub_cube)} and the image {hdr['NAXIS']}")
        print(sub_im)
        data = Cube[0].data[sub_im[0,0]:sub_im[0,1],sub_im[1,0]:sub_im[1,1]]

        hdr['NAXIS1'] = sub_cube[1,1]-sub_cube[1,0]
        hdr['NAXIS2'] = sub_cube[0,1]-sub_cube[0,0]
        hdr['CRPIX1'] = hdr['CRPIX1'] -sub_im[1,0]
        hdr['CRPIX2'] = hdr['CRPIX2'] -sub_im[0,0]

    Cube.close()
    fits.writeto(outname,data,hdr,overwrite = True)
    return outname

cutout_cube.__doc__ =f'''
 NAME:
    cutout_cube

 PURPOSE:
    Cut filename back to the size of subcube, update the header and write back to disk.

 CATEGORY:
    fits_functions

 INPUTS:
    filename = name of the cube to be cut
    outname = name of the output file
    sub_cube = array that contains the new size as
                [[z_min,z_max],[y_min,y_max], [x_min,x_max]]
                adhering to fits' idiotic way of reading fits files.

 OPTIONAL INPUTS:

 OUTPUTS:
    the cut cube is written to disk.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fit_gaussian(x,y, covariance = False,errors = None, \
    verbose= False):
    if verbose:
        print(f'''FIT_GAUSSIAN: Starting to fit a Gaussian.
{'':8s}x = {x}
{'':8s}y = {y}
''')
    # Make sure we have numpy arrays
    x= np.array(x,dtype=float)
    y= np.array(y,dtype=float)
    # First get some initial estimates
    est_peak = np.nanmax(y)

    if errors == None:
        errors = np.full(len(y),1.)
        absolute_sigma = False
    else:
        absolute_sigma = True
    peak_location = np.where(y == est_peak)[0]
    if peak_location.size > 1:
        peak_location = peak_location[0]
    est_center = float(x[peak_location])
    est_sigma = np.nansum(y*(x-est_center)**2)/np.nansum(y)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        succes = False
        maxfev= int(100*(len(x)))
        increase =False
        while not succes:
            if verbose:
                print(f'''FIT_GAUSSIAN: Starting the curve fit with {maxfev}
''')
            try:
                gauss_parameters, gauss_covariance = curve_fit(gaussian_function, \
                        x, y,p0=[est_peak,est_center,est_sigma],sigma= errors,\
                        absolute_sigma= absolute_sigma,maxfev=maxfev)

                succes = True
            except OptimizeWarning:
                maxfev =  2000*(len(x))
                print(f'maxfev = {maxfev} which should break the loop {1000*(len(x))}')
            except RuntimeError as e:
                split_error = str(e)
                if 'Optimal parameters not found: Number of calls to function has reached maxfev' in \
                    split_error:
                    maxfev += 100*int(len(x))
                    if verbose:
                        print(f'''FIT_GAUSSIAN: We failed to find an optimal fit due to the maximum number of evaluations. increasing maxfev to {maxfev}
''')
                else:
                    print(f'''FIT_GAUSSIAN: failed due to the following error {split_error}
''')
                    raise RuntimeError(split_error)
            if maxfev >  1000*(len(x)):
                gauss_parameters = np.array(['NaN','NaN','NaN'],dtype=float)
                gauss_covariance = np.array(['NaN','NaN','NaN'],dtype=float)
                succes = True
                #print(f'Failed to find a proper fit to the gaussian')
                #raise FittingError("FIT_GAUSSIAN: failed to find decent gaussian parameters")

    #gauss_parameters, gauss_covariance = curve_fit(gaussian_function, x, y,p0=[est_peak,est_center,est_sigma],sigma= errors,absolute_sigma= absolute_sigma)
    if covariance:
        return gauss_parameters, gauss_covariance
    else:
        return gauss_parameters

fit_gaussian.__doc__ =f'''
 NAME:
    fit_gaussian
 PURPOSE:
    Fit a gaussian to a profile, with initial estimates
 CATEGORY:
    supprt_functions

 INPUTS:
    x = x-axis of profile
    y = y-axis of profile
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    covariance = false
    return to covariance matrix of the fit or not



 OUTPUTS:
    gauss_parameters
    the parameters describing the fitted Gaussian

 OPTIONAL OUTPUTS:
    gauss_covariance
    The co-variance matrix of the fit

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def gaussian_function(axis,peak,center,sigma):
    return peak*np.exp(-(axis-center)**2/(2*sigma**2))


def get_model_DHI(filename):
    #Get the sbrs
    radi,sbr,sbr_2,systemic = load_tirific(filename,\
        Variables = ['RADI','SBR','SBR_2','VSYS'],array=True)
    #convert to solar_mass/pc^2
    sbr_msolar = columndensity(sbr*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    sbr_2_msolar = columndensity(sbr_2*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    # interpolate these to ~1" steps
    new_radii = np.linspace(0,radi[-1],int(radi[-1]))
    new_sbr_msolar = np.interp(new_radii,radi,sbr_msolar)
    new_sbr_2_msolar = np.interp(new_radii,radi,sbr_2_msolar)

    index_1 = np.where(new_sbr_msolar > 1.)[0]
    index_2 = np.where(new_sbr_2_msolar > 1.)[0]
    if np.sum(sbr_2) > 0.:
        if index_1.size > 0 and index_2.size > 0:
            DHI = float(new_radii[index_1[-1]]+new_radii[index_2[-1]])
        elif index_1.size > 0:
            DHI = float(new_radii[index_1[-1]])
        elif index_2.size > 0:
            DHI = float(new_radii[index_2[-1]])
        else:
            DHI = float('NaN')
    else:
        if index_1.size > 0:
            DHI = float(new_radii[index_1[-1]]*2.)
        else:
            DHI = float('NaN')
    return DHI
get_model_DHI.__doc__ =f'''
 NAME:
    get_DHI

 PURPOSE:
    get the DHI as determined by the SBR profiles in the fit from the Tirific Template

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


    Model = 'Finalmodel'
    location of the def file to get DHI from. it should be in the fitting dir in the {{Model}}/{{Model}}.def

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



def isiterable(variable):
    '''Check whether variable is iterable'''
    #First check it is not a string as those are iterable
    if isinstance(variable,str):
        return False
    try:
        iter(variable)
    except TypeError:
        return False

    return True
isiterable.__doc__ =f'''
 NAME:
    isiterable

 PURPOSE:
    Check whether variable is iterable

 CATEGORY:
    support_functions

 INPUTS:
    variable = variable to check

 OPTIONAL INPUTS:

 OUTPUTS:
    True if iterable False if not

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Convert Jansky/beam to Jansky bases on the header
def JB_to_Jy(array,header):
    beam_in_pixels = pixels_in_beam(header)
    return array/beam_in_pixels

def load_tirific(def_input,Variables = None,array = False,\
        ensure_rings = False ,dict=False):
    #Cause python is the dumbest and mutable objects in the FAT_defaults
    # such as lists transfer
    if Variables == None:
        Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI',\
                     'VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS',\
                     'SDIS','VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','XPOS_2',\
                     'YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']


    # if the input is a string we first load the template
    if isinstance(def_input,str):
        def_input = tirific_template(filename = def_input )

    out = []
    for key in Variables:

        try:
            out.append([float(x) for x  in def_input[key].split()])
        except KeyError:
            out.append([])
        except ValueError:
            out.append([x for x  in def_input[key].split()])

    #Because lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand

    if array:
        tmp = out
        #We can ensure that the output has the same number of values as there are rings
        if ensure_rings:
            length=int(def_input['NUR'])
        else:
            #or just take the longest input as the size
            length = max(map(len,out))
        #let's just order this in variable, values such that it unpacks properly into a list of variables
        out = np.zeros((len(Variables),length),dtype=float)
        for i,variable in enumerate(tmp):
            if len(variable) > 0.:
                out[i,0:len(variable)] = variable[0:len(variable)]

    if dict:
        tmp = {}
        for i,var in enumerate(Variables):
            tmp[var] = out[i]
        out = tmp
    elif len(Variables) == 1:
        out= out[0]
    #print(f'''LOAD_TIRIFIC: We extracted the following profiles from the Template.
#{'':8s}Requested Variables = {Variables}
#{'':8s}Extracted = {out}
#''')
    #Beware that lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand
    # However if you make a np. array from it make sure that you specify float  or have lists of the same length else you get an array of lists which behave just as dumb

    return out
load_tirific.__doc__ =f'''
 NAME:
    load_tirific

 PURPOSE:
    Load values from variables set in the tirific files

 CATEGORY:
    common_functions

 INPUTS:
    def_input = Path to the tirific def file or a FAT tirific template dictionary

 OPTIONAL INPUTS:
    Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']


    array = False
        Specify that the output should be an numpy array with all varables having the same length

    ensure_rings =false
        Specify that the output array should have the length of the NUR parameter in the def file

    dict = False
        Return the output as a dictionary with the variable names as handles
 OUTPUTS:
    outputarray list/array/dictionary with all the values of the parameters requested

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    This function has the added option of a dictionary compared to pyFAT
'''
def pixels_in_beam(hdr):
    beamarea=(np.pi*abs(hdr['BMAJ']*hdr['BMIN']))/(4.*np.log(2.))
    beam_in_pixels = beamarea/(abs(hdr['CDELT1'])*abs(hdr['CDELT2']))
    return beam_in_pixels

def reduce_header_axes(hdr,axes= 3):
    ax = ['CDELT','CTYPE','CUNIT','CRPIX','CRVAL','NAXIS']
    while hdr['NAXIS'] > axes:
        for par in ax:
            try:
                hdr.remove(f"{par}{hdr['NAXIS']}")
            except KeyError:
                pass

        hdr['NAXIS'] -= 1
    return hdr

def reduce_data_axes(data,axes= 3):
    while len(data.shape) > axes:
        data = data[0,:]
    return data


def rotateImage(image, angle, pivot,order=1):
    padX = [int(image.shape[1] - pivot[0]), int(pivot[0])]
    padY = [int(image.shape[0] - pivot[1]), int(pivot[1])]
    imgP = np.pad(image, [padY, padX], 'constant')
    imgR = rotate(imgP, angle, axes=(1, 0), reshape=False,order=order)
    return imgR[padY[0]: -padY[1], padX[0]: -padX[1]]

def rotateCube(Cube, angle, pivot,order=1):
    padX= [int(Cube.shape[2] - pivot[0]), int(pivot[0])]
    padY= [int(Cube.shape[1] - pivot[1]), int(pivot[1])]
    imgP= np.pad(Cube, [[0, 0], padY, padX], 'constant')
    #Use nearest neighbour as it is exact enough and doesn't mess up the 0. and is a lot faster
    imgR = rotate(imgP, angle, axes =(2, 1), reshape=False,order=order)
    return imgR[:, padY[0]: -padY[1], padX[0]: -padX[1]]
rotateCube.__doc__=f'''
 NAME:
    rotateCube(Cube, angle, pivot)

 PURPOSE:
    rotate a cube in the image plane

 CATEGORY:
    common_functions

 INPUTS:
    Cube = the cube data array
    angle = the angle to rotate under
    pivot = the point around which to rotate

 OPTIONAL INPUTS:

 OUTPUTS:
    the rotated cube

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def regrid_array(oldarray, Out_Shape):
    oldshape = np.array(oldarray.shape)
    newshape = np.array(Out_Shape, dtype=float)
    ratios = oldshape/newshape
        # calculate new dims
    nslices = [ slice(0,j) for j in list(newshape) ]
    #make a list with new coord
    new_coordinates = np.mgrid[nslices]
    #scale the new coordinates
    for i in range(len(ratios)):
        new_coordinates[i] *= ratios[i]
    #create our regridded array
    newarray = map_coordinates(oldarray, new_coordinates,order=1)
    if any([x != y for x,y in zip(newarray.shape,newshape)]):
        print("Something went wrong when regridding.")
    return newarray
regrid_array.__doc__ =f'''
 NAME:
    regridder
 PURPOSE:
    Regrid an array into a new shape through the ndimage module
 CATEGORY:
    fits_functions

 INPUTS:
    oldarray = the larger array
    newshape = the new shape that is requested

 OPTIONAL INPUTS:

 OUTPUTS:
    newarray = regridded array

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    scipy.ndimage.map_coordinates, np.array, np.mgrid

 NOTE:
'''
def update_disk_angles(Tirific_Template, verbose = False):
    extension = ['','_2']
    for ext in extension:
        PA = np.array(load_tirific(Tirific_Template,Variables = [f'PA{ext}']),dtype=float)
        inc = np.array(load_tirific(Tirific_Template,Variables = [f'INCL{ext}']),dtype=float)
        if verbose:
            print(f'''UPDATE_DISK_ANGLES: obtained  this from the template
{'':8s} inc{ext} = {inc}
{'':8s} PA{ext} = {PA}
''')
        angle_adjust=np.array(np.tan((PA[0]-PA)*np.cos(inc*np.pi/180.)*np.pi/180.)*180./np.pi,dtype = float)
        if ext == '_2':
            angle_adjust[:] +=180.
        if verbose:
            print(f'''UPDATE_DISK_ANGLES: adusting AZ1P{ext} with these angles
{'':8s}{angle_adjust}
''')
        Tirific_Template.insert(f'AZ1W{ext}',f'AZ1P{ext}',f"{' '.join([f'{x:.2f}' for x in angle_adjust])}")
update_disk_angles.__doc__ =f'''
 NAME:
    update_disk_angles

 PURPOSE:
    Update the AZ1W and AZ1P parameters to match the warp

 CATEGORY:
    modify_template

 INPUTS:
    Tirific_Template = Standard FAT Tirific Template, so a proper dictionary

 OPTIONAL INPUTS:


 OUTPUTS:
    Updated template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def tirific_template(filename = ''):
    if filename == '':
        raise InputError(f'Tirific_Template does not know a default')
    else:
        with open(filename, 'r') as tmp:
            template = tmp.readlines()
    result = Proper_Dictionary()
    counter = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip().upper())
        if key == '':
            result[f'EMPTY{counter}'] = line
            counter += 1
        else:
            result[key] = str(line.split('=')[1].strip())
    return result
tirific_template.__doc__ ='''
 NAME:
    tirific_template

 PURPOSE:
    Read a tirific def file into a dictionary to use as a template.
    The parameter ill be the dictionary key with the values stored in that key

 CATEGORY:
    read_functions

 INPUTS:
    filename = Name of the def file

 OPTIONAL INPUTS:
    filename = ''
    Name of the def file, if unset the def file in Templates is used



 OUTPUTS:
    result = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      split, strip, open

 NOTE:
'''

def write_tirific(Tirific_Template, name = 'tirific.def',\
                full_name = False  ):
    #IF we're writing we bump up the restart_ID and adjust the AZ1P angles to the current warping
    if int(Tirific_Template['NUR']) == 2:
        update_disk_angles(Tirific_Template )
    if 'RESTARTID' in Tirific_Template:
        Tirific_Template['RESTARTID'] = str(int(Tirific_Template['RESTARTID'])+1)
   
    if full_name:
        file_name = name
    else:
        current_dir = os.getcwd()
        file_name = f'{current_dir}/{name}'
    with open(file_name, 'w') as file:
        for key in Tirific_Template:
            if key[0:5] == 'EMPTY':
                file.write('\n')
            else:
                file.write((f"{key}= {Tirific_Template[key]} \n"))
write_tirific.__doc__ =f'''
 NAME:
    tirific

 PURPOSE:
    Write a tirific template to file

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:


    name = 'tirific.def'
    name of the file to write to

 OUTPUTS:
    Tirific def file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 '''