# pk_common_functions
There are many functions in pyFAT and pyHIARD that are useful. However we do not always want to dress programs up with all pyFAT machinery. This is to have a central repository for these functions to keep them up to date. !!! --- These Functions should not be used for code that is intended for distribution.

There is no pypi but downloading the repository and then installing from the directory should allow for a pip install. Documentation is limited to this README

i.e:

pip install <path_to_repository>

should install the module pk_common_functions.

function can then be accessed through:

from pk_common_functions.functions import functions

or

import pk_common_functions.functions


# Script functions
--------------------------------------

*add_bar_to_model*

add a bar to your favorite tirific model

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
  
  




