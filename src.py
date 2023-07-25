## Collection of scripts to calculate disk lifetimes based on X-ray photoevaporation.

import numpy as np

def scaling_radius(Mstar):
    """
    Calculates the "characteristic" radius of a disk of given stellar mass.
    scaling from Eq. 8 in Wilhelm & Zwart (2021), but choose 50au rather than 200au as "characteristic radius" 
    """
    return 50.*Mstar**0.45


def viscous_timescale(R1, nu0):
    """
    viscous timescale (tnu) at R1 in years
    """
    return R1**2./(3.*nu0)


def nu0(Mdot0, R1, Mdisc0):
    """
    initial viscosity at R1 in the disc
    """
    return 2.*Mdot0*R1**2./(3.*Mdisc0)


def viscous_Mdot_init(Mstar):
    """
    Viscous accretion rates in Lupus from Alcala+(2017); in solar masses/yr
    """
    if Mstar < 0.2:
        return 10.**(4.58*np.log10(Mstar)-6.11)
    else:
        return 10.**(1.37*np.log10(Mstar)-8.46)
    

def grav_radius(Mstar):
    """
    Gravitational radius in astronomical units; from Owen+2012, Eq. 1
    """
    Tgas = 1e4 # K;  temperature of the heated gas for EUV photoevaporation
    
    return 8.9 * 1./(Tgas/1e4) * Mstar


def MdotMstar(Mstar, Lx='variable'):
    """
    Scaling of wind mass-loss rate from XEUV-photoevaporation as a function of stellar mass.
    From Picogna et al. 2021
    """
    if Lx=='variable':
        return 3.93e-8 * Mstar
    elif Lx=='fixed':
        return 2.34e-8 * Mstar + 6.23e-9
    
 
def MdotLx(Lx, fit='P19'):
    """
    Scaling of wind mass-loss rate from XEUV-photoevaporation as a function of X-ray luminosity.
    From...
    'P19': Picogna et al. (2019) 
    'E21': Eq. 9 in Ercolano et al. (2021) --> as a function of soft X-ray luminosity only!!!!!
    """
    if fit=='P19':
        #AL = -2.7326
        #BL = 3.3307
        #CL = -2.9868e-3
        #DL = -7.2580
        #return 10. ** (AL * np.exp(((np.log(np.log10(Lx)) - BL) ** 2.) / CL) + DL)
        
        # updated fit, which is better for low Lx
        x0 = 29.11676386
        a = 3.17632301
        b = 1.83224246
        c = -10.43121048
        return 10.**(a / (1. + np.exp(-b*(np.log10(Lx)-x0)))+c)
    
    elif fit=='E21':
        AL = -1.947e17
        BL = -1.572e-4
        CL = -2.866e-1
        DL = -6.694
        return 10. ** (AL * np.exp(((np.log(np.log10(Lx)) - BL) ** 2.) / CL) + DL)
  

def Lxsoft(Lx):
    """
    Calculates the value of the "soft" X-ray luminosity for a given LX,
    based on a fit of the data points in Table 4 of Ercolano+2021.
    """
    a = 4.17087087e-01 
    b = 3.11111114e+28
    
    return a*Lx+b

def LxMstar(Mstar, region='Taurus'):
    """ Calculation of the mean X-ray luminosity for a given stellar mass.
    """
    if region == 'Taurus':
        # from Guedel et al. 2007
        return 10. ** (1.54 * np.log10(Mstar) + 30.31)
    elif region == 'ONC':
        # from Preibisch+2005
        return 10. ** (1.44 * np.log10(Mstar) + 30.37)
    else:
        print("Choose either 'Taurus' or 'ONC' as input regions.")
    
    
def tdisc_mass(logLx, Mstar):
    """
    total disc lifetime in Myr for the new mass-dependent XPE profiles by Picogna+2021
    """
    R1 = scaling_radius(Mstar) # scaling radius of the disk
    Mdisc0 = 0.14*Mstar # initial disk mass
    Mdot_acc = viscous_Mdot_init(Mstar) # viscous accretion rate
    Mdot0 = (0.5*Mdisc0/(5e6))**3./(Mdot_acc**2.) # INITIAL accretion rate at t=0
    Rg = grav_radius(Mstar) # gravitational radius

    # wind mass-loss rate
    Mdot_Lx = MdotMstar(Mstar, Lx='variable')*MdotLx(Lxsoft(10.**logLx), fit='E21')/MdotLx(Lxsoft(LxMstar(Mstar)), fit='E21')
    
    # viscous timescale
    tnu = viscous_timescale(Rg, nu0(Mdot0, R1, Mdisc0))
    
    tclear = (Mdisc0/(2.*Mdot_acc**(1./3.)*Mdot_Lx**(2./3.))) # clearing timescale
    tdisc = (tclear+tnu)/1e6 # disk lifetimes in Myr
    
    return tdisc


def tdisc_Picogna(logLx):
    """
    total disc lifetime in Myr for the XPE profile by Picogna+2019
    """
    x0 = 28.4869639
    a = 59.64474636
    b = -2.59269424
    c = 1.07840619

    return a / (1. + np.exp(-b*(logLx-x0)))+c


def tdisc_Owen(logLx):
    """
    total disc lifetime in Myr for the XPE profile by Owen+2012
    """
    a = 8.87988245e+19
    b = 1.47988442e+00
    c = -8.61924090e-02

    return a * np.exp(-b * logLx) + c
