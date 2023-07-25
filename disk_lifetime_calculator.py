import numpy as np
import argparse
import sys
from src import *

def main():
    """
    This script caculates a disk lifetime for given stellar mass and X-ray luminosity.
    You call it via
    
    python demo.py Mstar Lx
    
    where Mstar is the stellar mass in units of solar masses and Lx is the corresponding X-ray luminosity in units of erg/s. The resulting disk lifetime is then given in units of million years.
    """
    if len(sys.argv) != 3:
        print("Please provide exactly two command-line arguments: python disk_lifetime_calculator.py Mstar[Msol] and Lx[erg/s]. ")
        return

    Mstar = float(sys.argv[1])
    Lx = float(sys.argv[2])

    print("Mstar: ", Mstar, " Msun")
    print("Lx: ", Lx, " erg/s")
    
    tdisc = tdisc_mass(np.log10(Lx), Mstar)
    
    print("The resulting disk lifetime is ", np.round(tdisc,2), " Myr." )

if __name__ == "__main__":
    main()
