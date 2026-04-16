#
#region Titles and Header
# Nature: Optimization
# Methodology: Set trimming
##################################################################################################################
# VERSION        DATE            AUTHOR                    DESCRIPTION OF CHANGES MADE
#   0.0          01-Dec-2024     Mariana Mello              Proposed
#   0.2          12-May-2025     Mariana Mello             Changed name from 'pd' to 'm_p'
##################################################################################################################
#endregion


#region Import Library
from STHE.Calculations import Calculations_STHE_countingtable
from math import pi
#endregion

#region Calculations


# def STHE_Dist_Model(Ds, dte, Npt, rp, lay, L)     # , m_p):
    #
'm_p': {
            # Hot stream
            'mh': 20,           # Flow rate (kg*s**-1)
            'Cph': 2840,        # Heat capacity (J*(kg*K)**-1)
            'mih': 0.002,       # Viscosity (Pa*s)
            
           # Cold stream
            'mc': 60,           # Flow rate (kg*s**-1)
            'Cpc': 4187,        # Heat capacity (J*(kg*K)**-1)
            'mic': 0.0005,      # Viscosity (Pa*s)
 
            'Tci': 47,                    # Inlet temperature of the cold stream (oC)
            'Tco': 56,                    # Outlet temperature of the cold stream (oC)
            'Thi': 120,                   # Inlet temperature of the hot stream (oC)
            'Tho': 80,                    # Outlet temperature of the hot stream (oC)
            'dte': xx,                    
            'L': xx
            'Nb': 5,
             'U': 500
                 },
#    m_p['Thi'], m_p['Tci']; m_p['mh'], m_p['Cph']; m_p['mc], m_p['Cpc']
    # Calculate total area
       # Ntt = Calculations_STHE_countingtable.STHE_counting_table(Ds, dte, Npt, rp, lay, m_p)
    A = Ntt * pi * dte * L
    # Divide by the number of baffles to get 
    AB= A/ (Nb+1)
    # generqate the vectors for the internal temperatures
    # solve system of equations  to obtain Tho, Tco
    # 
print  Tho, Tco
