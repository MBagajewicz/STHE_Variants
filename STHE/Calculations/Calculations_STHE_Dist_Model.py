#
#region Titles and Header
# Nature: Distributed Model with Rita
# Methodology: Solve System of Equations
##################################################################################################################
# VERSION        DATE            AUTHOR                    DESCRIPTION OF CHANGES MADE
#   0.0          16-Apr-2026    Miguel Bagajewicz            First draft
##################################################################################################################
#endregion


#region Import Library
#from STHE.Calculations import Calculations_STHE_countingtable
from math import pi
import numpy as np
#from scipy.optimize import fsolve
from scipy.optimize import least_squares
import math
#endregion

#region Calculations
N = int(input("N:"))




def STHE_Dist_Model(N):#(Ds, dte, Npt, rp, lay, L):
    m_p =   {
                # Hot stream
                'mh': 20,          # Flow rate (kg*s**-1)
                'Cph': 2840,        # Heat capacity (J*(kg*K)**-1)
                #'mih': 0.002,       # Viscosity (Pa*s)
                
                # Cold stream
                'mc': 60,           # Flow rate (kg*s**-1)
                'Cpc': 4187,        # Heat capacity (J*(kg*K)**-1)
                #'mic': 0.0005,      # Viscosity (Pa*s)
    
                'Tci': 47,                    # Inlet temperature of the cold stream (oC)
                #'Tco': 56,                    # Outlet temperature of the cold stream (oC)
                'Thi': 120,                   # Inlet temperature of the hot stream (oC)
                #'Tho': 80,                    # Outlet temperature of the hot stream (oC)
                'dte': 100,                    
                'L': 50,
                'Nb': 5,
                'U': 500,
                'Ntt':100
            }
        #m_p['Thi'], m_p['Tci']; m_p['mh'], m_p['Cph']; m_p['mc], m_p['Cpc']
        # Calculate total area
        # Ntt = Calculations_STHE_countingtable.STHE_counting_table(Ds, dte, Npt, rp, lay, m_p)
    A = m_p['Ntt'] * pi *  m_p['dte'] *  m_p['L']
    # Divide by the number of baffles to get 
    AB= A/ (m_p['Nb']+1)
    # generqate the vectors for the internal temperatures
    # solve system of equations  to obtain Tho, Tco
    UAB=m_p['U']*AB
    

    def LMTD(Tc_out,Tc_in,Th_out,Th_in):
        deltaT1=Th_out-Tc_in
        deltaT2=Th_in-Tc_out
        if deltaT1>0 and deltaT2>0 and abs(deltaT1-deltaT2)>1e-6:
            lmtd=(deltaT1-deltaT2)/math.log(deltaT1/deltaT2)
            return lmtd
        return 0



    def equations(vars):
        vars_positive=np.abs(vars)
        Th = vars_positive[:N]      #Tho=ThN=Th[N-1]   ThN=Th[N-1]
        Tc = vars_positive[N:]      #Tco=Tc1=Tc[0]     TcN=Tc[N-1]

        eqs = []


       ########################################################################################
       # model #1
        #eq1 = m_p['mc']*m_p['Cpc']*(Tc[N-1]-m_p['Tci']) - m_p['mh']*m_p['Cph']*(Th[N-2]-Th[N-1])
        #eq2 = m_p['mc']*m_p['Cpc']*(Tc[N-1]-m_p['Tci']) - UAB*(Th[N-1]-Tc[N-1])
        #eq5 = m_p['mc']*m_p['Cpc']*(Tc[1]-Tc[0]) - m_p['mh']*m_p['Cph']*(m_p['Thi']-Th[0])
        #eq6 = m_p['mc']*m_p['Cpc']*(Tc[1]-Tc[0]) - UAB*(Th[0]-Tc[0])
       # end of model #1 
       #########################################################################################
      
        eq1 = m_p['mc']*m_p['Cpc']*(Tc[N-1]-m_p['Tci']) - m_p['mh']*m_p['Cph']*(Th[N-2]-Th[N-1])
        eq2 = m_p['mc']*m_p['Cpc']*(Tc[N-1]-m_p['Tci']) - UAB*LMTD(Tc[N-1],m_p['Tci'],Th[N-1],Th[N-2])
        eq5 = m_p['mc']*m_p['Cpc']*(Tc[0]-Tc[1]) - m_p['mh']*m_p['Cph']*(m_p['Thi']-Th[0])
        eq6 = m_p['mc']*m_p['Cpc']*(Tc[0]-Tc[1]) - UAB*LMTD(Tc[0],Tc[1],Th[0],m_p['Thi'])
        
        eqs.append(eq1)
        eqs.append(eq2)
        eqs.append(eq5)
        eqs.append(eq6)



        for i in range(2,N):
        ########################################################################################
        # model #1
            #eq3 = m_p['mc']*m_p['Cpc']*(Tc[i-1]-Tc[i]) - m_p['mh']*m_p['Cph']*(Th[i-2]-Th[i-1])
            #eq4 = m_p['mc']*m_p['Cpc']*(Tc[i-1]-Tc[i]) - UAB*(Th[i-1]-Tc[i-1])

        # end of model #1 
        ########################################################################################

            eq3 = m_p['mc']*m_p['Cpc']*(Tc[i-1]-Tc[i]) - m_p['mh']*m_p['Cph']*(Th[i-2]-Th[i-1])
            eq4 = m_p['mc']*m_p['Cpc']*(Tc[i-1]-Tc[i]) - UAB*LMTD(Tc[i-1],Tc[i],Th[i-1],Th[i-2])

            eqs.append(eq3)
            eqs.append(eq4)

            eqs.append(1e6 * max(Th[i-1] - Th[i-2] + 1e-6, 0)) #the temperature is decreasing
            eqs.append(1e6 * max(Th[N-1] - Th[N-2] + 1e-6, 0))
            eqs.append(1e6 * max(Tc[i-1] - Tc[i-2] + 1e-6, 0))
            eqs.append(1e6 * max(Tc[N-1] - Tc[N-2] + 1e-6, 0))

            #eqs.append(1e10 * max(Tc[i-1]- Th[i-2] + 1e-6, 0)) # i the same area, Tc must be lower than Th

         


        return eqs
    
    vars = np.full(2*N,66)
    #lb=[m_p['Tci']]*(2*N)
    #ub=[m_p['Thi']]*(2*N)
    #,bounds=(lb,ub)
    result = least_squares(equations,vars)
    return result.x

T=STHE_Dist_Model(N)   #T includes Th1 Th2 ... ThN Tc1 Tc2...TcN
Th_sol = T[:N]
Tc_sol = T[N:]

print("Th:", Th_sol)
print("Tc:", Tc_sol)

print("Tho:", Th_sol[N-1])
print("Tco:", Tc_sol[0])
