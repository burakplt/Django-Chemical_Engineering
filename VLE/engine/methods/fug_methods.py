from math import exp, sqrt, log
from numpy import roots as np_roots
from ..chemsep_operation import EosInterface as dbcall
from chemeasy.settings import BASE_DIR

MODELS_URL = BASE_DIR+"/VLE/engine/Models/"

class PR76():

    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """PENG-ROBINSON equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.37464+1.54226*w-0.26992*w**2 #PR kappa value
            c = 0.45724*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.07780*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None:
               
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        if i==j:
                            kijs[(i,j)] = 0
                        else:
                            if kij_tune.get((i,j),None)!=None:
                                kijs[(i,j)] = calculate_kij(cs[i],cs[j],kij_tune[(i,j)] )
                            else:
                                kijs[(i,j)] = kijs[(j,i)]
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if i==j:
                        kijs[(i,j)] = 0
                    else:
                        if kij_input.get((i,j),None):
                            if abs(kij_input.get((i,j))) < 0.3:
                                kijs[(i,j)] = kij_input[(i,j)]
                            else:
                                kijs[(i,j)] = 0
                        else:
                            kijs[(i,j)] = kijs[(j,i)]
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)
        
        fugacity_coefficients = []
        for i in range(0,len(cs)):
            fugacity_coefficients.append( calculate_phi(i,T))

        return fugacity_coefficients, kijs  

class PR78 ():
    
    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """PENG-ROBINSON 78 equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            if w <= 491:
                kappa = 0.37464 + 1.54226*w - 0.26992*w**2 #PR kappa value
            else:
                kappa = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
            
            c = 0.457235*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.077796*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None:
               
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        if i==j:
                            kijs[(i,j)] = 0
                        else:
                            if kij_tune.get((i,j),None)!=None:
                                kijs[(i,j)] = calculate_kij(cs[i],cs[j],kij_tune[(i,j)] )
                            else:
                                kijs[(i,j)] = kijs[(j,i)]
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if i==j:
                        kijs[(i,j)] = 0
                    else:
                        if kij_input.get((i,j),None):
                            if abs(kij_input.get((i,j))) < 0.3:
                                kijs[(i,j)] = kij_input[(i,j)]
                            else:
                                kijs[(i,j)] = 0
                        else:
                            kijs[(i,j)] = kijs[(j,i)]

        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)
            
        fugacity_coefficients = []
        for i in range(0,len(cs)):
            fugacity_coefficients.append( calculate_phi(i,T))

        return fugacity_coefficients, kijs

class RK ():
    
    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """Redlich-Kwong equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            a = 0.427480* (R**2) * (Tc**2.5) /Pc
            return a

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.086640*R*Tc)/Pc 
            return b

        kijs = {}
    
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None:
               
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        if i==j:
                            kijs[(i,j)] = 0
                        else:
                            if kij_tune.get((i,j),None)!=None:
                                kijs[(i,j)] = calculate_kij(cs[i],cs[j],kij_tune[(i,j)] )
                            else:
                                kijs[(i,j)] = kijs[(j,i)]
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if i==j:
                        kijs[(i,j)] = 0
                    else:
                        if kij_input.get((i,j),None):
                            if abs(kij_input.get((i,j))) < 0.3:
                                kijs[(i,j)] = kij_input[(i,j)]
                            else:
                                kijs[(i,j)] = 0
                        else:
                            kijs[(i,j)] = kijs[(j,i)]

        def calculate_amix(y):
            """a value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)]
                    ai = calculate_a(cs[i]) #ai value
                    aj = calculate_a(cs[j]) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2.5) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            root = np_roots(coefficients)
            return max(root)# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp)
            b = calculate_b(comp)
            Ai = calculate_A(a,T)
            Bi = calculate_B(b,T)
            
            phi = Bi/B*(Z-1) - log(Z-B)+ A/B*(Bi/B - 2*(Ai/A)**0.5)*log(1+B/Z)
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T

        phi = []
        for i in range(len(cs)):
            phi.append(calculate_phi(i,T))
        
        return phi, kijs

class SRK():
    
    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """Soave-Redlich-Kwong equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.48 + 1.574*w - 0.176*w**2 #SRK kappa value
            c = 0.42747*(R**2)*(Tc**2)/Pc #SRK multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #SRK alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.08664*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None:
               
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        if i==j:
                            kijs[(i,j)] = 0
                        else:
                            if kij_tune.get((i,j),None)!=None:
                                kijs[(i,j)] = calculate_kij(cs[i],cs[j],kij_tune[(i,j)] )
                            else:
                                kijs[(i,j)] = kijs[(j,i)]
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if i==j:
                        kijs[(i,j)] = 0
                    else:
                        if kij_input.get((i,j),None):
                            if abs(kij_input.get((i,j))) < 0.3:
                                kijs[(i,j)] = kij_input[(i,j)]
                            else:
                                kijs[(i,j)] = 0
                        else:
                            kijs[(i,j)] = kijs[(j,i)]
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/B*(2*ak/amix - b/bmix)*log((Z+B)/Z)
            return exp(phi)
        fug_phi = []
        for i in range(0,len(cs)):
            fug_phi.append( calculate_phi(i,T) )
        
        return fug_phi, kijs

    def phi_liquid(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """Soave-Redlich-Kwong equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.48 + 1.574*w - 0.176*w**2 #SRK kappa value
            c = 0.42747*(R**2)*(Tc**2)/Pc #SRK multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #SRK alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.08664*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None:
               
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        if i==j:
                            kijs[(i,j)] = 0
                        else:
                            if kij_tune.get((i,j),None)!=None:
                                kijs[(i,j)] = calculate_kij(cs[i],cs[j],kij_tune[(i,j)] )
                            else:
                                kijs[(i,j)] = kijs[(j,i)]
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if i==j:
                        kijs[(i,j)] = 0
                    else:
                        if kij_input.get((i,j),None):
                            if abs(kij_input.get((i,j))) < 0.3:
                                kijs[(i,j)] = kij_input[(i,j)]
                            else:
                                kijs[(i,j)] = 0
                        else:
                            kijs[(i,j)] = kijs[(j,i)]
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            roots = np_roots(coefficients)
            for root in roots:
                if root > 0 and root < max(roots):
                    min_root = root         
            return min_root # Return smallest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
                
            phi = b*(Z-1)/bmix - log(Z-B) - A/B*(2*ak/amix - b/bmix)*log((Z+B)/Z)
            return exp(phi)
        fug_phi = []
        for i in range(0,len(cs)):
            fug_phi.append( calculate_phi(i,T) )
        return fug_phi, kijs

class Ideal():
    """Ideal property method"""
    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):

        phi = []; kijs = {"Parameters":"No interaction parameters for Ideal system."}
        for i in range(0, len(components)):
            phi.append(1)
        return phi, kijs

    def gamma(components, temp, fractions):

        gammas = []; kijs = {"Parameters":"No interaction parameters for Ideal system."}
        for i in range(0, len(components)):
            gammas.append(1)
        return gammas, kijs

class Uniquac():
    """UNIQUAC model based activity coefficient calculations."""

    def gamma(components,temperature,fractions):
        
        cs = components 
        T = temperature
        x = fractions
        for item in x:
            if item == 0:
                item = 1E-05
              
        r = []; q = []; qp = []
        for k in range(0,len(cs)):
            r.append( float(cs[k].UniquacR) )
            q.append( float(cs[k].UniquacQ) )
            qp.append( float(cs[k].UniquacQP) )
        
        #---Calculate teta and fi values for each substance
        teta = []; fi = []; tetap = [] #tetap = teta' prime value for gR calc.
        for i in range(0,len(cs)):
            fi_nom = x[i]*r[i] #fi Nominator
            fi_denom = 0  #fi Denominator
            teta_nom = x[i]*q[i]
            teta_denom = 0
            tetap_nom = x[i]*qp[i]
            tetap_denom = 0
            for j in range(0,len(cs)):
                fi_denom += x[j]*r[j]
                teta_denom += x[j]*q[j]
                tetap_denom += x[j]*qp[j]
            fi.append(fi_nom/fi_denom) #Fi value of the i. component
            teta.append(teta_nom/teta_denom)  #teta value of the i. component
            tetap.append(tetap_nom/tetap_denom) #teta' prime value of the i. component
        
        def a_ij(id1, id2):
            file_path = MODELS_URL+"uniquac.txt"
            with open(file_path, 'r') as f:
                isFound = False # Is parameters found?
                for line in f.readlines():
                    aux = line.split(';')
                    if aux[0] == id1 and aux[1] == id2:
                        a12 = aux[2]
                        isFound = True
                    elif aux[0] == id2 and aux[1] == id1:
                        a12 = aux[3]
                        isFound = True
            if isFound:
                return float(a12) #units.mol_enthalpy(float(a12),"CGS","SI") #Convert to kJ/kmol
            else: 
                print('No parameters were found!')

        def tau(i,j):
            """Calculates tau_ij values"""
            if i == j:
                return 1
            else:
                id1 = cs[i].LibraryIndex
                id2 = cs[j].LibraryIndex
                return exp( -a_ij(id1,id2)/(1.9872*T)) #R = 1.9872 cal/mol.K 
        
        taus = {}
        for i in range(0,len(cs)):
            for j in range(0,len(cs)):
                taus[(i,j)] = tau(i,j)

        def unsymmetric():
            l = [0,0]
            l[0] = 5*(r[0]-q[0]) - (r[0]-1)
            l[1] = 5*(r[1]-q[1]) - (r[1]-1)
            C1 = log(fi[0]/x[0]) + 5*q[0]*log(teta[0]/fi[0]) + fi[1]*(l[0]-r[0]*l[1]/r[1]) 
            R1 = qp[0]*log(tetap[0]+tetap[1]*taus[(1,0)]) + tetap[1]*qp[0]*(taus[(1,0)]/(tetap[0]+tetap[1]*taus[(1,0)]) - taus[(0,1)]/(tetap[1]+tetap[0]*taus[(0,1)]))  
            return exp(C1+R1)
        
        def symmetric():
            C = []; R = []
            for i in range(0,len(cs)):
                C.append( 1 + log(fi[i]/x[i]) - fi[i]/x[i] -5*q[i]*( 1+ log(fi[i]/teta[i])- fi[i]/teta[i] ))
                for j in range(0,len(cs)):
                    if i != j :
                        R.append( qp[i]*( 1- log( tetap[j]*taus[(j,i)]+tetap[i] )- tetap[j]*taus[(i,j)]/(tetap[j]+tetap[i]*taus[(i,j)]) - tetap[i]/(tetap[j]*taus[(j,i)]+tetap[i]) ) )            
            return exp(C[0]+R[0]), exp(C[1]+R[1])
        
        return symmetric(), taus
        
class NRTL():
    """NRTL Activity coefficient calculations"""

    def gamma(components, temperature, fractions):

        cs = components
        T = temperature
        x = fractions

        def a_ij(id1, id2):
            file_path = MODELS_URL+"nrtl.txt"
            with open(file_path, 'r') as f:
                isFound = False # Is parameters found?
                for line in f.readlines():
                    aux = line.split(';')
                    if aux[0] == id1 and aux[1] == id2:
                        a12 = aux[2]
                        alfa = aux[4]
                        isFound = True
                    elif aux[0] == id2 and aux[1] == id1:
                        alfa = aux[4]
                        a12 = aux[3]
                        isFound = True
            if isFound:
                return float(a12), float(alfa) #units.mol_enthalpy(float(a12),"CGS","SI") #Convert to kJ/kmol
            else: 
                print('WARNING!: No parameters were found for a_ij! Default parameters were used')
                return 100, 0.5 #Default parameters 

        aij = {}
        for i in range(0,len(cs)):
            for j in range(0,len(cs)):
                if i != j:
                    aij[(i,j)] = a_ij( cs[i].LibraryIndex, cs[j].LibraryIndex )
        
        def tau(i,j):
            """Calculates tau_ij values"""
            if i == j:
                return 0, 1
            else:
                aux_tau = aij[(i,j)]
                return aux_tau[0]/(1.9872*T), aux_tau[1]  #R = 1.9872 cal/mol.K
            
        def G(i,j):
            """Calculates Gij value"""
            if i == j:
                return 1
            else:
                aux_G = tau(i,j)
                return exp(-aux_G[1] * aux_G[0] )

        S = []; C = []
        for i in range (0, len(cs)):
            aux1 = 0; aux2 = 0
            for j in range(0,len(cs)):
                aux1 += x[j]*G(j,i) 
                aux2 += x[j]*G(j,i) * tau(j,i)[0]
            S.append(aux1)
            C.append(aux2)
        
        gamma = []
        for i in range(0,len(cs)):
            aux_k = 0
            for k in range(0,len(cs)):
                aux_k += x[k]*G(i,k)*(tau(i,k)[0] - C[k]/S[k])/S[k]
            gamma.append( exp( C[i]/S[i] + aux_k) )

        return gamma, aij
    
class Dortmund():
    """Modified Unifac Dortmund model"""
    
    def gamma(components,temperature,fractions):
        
        cs = components 
        T = temperature
        x = fractions
        for item in x:
            if item == 0:
                item = 1E-05
        
        # Get Q and R values for groups
        groupi = []; groupk = {}; ip = {}
        file_path = MODELS_URL+"modfac.txt"
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for i in range(0,len(cs)):
                groups = cs[i].ModifiedUnifac
                rk_data = []
                for pair in groups:
                    for line in lines:
                        aux = line.split(';')
                        if aux[3] == str(pair[0]):
                            ip[pair[0]] = int(aux[0]) 
                            if pair[0] in groupk.keys():
                                groupk[pair[0]][0].append((i,pair[1]))
                            else:
                                groupk[pair[0]] = ([(i, pair[1])], float(aux[4]), float(aux[5]))
                            rk_data.append( (pair[0], pair[1], float(aux[4]), float(aux[5])) )
                            break 
                groupi.append(rk_data)               
        #groupk= {17: ([(0, 1)], 0.92, 1.4), 1: ([(1, 1)], 0.9011, 0.848), 2: ([(1, 1)], 0.6744, 0.54), 15: ([(1, 1)], 1.0, 1.2)}
        
        #Calculate r and q values for components
        r = []; q = []
        for i in range(0,len(cs)):
            ri = 0; qi = 0
            for data in groupi[i]:
                ri += data[1]*data[2]
                qi += data[1]*data[3]
            r.append(ri)
            q.append(qi)        
        
        # Calculation of residual and combinatorial parts
        # ln gamma_k = Qk*[ 1-log(sum(tetai*taui,k)) - sum [ (tetai*taui,m)/sum(tetaj*tauj,m)]
        # Calculate activity coefficients for each group
        
        group_names = [] # Get group numbers 
        for key in groupk.keys():
            group_names.append(key)
        
        def X(k):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*x[i]
            
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += x[itm[0]]*itm[1]
            return aux1/aux2
        
        def tau(m,n):
            if m == n:
                return 1
            else:
                file_name = MODELS_URL+"modfac_ip.txt"
                found = False
                m = ip[m]; n = ip[n]
                with open(file_name, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.split()

                        if m == n:
                            aij = 0; bij = 0; cij = 0
                            found = True
                            break
                        elif int(line[0]) == m and int(line[1]) == n:
                            aij = float(line[2])
                            found = True
                            bij = float(line[3])
                            cij = float(line[4])
                            break
                        elif int(line[0]) == n and int(line[1]) == m:
                            aij = float(line[5])
                            found = True
                            bij = float(line[6])
                            cij = float(line[7])
                            break
                if found:
                    return exp( -(aij + bij*T + cij*T**2)/T)
                else:
                    print("WARNING! No MODFAC interaction parameters were found for groups",m,n)
                    return exp(-50/T) #default value
        
        taus = {}
        for m in group_names:
            for n in group_names:
                taus[(m,n)] = tau(m,n)

        Xk = [] #Calculate and store Xk values
        for k in group_names:
            Xk.append(X(k))
        
        Xi = [] #Calculate and store Xk values for pure components
        def X2(k, xi):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*xi[i]
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += xi[itm[0]]*itm[1]
            return aux1/aux2

        def teta(k):
            """Teta value for group m"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xk[nk]
            tet = groupk[k][2]*Xk[kk]/aux
            return tet

        for i in range(0,len(cs)):
            ki = []
            for k in group_names:# TODO loop sırasını değiştir i dışa k içe
                xi = x.copy()
                for j in range(0,len(xi)): 
                    if i==j:
                        xi[j] = 1
                    else:
                        xi[j] = 0
                ki.append( X2(k,xi) ) 
            Xi.append(ki)
        
        def tetai(k, i):
            """Teta value for group m in pure component"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xi[i][nk]
            teti = groupk[k][2]*Xi[i][kk]/aux
            return teti
        
        teta_k = []; teta_ki = []
        for i in range(0,len(cs)):
            pure_k = []
            for k in group_names:
                pure_k.append(tetai(k,i))
            teta_ki.append(pure_k)
        
        for k in group_names:
                teta_k.append(teta(k))
        
        activity_R = [] #Residual part for activity coefficient ln gammaR
        for i in range(0,len(cs)):
            ln_gamma_R = 0
            
            for k in group_names:
                vk = 0
                for t in groupk[k][0]:
                    if t[0] == i:
                        vk = t[1]
                
                Qk = groupk[k][2]
                kk = group_names.index(k)
                nom = 0; aux = 0
                nom_i = 0; aux_i = 0
                
                for m in group_names:
                    denom_i = 0; denom = 0
                    mm = group_names.index(m)
                    for n in group_names:
                        nn = group_names.index(n)
                        denom += teta_k[nn]*taus[(n,m)]
                        denom_i += teta_ki[i][nn]*taus[(n,m)]
                        
                    nom += teta_k[mm]*taus[(k,m)]/denom
                    aux += teta_k[mm]*taus[(m,k)]
                    nom_i += teta_ki[i][mm]*taus[(k,m)]/denom_i
                    aux_i += teta_ki[i][mm]*taus[(m,k)]
                
                ln_gamma_k = Qk*(1- log(aux) - nom )
                ln_gamma_ki = Qk*(1- log(aux_i) - nom_i )
                ln_gamma_R += vk*(ln_gamma_k - ln_gamma_ki)
                #print(k, ln_gamma_k)
            activity_R.append(ln_gamma_R)
        
        activity_C = []
        #Gamma combinatorial for components
        V = []; F = []; Vp = [] # V' modified dortmund
        for i in range (0,len(cs)):
            aux_r = 0; aux_q = 0; aux_rp = 0
            for j in range(0,len(cs)):
                aux_r += r[j]*x[j]
                aux_rp += (r[j]**0.75)*x[j]
                aux_q += q[j]*x[j]
            V.append(r[i]/aux_r)
            Vp.append((r[i]**0.75)/aux_rp)
            F.append(q[i]/aux_q)
        
        for i in range(0, len(cs)):
            aux = 1 - Vp[i]+ log(Vp[i]) - 5*q[i]*( 1- V[i]/F[i]+ log(V[i]/F[i]) )
            activity_C.append(aux)
        
        activity_coefficients = []
        for i in range(0,len(cs)):
            activity_coefficients.append( exp(activity_C[i] + activity_R[i]) )
        
        return activity_coefficients, taus

class Unifac():
    """Unifac model activity coefficient"""
    def gamma(components,temperature,fractions):
        
        cs = components 
        T = temperature
        x = fractions
        for item in x:
            if item == 0:
                item = 1E-05
        
        # Get Q and R values for groups
        groupi = []; groupk = {}; ip = {}
        file_path = MODELS_URL+"unifac.txt"
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for i in range(0,len(cs)):
                groups = cs[i].UnifacVLE
                rk_data = []
                for pair in groups:
                    for line in lines:
                        aux = line.split(',')
                        if aux[1] == str(pair[0]):
                            ip[pair[0]] = int(aux[0]) 
                            if pair[0] in groupk.keys():
                                groupk[pair[0]][0].append((i,pair[1]))
                            else:
                                groupk[pair[0]] = ([(i, pair[1])], float(aux[4]), float(aux[5]))
                            rk_data.append( (pair[0], pair[1], float(aux[4]), float(aux[5])) )
                            break 
                groupi.append(rk_data)               
        #groupk= {17: ([(0, 1)], 0.92, 1.4), 1: ([(1, 1)], 0.9011, 0.848), 2: ([(1, 1)], 0.6744, 0.54), 15: ([(1, 1)], 1.0, 1.2)}
        
        #Calculate r and q values for components
        r = []; q = []
        for i in range(0,len(cs)):
            ri = 0; qi = 0
            for data in groupi[i]:
                ri += data[1]*data[2]
                qi += data[1]*data[3]
            r.append(ri)
            q.append(qi)        
       
        # Calculation of residual and combinatorial parts
        # ln gamma_k = Qk*[ 1-log(sum(tetai*taui,k)) - sum [ (tetai*taui,m)/sum(tetaj*tauj,m)]
        # Calculate activity coefficients for each group
        
        group_names = [] # Get group numbers 
        for key in groupk.keys():
            group_names.append(key)
        
        def X(k):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*x[i]
            
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += x[itm[0]]*itm[1]
            return aux1/aux2
        
        def tau(m,n):
            if m == n:
                return 1
            else:
                file_name = MODELS_URL+"unifac_ip.txt"
                found = False
                m = ip[m]; n = ip[n]
                with open(file_name, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.split("\t")
                        if int(line[0]) == m and int(line[2]) == n:
                            aij = float(line[4]); found = True
                        elif int(line[0]) == m and int(line[2]) == n:
                            aij = float(line[5]); found = True
                if found:
                    return exp(-aij/T)
                else:
                    print("WARNING! No UNIFAC interaction parameters were found for groups",m,n)
                    return exp(-50/T) #default value
        
        taus = {}
        for m in group_names:
            for n in group_names:
                taus[(m,n)] = tau(m,n)

        Xk = [] #Calculate and store Xk values
        for k in group_names:
            Xk.append(X(k))
        
        Xi = [] #Calculate and store Xk values for pure components
        def X2(k, xi):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*xi[i]
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += xi[itm[0]]*itm[1]
            return aux1/aux2

        def teta(k):
            """Teta value for group m"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xk[nk]
            tet = groupk[k][2]*Xk[kk]/aux
            return tet

        for i in range(0,len(cs)):
            ki = []
            for k in group_names:
                xi = x.copy()
                for j in range(0,len(xi)): 
                    if i==j:
                        xi[j] = 1
                    else:
                        xi[j] = 0
                ki.append( X2(k,xi) ) 
            Xi.append(ki)

        def tetai(k, i):
            """Teta value for group m in pure component"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xi[i][nk]
            teti = groupk[k][2]*Xi[i][kk]/aux
            return teti
        
        teta_k = []; teta_ki = []
        for i in range(0,len(cs)):
            pure_k = []
            for k in group_names:
                pure_k.append(tetai(k,i))
            teta_ki.append(pure_k)
        
        for k in group_names:
                teta_k.append(teta(k))
        
        activity_R = [] #Residual part for activity coefficient ln gammaR
        for i in range(0,len(cs)):
            ln_gamma_R = 0
            
            for k in group_names:
                vk = 0
                for t in groupk[k][0]:
                    if t[0] == i:
                        vk = t[1]
                
                Qk = groupk[k][2]
                kk = group_names.index(k)
                nom = 0; aux = 0
                nom_i = 0; aux_i = 0
                
                for m in group_names:
                    denom_i = 0; denom = 0
                    mm = group_names.index(m)
                    for n in group_names:
                        nn = group_names.index(n)
                        denom += teta_k[nn]*taus[(n,m)]
                        denom_i += teta_ki[i][nn]*taus[(n,m)]
                        
                    nom += teta_k[mm]*taus[(k,m)]/denom
                    aux += teta_k[mm]*taus[(m,k)]
                    nom_i += teta_ki[i][mm]*taus[(k,m)]/denom_i
                    aux_i += teta_ki[i][mm]*taus[(m,k)]
                
                ln_gamma_k = Qk*(1- log(aux) - nom )
                ln_gamma_ki = Qk*(1- log(aux_i) - nom_i )
                ln_gamma_R += vk*(ln_gamma_k - ln_gamma_ki)
                
            activity_R.append(ln_gamma_R)
        
        activity_C = []
        #Gamma combinatorial for components
        V = []; F = []
        for i in range (0,len(cs)):
            aux_r = 0; aux_q = 0
            for j in range(0,len(cs)):
                aux_r += r[j]*x[j]
                aux_q += q[j]*x[j]
            V.append(r[i]/aux_r)
            F.append(q[i]/aux_q)
        
        for i in range(0, len(cs)):
            aux = 1 - V[i]+ log(V[i]) - 5*q[i]*( 1- V[i]/F[i]+ log(V[i]/F[i]) )
            activity_C.append(aux)

        activity_coefficients = []
        for i in range(0,len(cs)):
            activity_coefficients.append( exp(activity_C[i] + activity_R[i]) )
        print(taus)
        return activity_coefficients, taus