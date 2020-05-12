from .methods import fug_methods as fm

def fugacity_vapor(components, temp, pressure, fractions, method=None, kij_input = None, kij_tune=None):
    """ Equation of state solver for vapor phase.
    :param components: Array that contains chemicals.
    :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
    :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
    """
    methods = {"SRK":fm.SRK.phi_vapor,"Ideal":fm.Ideal.phi_vapor,"PR76":fm.PR76.phi_vapor,"PR78":fm.PR76.phi_vapor,"RK":fm.RK.phi_vapor}
    if method != None:
        result = methods[method](components, temp, pressure, fractions, kij_input, kij_tune)
    else:
        result = methods["PR78"](components, temp, pressure, fractions, kij_input, kij_tune)

    return result

def fugacity_liquid(components, temperature, fractions, method=None, kij_input = None, kij_tune=None):
    
    methods = {"SRK":fm.SRK.phi_liquid,"Ideal":fm.Ideal.gamma, "NRTL":fm.NRTL.gamma,"Uniquac":fm.Uniquac.gamma,"Unifac":fm.Unifac.gamma,"Dortmund":fm.Dortmund.gamma}
    if method != None and method != "SRK":
        gamma = methods[method](components, temperature, fractions)
    else:
        gamma = methods["SRK"](components, temperature, 101325, fractions, kij_input, kij_tune)

    return gamma

def activities(components, temperature, fractions, method=None):
    
    methods = {"Ideal":fm.Ideal.gamma, "NRTL":fm.NRTL.gamma,"Uniquac":fm.Uniquac.gamma,"Unifac":fm.Unifac.gamma,"Dortmund":fm.Dortmund.gamma}
    gamma = methods[method](components, temperature, fractions)

    return gamma