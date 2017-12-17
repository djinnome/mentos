import pyOpt
    
from scipy.stats import entropy
import numpy as np
import cvxpy as cvx
import pandas as pd
import numpy as np
import cobra, re, os


pd.options.display.float_format = '{:.3g}'.format

def find_equilibrium(met_bounds,efflux_mets, uptake_mets, fullS, mu0, R, T ):
    rxns, mets = fullS.columns, fullS.index
    efflux = [mets.get_loc( e ) for e in efflux_mets ]
    uptake = [mets.get_loc( u ) for u in uptake_mets ]
    mu_uptake = mu0[uptake_mets] + R*T*met_bounds[uptake_mets].apply(np.log)
    log_c_efflux = cvx.Variable(len(efflux))
    mu_efflux = log_c_efflux*R*T + mu0[efflux_mets].values
    p = cvx.Problem( cvx.Minimize(cvx.norm2( mu_efflux )),
                     [cvx.sum_entries(fullS.T[efflux_mets].as_matrix()*mu_efflux) == -fullS.T[uptake_mets].dot(mu_uptake).sum()] )
    p.solve()
    Reactant_potential = -fullS.T[uptake_mets].dot(mu_uptake)
    Product_potential=   cvx2a(fullS.T[efflux_mets].as_matrix()*mu_efflux.value)
    return pd.DataFrame(dict(log_c_efflux=cvx2a(log_c_efflux.value),
                             mu_efflux=cvx2a(mu_efflux.value), \
                             mu0_efflux=mu0[efflux_mets]), index=efflux_mets), \
          pd.DataFrame(dict(Reactant_potential=Reactant_potential,
                             Product_potential=Product_potential,
                             deltaG=Product_potential - Reactant_potential), index=rxns), \
           pd.DataFrame(dict(log_c_uptake=met_bounds[uptake_mets].apply(np.log), \
                             mu_uptake = mu_uptake,
                             mu0_uptake=mu0[uptake_mets]),index=uptake_mets)

def find_equilibrium2(log_met_bounds, fullS, mu0, R, T, efflux_mets, uptake_mets ):
    """
    Find the equilibrium given the fixed metabolites in met_bounds
    """
    rxns, mets = fullS.columns, fullS.index
    log_c = cvx.Variable(len(mets))
    mu = log_c*R*T + mu0.values
    deltaG = fullS.T.as_matrix()*mu
    efflux_idx = [mets.get_loc(e) for e in efflux_mets]
    uptake_idx = [mets.get_loc(u) for u in uptake_mets]

    
    p = cvx.Problem( cvx.Minimize(cvx.norm2( deltaG )),
                     [log_c[mets.get_loc(met)] == log_met_bounds[met] for met in log_met_bounds.index])
    p.solve(verbose=True)
    mu_efflux = pd.Series(cvx2a(mu[efflux_idx].value),index=efflux_mets)
    mu_uptake = pd.Series(cvx2a(mu[uptake_idx].value),index=uptake_mets)
    G_products = fullS.T[efflux_mets].dot(mu_efflux)
    G_reactants = fullS.T[uptake_mets].dot(mu_uptake)
    return pd.DataFrame({'c': np.exp(cvx2a(log_c.value)),'$\log c$': cvx2a(log_c.value), '$RT\log c$':R*T*cvx2a(log_c.value), '$\mu^0$': mu0, '$\mu$':cvx2a(mu.value)},index=mets), pd.DataFrame({'$\Delta G$':cvx2a(deltaG.value), '$G_{products}$': G_products , '$G_{reactants}$': G_reactants},index=rxns), fullS.T.dot(cvx2a(mu.value))

def find_equilibrium3(log_met_bounds, fullS, mu0, R, T, efflux_mets, uptake_mets, ref_mols=1e-3 ):
    """
    Find the equilibrium given the fixed metabolites in met_bounds
    """
    rxns, mets = fullS.columns, fullS.index
    log_c = cvx.Variable(len(mets))
    #r = cvx.Variable(len(rxns))
    mu = log_c*R*T + mu0.values
    efflux_idx = [mets.get_loc(e) for e in efflux_mets]
    uptake_idx = [mets.get_loc(u) for u in uptake_mets]

    deltaG = fullS.T.as_matrix()*mu
    
    min_ref_mols_obj =  cvx.Minimize(cvx.norm2(log_c -np.log(ref_mols)))
    metabolite_constraint = [log_c[mets.get_loc(met)] == log_met_bounds[met] for met in log_met_bounds.index] 
    equilibrium_constraint = [deltaG == 0]
    p = cvx.Problem(min_ref_mols_obj,   equilibrium_constraint + metabolite_constraint) # + equilibrium_constraint)
    p.solve()
    
    #p = cvx.Problem( cvx.Minimize(cvx.norm2( deltaG )),
    #                 metabolite_constraint)
    #p.solve(verbose=True)
    
    
    mu_efflux = pd.Series(cvx2a(mu[efflux_idx].value),index=efflux_mets)
    mu_uptake = pd.Series(cvx2a(mu[uptake_idx].value),index=uptake_mets)

    G_products = fullS.T[efflux_mets].dot(mu_efflux)
    G_reactants = fullS.T[uptake_mets].dot(mu_uptake)
       
   
    potential=pd.DataFrame({'c': np.exp(cvx2a(log_c.value)),
                            '$\log c$': cvx2a(log_c.value), 
                            '$RT\log c$':R*T*cvx2a(log_c.value), 
                            '$\mu^0$': mu0, 
                            '$\mu$':cvx2a(mu.value)},index=mets) 
    deltaG = pd.DataFrame({'$\Delta G$':cvx2a(deltaG.value), 
                            '$G_{products}$': G_products , 
                            '$G_{reactants}$': G_reactants},index=rxns)
    
    total_reactant_potential = (fullS.T[uptake_mets]*cvx2a(mu[uptake_idx].value)).T
    total_reactant_potential['stoichiometry'] = fullS.T[uptake_mets].squeeze()
    total_reactant_potential['$\mu$'] = cvx2a(mu[uptake_idx].value)
    total_product_potential = (fullS.T[efflux_mets]*cvx2a(mu[efflux_idx].value)).T
    total_product_potential['stoichiometry'] = fullS.T[efflux_mets].squeeze()
    total_product_potential['$\mu$'] = cvx2a(mu[efflux_idx].value)
    return dict(potential= potential, 
                deltaG=deltaG, 
                total_reactant_potential=total_reactant_potential,
                total_product_potential = total_product_potential,
                metabolite_constraint=pd.DataFrame({'metabolite_constraint':[m.dual_value 
                                                            for m in metabolite_constraint]}
                                        ,index=log_met_bounds.index), 
                equilibrium_constraint = pd.DataFrame({'equilibrium_constraint':[e.dual_value 
                                                            for e in equilibrium_constraint]},
                                                        index=rxns), 
                status=p.status)
                                 

def cvx2a( matrix ):
    """
    cvx variable or value to ndarray
    """
    if type(matrix) is cvx.expressions.variables.variable.Variable:
        return np.squeeze(np.asarray(matrix.value))
    elif type(matrix) is np.matrixlib.defmatrix.matrix:
        return np.squeeze(np.asarray(matrix))
    else: # if array-like
        return np.squeeze(matrix)

def get_rates_from_log_likelihood_and_flux( log_likelihood, v):
    """
From the definition of forward likelihood and net flux:

.. math::
    \begin{eqnarray}
    L_+ &=& \frac{r_{+}}{r_{-}} \\
    v &=& r_{+} - r_{-} \\
    \end{eqnarray}$$

We can solve for the forward reaction rates :math:`r_{+}$ and $r_{-}`:

.. math::
    \begin{eqnarray}
    r_{-} & = & \frac{r_{+}}{L_+} \\
    r_{+}  - r_{-} & = & v\\
    r_{+} - \frac{r_{+}}{L_+} & = & v \\
    r_{+}(1 -\frac{1}{L_+}) & = & v \\
    r_{+} & = & \frac{v}{1-\frac{1}{L_+}} \\
       & = & \frac{L_+v}{L_+ - 1} \\
    r_{-} & = & \frac{v}{L_+ - 1}  \\
    \end{eqnarray}

"""
    net_flux = cvx2a(v)
    forward_likelihood = cvx2a(log_likelihood )
    forward_likelihood = np.exp(forward_likelihood)
    forward_rate = net_flux*forward_likelihood/(forward_likelihood - 1)
    backward_rate = net_flux/(forward_likelihood - 1)
    return forward_rate, backward_rate

def getMetBounds( S,fluxes, c_L, c_U, ext=r'_ext?'):
    """Set external uptake metabolites to upper bound concentration.
    Set external eflux metabolites to lower bound concentration
    if flux * stoichiometric coefficient is negative, then it is an uptake metabolite
    if flux * stoichiometric coefficient is positive, then it is an efflux metabolite"""
    external_mets = S.loc[S.index.str.contains(ext)].dot(fluxes).squeeze()
    met_bounds = pd.Series(index=external_mets.index)

    met_bounds[external_mets < 0] = c_U
    met_bounds[external_mets >  0] = c_L
    return met_bounds

def predict_log_likelihoods(S, R, T, mu0, uptake_met_concs, excreted_met_concs ):
    m,n = S.shape
    rxns = S.columns
    mets = S.index
    log_c = cvx.Variable(m)
    log_likelihood = cvx.Variable(n)
    deltaG0 = S.T.dot(mu0)
    log_Q = S.as_matrix().T*log_c
    log_K = cvx.Constant(-1.0/(R*T))*deltaG0
    mu = R*T*log_c + mu0.values

    uptake_met_indices = S.index.get_loc(uptake_met_concs.index)
    excreted_met_indices = S.index.get_loc(excreted_met_concs.index)
    
    obj = cvx.Maximize(cvx.sum_entries( cvx.entr( log_likelihood )))
    constraints = [log_c[uptake_met_indices] == uptake_met_concs.values,
                   log_c[excreted_met_indices] == excreted_met_concs.values,

                   log_likelihood == log_K - log_Q,

                   mu >= np.sum(mu[excreted_met_indices]),
                   mu <= np.sum(mu[uptake_met_indicies])
                   ]
    prob = cvx.Problem(obj, constraints)
    prob.solve(verbose=True)
    return pd.DataFrame(log_c.value, index=mets), pd.DataFrame(log_likelihood.value, index=rxns)

def get_rxn_bounds_from_log_likelihood( log_likelihood ):
    """
.. math::
    \begin{eqnarray}
    v_{lower}(\log L) & = & \min\left(-1,\text{sign}(\log L)\cdot L^{\text{sign}(\log L)}\right) + 1 \\
    v_{upper}(\log L) & = &  \max\left(1,\text{sign}(\log L)\cdot L^{\text{sign}(\log L)}\right) - 1
    \end{eqnarray}

    """
    log_likelihood = cvx2a(log_likelihood)
    n = len(log_likelihood)
    thermodynamic_driving_force =  np.sign(log_likelihood)*np.power(np.exp(log_likelihood), np.sign(log_likelihood))
    return np.minimum(thermodynamic_driving_force, -1) + 1, np.maximum(thermodynamic_driving_force, 1) - 1
 


def predict_fluxes( S, rxn_bounds, biomass ):
    m,n = S.shape
    v = cvx.Variable(n)
    obj = cvx.Maximize( v[S.columns.get_loc(biomass)] )
    constraints = [S.as_matrix()*v == 0,
                   rxn_bounds['lower'].values <= v,
                   v <= rxn_bounds['upper'].values ]
    prob = cvx.Problem(obj, constraints)
    prob.solve( verbose=True)
    return pd.DataFrame(v.value, index=S.columns)
def make_variables( x, fullS, mu0, deltaG0,R = 8.3144598/1000.0,  T = 298.15 ): # ideal gas constant in cals/mol ):
    m,n = fullS.shape
    log_c = x[:m]
    forward_rate = np.abs(x[m:m+n])
    backward_rate = np.abs(x[m+n:m+2*n])
    log_Q = np.dot(fullS.T,log_c)   # log of the Reaction quotient
    log_K = -1.0/(R*T)*deltaG0.as_matrix() 
    forward_likelihood = forward_rate/backward_rate
    backward_likelihood = backward_rate/forward_rate
    total_likelihood = forward_likelihood.sum()+ backward_likelihood.sum()
    forward_probability = forward_likelihood/total_likelihood
    backward_probability= backward_likelihood/total_likelihood
    sign = np.sign(np.log(forward_likelihood))
    thermodynamic_driving_force = sign*np.power(forward_likelihood, sign)
    net_flux = forward_rate - backward_rate
    mu = mu0 + R*T*log_c
    return log_c, forward_rate, backward_rate, log_Q, log_K, forward_likelihood, backward_likelihood, \
    forward_probability, backward_probability, mu, thermodynamic_driving_force, net_flux



def generate_metabolite_report( log_c, forward_rate, backward_rate, S, metabolites, internal_mets, rxns, fullS, mu0,     T = 298.15,     V = 1e-15,     R = 8.3144598/1000.0  ) :
    n_A = 6.022e23       # Avogadros number
    forward_likelihood = forward_rate/backward_rate
    backward_likelihood = backward_rate/forward_rate
    forward_probabilities = pd.DataFrame(forward_likelihood/np.sum(forward_likelihood + backward_likelihood), index=rxns)
    backward_probabilities = pd.DataFrame(backward_likelihood/np.sum(forward_likelihood + backward_likelihood), index=rxns)
    net_probabilities = forward_probabilities - backward_probabilities
    net_flux = forward_rate - backward_rate
    mets = pd.DataFrame(n_A*V*np.exp(log_c), index=metabolites, columns=['Counts'], dtype=int)
    mets['Concentrations'] = pd.DataFrame(np.exp(log_c), index=metabolites)
    mets['Log Concentrations'] = pd.DataFrame(log_c, index=metabolites)
    mets['Standard Chemical potential'] = mu0
    mets['Chemical potential'] = mu0 + R*T*log_c
    mets['Normalized Concentrations'] = pd.DataFrame(mets['Concentrations']/mets['Concentrations'].sum(),index=metabolites)
    mets['Absolute activities'] = np.exp(mets['Chemical potential'])
    mets['Normalized absolute activities'] = mets['Absolute activities']/mets['Absolute activities'].sum()
    mets['Entropy'] = -mets['Normalized absolute activities']*mets['Normalized absolute activities'].apply(np.log)
    mets['Entropy of Concentrations'] = -mets['Normalized Concentrations']*mets['Normalized Concentrations'].apply(np.log)
    mets['fullS*forward_rate'] = pd.DataFrame(np.dot(S,forward_rate), index= internal_mets)
    mets['fullS*backward_rate'] = pd.DataFrame(np.dot(S,backward_rate), index=internal_mets)
    mets['fullS*net_flux'] = pd.DataFrame(np.dot(S,net_flux), index=internal_mets)
    mets['fullS*net_likelihood'] = pd.DataFrame(np.dot(S,forward_likelihood - backward_likelihood), index=internal_mets)
    mets['fullS*net_probabilities'] = fullS.dot(net_probabilities)
    #mets['Steady state constraints'] = pd.DataFrame(constraints[-1].dual_value, index=internal_metabolites)
    return mets.astype(np.float64)

def print_report( report_dir, out_template, df ):
    nonalphanumRE = re.compile(r'[^A-Za-z0-9_]')
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
    for c in df.columns:
        df[c].to_csv(os.path.join(report_dir,out_template.format(nonalphanumRE.sub('_', c))),header=True)
    df.to_csv(os.path.join(report_dir, out_template.format(os.path.basename(report_dir))))
def generate_rxn_report(metabolites, log_c, log_Q, log_K,forward_rate, backward_rate, rxns, deltaG0, biomass_rxn, T=298.15, V=1e-15,     R = 8.3144598/1000.0):
 # ideal gas constant
    n_A = 6.022e23       # Avogadros number

    forward_likelihood = forward_rate/backward_rate
    backward_likelihood = backward_rate/forward_rate
    df = pd.DataFrame(forward_likelihood,index=rxns, columns=['Forward likelihoods'])
    df['Backward likelihoods'] = pd.DataFrame(backward_likelihood, index=rxns) #/np.log(df['Rxn likelihoods']).sum()
    df['Delta G'] = pd.DataFrame((log_K - log_Q)*(-R*T), index=rxns)
    #df['Q_r'] = pd.DataFrame(np.exp(log_Q.value),index=rxns)
    #df['K_eq'] = pd.DataFrame(np.exp(log_K.value),index=rxns)
    df['log Q_r'] = pd.DataFrame(log_Q, index=rxns)
    df['log K_eq'] = pd.DataFrame(log_K, index=rxns)
    df['Delta G0'] = pd.DataFrame(deltaG0, index=rxns)
    df['logK - logQ'] = log_K - log_Q #/np.log(df['Rxn likelihoods']).sum()
    df['logQ - logK'] = log_Q - log_K #/np.log(df['Rxn likelihoods']).sum()
    df['Forward Log likelihoods'] = np.log(forward_likelihood) #/np.log(df['Rxn likelihoods']).sum()
    df['Backward Log likelihoods'] = np.log(backward_likelihood)#/np.log(df['Rxn likelihoods']).sum()

    df['Forward probabilities'] = pd.DataFrame(forward_likelihood/np.sum(forward_likelihood + backward_likelihood), index=rxns)
    df['Backward probabilities'] = pd.DataFrame(backward_likelihood/np.sum(forward_likelihood + backward_likelihood), index=rxns)
    df['Net probabilities'] = df['Forward probabilities'] - df['Backward probabilities']
    df['Net likelihoods'] = df['Forward likelihoods'] - df['Backward likelihoods']
    df['Forward rate'] = pd.Series(forward_rate, index=rxns)
    df['Backward rate'] = pd.Series( backward_rate, index=rxns)
    sgn =  np.sign(np.log(forward_likelihood))
    df['Thermodynamic driving force'] = -df['Delta G']
    df['Net flux'] = df['Forward rate'] - df['Backward rate']
    #df['Forward likelihood constraints'] = pd.DataFrame(constraints[0].dual_value, index=rxns)
    #df['Backward likelihood constraints'] = pd.DataFrame(constraints[1].dual_value, index=rxns)
    
    df['Forward Entropy Production'] = -df['Forward probabilities']*df['Forward probabilities'].apply(np.log)
    df['Backward Entropy Production'] = -df['Backward probabilities']*df['Backward probabilities'].apply(np.log)
    df['Reaction Entropy Production'] = df['Forward Entropy Production'] + df['Backward Entropy Production']
    df['Net Reaction Entropy Production'] = df['Forward Entropy Production'] - df['Backward Entropy Production']
    df['Microscopic Reaction Entropy Production Rate'] = df['Forward Entropy Production']*df['Forward rate'] + df['Backward Entropy Production']*df['Backward rate']
    df['Microscopic Reaction Entropy Production Net Flux'] = df['Forward Entropy Production']*df['Net flux'] + df['Backward Entropy Production']*df['Net flux']
    df['Microscopic Reaction Entropy Production Net Flux Difference'] = df['Forward Entropy Production']*df['Net flux'] - df['Backward Entropy Production']*df['Net flux']
    df['Total Entropy Production'] = df['Reaction Entropy Production'].sum()
    df['Total Microscopic Entropy Production Rate'] = df['Microscopic Reaction Entropy Production Net Flux'].sum()
    df['Total Microscopic Entropy Production Rate'] = df['Microscopic Reaction Entropy Production Net Flux'].sum()
    df['Macroscopic Reaction Entropy Production Net Flux'] = df['Forward Entropy Production']*df.loc[biomass_rxn,'Net flux'] + df['Backward Entropy Production']*df.loc[biomass_rxn,'Net flux']
    df['Macroscopic Reaction Entropy Production Net Flux Difference'] = df['Forward Entropy Production']*df.loc[biomass_rxn,'Net flux'] - df['Backward Entropy Production']*df.loc[biomass_rxn,'Net flux']
    df['Total Macroscopic Entropy Production Net Flux'] = df['Macroscopic Reaction Entropy Production Net Flux'].sum()
    df['Total Macroscopic Entropy Production Net Flux Difference'] = df['Macroscopic Reaction Entropy Production Net Flux Difference'].sum()

    
    return df


def run_mentos(fullS, S, internal_mets, deltaG0, mu0, c_L, c_U, v_L, v_U,initial_log_c,initial_forward_rate, initial_backward_rate, obj, biomass_rxn, R = 8.3144598/1000.0,T = 298.15):
    """
    run_mentos computes the nonconvex optimization using pyOpt.  
.. math::
    \underset{x}{\min} & f(x) & 

     & g_j(x) = 0, & j = 1, \ldots, m_e  

     & g_j(x) \leq 0, & j = m_e + 1, \ldots, m  

     & x_{iL} \leq x_i \leq x_{iU}, & i = 1,\ldots, n 


where:

    * x is the vector of design variables
    * f(x) is a nonlinear function
    * g(x) is a linear or nonlinear function
    * n is the number of design variables
    * m_e is the number of equality constraints
    * m is the total number of constraints (number of equality constraints: m_i = m - m_e)

    The arguments to this function are as follows

:param fullS:  The stoichometric matrix data frame object containing all metabolites
:type fullS: pandas.DataFrame with reactions as columns and metabolites as row index
:param S:  The stoichiometric matrix numpy matrix containing only those metabolites in steady state
:type S: numpy.matrix
:param deltaG0: The standard change in Gibbs free energy for all reactions
:type deltaG0: numpy.array
:param mu0: The standard chemical potential numpy array for all metabolites
:type mu0: numpy.array
:param c_L:  The lower bound numpy array on concentrations
:type c_L: numpy.array
:param c_U:  The upper bound numpy array on concentrations
:type c_U: numpy.array
:param v_L:  The lower bound numpy array on fluxes
:type v_L: numpy.array
:param v_U:  The upper bound numpy array on fluxes
:type v_U: numpy.array
:param initial_log_c: The initial log concentrations (perhaps computed from convex mentos) numpy array
:type initial_log_c: numpy.array
:param initial_forward_rate: The initial forward rate (perhaps computed from convex mentos) numpy array
:type initial_forward_rate: numpy.array
:param initial_backward_rate: The initial backward rate (perhaps computed from convex mentos) numpy array
:type initial_backward_rate: numpy.array
:param obj: The objective function. This is a python function of one parameter x.  Returns the value of the objective function f, the constraint function g, and a fail flag (0 is pass)
:type obj: function of 1 argument that returns the value of f(x), g(x), fail
:param biomass_rxn: The name of the biomass reaction in the fullS dataframe.
:type biomass_rxn: string
:return: reports for metabolite-sized data and reaction-sized data
:rtype:  tuple of pandas.DataFrame
    """
    m,n = fullS.shape
    mets = fullS.index
    rxns = fullS.columns
    i, n = S.shape
    mentos = pyOpt.Optimization(obj.__name__,obj)
    mentos.addObj('f')
    mentos.addVarGroup('log_c', m, type='c', value=initial_log_c, 
                       lower=np.log(c_L), upper=np.log(c_U))
    mentos.addVarGroup('forward_rates', n, type='c', value=initial_forward_rate, lower=0, upper=np.inf)
    mentos.addVarGroup('backward_rates', n, type='c', value=initial_backward_rate, lower=0, upper=np.inf)

    mentos.addConGroup(name='steady_state_flux',ncons=i,type='e')
    mentos.addConGroup(name='thermodynamics', ncons=n, type='e')

    mentos.addConGroup(name='boundary_conditions', ncons=m-i, type='e')
    mentos.addConGroup(name='energy_barrier',ncons=i, type='i')
    mentos.addConGroup(name='energy_sink',ncons=i, type='i')
    mentos.addConGroup(name='flux_lower_bounds',ncons=n, type='i',lower=v_L, upper=v_U)
    mentos.addConGroup(name='flux_upper_bounds',ncons=n, type='i',lower=v_L, upper=v_U)

    opt = pyOpt.PSQP()
    f_star, x_star, message = opt(mentos, sens_type='FD', disp_opts=True, sens_mode='')
    log_c = x_star[:m]
    forward_rate = np.abs(x_star[m:m+n])
    backward_rate = np.abs(x_star[m+n:m+2*n])
    log_Q = np.dot(fullS.T,log_c)   # log of the Reaction quotient
    log_K = -1.0/(R*T)*deltaG0.as_matrix() 
    metab = generate_metabolite_report(log_c, forward_rate, backward_rate, rxns, S, mets, internal_mets, rxns, fullS, mu0 )
    reactions = generate_rxn_report(mets, log_c, log_Q, log_K,forward_rate, 
                                                backward_rate, rxns, deltaG0, biomass_rxn)
    return metab, reactions


def run_mentos_ext_free_energy(fullS, S, internal_mets, deltaG0, mu0,  variable_group, constraint_group,  obj, biomass_rxn, R = 8.3144598/1000.0,T = 298.15):
    """
    run_mentos computes the nonconvex optimization using pyOpt.  
.. math::
    \underset{x}{\min} & f(x) & 

     & g_j(x) = 0, & j = 1, \ldots, m_e  

     & g_j(x) \leq 0, & j = m_e + 1, \ldots, m  

     & x_{iL} \leq x_i \leq x_{iU}, & i = 1,\ldots, n 


where:

    * x is the vector of design variables
    * f(x) is a nonlinear function
    * g(x) is a linear or nonlinear function
    * n is the number of design variables
    * m_e is the number of equality constraints
    * m is the total number of constraints (number of equality constraints: m_i = m - m_e)

    The arguments to this function are as follows

:param fullS:  The stoichometric matrix data frame object containing all metabolites
:type fullS: pandas.DataFrame with reactions as columns and metabolites as row index
:param S:  The stoichiometric matrix numpy matrix containing only those metabolites in steady state
:type S: numpy.matrix
:param deltaG0: The standard change in Gibbs free energy for all reactions
:type deltaG0: numpy.array
:param mu0: The standard chemical potential numpy array for all metabolites
:type mu0: numpy.array
:param variable_group: variable names, variable sizes, and intial values for each variable group
:type variable_group: list of dict
:param constraint_goup: Dictionary of constraint names, constraint sizes, constraint types for each constraint group
:type constraint_group: list of dict
:param obj: The objective function. This is a python function of one parameter x.  Returns the value of the objective function f, the constraint function g, and a fail flag (0 is pass)
:type obj: function of 1 argument that returns the value of f(x), g(x), fail
:param biomass_rxn: The name of the biomass reaction in the fullS dataframe.
:type biomass_rxn: string
:return: reports for metabolite-sized data and reaction-sized data
:rtype:  tuple of pandas.DataFrame
    """
    mets = fullS.index
    rxns = fullS.columns
    m,n = fullS.shape
    i, n = S.shape
    mentos = pyOpt.Optimization(obj.__name__,obj)
    mentos.addObj('f')
    for vg in variable_group:
        mentos.addVarGroup(vg['name'], nvars=vg['nvars'], value=vg['initial_value'], type=vg['type'] )
    for cg in constraint_group:
        mentos.addConGroup(name=cg['name'], ncons=cg['ncons'], type=cg['type'] )
    opt = pyOpt.PSQP()
    f_star, x_star, message = opt(mentos, sens_type='FD', disp_opts=True, sens_mode='')
    log_c = x_star[:m]
    forward_rate = np.abs(x_star[m:m+n])
    backward_rate = np.abs(x_star[m+n:m+2*n])
    log_Q = np.dot(fullS.T,log_c)   # log of the Reaction quotient
    log_K = -1.0/(R*T)*deltaG0.as_matrix() 
    metab = generate_metabolite_report(log_c, forward_rate, backward_rate, S, mets, internal_mets, rxns, fullS, mu0 )
    reactions = generate_rxn_report(mets, log_c, log_Q, log_K,forward_rate, 
                                                backward_rate, rxns, deltaG0, biomass_rxn)
    return metab, reactions

    
def compare_frames(**frames):
    return  pd.Panel.from_dict(frames,orient='minor').swapaxes(1,0).to_frame(filter_observations=False)

def frame_differences( left, right ):
    return left[pd.DataFrame(
        np.isclose(left, right ),
        index=left.index,
        columns=left.columns)].T
