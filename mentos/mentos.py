import pyOpt
from scipy.stats import entropy
import numpy as np
import cvxpy as cvx
import pandas as pd
import numpy as np
import cobra


pd.options.display.float_format = '{:.3g}'.format

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
    n = len(log_likelihood.index)
    zero = np.squeeze(np.zeros((n,1)))
    ll = np.sign(log_likelihood.values)*np.exp(log_likelihood.values)
    return pd.DataFrame({'lower_bound':np.minimum(ll, zero), 'upper_bound':np.maximum(ll, zero)}, index=log_likelihood.index)


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

                      

                   
def make_variables( x, fullS, mu0, deltaG0, R, T ):
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



def generate_metabolite_report( log_c, forward_rate, backward_rate, S, metabolites, internal_mets, mu0,     T = 298.15,     V = 1e-15  ) :
    n_A = 6.022e23       # Avogadros number

    R = 8.3144598/1000.0 # ideal gas constant

    forward_likelihood = forward_rate/backward_rate
    backward_likelihood = backward_rate/forward_rate
    net_flux = forward_rate - backward_rate
    mets = pd.DataFrame(n_A*V*np.exp(log_c), index=metabolites, columns=['Counts'], dtype=int)
    mets['Concentrations'] = pd.DataFrame(np.exp(log_c), index=metabolites)
    mets['Log Concentrations'] = pd.DataFrame(log_c, index=metabolites)
    mets['Standard Chemical potential'] = mu0
    mets['Chemical potential'] = mu0 + R*T*log_c
    mets['S*forward_rate'] = pd.DataFrame(np.dot(S,forward_rate), index= internal_mets)
    mets['S*backward_rate'] = pd.DataFrame(np.dot(S,backward_rate), index=internal_mets)
    mets['S*net_flux'] = pd.DataFrame(np.dot(S,net_flux), index=internal_mets)
    #mets['Steady state constraints'] = pd.DataFrame(constraints[-1].dual_value, index=internal_metabolites)
    return mets

def generate_rxn_report(metabolites, log_c, log_Q, log_K,forward_rate, backward_rate, rxns, deltaG0, biomass_rxn, T=298.15, V=1e-15):
    R = 8.3144598/1000.0 # ideal gas constant
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
    df['Net flux'] = df['Forward rate'] - df['Backward rate']
    #df['Forward likelihood constraints'] = pd.DataFrame(constraints[0].dual_value, index=rxns)
    #df['Backward likelihood constraints'] = pd.DataFrame(constraints[1].dual_value, index=rxns)
    
    df['Forward Entropy Production'] = -df['Forward probabilities']*df['Forward probabilities'].apply(np.log)
    df['Backward Entropy Production'] = -df['Backward probabilities']*df['Backward probabilities'].apply(np.log)
    df['Reaction Entropy Production'] = df['Forward Entropy Production'] + df['Backward Entropy Production']
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
    metab = generate_metabolite_report(log_c, forward_rate, backward_rate, S, mets, internal_mets, mu0 )
    reactions = generate_rxn_report(mets, log_c, log_Q, log_K,forward_rate, 
                                                backward_rate, rxns, deltaG0, biomass_rxn)
    return metab, reactions
