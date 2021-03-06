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
    return pd.DataFrame({'c': np.exp(cvx2a(log_c.value)),'$\log c$': cvx2a(log_c.value), '$RT\log c$':R*T*cvx2a(log_c.value), '$\mu^0$': mu0, '$\mu$':cvx2a(mu.value)},index=mets), pd.DataFrame({'$\Delta G$':cvx2a(deltaG.value), '$G_{products}$': G_products , '$G_{reactants}$': G_reactants},index=rxns), fullS.T*cvx2a(mu.value)

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
def make_variables_from_rates( x, fullS, mu0, deltaG0,R = 8.3144598/1000.0,  T = 298.15 ): # ideal gas constant in cals/mol ):
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

make_variables = make_variables_from_rates

def make_variables_from_likelihoods( x, fullS, mu0, deltaG0,R = 8.3144598/1000.0,  T = 298.15 ): # ideal gas constant in cals/mol ):
    m,n = fullS.shape
    log_c = x[:m]
    forward_likelihood = np.abs(x[m:m+n])
    backward_likelihood = np.reciprocal(forward_likelihood)
    log_Q = np.dot(fullS.T,log_c)   # log of the Reaction quotient
    log_K = -1.0/(R*T)*deltaG0.as_matrix()
    net_likelihood = forward_likelihood - backward_likelihood
    
    forward_rate = forward_likelihood*net_likelihood/(forward_likelihood - 1)
    backward_rate = net_likelihood/(forward_likelihood - 1 )
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
    mets['fullS*forward_rate'] = fullS.dot(forward_rate)
    mets['fullS*backward_rate'] =fullS.dot(backward_rate)
    mets['fullS*net_flux'] = fullS.dot(net_flux)
    mets['fullS*net_likelihood'] = fullS.dot(forward_likelihood - backward_likelihood)
    mets['fullS*net_probabilities'] = fullS.dot(net_probabilities)
    #mets['Steady state constraints'] = pd.DataFrame(constraints[-1].dual_value, index=internal_metabolites)
    return mets.astype(np.float64)

def print_report( report_dir, out_template, df ):
    nonalphanumRE = re.compile(r'[^A-Za-z0-9_]')
    if not os.path.isdir(report_dir):
        os.makedirs(report_dir)
    for c in df.columns:
        df[c].to_csv(os.path.join(report_dir,out_template.format(nonalphanumRE.sub('_', c))),header=True)
    df.to_csv(os.path.join(report_dir, out_template.format(nonalphanumRE.sub('_', os.path.basename(report_dir)))), sep='\t')
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
    df['Reduced rate'] = df['Forward rate']*df['Backward rate']/(df['Forward rate'] + df['Backward rate'])
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


    
def compare_frames(**frames):
    return  pd.Panel.from_dict(frames,orient='minor').swapaxes(1,0).to_frame(filter_observations=False)

def frame_differences( left, right ):
    return left[pd.DataFrame(
        np.isclose(left, right ),
        index=left.index,
        columns=left.columns)].T


def check_constraints( x0, constraints,rxns, mets ):
        internal_mets = [m for m in mets if '_e' not in m]
        for constraint in constraints: 
            constraint_eval = constraint['fun'](x0)
            if constraint['type'] == 'eq':
                if np.allclose(constraint_eval, 0):
                    display(Latex('{}'.format(constraint['fun'].__doc__)))
                    if 'jac' in constraint:
                        display(Latex('{}'.format(constraint['jac'].__doc__)))
                else:
                    display('{} violated'.format(constraint['fun'].__name__))
                    if len(constraint_eval) == len(internal_mets):
                        constraint_s = pd.Series(constraint_eval,index=internal_mets)
                        
                        display(constraint_s[~constraint_s.apply(lambda x: np.isclose(x,0))].to_frame(constraint['fun'].__doc__))
                    elif len(constraint_eval) == len(rxns):
                        display(pd.DataFrame(constraint_eval, index=rxns))
                    
            elif constraint['type'] == 'ineq':
                if np.all(constraint_eval >= 0):
                    display(Latex('{}'.format(constraint['fun'].__doc__)))
                else:
                    display('{} violated'.format(constraint['fun'].__name__))
                    if len(constraint_eval) == len(rxns):
                        display(pd.DataFrame(constraint_eval,index=rxns))
                
            
def make_logc_bounds(mets, met_bounds,epsilon=1, log_c_L = np.log(1e-8), log_c_U=np.log(1e-3)):
        logc_bounds = []
        for met in mets:
            if met in met_bounds:
                logc_bounds.append((np.log(met_bounds[met]), np.log(met_bounds[met])))
            else:
                logc_bounds.append((log_c_L, log_c_U))
        return logc_bounds
    
def make_rate_bounds( rxns, r_L=0, r_U=None):
        return [(r_L, r_U) for rxn in rxns]

def entropy_production_thermo_net_likelihood_ss_obj( x ):
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = entropy(np.concatenate( (forward_probability, backward_probability ) ) ) \
            - 1000*np.sum( np.square( np.dot( S, forward_likelihood - backward_likelihood ))) \
            - 1000*np.sum( np.square( np.log( forward_likelihood ) - log_K + log_Q )) 
        return -f
def entropy_production_obj( x ):
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = entropy(np.concatenate( (forward_probability, backward_probability ) ) )
        return -f
def micro_entropy_production_rate_obj(x):
        """Sum of Elementwise product of entropy and rate for forward and backward probabilities"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = np.sum( entr(forward_probability)*forward_rate)  + np.sum( entr(backward_probability)*backward_rate )
        return -f
    
def micro_entropy_production_net_flux_obj( x ):
        """Sum of Elementwise product of entropy and net flux for forward and backward probabilities"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = np.sum( entr( forward_probability) * net_flux ) + np.sum( entr( backward_probability ) * net_flux )
        return -f
    
def macro_entropy_production_rate_obj(x):
        """Maximize product of  entropy production difference and macroscopic biomass growth rate of the ABC model"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = entropy(np.concatenate( (forward_probability, backward_probability ) ) )*net_flux[biomass] #- 1000*np.sum(np.abs(slack))
        return -f
    
def macro_entropy_production_rate_and_net_likelihood_ss_obj(x):
        """Maximize product of  entropy production difference and macroscopic biomass growth rate of the ABC model"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = entropy(np.concatenate( (forward_probability, backward_probability ) ) )*net_flux[biomass] - 1000*np.sum(np.square(np.dot(S, forward_likelihood - backward_likelihood)))
        return -f
        
def macro_entropy_production_rate_thermo_and_net_likelihood_ss_obj(x):
        """Maximize product of  entropy production difference and macroscopic biomass growth rate of the ABC model"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = entropy( np.concatenate(( forward_probability, backward_probability ) ) )*net_flux[biomass] \
            - 1000*np.sum( np.square( np.dot( S, forward_likelihood - backward_likelihood ))) \
            - 1000*np.sum( np.square( np.log( forward_likelihood ) - log_K + log_Q ))
        return -f
        
def macro_entropy_production_rate_and_thermo_obj(x):
        """Maximize product of  entropy production difference and macroscopic biomass growth rate of the ABC model"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        f = entropy(np.concatenate( (forward_probability, backward_probability ) ) )*net_flux[biomass] - 1000*np.sum(np.square(np.log(forward_likelihood) - log_K + log_Q))
        return -f
def steady_state_net_likelihood_constraint(x):
        """$$S\cdot L_{net} = 0$$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return np.dot(S,forward_likelihood-backward_likelihood)
def steady_state_constraint(x):
        """$$S\cdot v = 0$$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return np.dot(S,net_flux)
def thermodynamic_constraint(x):
        """$$\log L_+ + \log Q - \log K = 0$$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return np.log(forward_likelihood) + log_Q - log_K
    
def thermodynamic_jacobian(x):
        """$$J\left[\log L_+ - \log K + \log Q\right]: {\mathscr R}^{2n+m}\rightarrow {\mathscr R}^{n\times 2n+m}=\left[\text{diag }(L_+^{-1})  & S^T\right] $$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables( x, fullS, mu0, deltaG0, R, T )
        return np.append( np.diag(np.reciprocal(forward_likelihood)),      
                                    fullS.T.as_matrix(), axis=1)
def energy_barrier_constraint( x ):
        """$$ \|\Delta G_{ext}\| - S^T\mu \geq 0$$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return   np.abs(external_free_energy) - np.dot(fullS.T.as_matrix(),mu)
def energy_sink_constraint( x ):
        """$$S^T\mu + \|\Delta G_{ext}\| \geq 0$$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux = make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return  np.dot(fullS.T.as_matrix(),mu) + np.abs(external_free_energy) 
def second_law_constraint( x ):
        """$$\log L_+ \geq 0 \iff v \geq 0$$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux= make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return np.log(forward_likelihood)*net_flux
def flux_upper_constraint( x ):
        """$$ v_U - r_+ + r_- \geq 0 $$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux= make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return v_U - forward_rate + backward_rate
def flux_lower_constraint( x ):
        """$$ r_+ - r_- - v_L \geq 0  $$"""
        log_c, \
        forward_rate, backward_rate, \
        log_Q, log_K, \
        forward_likelihood, backward_likelihood, \
        forward_probability, backward_probability, mu, \
        thermodynamic_driving_force, \
        net_flux= make_variables_from_likelihoods( x, fullS, mu0, deltaG0, R, T )
        return  forward_rate - backward_rate - v_L
def cobra_model_to_mentos( cobra_model, 
                    external_compartment='e', 
                    std_rxn_free_energy='deltaG0', 
                    std_chemical_potential='mu0',
                    boundary_concentration = 'boundary_concentration',
                    c_L = 1e-8,
                    c_U = 1e-3,
                    v_L = 0,
                    v_U = 1000):
    fullS = pd.DataFrame(cobra.util.array.create_stoichiometric_matrix(cobra_model), 
                    index=[m.id for m in cobra_model.metabolites],
                    columns = [r.id for r in cobra_model.reactions]),
                
    rxns, mets = fullS.columns, fullS.index
    external_mets = [met.id for met in cobra_model.metabolites if met.compartment == external_compartment ]
    internal_mets = [met.id for met in cobra_model.metabolites if met.compartment != external_compartment ]
    S = fullS.loc[internal_mets]
    deltaG0 = pd.Series(dict([(rxn.id, rxn.notes[std_rxn_free_energy] )
                             for rxn in cobra_model.reactions 
                                 if std_rxn_free_energy in rxn.notes]), 
                        index=rxns)
    mu0     = pd.Series(dict([(met.id, met.notes[std_chemical_potential] )
                             for met in cobra_model.metabolites
                                 if std_chemical_potential in met.notes]), 
                        index=mets)
    met_bounds = pd.Series(dict([(met.id, met.notes[boundary_concentration]) 
                                     for met in cobra_model.metabolites 
                                      if boundary_condition in met.notes]),
                            index=external_mets)
    biomass_rxn = [rxn.id 
               for rxn in cobra_model.reactions 
               if rxn.objective_coefficient != 0][0]
    biomass = rxns.get_loc(biomass_rxn)
    log_c_L = np.log(c_L)
    log_c_U = np.log(c_U)
    return dict(fullS=fullS, 
                S = S,
                rxns=rxns,
                mets=mets,
                external_mets = external_mets,
                internal_mets = internal_mets,
                met_bounds = met_bounds,
                mu0 = mu0,
                deltaG0 = deltaG0,
                log_c_L = log_c_L,
                log_c_U = log_c_U
                )
