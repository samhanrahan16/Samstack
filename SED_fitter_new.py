import numpy as np
import emcee
import corner
from astropy import constants
import matplotlib.pyplot as plt

_h = constants.h.value
_c = constants.c.value
_kb = constants.k_B.value
_nu0 = _c/100e-6

def bbf(T,nu):
    return 2.3e-11 * (2 * _h * nu**3)/(_c**2) * (1/(np.exp((_h * nu)/(_kb * T))))

def bb(wave,T,beta,A):
    nu = _c/wave
    B = bbf(T,nu)
    flux = 1e26 * A * B * ((nu/_nu0)**beta)

    return flux

def lnlike(theta,wave,flux,error):
    T,beta,log_A = theta
    ll = 0
    diff = flux - bb(wave,T,beta,10**log_A)
    sig = (diff**2)/(2*error**2)
    ll = -0.5 * np.sum(sig)

    return ll

def lnprior(theta):
    T,beta,log_A = theta
    if 3 < T < 100 and 0 < beta < 3 and -3 < log_A < 10:
        return 0.0
    return -np.inf


def lnprob(theta,wave,flux,error):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,wave,flux,error)


class SED:
    '''Fits a modified black body SED to a set of fluxes, takes args of the fluxes,errors,wavelengths and initial guesses to the parameters'''

    def __init__(self,waves,fluxes,errors,theta0,nwalkers):

        self.wvs = np.array(waves)
        self.waves =np.array(waves) * 1e-6
        self.fluxes = np.array(fluxes) 
        self.errors = np.array(errors) 
        T0,beta0,log_A0 = theta0
        self.ndim = np.shape(theta0)[0]
        self.nwalkers = nwalkers

        self.Tinit = T0 + 10*np.random.random(nwalkers)
        self.betainit = beta0 + 0.3*np.random.random(nwalkers)
        self.log_Ainit = log_A0 + 2*np.random.random(nwalkers)

    def run_MCMC(self):

        pos = [np.array(x) for x in zip(self.Tinit,self.betainit,self.log_Ainit)]

        sampler = emcee.EnsembleSampler(self.nwalkers,self.ndim,lnprob,args=(self.waves,self.fluxes,self.errors))

        sampler.run_mcmc(pos,600)

        samples = sampler.chain[:,100:,:].reshape(-1,self.ndim)

        return sampler,samples

    def MCMC_plots(self):

        sampler,samples = self.run_MCMC()
        labels = {}
        labels[1] = 'T'
        labels[2] = r'$\beta$'
        labels[3] = r'$log{A}$'

        for x in range(self.ndim):
            plt.figure()
            for i in range(self.nwalkers):
                vals = []
                for j in range(600):
                    val = sampler.chain[i][j][x]
                    vals.append(val)
                plt.plot(range(1,601,1),vals)
                plt.xlabel('Step')
                plt.ylabel(labels[x+1])
            plt.show()

        fig = corner.corner(samples,labels=['$T$',r'$\beta$',r'$log{A}$'])
        fig.show()

    def param_vals(self):

        sampler,samples = self.run_MCMC()

        T_mcmc,beta_mcmc,log_A_mcmc = map(lambda v: (v[1],v[2] - v[1],v[1]-v[0]),zip(*np.percentile(samples,[32,50,68],axis=0)))

        return T_mcmc,beta_mcmc,log_A_mcmc

    def plot_SED(self):

        T_mcmc,beta_mcmc,log_A_mcmc = self.param_vals()

        plt_wvs = np.arange(75,600,1)
        plt_nus = [_c/(wv * 1e-6) for wv in plt_wvs]
        plt_fluxes = list(bb(plt_wvs,T_mcmc[0],beta_mcmc[0],10**(log_A_mcmc[0])))

        plt.plot(plt_wvs,plt_fluxes)
        plt.scatter(self.wvs,self.fluxes)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux/Jy')
        plt.show()
                    

        
    
