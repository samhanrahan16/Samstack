import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
import scipy.optimize as opt

def mod_bb_func(T,beta,A,freq,h,c,kb):
    pre = 2*h*(freq**3)/(c**2)
    exp = h*freq/(kb*T)
    flux = A * pre * (freq**beta) * (1/(np.exp(exp) - 1))

    return flux

def loglike((T,beta,log_A),freqs,fluxes,errors,h,c,kb):
    A = np.exp(log_A)
    model_vals = np.array([mod_bb_func(T,beta,A,freq,h,c,kb) for freq in freqs])
    flux = np.array(fluxes)
    errs = np.array(errors)
    inv_sigma2 = 1/errs**2
    like = -0.5*np.sum((flux-model_vals)**2 * inv_sigma2 - np.log(inv_sigma2))

    return like

def loglike_opt((T,beta,log_A),freqs,fluxes,errors,h,c,kb):
    like = loglike((T,beta,log_A),freqs,fluxes,errors,h,c,kb):
    return -like

def lnprior(T,beta,log_A):
    if T>0 and -3<log_A<12 and 0.0<beta<2.0:
        return 0.0
    return -np.inf

def lnprob((T,beta,log_A),freqs,fluxes,errors,h,c,kb):
    lp = lnprior(T,beta,log_A)
    if not np.isfinite(lp):
        return -np.inf
    return lp + loglike((T,beta,log_A),freqs,fluxes,errors,h,c,kb)
    

class modified_bb:
    '''Fits a modified blackbody to the stacked flux at each wavelength,wavelengths in micrometers'''

    def __init__(self,fluxes,wavelengths,errors,T0,beta0,log_A0):

        self.fluxes = fluxes
        self.wavelengths = wavelengths
        self.errors = errors
        self.c = 3e8
        self.h = 6.63e-34
        self.kb = 1.38e-23
        self.frequencies = [c/(wavelength * 1e-6) for wavelength in self.wavelengths]
        self.T0 = T0
        self.beta0 = beta0
        self.log_A0 = log_A0

    def fit_SED(self):

        T0 = np.array([self.T0])
        beta0 = np.array([self.beta0])
        log_A0 = np.array([self.log_A0])

        result = opt.minimize(loglike_opt,(T0,beta0,log_A0),method='Nelder-mead',args = (self.frequencies,self.fluxes,self.errors,self.h,self.c,self.kb))

        ndim,nwalkers = 3,100
        pos = [result.x * + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(self.frequencies,self.fluxes,self.errors,self.h,self.c,self.kb))

        sampler.rum_mcmc(pos,500)

        self.samples = sampler.chain[:,50,:].reshape((-1,ndim))

        T_mcmc,beta_mcmc,log_A_mcmc = map(lambda v: (v[1],v[2]-v[1],v[1]-v[0]),zip(*np.percentile(self.samples,[16,50,84],axis=0)))

        return T_mcmc,beta_mcmc,log_A_mcmc

    def plot_corner(self,save_path):

        T_mcmc,beta_mcmc,log_A_mcmc = self.fit_SED()

        fig = corner.corner(self.samples,labels=['$T/K$','$beta$','$log_A$'])
        file_path = save_path + '/corner_plots.png'
        fig.savefig(file_path)
        fig.show()

    def plot_modified_bb(self,save_path):

        T_mcmc,beta_mcmc_log_A_mcmc = self.fit_SED()

        plot_wvs = np.arange(15.,600.,0.1)
        plot_freqs = [self.c/(wv*1e-6) for wv in plot_wvs]
        A_mcmc = np.exp(log_A_mcmc[0])

        plot_fluxes = [mod_bb_func(T_mcmc[0],beta_mcmc[0],A_mcmc,freq,self.h,self.c,self.kb) for freq in plot_freqs]

        image_name = 'modified_bb_' + str(np.round(T_mcmc[0]))
        save_name = image_name + '.png'
        file_path = save_path + '/' + save_name

        plt.plot(plot_wvs,plot_fluxes)
        plt.scatter(self.wavelengths,self.fluxes)
        plt.xlabel('Wavelength/' + r'$\mu$' + 'm')
        plt.ylabel('Flux/Jy')
        plt.title(image_name)
        plt.savefig(file_path)
        plt.show()

        
