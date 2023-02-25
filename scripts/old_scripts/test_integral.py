import scipy.integrate as integrate
import scipy.special as special
import numpy
import sympy as smp


def prob_n_reads_old(n, N, mean_, beta_):

    return (special.gamma(beta_ + n)/(special.factorial(n)*special.gamma(beta_))) * (((mean_*N)/(beta_ + mean_*N))**n) * ((beta_/(beta_ + mean_*N))**beta_)


def prob_n_reads(n, N, mean_, beta_):

    return numpy.exp( special.gammaln(beta_+n) - special.gammaln(n+1) - special.gammaln(beta_) )   * (((mean_*N)/(beta_ + mean_*N))**n) * ((beta_/(beta_ + mean_*N))**beta_)


def integrand(n, N, mean_, beta_):
    return (n/N)*numpy.log(n/N)* prob_n_reads(n, N, mean_, beta_)







#(special.gamma(beta_ + n)/(special.factorial(n)*special.gamma(beta_))) * 

I = integrate.quad(integrand, a=0, b=10000,  args=(10000, 0.001, 1.4))

first_moment_x_logx = I[0]

print(first_moment_x_logx)

#print(first_moment_x_logx)

#beta_ = 5

#for n_i in range(200):

#    term = (special.gamma(beta_ + n_i)/(special.factorial(n_i)*special.gamma(beta_)))

#    term_approx = (beta_**n_i) / special.factorial(n_i)
#    print(term, term_approx )


