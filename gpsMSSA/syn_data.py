'''
Generate synthetic data for testing mssa code

library for getting building functions
pass the objects following into the function object:
	heavy
	invHeavy
	trend
	periodAnn
	periodSemiAnn
	jump
	exponential
	sse
'''
import numpy as np
import pandas as pd
import inspect
import scipy


class function(object):
	'''
	function is the function that describes the gps time series

	Attributes:
		equation : equation describing function
	'''

	def __init__(self, equation):
		self.equation = equation

	def get_equation(self):
		return self.equation

	def get_params(self):
		params = [eq.__dict__ for eq in self.equation]
		return params

	def eq2latex(self):
		'''
		print a latex string using matplotlib representing the function
		'''
		return

def model_tseries(times,funList):
      '''
      takes in measurement times and computes the value of the function at those times
      '''
      modeledData = []
      for t in times:
		  for fun in funList:
			  fun.t = t
			  if fun.name == 'trend' or fun.name == 'periodAnn' or fun.name == 'periodSemiAnn':
				  fun.t0 = min(times)
			  if fun.name =='periodAnnVar':
				if (t-times[0])%1.0 < 0.002:
				   fun.A1  = (np.random.sample())
				   fun.A2 = (np.random.sample())
		  fun = function(funList)
		  modeledData.append(compute(fun))
      return modeledData

class heavy(object):
	'''
	the heaviside function

	Atributes:
		t = time
		t0 = the time at which to apply the heaviside function
	'''
	def __init__(self,t,t0):
	  self.t = t
	  self.t0 = t0
	  self.name = 'heavy'

	def heavyside(self):
	  if self.t<self.t0:
	    x = 0
	  else:
	    x = 1
	  return x

class inverseHeavy(object):
	'''
	the inverse heaviside function

	Atributes:
		t = time
		t0 = the time at which to remove the heaviside function
	'''
	def __init__(self,t,t0):
		self.t = t
		self.t0 = t0
		self.name = 'invHeavy'

	def invHeavy(self):
		if self.t>self.t0:
			x = 0
		else:
			x = 1
		return x

class trend(object):
	'''
	trend : y = mx+b

	Attributes
		t : time
		m : the slope
		b : the bias
	'''
	def __init__(self,t0,t,m,b):
		self.t0 = t0
		self.t = t
		self.m = m
		self.b = b
		self.name = 'trend'

	def trend(self):
		y = self.m*(self.t-self.t0) + self.b
		return y

class periodAnn(object):
	'''
	annual period

	Attributes:
		t = time
		A1 = Annual Cosine Amplitude
		A2 = Annual Sin Amplitude
	'''
	def __init__(self,t,t0,A1,A2):
		self.t0 = t0
		self.t 	= t
		self.A1 = A1
		self.A2 = A2
		self.name = 'periodAnn'

	def annual(self):
		y = self.A1*np.cos(2.0*np.pi*(self.t-self.t0)) + self.A2*np.sin(2.0*np.pi*(self.t-self.t0))
		return y

class periodAnnVar(object):
	'''
	annual period

	Attributes:
		t = time
		A1 = Annual Cosine Amplitude
		A2 = Annual Sin Amplitude
	'''
	def __init__(self,t,t0,A1,A2):
		self.t0 = t0
		self.t 	= t
		self.A1 = A1
		self.A2 = A2
		self.name = 'periodAnnVar'

	def annual(self):
		y = self.A1*np.cos(2.0*np.pi*(self.t-self.t0)) + self.A2*np.sin(2.0*np.pi*(self.t-self.t0))
		return y

class periodSemiAnn(object):
	'''
	semi-annual period

	Attributes:
		t = time
		A3 = Annual Cosine Amplitude
		A4 = Annual Sin Amplitude
	'''
	def __init__(self,t,t0,A3,A4):
		self.t = t
		self.t0 = t0
		self.A3 = A3
		self.A4 = A4
		self.name = 'periodSemiAnn'

	def semiannual(self):
		y = self.A3*np.cos(4.0*np.pi*(self.t-self.t0)) + self.A4*np.sin(4.0*np.pi*(self.t-self.t0))
		return y

class jump(object):
	'''
	jump (either antenna change or earthquake

	Attributes:
		t : time
		t0 : time of jump
		G : Amplutude of jump
	'''
	def __init__(self,t,t0,G):
		self.t = t
		self.t0 = t0
		self.G = G
		self.name = 'jump'

	def jump(self):
		y = self.G*heavy(self.t,self.t0).heavyside()
		return y

class exponential(object):
	'''
	The postsiesmic exponential deformation

	Attributes:
		t = time
		C = Amplitde of the decay
		eqt = Time of the start of the decay
		tau = relaction time
	'''
	def __init__(self,t,C,eqt,tau):
		self.t = t
		self.C = C
		self.eqt = eqt
		self.tau = tau
		self.name = 'exponential'

	def exponential(self):
		h = heavy(self.t,self.eqt)
		y = h.heavyside()*self.C*(1.0-np.exp(-(self.t-self.eqt)/(self.tau/365.0)))
		return y

class sse(object):
	'''
	arctan as approximation for sse from Holtkamp and Brudzinski 2008
	Attributes:
	t = time
		t01 =  median time of the sse
		tau = period over which the sse takes place.
		U = the Amplitude
	'''
	def __init__(self,t,t01,tau,U):
		self.t = t
		self.t01 = t01
		self.tau = tau
		self.U =  U
		self.name = 'sse'

	def sse(self):
		y = 0.5*self.U*(np.tanh((self.t-self.t01)/(self.tau/365.0))-1)
		return y

def compute(fun):
  '''
  computes a function at the value t specified in the function class
  '''
  value = 0
  for i in range(len(fun.equation)):
    method = getattr(fun.equation[i],inspect.getmembers(fun.equation[i],predicate=inspect.ismethod)[1][0])
    value += method()
  return value


def make_syn_df(times,funlist,numChannels):
    '''
    generate synthetic time series
    '''
    x = model_tseries(times,funlist)
    columns = ['C'+str(i) for i in range(numChannels)]
    data = np.zeros((len(times),numChannels))
    for i in range(numChannels):
        #noise = np.random.normal(0,1,len(times))
        data[:,i] = x #np.add(x,noise)
    df = pd.DataFrame(data,index=times,columns=columns)
    return df

def f_alpha_gaussian(n, q_d, alpha):

#*****************************************************************************80
#
## F_ALPHA_GAUSSIAN generates noise using a Gaussian distribuion.
#
#  Discussion:
#
#    This function generates a discrete colored noise vector X of size N
#    with a power spectrum distribution that behaves like 1 / f^ALPHA.
#    The white noise used in the generation is sampled from a Gaussian
#    distribution with zero mean and a variance of Q_D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    30 March 2011
#
#  Author:
#
#    Miroslav Stoyanov
#
#  Reference:
#
#    Jeremy Kasdin,
#    Discrete Simulation of Colored Noise and Stochastic Processes
#    and 1/f^a Power Law Noise Generation,
#    Proceedings of the IEEE,
#    Volume 83, Number 5, 1995, pages 802-827.
#
#    Miroslav Stoyanov, Max Gunzburger, John Burkardt,
#    Pink Noise, 1/F^Alpha Noise, and Their Effect on Solution
#    of Differential Equations,
#    submitted, International Journal for Uncertainty Quantification.
#
#  Parameters:
#
#    Input, integer N, the number of elements in the discrete noise vector.
#
#    Input, real Q_D, the variance of the Gaussian distribution.  A standard
#    Gaussian distribution has a variance of 1.  The variance is the square
#    of the standard deviation.
#
#    Input, real ALPHA, specifies that the computed noise is to have a
#    power spectrum that behaves like 1/f^alpha.  Normally 0 <= ALPHA <= 2.
#
#    Output, real X(N), the computed discrete noise vector.
#

#
#  Set the standard deviation of the white noise.
#
    stdev = np.sqrt(np.abs(q_d))
#
#  Generate the coefficients.
#
    hfa = np.zeros(2 * n)
    hfa[0] = 1.0
    for i in range(1,n):
        hfa[i] = hfa[i-1] * (0.5*alpha + (i-1)) / (i)
    hfa[n:2*n] = 0.0
#
#  Sample the white noise.
#
    wfa = stdev*np.random.randn(n)
#
#  Pad the array with zeros in anticipation of the FFT process.
#
    z = np.zeros(n)
    wfa = np.append(wfa,z)
#
#  Perform the discrete Fourier transforms.
#
    fh = np.fft.fft(hfa)
    fw = np.fft.fft(wfa)
#
#  Multiply the two complex vectors.
#
  #fw = fw(1:n+1)
    fh = fh[0:n+1]
    fw = fw[0:n+1]

    fw = fh*fw
    #fw = fh.*fw
#
#  Scale to match the conventions of the Numerical Recipes FFT code.
#
    fw[0]   = fw[0]/2.0
    fw[:-1] = fw[:-1]/2.0
#
#  Pad the array with zeros in anticipation of the FFT process.
#
    z = np.zeros(n - 1)
    fw = np.append(fw,z)
    #fw = [ fw; z ];
#
#  Take the inverse Fourier transform.
#
    x = np.fft.ifft(fw)
#
#  Only the first half of the inverse Fourier transform is useful.
#
    x = 2.0*np.real(x[0:n])
    return x
