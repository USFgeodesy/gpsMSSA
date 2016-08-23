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
		self.t0 = t0import sys
from distutils.core import setup, import sys
from distutils.core import setup, Extension
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

times = np.arange(2012,2016,0.01)
fun  = 'sse'
args = [0,2014.0,40.0,3.0]
funList = []
funList.append(eval(fun)(*args))
df = make_syn_df(times,funList,10)
