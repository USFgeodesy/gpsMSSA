'''
Utilities for reading in GPS data for gpsMssa
'''
import numpy as np
import scipy.stats.stats as stats
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pylab as plt
from scipy.signal import detrend
import matplotlib.mlab as mlab

def readData(filename):
    '''
    read in data in format:
    date location uncertainty
    '''
    #read in filename splitting into useful python array
    date,loc,uncert = np.genfromtxt(filename,dtype=float,usecols = (0,1,2), unpack = True)
    dateNum = [float(d) for d in date]
    loc = [float(l) for l in loc]
    uncert = [float(x) for x in uncert]
    comp = filename.split('.')[1]
    return timeseries(comp,dateNum,loc,uncert)

def decimal2datetime(number):
    '''
    convert decimal year to date time object
    '''
    year = int(number)
    d = timedelta(days=(number - year)*365)
    day_one = datetime(year,1,1)
    date = d + day_one
    return date.date()

class timeseries(object):
      '''
      timeseries is an object for holding individual components of GPS timeseries
      Attributes:
	  component : component name i.e lat,lon, rad
	  times  :times of measurement
	  location : displacment
	  uncertainty : uncertainty
      '''
      def __init__(self,component,times,location,uncertainty):
	        self.component = component
	        self.times = times
	        self.location = location
	        self.uncertainty = uncertainty

      def get_component(self):
	        return self.component

      def get_times(self):
	        return self.times

      def get_location(self):
	        return self.location

      def get_uncertainty(self):
	        return self.uncertianty
      def ts2pandas(self,start,end):
          '''
          convert time series object to pandas series
          '''
          #self.location = stats.zscore(self.location)
          sr = pd.Series((self.location),[decimal2datetime(time) for time in self.times], name = self.component)
          sr = sr.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')
          sr = sr[pd.to_datetime(decimal2datetime(start)).date():pd.to_datetime(decimal2datetime(end)).date()]
          sr = pd.Series(mlab.detrend_linear(sr.values.T[0]),index = sr.index)
          times = pd.date_range(pd.to_datetime(decimal2datetime(start)), pd.to_datetime(decimal2datetime(end)))
          sr = sr.reindex(index = times)
          #sr = sr.dropna()
          #sr = pd.Series(sr.get_values(),sr.index)

          #times = pd.date_range(pd.to_datetime(decimal2datetime(start)), pd.to_datetime(decimal2datetime(end)))
          #sr = sr.reindex(index = times)
          #sr = sr - sr.mean()
          #sr = pd.Series(mlab.detrend_linear(sr.values.T[0]),index = sr.index)
          sr2 = sr 
          for i,s in enumerate(sr.values):
              if np.isnan(s):
                  sr.values[i] = np.random.normal(scale = np.std(sr))
          
          return sr,sr2 #,times

      def plot(self):
          plt.errorbar(self.times, self.location, yerr= self.uncertainty)
          plt.title = self.component
          plt.show()
          return

      def fillGaps(self,start,end):
            '''
            assume first and last epoch have a position of zero
            fill all gaps in between in the data assuming  a linear trend
            and gausian noise with std 1
            not implemented
            '''
            #calculate number of days in observation period
            numdays = (decimal2datetime(end)-decimal2datetime(start)).days



      def detrend(self):
          return
      def zscoreData(self):
           '''
           remove the mean and divide by the standard deviation
           '''
           self.location = stats.zscore(self.location)
           return
