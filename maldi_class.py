import numpy as np

from scipy.integrate import cumtrapz
from scipy.stats import binom
from itertools import cycle
import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure
import sys
import re

class Maldi(object):
	au = [(196.966543,1.00),]
	pd = [(101.905634,0.0102),(103,0),(103.904029,0.1114),(104.905079,0.2233),(105.903478,0.2733),(107,0),(107.903895,0.2646),(109,0),(109.905167,0.1172)]
	ag = [(106.905092,0.51839),(108,0),(108.904756,0.48161)]
	s = [(31.97207070,0.9493),(32.97145843,0.0076),(33.96786665,0.0429),(35,0),(35.96708062,0.0002)]
	c = [(12.000,0.9893),(13.003,0.0107)]
	h = [( 1.007825032,0.999885),(2.014101778,0.000115)]
	o = [(15.99491463,0.99757),(16.9991312,0.00038),(17.9991603,0.00205)]

	element_dict = {'Ag':ag,'Au':au,'Pd':pd,'S':s,'C':c,'O':o,'H':h}

	def __init__(self,maldi_filename,metal,ligand1,ligand2,ligand_number,offset=0):
		self.maldi_filename = maldi_filename
		self.maldi_data = np.loadtxt(self.maldi_filename,skiprows=1)
		self.maldi_data[:,0]=self.maldi_data[:,0]+offset
		self.metal = metal
		self.ligand1 = ligand1
		self.ligand2 = ligand2
		self.ligand_number = ligand_number

	def find_peak(self,molecule_string):
		element_matches = re.finditer(r"[a-zA-Z]+",molecule_string)
		count_matches = re.finditer(r"[0-9]+",molecule_string)
		molecule = []
		for count_string,element_string in zip(count_matches,element_matches):
			element = self.element_dict[element_string.group(0)]
			element_probs = [isotope[1] for isotope in element]
			molecule.append((element,int(count_string.group(0))))
		peak = np.array([1])
		peakx_start = 0
		for element in molecule:
			element_probs = [isotope[1] for isotope in element[0]]
			peakx_start+=element[0][0][0]*element[1]
			for count in range(element[1]):
				peak = np.convolve(peak,element_probs)
		peakx = np.arange(peak.shape[0],dtype=float)
		peakx+=peakx_start
		return (peakx,peak)

	def get_peak_family(self):
		peaks = np.empty((self.ligand_number+1,3))
		for i in range(self.ligand_number+1):
			peakx,peak=self.find_peak(self.metal+self.ligand1*i+self.ligand2*(self.ligand_number-i))
			peaks[i,0]= peakx[max(np.where(np.cumsum(peak)<0.005)[0])] if np.any(np.cumsum(peak)<0.005) else peakx[0]
			peaks[i,1]=np.average(peakx,weights=peak)
			peaks[i,2]= peakx[min(np.where(np.cumsum(peak)>0.995)[0])] if np.any(np.cumsum(peak)>0.995) else peakx[-1]
		return peaks

	def print_peak_family(self):
		peaks = np.empty((self.ligand_number+1,3))
		for i in range(self.ligand_number+1):
			peakx,peak=self.find_peak(self.metal+self.ligand1*i+self.ligand2*(self.ligand_number-i))
			peaks[i,0]= peakx[max(np.where(np.cumsum(peak)<0.005)[0])] if np.any(np.cumsum(peak)<0.005) else peakx[0]
			peaks[i,1]=np.average(peakx,weights=peak)
			peaks[i,2]= peakx[min(np.where(np.cumsum(peak)>0.995)[0])] if np.any(np.cumsum(peak)>0.995) else peakx[-1]
			print("Peak "+str(i)+":")
			print(str(peaks[i,0])+" "+str(peaks[i,1])+" "+str(peaks[i,2]))

	def get_noise(self,noise_start,noise_end):
		noise_data = self.maldi_data[(self.maldi_data[:,0]>=noise_start) & (self.maldi_data[:,0]<=noise_end)]
		return(np.average(noise_data),np.std(noise_data))

	def integrate_peak(self,molecule_string):
		(peakxs,peaks) = self.find_peak(molecule_string)
		#maldi_data = np.loadtxt(self.maldi_filename,skiprows=1)
		lowerbound = peakxs[max(np.where(np.cumsum(peaks)<0.005)[0])] if np.any(np.cumsum(peaks)<0.005) else peakxs[0]
		upperbound = peakxs[min(np.where(np.cumsum(peaks)>0.995)[0])] if np.any(np.cumsum(peaks)>0.995) else peakxs[-1]
		peak_data = self.maldi_data[(self.maldi_data[:,0]>lowerbound) & (self.maldi_data[:,0]<upperbound)]
		area = cumtrapz(peak_data[:,1],peak_data[:,0])
		return area[-1]

	def integrate_peak_with_noise(self,molecule_string,noise):
		(peakxs,peaks) = self.find_peak(molecule_string)
		#maldi_data = np.loadtxt(self.maldi_filename,skiprows=1)
		lowerbound = peakxs[max(np.where(np.cumsum(peaks)<0.005)[0])] if np.any(np.cumsum(peaks)<0.005) else peakxs[0]
		upperbound = peakxs[min(np.where(np.cumsum(peaks)>0.995)[0])] if np.any(np.cumsum(peaks)>0.995) else peakxs[-1]
		peak_data = self.maldi_data[(self.maldi_data[:,0]>lowerbound) & (self.maldi_data[:,0]<upperbound)]
		peak_data[:,1] = peak_data[:,1]-np.random.normal(loc=noise[0]-2*noise[1],scale=noise[1],size=len(peak_data[:,1]))
		area = cumtrapz(peak_data[:,1],peak_data[:,0])
		return area[-1]


	def integrate_peak_family(self):
		peaks = np.zeros(self.ligand_number+1)
		for i in range(self.ligand_number+1):
			peaks[i]=self.integrate_peak(self.metal+self.ligand1*i+self.ligand2*(self.ligand_number-i))
		return peaks

	def integrate_peak_family_with_noise(self,noise):
		peaks = np.zeros(self.ligand_number+1)
		for i in range(self.ligand_number+1):
			peaks[i]=self.integrate_peak_with_noise(self.metal+self.ligand1*i+self.ligand2*(self.ligand_number-i),noise)
		return peaks

	def get_total_integrated_area(self):
		peaks = self.integrate_peak_family()
		return sum(peaks)

	def get_ligand1_fraction(self):
		peaks = self.integrate_peak_family()
		normed_peaks = peaks/np.sum(peaks)
		bins = np.arange(self.ligand_number+1)
		ligand1_fraction = np.sum(bins*normed_peaks/float(self.ligand_number))
		return ligand1_fraction

	def get_ssr(self):
		peaks = self.integrate_peak_family()
		normed_peaks = peaks/np.sum(peaks)
		bins = np.arange(self.ligand_number+1)
		ligand1_fraction = np.sum(bins*normed_peaks/float(self.ligand_number))
		binomial = map(lambda k: binom.pmf(k,self.ligand_number,ligand1_fraction),np.arange(self.ligand_number+1))
		ssr = np.sum((normed_peaks-binomial)**2)
		return ssr
	
	def plot_peak_family(self,filename):
		peaks = self.get_peak_family()	
		minimum_mz = peaks[0,0]
		maximum_mz = peaks[-1,2]
		data_within_range = self.maldi_data[np.where((self.maldi_data[:,0]>minimum_mz) & (self.maldi_data[:,0]<maximum_mz))[0],:]
		plt.plot(data_within_range[:,0],data_within_range[:,1])
		plt.savefig(filename,format='png')
	
	def get_ssr_noisiness(self,noise,numsamples):
		ssrs = np.empty(numsamples)
		for i in range(numsamples):
			#if ((i+1)%10)==0:
			#	print("On sample "+str(i+1))
			peaks = self.integrate_peak_family_with_noise(noise)
			normed_peaks = peaks/np.sum(peaks)
			bins = np.arange(self.ligand_number+1)
			ligand1_fraction = np.sum(bins*normed_peaks/float(self.ligand_number))
			binomial = map(lambda k: binom.pmf(k,self.ligand_number,ligand1_fraction),np.arange(self.ligand_number+1))
			ssrs[i] = np.sum((normed_peaks-binomial)**2)
		return ssrs
