import os,sys
import glob
from maldi_class import Maldi
import numpy as np
import argparse as ap
import subprocess as sb
import matplotlib.pyplot as plt
#folder = os.path.abspath(sys.argv[1])

parser = ap.ArgumentParser(description='Reads and integrates MALDI spectrum for a specified fragment family, in addition calculation of SSR and SSR sensitivity is also performed.')
parser.add_argument('--heavy',dest='heavy_ligand',help='The empirical formula of the heavier ligand with the number of atoms delimiting each element (e.g. DDT = S1C12H25)')
parser.add_argument('--light',dest='light_ligand',help='The empirical formula of the lighter ligand with the number of atoms delimiting each element (e.g. DDT = S1C12H25)')
parser.add_argument('--fnum',dest='fragment_number',help='The number of total ligands in the fragment family')
parser.add_argument('--folder',dest='folder',help='The folder containing the raw MALDI-MS spectra')
parser.add_argument('--offset',dest='offset',default='0',help='Offset MALDI spectrum by constant to account for calibration')
parser.add_argument('--fmetal',dest='fmetal',default='Ag',help='Nanoparticle metal found in fragment (e.g. Au, Ag, Pd ...)')

args = parser.parse_args()

folder = args.folder
num_files = len([mld_file for mld_file in glob.iglob(folder+'/*.txt')])
fragment_number = int(args.fragment_number)
heavy_ligand = str(args.heavy_ligand)
light_ligand = str(args.light_ligand)
offset = int(args.offset)
results = np.empty((num_files,fragment_number+1+1+1+1+1))
metal_fragment = args.fmetal+str(fragment_number+1) if args.fmetal == 'Ag' else args.fmetal+str(fragment_number)
print("Metal fragment is: "+metal_fragment)

file_list = np.empty(num_files,dtype="S100")

for i,maldi_file in enumerate(glob.iglob(folder+'/*.txt')):
	mld = Maldi(maldi_file,metal_fragment,heavy_ligand,light_ligand,fragment_number,offset=offset)
	#mld.print_peak_family()
	peaks = mld.integrate_peak_family()
	noise = mld.get_noise(1700,1750)
	normed_peaks = peaks/np.sum(peaks)
	results[i,0:(fragment_number+1)] = normed_peaks
	results[i,(fragment_number+1)] = mld.get_ssr()
	results[i,fragment_number+2] = mld.get_ligand1_fraction()
	#ssrs = mld.get_ssr_noisiness(noise,200) not analyzing noise for speed
	results[i,fragment_number+3] = 0
	results[i,fragment_number+4] = mld.get_total_integrated_area()
	file_list[i] = maldi_file
	percentage = (float(i+1)/float(num_files))*100
	hashes = "#"*int(percentage)
	spaces = " "*(100-int(percentage))
	sb.call(["printf",'\r'+hashes+spaces+'(%2.0f%%)',str(percentage)])

sb.call(["printf","\n"])
np.savetxt('ssrs.txt',results,header="peak0\tpeak1\tpeak2\tpeak3\tpeak4\tpeak5\tSSR\txHeavy\tstdev\tpeak_area")
np.savetxt('file_list.txt',file_list,fmt="%s")

