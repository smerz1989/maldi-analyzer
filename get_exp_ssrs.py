import subprocess as sb
import numpy as np
from functools import reduce

sb.call(["printf","Analyzing Mercaptohexanol and Dodecanethiol Spectrum:\n"])
sb.call(["python2","analyze_folder.py","--fnum","5","--light","S1C6H13O1","--heavy","S1C12H25","--folder","../AgNP_MHA_DDT/","--offset","0"])
sb.call(["cp","ssrs.txt","mha_ddt_ssrs.txt"])

sb.call(["printf","Analyzing Mercaptoundecanol and Dodecanethiol-d25 Spectrum:\n"])
sb.call(["python2","analyze_folder.py","--fnum","5","--light","S1C12H25O1","--heavy","S1C12H50","--folder","../AgNP_MUDA_DDTd25/","--offset","0"])
sb.call(["cp","ssrs.txt","muda_ddt_ssrs.txt"])

sb.call(["printf","Analyzing Mercaptoundecanol and Hexanethiol Spectrum:\n"])
sb.call(["python2","analyze_folder.py","--fnum","5","--heavy","S1C12H25O1","--light","S1C6H13","--folder","../AgNP_MUDA_HT/","--offset","0"])
sb.call(["cp","ssrs.txt","muda_ht_ssrs.txt"])

mha_ddt_ssrs = np.loadtxt('mha_ddt_ssrs.txt')
muda_ddt_ssrs = np.loadtxt('muda_ddt_ssrs.txt')
muda_ht_ssrs = np.loadtxt('muda_ht_ssrs.txt')

mha_ddt_ssrs_50 = mha_ddt_ssrs[np.where((mha_ddt_ssrs[:,7]<0.55) & (mha_ddt_ssrs[:,7]>0.45))[0],:]
muda_ddt_ssrs_50 = muda_ddt_ssrs[np.where((muda_ddt_ssrs[:,7]<0.55) & (muda_ddt_ssrs[:,7]>0.45))[0],:]
muda_ht_ssrs_50 = muda_ht_ssrs[np.where((muda_ht_ssrs[:,7]<0.55) & (muda_ht_ssrs[:,7]>0.45))[0],:]

mha_ddt_ssr_mean,mha_ddt_ssr_std = (np.mean(mha_ddt_ssrs_50[:,6]),np.std(mha_ddt_ssrs_50[:,6]))
muda_ddt_ssr_mean,muda_ddt_ssr_std = (np.mean(muda_ddt_ssrs_50[:,6]),np.std(muda_ddt_ssrs_50[:,6]))
muda_ht_ssr_mean,muda_ht_ssr_std = (np.mean(muda_ht_ssrs_50[:,6]),np.std(muda_ht_ssrs_50[:,6]))

print("MHA/DDT:\n")
print(len(mha_ddt_ssrs_50))
print(mha_ddt_ssr_mean)
print(mha_ddt_ssr_std)
print("MUDA/DDT:\n")
print(len(muda_ddt_ssrs_50))
print(muda_ddt_ssr_mean)
print(muda_ddt_ssr_std)
print("MUDA/HT:\n")
print(len(muda_ht_ssrs_50))
print(muda_ht_ssr_mean)
print(muda_ht_ssr_std)

with open('ssr_data_experimental.csv','w') as f:
    f.write("Category,SSR,stdev\n")
    f.write("DDT_MHOH,"+str(mha_ddt_ssr_mean)+","+str(mha_ddt_ssr_std)+'\n')
    f.write("DDT_MUDA,"+str(muda_ddt_ssr_mean)+","+str(muda_ddt_ssr_std)+'\n')
    f.write("HT_MUDA,"+str(muda_ht_ssr_mean)+","+str(muda_ht_ssr_std)+'\n')

with open('maldi_spectrum_experimental.csv','w') as f:
	f.write("Category,Peak_Num,Mean,Stdev\n")
	f.write('\n'.join(["DDT_MHA,"+str(i)+","+str(np.mean(mha_ddt_ssrs_50[:,i]))+","+str(np.std(mha_ddt_ssrs_50[:,i])) for i in range(6)]))
	f.write('\n')
	f.write('\n'.join(["DDT_MUDA,"+str(i)+","+str(np.mean(muda_ddt_ssrs_50[:,i]))+","+str(np.std(muda_ddt_ssrs_50[:,i])) for i in range(6)]))
	f.write('\n')
	f.write('\n'.join(["HT_MUDA,"+str(i)+","+str(np.mean(muda_ht_ssrs_50[:,i]))+","+str(np.std(muda_ht_ssrs_50[:,i])) for i in range(6)]))
	#f.write("Category,Peak0,Stdev,Peak1,Stdev,Peak2,Stdev,Peak3,Stdev,Peak4,Stdev,Peak5,Stdev\n")
	#f.write("DDT_MHA,"+','.join([str(np.mean(mha_ddt_ssrs_50[:,i]))+","+str(np.std(mha_ddt_ssrs_50[:,i])) for i in range(6)])+'\n')
	#f.write("DDT_MUDA,"+','.join([str(np.mean(muda_ddt_ssrs_50[:,i]))+","+str(np.std(muda_ddt_ssrs_50[:,i])) for i in range(6)])+'\n')
	#f.write("HT_MUDA,"+','.join([str(np.mean(muda_ht_ssrs_50[:,i]))+","+str(np.std(muda_ht_ssrs_50[:,i])) for i in range(6)])+'\n')







