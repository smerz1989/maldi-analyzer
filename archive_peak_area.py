import numpy as np

muda_ddt_data = np.loadtxt('muda_ddt_ssrs.txt')
ht_muda_data = np.loadtxt('muda_ht_ssrs.txt')
mha_ddt_data = np.loadtxt('mha_ddt_ssrs.txt')

print("MUDA_DDT Peak Sum: "+str(np.sum(muda_ddt_data[:-1]))+" ;Count is "+str(len(muda_ddt_data)))
print("MUDA_HT Peak Sum: "+str(np.sum(ht_muda_data[:-1])))+" ;Count is "+str(len(ht_muda_data))
print("MHA_DDT Peak Sum: "+str(np.sum(mha_ddt_data[:-1])))+" ;Count is "+str(len(mha_ddt_data))

print(str(np.sum(muda_ddt_data[:-1])))
print(str(np.sum(ht_muda_data[:-1])))
print(str(np.sum(mha_ddt_data[:-1])))
