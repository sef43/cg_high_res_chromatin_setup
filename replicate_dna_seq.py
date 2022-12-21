import numpy as np
import sys

N=int(sys.argv[1])

seq=np.loadtxt("DNA_sequence.txt",dtype=str)

ids=np.array(seq[:,0],dtype=int)
bps=seq[:,1]

M = ids[-1]+2


new_ids_list=[ids]
new_bps_list=[bps]

print("replicate 1 1 "+str(N))

for n in range(1,N):
	new_ids=ids+M*n
	new_bps=bps[:]
	
	d1 = new_ids_list[-1][-1]
	d2 = new_ids[0]
	
	p11=str(d1+1)
	p12=str(d1+2)
	p21=str(d2+1)
	p22=str(d2+2)
	
	d1=str(d1)
	d2=str(d2)

	print("create_bonds single/bond 2 "+d1+" "+d2)
	print("create_bonds single/bond 3 "+p11+" "+p21)
	print("create_bonds single/bond 3 "+p11+" "+p22)
	print("create_bonds single/bond 3 "+p12+" "+p21)
	print("create_bonds single/bond 3 "+p12+" "+p22)

	new_ids_list.append(new_ids)
	new_bps_list.append(new_bps)
	

out_ids = np.array(np.concatenate(new_ids_list),dtype=str)
out_bps = np.concatenate(new_bps_list)
new_len = len(out_bps)
#print(out_ids.shape)
#print(out_bps.shape)

out=np.vstack((out_ids, out_bps)).T

np.savetxt("new_DNA_sequence.txt", out, fmt="%s", header = str(new_len))



