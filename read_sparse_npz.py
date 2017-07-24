import numpy as np
from scipy.sparse import csr_matrix

fpath = 'I:\\BU4D_2_mesh\\coords\\F001_1_066.npz'
t = np.load(fpath)
col = t['col'].astype(np.int32)
ptr = t['ptr'].astype(np.int32)
x_local = t['x_local'].astype(np.float32)
y_local = t['y_local'].astype(np.float32)

dim = len(ptr)-1
x_mat = csr_matrix((x_local, col, ptr),shape=(dim,dim))
y_mat = csr_matrix((y_local, col, ptr),shape=(dim,dim))

for i in range(0,dim,100):
	res = 0
	for j in range(dim):
		if x_mat[i,j] != 0:
			res += 1
	print(i,res)