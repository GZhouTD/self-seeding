import numpy as np
def moddist(f2o,gamrange):
	f = open(f2o) 
	alllines = f.readlines()
	del alllines[0:5]
	m = len(alllines)
	f.close
	data = ''.join(alllines)
	data = data.split()
	data = np.array(data)
	data = data.reshape(m,6)
	gamma = data[:,-1]
