import numpy as np
import pandas as pd
from matplotlib import pylab
import time
from IM_Model_precise import IM_Model_precise

#parameters for PS133k in ref: Costanzo, S. et al. Macromolecules 2016, 49, 3925-3935.
GN0 = 2.79E5
τc = 1.24E-3
τR = 0.124
λmax = 3.34
CQ = 6
β = 0.25
Gi = np.array([1.14E2, 2.52, 8.80E-1, 3.81E-1, 2.07E-1, 1.70E-1, 1.76E-1, 1.37E-1, 2.03E-1]) * GN0
τi_eq = np.array([2.16E-3, 5.90E-2, 3.33E-1, 1.53, 6.87, 2.9E1, 1.11E2, 3.63E2, 8.92E2]) * τc
shear_rate = 31.6 # 1/s
δt = 0.0001 # seconds
finish_time = 10 # senconds

t0 = time.time()
model = IM_Model_precise(λmax=λmax,
                τR=τR,
                shear_rate=shear_rate,
                β = β,
                δt=δt,
                finish_time=finish_time,
                Gi=Gi,
                τi_eq=τi_eq,
                CQ=CQ)
model.IM_Tumbling_Multimode_integral()
t_span, σ_span = model.calculate_stress()
print('Cost time: ', time.time() - t0, ' s')

pylab.loglog(model.t_span,σ_span[:,0,1]/shear_rate,'bo')
pylab.loglog(model.t_span,(σ_span[:,0,0]-σ_span[:,1,1])/shear_rate,'r-')
pylab.loglog(model.t_span,(σ_span[:,2,2]-σ_span[:,1,1])/shear_rate,'g.')