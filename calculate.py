import numpy as np
import pandas as pd
from matplotlib import pylab
import time
from IM_Model import IM_Tumbling_Multimode_integral
import csv

#parameters
τR = 17.236 #125oC
λmax = 4.8
CQ = 6
β = 1.0

τi_eq = np.array([0.0001,0.000630957,0.00158489,0.01,0.025118864,0.063095734,0.158489319,0.398107171,1,2.511886432,6.309573445,15.84893192,100])
Gi = np.array([292271168.6,1252840.227,4846363.131,987743.419,7718.216398,258428.423,233109.8931,56871.30167,107624.516,25889.62968,75359.12278,67691.40492,1681.8998])
#shear_rate = 1.0 # 1/s
total_strain = 200
total_points_n = 1e5
shear_rates = [10.0, 5.62, 3.16, 1.78, 1.0, 0.562, 0.316, 0.178, 0.1, 0.0562, 0.0316, 0.0178, 0.01, 0.00562, 0.00316, 0.00178, 0.001]
print("shear rates: {}".format(shear_rates))
for shear_rate in shear_rates:
    finish_time = total_strain/shear_rate # seconds
    if finish_time > 1000:
        finish_time = 1000
    elif finish_time < 120:
        finish_time = 120
    δt = finish_time/total_points_n # seconds
    t0 = time.time()
    t_span, τd_span, σ_span, λ_span, S_span = IM_Tumbling_Multimode_integral(λmax=λmax,
                                                                            τR=τR,
                                                                            shear_rate=shear_rate,
                                                                            β = β,
                                                                            δt=δt,
                                                                            finish_time=finish_time,
                                                                            Gi=Gi,
                                                                            τi_eq=τi_eq,
                                                                            CQ=CQ)
    print("Cost time:{:.1f}s".format(time.time() - t0))
    #write in file
    with open(f"{shear_rate:.5f}s-1.csv",mode="w",newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(["strain","time","stress","viscosity","N1","-N2","N1appCPP4R8","N1appCPP4R15"])
        csv_writer.writerow(["--", "s","Pa","Pa-s","Pa","Pa","Pa","Pa"])
        csv_writer.writerow([f"{shear_rate:.5f}s\\+(-1)"]*8)
        stress = σ_span[:,0,1]
        Viscosity = σ_span[:,0,1]/shear_rate
        N1 = σ_span[:,0,0]-σ_span[:,1,1]
        N2 = -σ_span[:,2,2]+σ_span[:,1,1]
        N1appCPP4R8 = N1+2*(N1+2*N2)*0.693147
        N1appCPP4R15 = N1+2*(N1+2*N2)*1.3217558
        for ind, row in enumerate(zip(t_span * shear_rate, t_span, stress, Viscosity, N1, -N2, N1appCPP4R8, N1appCPP4R15)):
            if ind > 100 and ind < 1000 and ind % 3 != 0:
                continue
            elif ind > 1000 and ind < 1E4 and ind % 10 != 0:
                continue
            elif ind >1E4 and ind<1E5 and ind % 30 !=0:
                continue
            elif ind >1E5 and ind % 100 !=0:
                continue
            else:
                csv_writer.writerow(row)
    max_index = np.argmax(stress)
    max_stress = stress[max_index]
    overshoot_strain = t_span[max_index] * shear_rate
    with open("steady.csv",mode="a",newline='') as f:
       csv_writer = csv.writer(f)
       csv_writer.writerow([shear_rate, Viscosity[-1], N1[-1], -N2[-1], overshoot_strain, max_stress]) 
