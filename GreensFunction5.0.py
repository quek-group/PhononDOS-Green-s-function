#1. Compute surface Green's function
#2. Get phonon sel-energy
#3. Get device's Green's function
#4. Compyte spectra

import numpy as np
from numpy import ndarray
from numpy.linalg import inv
import math

n_atm = 12
n_atm_dvc = 14
p = 3.14159265359

#1. Compute surface Green's function
# Read Dynamic matrix of SiC unitcell.(Computed from a 1x1x2 supercell)

fr_D_bulk = open('dyn_3Nx3N-bulk.dat','r')



i = 0

fr_D_bulk.seek(0)
D_bulk = [[complex(0,0) for i in range(3*n_atm)] for i in range(3*n_atm)]
for line in fr_D_bulk.readlines():
    line = line.strip()
    num = line.split()
    j = 0
    while j < 3*n_atm:
        D_bulk[i][j] = complex(float(num[2*j]),float(num[2*j+1]))
        j = j+1
    i = i+1
fr_D_bulk.close()
D_bulk=np.array(D_bulk)

D_bulk_dag = [[complex(0,0) for i in range(3*n_atm)] for i in range(3*n_atm)]
D_bulk_dag = D_bulk.conjugate().T
D_bulk = (D_bulk + D_bulk_dag)/2

#Read Connection matrix between two SiC unitcell.(Computed from a 1x1x2 supercell)
fr_D_cnct = open('dyn_3Nx3N-connection.dat','r')
i = 0 

fr_D_cnct.seek(0)
D_cnct = [[complex(0,0) for i in range(3*n_atm)] for i in range(3*n_atm)]
for line in fr_D_cnct.readlines():
    line = line.strip()
    num = line.split()
    j = 0
    while j < 3*n_atm:
        D_cnct[i][j] = complex(float(num[2*j]),float(num[2*j+1]))
        j = j+1
    i = i+1
fr_D_cnct.close()
D_cnct=np.array(D_cnct)


D_cnct_dag = D_cnct.conjugate().T

#i = 0
#while i < 35:
#    print(D_cnct_dag[i][0])
#    i = i+1


# Read device dynamic matrix

fr_D_dvc = open('dyn_3Rx3R-device.dat','r')
i = 0 

fr_D_dvc.seek(0)
D_dvc = [[complex(0,0) for i in range(3*n_atm_dvc)] for i in range(3*n_atm_dvc)]

for line in fr_D_dvc.readlines():
    line = line.strip()
    num = line.split()
    j = 0
    while j < 3*n_atm_dvc:
        D_dvc[i][j] = complex(float(num[2*j]),float(num[2*j+1]))
        j = j+1
    i = i+1
fr_D_dvc.close()
D_dvc=np.array(D_dvc)


D_dvc_dag = [[complex(0,0) for i in range(3*n_atm_dvc)] for i in range(3*n_atm_dvc)]
D_dvc_dag = D_dvc.conjugate().T
D_dvc = (D_dvc + D_dvc_dag)/2

# Construct Tdc : connection matrix between device and contact

#A = [[complex(0,0) for i in range(3*n_atm)] for i in range(3*(n_atm_dvc-n_atm))]
# Elements corresponding to interaction between Ga and bottom SiC sub is 0.

#Tdc = np.vstack((A,D_cnct_dag))

#Tdc_dag = Tdc.conjugate().T


n_diff = n_atm_dvc-n_atm
if n_diff > 0:
    A = [[complex(0,0) for i in range(3*n_atm)] for i in range(3*(n_atm_dvc-n_atm))]
    Tdc = np.vstack((A,D_cnct_dag))
elif n_diff < 0:
    Tdc = np.delete(D_cnct_dag ,i in range (n_diff+1) ,axis = 1)
else :
    Tdc = D_cnct_dag

Tdc_dag = Tdc.conjugate().T


#write Tnb
##fw_Tdc = open('Tdc.dat','w')
##
##i = 0
##j = 0
##while i < 42:
##    j = 0
##    while j < 36:
##        fw_Tdc.write("%20.6e %12.6ei" % (Tdc[i][j].real,Tdc[i][j].imag) )
##        j = j+1
##    fw_Tdc.write('\n')
##    i = i+1
##    #print(i)
##fw_Tdc.close()


#Decimation method: compute surface Green's function matrix g0,0

omega = 0.00024622101664683997 # 27cm-1 -> Ry
I = np.identity(3*n_atm)
I2 = np.identity(3*n_atm_dvc)
zero = 1E-11
fw_spect = open('Spectra.dat','w')
fw_dos = open('Dos.dat','w')
#fw_spect_ori = open('Spectra-ori.dat','w')
Ry_cm = 13.605662285137*8065.54429

omega = 0

##D_cnct_dag: H01; D_cnct:H10
H01 = D_cnct_dag
H01_d = D_cnct
#D_bulk is H00
H00 = D_bulk


while omega < 0.00118464979: #130 cm-1
#while omega < 0.0018225381523347935: # 200 cm-1
#while omega < 0.005467614457: # 600 cm-1
#while omega < 0.009112690761673967:
    invg0 = ((omega+zero*complex(0,1))**2)*I-H00
    g0 = inv(invg0)
    t0 = np.dot(g0,H01_d)
    t0_t = np.dot(g0,H01)
    ti = t0
    ti_t = t0_t

    Tn_factor = I
    Tmat = t0
    Tnt_factor = I
    Tmat_t = t0_t


    i = 0
    while i < 40:
        ti_sq = np.dot(ti,ti)
        ti_t_sq = np.dot(ti_t,ti_t)
        titi_t = np.dot(ti,ti_t)
        ti_tti = np.dot(ti_t,ti)
        invfactor = inv(I-titi_t-ti_tti)

        Tn_factor = np.dot(Tn_factor,ti_t)
        Tnt_factor = np.dot(Tnt_factor,ti)
        ti = np.dot(invfactor,ti_sq)
        ti_t = np.dot(invfactor,ti_t_sq)
        
        Tmat_n =np.dot(Tn_factor,ti)
        Tmat_t_n = np.dot(Tnt_factor,ti_t)

        Tmat = Tmat + Tmat_n
        Tmat_t = Tmat_t + Tmat_t_n

        i = i+1
    
    #Surface GF
    g00=inv(invg0-np.dot(H01,Tmat))

    #Bulk GF
    #gb = inv(invg0 - np.dot(H01,Tmat)-np.dot(H01_d,Tmat_t))
    #gb = inv(invg0)
    #diff = gb - gb.conjugate().T
    #A_omega = np.dot(complex(0,1), diff) #*omega
    #dos = A_omega*omega

    #self energy
    Sigma = np.dot(np.dot(Tdc,g00),Tdc_dag)
    #GF of device
    Gd = inv(omega**2*I2-D_dvc-Sigma)
    Gd_dag = Gd.conjugate().T
    #Spectra
    Diff_Gd = Gd-Gd_dag
    A_omega = np.dot(complex(0,1), Diff_Gd) #*omega
    dos = np.dot(complex(0,1), Diff_Gd)*omega
    Spectra = np.trace(A_omega).real
    Dos_spect = np.trace(dos).real/p

    omega_cm = omega*Ry_cm
    fw_spect.write("%12.6f %20.6e" % (omega_cm,Spectra))
    fw_spect.write('\n')
    fw_dos.write("%12.6f %20.6e" % (omega_cm,Dos_spect))
    fw_dos.write('\n')

###No self-energy: spectra should be zero energywhere as the D_dvc matrix elements are real.
##    Gd_ori = inv(omega**2*I2-D_dvc)
##    print(Gd_ori)
##    Gd_ori_dag = Gd_ori.conjugate().T
##    Diff_Gd_ori = Gd_ori-Gd_ori_dag
##    A_omega_ori = np.dot(complex(0,1), Diff_Gd_ori)
##    Spectra_ori = np.trace(A_omega_ori).real
##
##    fw_spect_ori.write("%12.6f %20.6e" % (omega_cm,Spectra_ori))
##    fw_spect_ori.write('\n')
    
    
    #omega = omega + 0.000005
    #omega = omega + 1
    omega = omega + 0.0000001


fw_spect.close()
fw_dos.close()
























