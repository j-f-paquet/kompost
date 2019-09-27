import numpy as np
import os
import re
from eos import T_qcd_fct
import configparser

hbarc=0.1973

# Fetch all directories corresponding to parameter points
pre_kompost_output=os.listdir(path="./")
subdir_regex = re.compile('tIn02.tOut([0-9\.]+).music_init_flowNonLinear_pimunuTransverse.txt')
#kompost_output=[(subdir,match.group(1)) for subdir in pre_input_subdirs if match := subdir_regex.match(subdir)]

kompost_output={}
for tmp_file in pre_kompost_output:

    match=subdir_regex.match(tmp_file)
    if (match != None):
        tau=float(match.group(1))
        kompost_output[tau]=tmp_file
        

tau_list=np.sort(np.array(list(kompost_output.keys())))
tau_min=tau_list[0]
tau_max=tau_list[-1]
dtau=tau_list[1]-tau_list[0]

#print(tau_list)
#print(tau_min)
#print(tau_max)
#print(dtau)

############################################
############ Equations of state ############
############################################

# T_{ideal}(e) with e in GeV/fm3 and T in GeV? 
def T_ideal(e_in_GeV_over_fm3,Nf=0):

    Nc=3

    pi=np.pi

    e_in_one_over_fm4=e_in_GeV_over_fm3/hbarc

    T_in_fm=np.power(90.0/(pi*pi)*(e_in_one_over_fm4/3.0)/(2.*(Nc*Nc-1)+7./2*Nc*Nf), .25)

    return T_in_fm*hbarc



# T(e) with e in GeV/fm3 and T in GeV? 
def T_qcd(e_in_GeV_over_fm3): 

    #return T_qcd_fct(e_in_GeV_over_fm3)
    return T_qcd_fct(e_in_GeV_over_fm3)

eos_list={
    'qcd':T_qcd, 
    'ideal_gluon':lambda e: T_ideal(e,Nf=0),
    'ideal_Nf3':lambda e: T_ideal(e,Nf=3)
}



########################################################################################
############ Load KoMPoST output and save in format used to compute photons ############ 
########################################################################################


final_array_to_save={key:None for key in eos_list.keys()}

# Loop over time
for tau in tau_list:

    filename=kompost_output[tau]

    #//profile >> dummy1 >> dummy2 >> dummy3 >> density >> utau >> ux >> uy  >> dummy  >> dummy  >> dummy  >> dummy;
    #//energy density should be in GeV/fm^3
    #//utau should be in fm^-1 (T^{\tau/x/y} in fm^-4, T^{\tau/x/y \eta} in fm^-5, T^{\eta\eta} in fm^-6)
    #outfile << "0 " << xpos << " " << ypos << " " << energy_density/(M_HBARC*M_HBARC*M_HBARC) << " " << flow[0] << " " << flow[1] << " " << flow[2] << " " << flow[3]/M_HBARC << " ";
    #//\pi^{\mu\nu} should be in fermi
    #//order: pitautau >> pitaux >> pitauy >> pitaueta >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta
    tmp_data=np.loadtxt(filename)
    e_dens=tmp_data[:,3]
    utau=tmp_data[:,4]
    ux=tmp_data[:,5]
    uy=tmp_data[:,6]
    ueta=tmp_data[:,7]

    vx=np.divide(ux,utau)
    vy=np.divide(uy,utau)

    zeros=np.zeros_like(vx)

    for key, item in final_array_to_save.items():

        tmp_T_fct=eos_list[key]

        # 
        T=tmp_T_fct(e_dens)

        #final_array_to_save=np.concatenate(final_array_to_save,array_to_save)
        array_to_save=np.transpose([T,zeros,vx,vy,zeros])

        if (item is None):
            final_array_to_save[key]=array_to_save
        else:
            final_array_to_save[key]=np.concatenate((final_array_to_save[key],array_to_save))

# Output to file
for key in final_array_to_save.keys():

    filename="evolution_xyeta_eos_"+key+".dat"

    np.savetxt(filename, final_array_to_save[key])


###########################################################
############ Write header file with useful info ###########
###########################################################

# Read header of kompost output to figure out what is the grid size et al
filename=kompost_output[tau_min]
f = open(filename, "r")
kompost_output_header = f.readline()
f.close()

subdir_regex = re.compile('# tau_in_fm ([0-9\.]+) etamax= 1 xmax= ([0-9]+) ymax= ([0-9]+) deta= 0 dx= ([0-9\.]+) dy= ([0-9\.]+)')
match=subdir_regex.match(kompost_output_header)
tau0=match.group(1)
Nx=match.group(2)
Ny=match.group(3)
dx=match.group(4)
dy=match.group(5)

filename="hydro_info_header_h"
f = open(filename, "w")
f.write("const int MUSIC_real_nx ="+str(Nx)+";\n")
f.write("const int MUSIC_real_ny ="+str(Ny)+";\n")
f.write("const int MUSIC_real_neta =1;\n")
f.write("const double MUSIC_tau0 = "+str(tau0)+";\n")
f.write("const double MUSIC_dx = "+str(dx)+";\n")
f.write("const double MUSIC_dy = "+str(dy)+";\n")
f.write("const double MUSIC_deta = "+str(0)+";\n")
f.write("const double MUSIC_dtau = "+str(dtau)+";\n")
f.write("const bool MUSIC_with_shear_viscosity = false;\n")
f.write("const bool MUSIC_with_bulk_viscosity = false;\n")
f.write("const bool MUSIC_with_rhob = false;\n")
f.write("const bool MUSIC_with_diffusion = false;\n")
f.write("const bool MUSIC_outputBinaryEvolution=false;\n")
f.close()
