import numpy as np
import os
import re
from eos import T_qcd_fct, p_in_GeVfm3_from_e_in_GeVfm3, cs2_from_e_in_GeVfm3
import configparser

hbarc=0.1973

############################################
############     Load files     ############
############################################

# Fetch all directories corresponding to parameter points
pre_kompost_output=os.listdir(path="./")
subdir_regex = re.compile('.*tOut([0-9\.]+).music_init_flowNonLinear_pimunuTransverse.txt')
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


final_array_to_save_old_format={key:None for key in eos_list.keys()}
final_array_to_save_new_format={key:None for key in eos_list.keys()}

# Loop over time
# Skip the last time step, which is the same as the hydro initial conditions (avoid doublecounting)
for tau in tau_list[:-1]:

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

    # New format
    deta_dummy=0.1
    volume=float(dx)*float(dy)*deta_dummy*dtau*tau*(zeros+1)
    eta=zeros
    pitautau=tmp_data[:,8]
    pitaux=tmp_data[:,9]
    pitauy=tmp_data[:,10]
    pitaueta=tmp_data[:,11]
    pixx=tmp_data[:,12]
    pixy=tmp_data[:,13]
    pixeta=tmp_data[:,14]
    piyy=tmp_data[:,15]
    piyeta=tmp_data[:,16]
    pietaeta=tmp_data[:,17]
    # Define Pi_b as (epsilon/3-pressure)
    tmp_pressure=p_in_GeVfm3_from_e_in_GeVfm3(e_dens)
    Pi_b=(e_dens/3.-tmp_pressure)/hbarc

    eps_plus_P_in_fm=(e_dens+tmp_pressure)/hbarc
    cs2=cs2_from_e_in_GeVfm3(e_dens)


    # For old format
    pitt = pitautau #pow(cosh(eta),2)*Wtautau + pow(tau*sinh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
    pitx = pitaux #cosh(eta)*Wtaux + tau*sinh(eta)*Wxeta*pow(tau,-1);
    pity = pitauy #cosh(eta)*Wtauy + tau*sinh(eta)*Wyeta*pow(tau,-1);
    pitz = tau*pitaueta #cosh(eta)*sinh(eta)*Wtautau + tau*( pow(cosh(eta),2) + pow(sinh(eta),2) )*Wtaueta*pow(tau,-1) + tau*tau*cosh(eta)*sinh(eta)*Wetaeta*pow(tau,-2);
    pixz = tau*pixeta #sinh(eta)*Wtaux + tau*cosh(eta)*Wxeta*pow(tau,-1);
    piyz = tau*piyeta #sinh(eta)*Wtauy + tau*cosh(eta)*Wyeta*pow(tau,-1);
    pizz = tau*tau*pietaeta #pow(sinh(eta),2)*Wtautau + pow(tau*cosh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);


    for key, item in final_array_to_save_old_format.items():

        tmp_T_fct=eos_list[key]

        # 
        T=tmp_T_fct(e_dens)

        # Old format
        array_to_save_old_ideal=np.transpose([T,zeros,vx,vy,zeros])
        array_to_save_old_shear=np.transpose(np.divide([pitt,pitx,pity,pitz,pixx,pixy,pixz,piyy,piyz,pizz],eps_plus_P_in_fm))
        array_to_save_old_bulk=np.transpose([Pi_b,eps_plus_P_in_fm,cs2])

        # New format
        # For new file format, save only cells with T>Tref
        Tref=0.105
        indices_to_save=T>Tref
        #volume=dx*dy*deta*dtau*tau, eta, T, ux, uy, ueta, pixx, pixy, pixeta, piyy,piyeta, pi_b
        array_to_save_new=np.transpose([volume, eta, T,ux,uy, tau*ueta, pixx/eps_plus_P_in_fm, pixy/eps_plus_P_in_fm, tau*pixeta/eps_plus_P_in_fm, piyy/eps_plus_P_in_fm, tau*piyeta/eps_plus_P_in_fm, Pi_b])
        array_to_save_new_trimmed=array_to_save_new[indices_to_save,:]

        if (item is None):
            final_array_to_save_old_format[key]={}
            final_array_to_save_old_format[key]['ideal']=array_to_save_old_ideal
            final_array_to_save_old_format[key]['shear']=array_to_save_old_shear
            final_array_to_save_old_format[key]['bulk']=array_to_save_old_bulk

            final_array_to_save_new_format[key]=array_to_save_new_trimmed
        else:
            final_array_to_save_old_format[key]['ideal']=np.concatenate((final_array_to_save_old_format[key]['ideal'],array_to_save_old_ideal))
            final_array_to_save_old_format[key]['shear']=np.concatenate((final_array_to_save_old_format[key]['shear'],array_to_save_old_shear))
            final_array_to_save_old_format[key]['bulk']=np.concatenate((final_array_to_save_old_format[key]['bulk'],array_to_save_old_bulk))

            final_array_to_save_new_format[key]=np.concatenate((final_array_to_save_new_format[key],array_to_save_new_trimmed))

# Output to file
for key in final_array_to_save_old_format.keys():


    # Old format
    to_save=final_array_to_save_old_format[key]['ideal']
    filename="evolution_xyeta_eos_"+key+"_old_format.dat"
    np.savetxt(filename, to_save)

    to_save=final_array_to_save_old_format[key]['shear']
    filename="evolution_Wmunu_over_epsilon_plus_P_xyeta_eos_"+key+"_old_format.dat"
    np.savetxt(filename, to_save)

    to_save=final_array_to_save_old_format[key]['bulk']
    filename="evolution_bulk_pressure_xyeta_eos_"+key+"_old_format.dat"
    np.savetxt(filename, to_save)

    # New format
    to_save=final_array_to_save_new_format[key]

    filename="evolution_xyeta_eos_"+key+"_new_format.dat"

    np.savetxt(filename, final_array_to_save_new_format[key])


###########################################################
############ Write header file with useful info ###########
###########################################################


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
f.write("const bool MUSIC_with_shear_viscosity = true;\n")
f.write("const bool MUSIC_with_bulk_viscosity = true;\n")
f.write("const bool MUSIC_with_rhob = false;\n")
f.write("const bool MUSIC_with_diffusion = false;\n")
f.write("const bool MUSIC_outputBinaryEvolution=false;\n")
f.close()
