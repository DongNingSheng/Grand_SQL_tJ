
The numerical code calculates the ground state of square lattice t1-t2-J1-J2 model with fixed chemical potential, see F Chen, FDM Haldane, DN Sheng, arXiv e-prints, arXiv: 2311.15092 for details.

The lattice size Num_site=Nx*Ny, and lengths parameters (Nx, Ny) can be found in common.f90.  The lattice geometry can be determined as square=1,  or triangle=1.   For different system sizes, one needs to change Nx, Ny accordingly.
Notice that Ny=6 (6-leg) can be down in 1-2 weeks if keeps 3000-6000 states.  Ny=8 usually takes much longer time.

The code main.f90 has more parameters you can adjust:

From line 23,  kept_min, and kept_max are the range of bond dimensions.
Set kept_max as initial disired bond dimensions.   kept_maxf as final bond
dimension. Sweep_num usually needs to be around 6-10+ (more sweeps give better convergence). For example, one can set:

kept_min=1600
kept_max=2000
Kept_maxf=4000
sweep_num=8
you job will be done with bond dimension 4000, after eight complete sweeps.

For triangular lattice:
jz(1:6, 1:num_site)=J1     nearest neighbor (NN) Heisenberg coupling
jz(7:12, 1:num_site)=J2   next nearest neighbor (NNN) Heisenberg coupling
jt(1:6, 1:num_site)=-3*dsqrt(2.0) NN hopping, the  factor sqrt(2.0) is due to 
                                    SU2 code (which adjusts coupling constant). 
jt(7:12, 1:num_site)=-3*dsqrt(2.0)*t2   NNN hopping,  the factor sqrt(2.0) is due to SU2 code.
jn is the coupling coef. for n_in_j terms,  we have jn=-0.25*jz
Here t2 is the parameter as the ratio of  NN_hopping/NNN_hopping.

For square lattice
jz(1:4, 1:num_site)=J1     nearest neighbor (NN) Heisenberg coupling
jz(5:8, 1:num_site)=J2   next nearest neighbor (NNN) Heisenberg coupling
jt(1:4, 1:num_site)=-3*dsqrt(2.0) NN hopping, the  factor sqrt(2.0) is due to 
                                    SU2 code (which adjusts coupling constant). 
jt(5:8, 1:num_site)=-3*dsqrt(2.0)*t2   NNN hopping,  the factor sqrt(2.0) is due to SU2 code.
jn is the coupling coef. for n_in_j terms,  we have jn=-0.25*jz

For this code with fixing chemical potential (esite),  we usually needs to start from a smaller chemical potential, 
to reduce particle numbers.   If we find the num_elec.dat with electron density around 1.00,  then we need to reduce
mu to smaller number. 

In the code esite=-\mu (chemical potential), which need to be justified therough testing.

To start the job, you do 
make 
which generates the compiled program tjsq.out
The intel new compiler ifort version 2021.5.0 is used.   For other version, one needs to do adjustment for the makefile.

For initial job,  you need to edit restart.dat to be in the form
0 0
0

The useful files and measurement results
 density_cor.dat   density correlations 
      3rd   4th      5th      6th  (columns)
    r_ji   <n_in_j>  <n_i>  <n_j>
 elec_cor.dat       hopping corr.      use 3rd and 4th columns
 spin_cor.dat       spin correlations,  use 3rd and 4th columns  
 num_elec.dat    onsite electron density,    x, y, n(x,y) format      
 Energy_All.dat      total energy
 entropy1.dat        X  entropy ,    X is subsystem length
 Kept_State.dat   actual bond dimension during the DMRG
 restart.dat       for initial of job, and continuation of the job.
To continue,   you may need to set kept_max  and kept_maxf to larger values.

 Sc_cor_x.dat
 Sc_cor_y.dat
	above are pairing correlations for measurements of different bonds.
	The first bond direction is indicated by the name of the files.
	data format  for SC_cor_x.dat
	4th   5th     6th    
         r    |Pxy|, |Pxx| 
	r is distance,  P is pairing corre .
	For other files, replacing x by y
 Sc_cor_y.dat
	4th   5th     6th    
         r    |Pyy|, |Pyx| 
For SC order, use the 3rd and 4th columns for x-bond pair and y-bond pair.
The position of the pair is shown in 1st column with an integer
which equals to  y+(x-1)*Ny

common.f90  have commonly used parameters and data structures 

main.f90 have parameters for the job and perform simulations for ground state and measurements.

dmrg_sub.f90  have subroutines for setting up the dmrg, set up neighboring coupings. It also contains the dmrg codes for doping infinite process (warm_up), left_to_right and  right_to_left sweepings.

den_mat find the eigenstates of reduced density matrix.

eigen.f90 perform the lanczos calculations for ground state.

get_sys update the dmrg blocks (including operators), perform truncations
 and also operator*vectors for lanczos.

measure.f90 perform the measurements

The DMRG method follows Steve White's original algorithm with including SU2 spin rotational symmetry.    The code was built based on original U1 code of Hong Chen Jiang (2008) and detailed note of Shoushu Gong (2020) for SU2 implimentation.   Many thanks to Hong Chen and Shoushu!

Please contact donna.sheng44@gmail.com (or donna.sheng@csun.edu) for more information.
