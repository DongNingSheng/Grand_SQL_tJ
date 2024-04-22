using ITensors
using LinearAlgebra
BLAS.set_num_threads(1)
using ITensors.HDF5
ITensors.Strided.disable_threads()
ITensors.enable_threaded_blocksparse()
using StatsBase: sample

Base.Sys.set_process_title(string("#CORES=",Threads.nthreads()))
M=2000
Ny = 6
Nx = 24
inv_delta=8
mu=5.4
sys=string("N",Ny,"x",Nx,"_mu",mu,"_M",M)
text=string(sys,"_node05")

delta=1/inv_delta
J1=1
t1=3
ratio=0.2
t2=t1*ratio
J2=J1*ratio^2

N = Nx*Ny
Ne=convert(Int64,N*(1-delta))
println()
@show sys
@show text
@show Threads.nthreads()

let
  sites = siteinds("tJ", N;
                   conserve_sz = true, conserve_nf=false)

  lattice = square_lattice_nnn(Nx, Ny; yperiodic = true)

  # Define the Heisenberg spin Hamiltonian on this lattice
  ampo = OpSum()

	for i=1:N
    ampo += -mu,"Ntot",i
	end


  for b in lattice
    if b.type=="1"
      #nn-bond
      ampo .+= J1/2, "S+", b.s1, "S-", b.s2
      ampo .+= J1/2, "S-", b.s1, "S+", b.s2
      ampo .+= J1,   "Sz", b.s1, "Sz", b.s2
      ampo .+= -0.25*J1, "Ntot", b.s1, "Ntot", b.s2
      ampo .+= -t1, "Cdagup", b.s1, "Cup", b.s2
      ampo .+= -t1, "Cdagup", b.s2, "Cup", b.s1
      ampo .+= -t1, "Cdagdn", b.s1, "Cdn", b.s2
      ampo .+= -t1, "Cdagdn", b.s2, "Cdn", b.s1
    else
      #nnn-bond
      ampo .+= J2/2, "S+", b.s1, "S-", b.s2
      ampo .+= J2/2, "S-", b.s1, "S+", b.s2
      ampo .+= J2,   "Sz", b.s1, "Sz", b.s2
      ampo .+= -0.25*J2, "Ntot", b.s1, "Ntot", b.s2
      ampo .+= -t2, "Cdagup", b.s1, "Cup", b.s2
      ampo .+= -t2, "Cdagup", b.s2, "Cup", b.s1
      ampo .+= -t2, "Cdagdn", b.s1, "Cdn", b.s2
      ampo .+= -t2, "Cdagdn", b.s2, "Cdn", b.s1
    end
  end

  H = MPO(ampo,sites)
  H = splitblocks(linkinds, H)

  f=h5open("H_MPO.h5","w")
  write(f,"H",H)
  close(f)

  ### Choose a initial state with randomly distributed electron

  sites_emp=inv_delta:inv_delta:N
	sites_Ne=setdiff(1:N,sites_emp)
	sites_up=sites_Ne[1:2:end]
	sites_dn=sites_Ne[2:2:end]

  state=String[]
  for i in 1:N
  	if i in sites_up
			state = vcat(state,"Up")
		elseif i in sites_dn
			state = vcat(state,"Dn")
		else
			state = vcat(state,"Emp")
		end
	end

  psi0 = randomMPS(sites,state,1000)

  sweeps = Sweeps(30)
  maxdim!(sweeps,1000,1000,1000,1000,1000,1000,1000,1000,1000,2000,2000,2000,2000,2000,2000,2000)
  cutoff!(sweeps,1E-8)
  noise!(sweeps,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-9,1E-9,1E-9,1E-9,0)

  @show sweeps

  file=string("psi_",sys)
  energy,psi = dmrg(H,psi0,sweeps; write_when_maxdim_exceeds=4000, write_step=file)
  mv(string(file,"_tmp.h5"),string(file,".h5");force=true)

  include("Measure_Grand.jl")
	Measure_func(Nx,Ny,sys;Occu=true, Drdu=true, Drud=true, Dadu=true, Dbdu=true)

  return
end
