function Measure_func(Nx::Int, Ny::Int, sys::String ; kwargs...)

	N=Nx*Ny
	ini=convert(Int64,Nx/4-1)
  Occu = get(kwargs, :Occu, false)
  Occu_M = get(kwargs, :Occu_M, false)
  Sz = get(kwargs, :Sz, false)
  Sz_M = get(kwargs, :Sz_M, false)
  Szz = get(kwargs, :Szz, false)
  Dbud = get(kwargs, :Dbud, false)
  Dbdu = get(kwargs, :Dbud, false)
  Drud = get(kwargs, :Drud, false)
  Drdu = get(kwargs, :Drdu, false)

	shift = get(kwargs, :shift, false)
  Dbud_shift = get(kwargs, :Dbud_shift, false)
  Dbdu_shift = get(kwargs, :Dbud_shift, false)
  Daud_shift = get(kwargs, :Daud_shift, false)
  Dadu_shift = get(kwargs, :Daud_shift, false)

  Dcud = get(kwargs, :Dcud, false)
  Dcdu = get(kwargs, :Dcdu, false)
  Dadu = get(kwargs, :Dadu, false)
  Daud = get(kwargs, :Daud, false)
  Pxx = get(kwargs, :Pxx, false)
  Pyy = get(kwargs, :Pyy, false)

  Phi_4e_p = get(kwargs, :Phi_4e_p, false)
  Phi_4e_m = get(kwargs, :Phi_4e_m, false)

  f = h5open(string("psi_",sys,".h5"),"r")
  psi=read(f,"psi",MPS)
  close(f)

  """
  ##############################################################################################################################
  Occu_M
  """
  if Occu_M
		Ne = expect(psi,"Ntot";sites = 1:N)
		open(string("Ne_M_",sys,".dat"), "w") do io
		    for i in 1:Ny
		    	for j in 1:Nx
		        println(io, i,"\t",j,"\t",real(Ne[(j-1)*Ny+i]))
		      end
		    end
		end
  end
  	
  """
  ##############################################################################################################################
  Occu
  """
  if Occu
		Ne = expect(psi,"Ntot";sites = 1:Ny:N)
		open(string("Ne_",sys,".dat"), "w") do io
		    for i in 1:1
		    	for j in 1:Nx
		        println(io, i,"\t",j,"\t",real(Ne[j]))
		      end
		    end
		end
#		FIG=plot(1:Nx,Ne);
#		savefig(FIG,string("Ne_",sys,".pdf"))
  end
  	
  """
  ##############################################################################################################################
  Sz
  """
	if Sz
  	Sz = expect(psi,"Sz";sites = 1:Ny:N)
  	open(string("Sz_",sys,".dat"), "w") do io
  	    for i in 1:1
  	    	for j in 1:Nx
  	        println(io, i,"\t",j,"\t",real(Sz[j]))
  	      end
  	    end
  	end
	end

  """
  ##############################################################################################################################
  Sz_M
  """
	if Sz_M
  	Sz_M = expect(psi,"Sz";sites = 1:N)
  	open(string("Sz_M_",sys,".dat"), "w") do io
  	    println(io, "y","\t","x","\t","Sz")
  	    for i in 1:Ny
  	    	for j in 1:Nx
  	        println(io, i,"\t",j,"\t",real(Sz_M[(j-1)*Ny+i]))
  	      end
  	    end
  	end
	end

  """
  ##############################################################################################################################
  Dcdu
  """
	if Dcdu
  	open(string("Dcdu_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    range_x= collect(0:(Nx-2)) .* Ny .+ 1
  	    Dc = Op_nn(psi,"Cdn","Cup";sites=range_x, x=Ny+1)
				for (ni,i) in enumerate(range_x)
  	    	println(io, i,"\t",i+Ny,"\t",real(Dc[ni]),"\t",imag(Dc[ni]))
				end
  	end
	end

  """
  ##############################################################################################################################
  Dcud
  """
	if Dcud
  	open(string("Dcud_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    range_x= collect(0:(Nx-2)) .* Ny .+ 1
  	    Dc = Op_nn(psi,"Cup","Cdn";sites=range_x, x=Ny+1)
				for (ni,i) in enumerate(range_x)
  	    	println(io, i,"\t",i+Ny,"\t",real(Dc[ni]),"\t",imag(Dc[ni]))
				end
  	end
	end

  """
  ##############################################################################################################################
  Dbdu_shift
  """
	if Dbdu_shift
		for si in shift
  		open(string("Dbdu_",sys,"_shift_",si,".dat"), "w") do io
  		    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  		    range_x= collect(0:(Nx-1)) .* Ny .+ si
  		    Db = Op_nn(psi,"Cdn","Cup";sites=range_x, x=1)
					for (ni,i) in enumerate(range_x)
  		    	println(io, i,"\t",i+Ny,"\t",real(Db[ni]),"\t",imag(Db[ni]))
					end
  		end
		end
	end

  """
  ##############################################################################################################################
  Dbud_shift
  """
	if Dbud_shift
		for si in shift
  		open(string("Dbud_",sys,"_shift_",si,".dat"), "w") do io
  		    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  		    range_x= collect(0:(Nx-1)) .* Ny .+ si
  		    Db = Op_nn(psi,"Cup","Cdn";sites=range_x, x=1)
					for (ni,i) in enumerate(range_x)
  		    	println(io, i,"\t",i+Ny,"\t",real(Db[ni]),"\t",imag(Db[ni]))
					end
  		end
		end
	end

  """
  ##############################################################################################################################
  Dbdu
  """
	if Dbdu
  	open(string("Dbdu_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    range_x= collect(0:(Nx-1)) .* Ny .+ 1
  	    Db = Op_nn(psi,"Cdn","Cup";sites=range_x, x=1)
				for (ni,i) in enumerate(range_x)
  	    	println(io, i,"\t",i+Ny,"\t",real(Db[ni]),"\t",imag(Db[ni]))
				end
  	end
	end

  """
  ##############################################################################################################################
  Dbud
  """
	if Dbud
  	open(string("Dbud_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    range_x= collect(0:(Nx-1)) .* Ny .+ 1
  	    Db = Op_nn(psi,"Cup","Cdn";sites=range_x, x=1)
				for (ni,i) in enumerate(range_x)
  	    	println(io, i,"\t",i+Ny,"\t",real(Db[ni]),"\t",imag(Db[ni]))
				end
  	end
	end

  """
  ##############################################################################################################################
  Dadu
  """
	if Dadu
  	open(string("Dadu_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    range_x= collect(0:(Nx-2)) .* Ny .+ 1
  	    Da = Op_nn(psi,"Cdn","Cup";sites=range_x, x=Ny)
				for (ni,i) in enumerate(range_x)
  	    	println(io, i,"\t",i+Ny,"\t",real(Da[ni]),"\t",imag(Da[ni]))
				end
  	end
	end

  """
  ##############################################################################################################################
  Daud
  """
	if Daud
  	open(string("Daud_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    range_x= collect(0:(Nx-2)) .* Ny .+ 1
  	    Da = Op_nn(psi,"Cup","Cdn";sites=range_x, x=Ny)
				for (ni,i) in enumerate(range_x)
  	    	println(io, i,"\t",i+Ny,"\t",real(Da[ni]),"\t",imag(Da[ni]))
				end
  	end
	end
	
  """
  ##############################################################################################################################
  Dadu_shift
  """
	if Dadu_shift
		for si in shift
  		open(string("Dadu_",sys,"_shift=",si,".dat"), "w") do io
  		    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  		    range_x= collect(0:(Nx-2)) .* Ny .+ si
  		    Da = Op_nn(psi,"Cdn","Cup";sites=range_x, x=Ny)
					for (ni,i) in enumerate(range_x)
  		    	println(io, i,"\t",i+Ny,"\t",real(Da[ni]),"\t",imag(Da[ni]))
					end
  		end
		end
	end

  """
  ##############################################################################################################################
  Daud_shift
  """
	if Daud_shift
		for si in shift
  		open(string("Daud_",sys,"_shift=",si,".dat"), "w") do io
  		    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  		    range_x= collect(0:(Nx-2)) .* Ny .+ si
  		    Da = Op_nn(psi,"Cup","Cdn";sites=range_x, x=Ny)
					for (ni,i) in enumerate(range_x)
  		    	println(io, i,"\t",i+Ny,"\t",real(Da[ni]),"\t",imag(Da[ni]))
					end
  		end
		end
	end

  """
  ##############################################################################################################################
  Pyy
  """
	if Pyy
  	open(string("Pyy_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-4)).*Ny .+ j
  	        bl=1
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = Pr1 - Pr2 + Pr3 - Pr4
  	    
  	       for i in 1:length(range_y)
  	           println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	       end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Pxx
  """
	if Pxx
  	open(string("Pxx_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+2):(Nx-4)).*Ny .+ j
  	        bl=Ny
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = Pr1 - Pr2 + Pr3 - Pr4
  	    
  	       for i in 1:length(range_y)
  	           println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	       end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Szz
  """
	if Szz
  	open(string("Szz_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-1)).*Ny .+ j
  	        Sr = correlation_vector_1N(psi,"Sz","Sz";site_x = range_x, sites_y = range_y)
  	
  	        for i in 1:length(range_y)
  	            println(io, range_x,"\t",range_y[i],"\t",real(Sr[i]),"\t",imag(Sr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Phi_4e_p
  """
	if Phi_4e_p
  	open(string("Phi_4e_p_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=ini:Nx-4
  	        range_x= j*Ny+1
  	        range_y= (j+1)*Ny + 1
  	        bl=1
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagup","Cdagdn","Cdagdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagdn","Cdagdn","Cdagup","Cdagup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagup","Cdagdn","Cdagup","Cdagdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdagdn","Cdagup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = Pr1 - Pr2 - Pr3 - Pr4
  	    
  	       for i in 1:length(range_y)
  	           println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	       end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Phi_4e_m
  """
	if Phi_4e_m
  	open(string("Phi_4e_m_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=ini:Nx-4
  	        range_x= j*Ny+1
  	        range_y= (j+1)*Ny + 1
  	        bl=1
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdagdn","Cdagup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagdn","Cdagup","Cdagup","Cdagdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagup","Cdagup","Cdagdn","Cdagdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagdn","Cdagup","Cdagup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr5 = correlation_4p(psi,"Cdagup","Cdagdn","Cdagup","Cdagdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr6 = correlation_4p(psi,"Cdagdn","Cdagup","Cdagdn","Cdagup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = 2*Pr1 + 2*Pr2 - Pr3 - Pr4 - Pr5 - Pr6
  	    
  	       for i in 1:length(range_y)
  	           println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	       end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Drdu
  """
	if Drdu
  	open(string("Drdu_",sys,".dat"), "w") do io
  	    range_x= 1:N
  	    Dr = correlation_matrix(psi,"Cdn","Cup";sites=range_x)
				writedlm(io,Dr)
  	end
	end

  """
  ##############################################################################################################################
  Drud
  """
	if Drud
  	open(string("Drud_",sys,".dat"), "w") do io
  	    range_x= 1:N
  	    Dr = correlation_matrix(psi,"Cup","Cdn";sites=range_x)
				writedlm(io,Dr)
  	end
	end
	
end
