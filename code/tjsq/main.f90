!=================================================
!For the t-J model on the square/triangle Lattice
!=================================================
program Quantum
	use pubdata
	implicit none
	real(8),parameter :: trun_err=1.0E-18
	integer :: i,j,x,y,trun_idx,levels
	real(8) :: bg_time,end_time
        integer delta
        open(1212,file='restart.dat')
        read(1212,*)restart1,sys_len_start
        read(1212,*)restart2
        close(1212)        
        if(sys_len_start==0)restart2=0
!!! change parameters here
        t2=0.20d00   !this is t2/t1       
                iter1=0
                if(restart1.ne.0)iter1=1  
                                                 xi1=ny*4+1
                                                xj1=ny*4+1
                        !! adjustable, for bulk measurements, cut-off boundary
        kept_max=600 
		sweep_num=6

		kept_min=min(3*kept_max/4, 3600)

                kept_maxf=1200

        if(restart1==0)then
                kept_min=min(1000, kept_max*3/4)
                kept_max=min(1200, kept_max)
                        endif

        tot_num_up=0
                jd=0
                jz=0
                jn=0
                jt=0
        jz(1:4, 1:num_site)=1.0d0
        jz(5:8, 1:num_site)=t2**2
        jd(1:18,1:num_site)=jz(1:18,1:num_site)*(-dsqrt(3.0d0))
        jt(1:4,1:num_site )=-3.0d0*dsqrt(2.0d0)
        jt(5:8, 1:num_site)=-3.0d0*dsqrt(2.0d0)*t2
        jn(1:12)=-jz(1:12,1)/4.0d00
!!! these are parameters for t-J model
        esite=-4.0d0  !! chemical potential


        call dmrg_setup0()
        call dmrg_setup(num_site/2-1, num_site/2-1)
        call dmrg_setup0()
	call Get_Lattice
	call Get_site_operator_su2_U0
	call Get_single_site_basis
	call Get_model_name(sysname,envname)
	call Get_wave_name
!	call Read_parameter
        iw6j1=0
        w6j1=0.d0
        iw6j2=0
        w6j2=0.d0
        iw6j3=0
        w6j3=0.d0
                                  call factrl
                                 !! for Wigner 3j 6j and 9j coefs^M
        infi_del=0
        cone=1.0d0
        czero=0.0d0
        openy=0

	!Start program
        pinz=1
        pind=1

        densityc=.true.
        superc=.true.
        meas_line=.true.

        trun_idx=1
		tot_num_down=0  

		!Get ground state
        if(restart1.ne.4)then
		levels=1
		call Get_Ground_State(trun_idx,trun_err,levels)
                endif

		!<2>: Measurement
		call Measurement(trun_idx)
             !	end do
end program Quantum


!Get ground_state wavefunction
subroutine Get_Ground_State(trun_idx,trun_err,levels)
	use pubdata
	implicit none
	real(8),intent(in) :: trun_err
	integer,intent(in) :: trun_idx,levels
        logical :: newflag
	integer :: i,iter,tmp_keep,point, idx1
	real(8) :: bg_time,end_time,tmp_err
	real(8),parameter :: lanerr=20.0e-6

	!<1>: For the Warm_up process
        iter=iter1
	call cpu_time(bg_time)
	open(7,file="Energy_All.dat",position='append')
	write(7,*) "Warm_up ", 't2=', t2, 'esite', esite(1),esite(Ny+1)
	write(7,*)"jz, jt",jz(1,1),jz(5,1),jt(1,1),jt(5,1)
	write(7,*) "kept",kept_min, kept_max, iter
	close(7)

	spectrum=.false.
	lanzero=lanerr*(10**3)
	tmp_err=trun_err*(200)
        if(restart1.eq.0)then
	call warm_up_point(trun_idx,Num_site/2,tmp_err,levels)
                endif

        if(restart1.ne.0)then
	lanzero=lanerr*2 
	tmp_err=trun_err*2 
                endif

        if(restart2.ne.0)then
	lanzero=lanerr 
	tmp_err=trun_err 
                endif
        if(iter.eq.0)then
	lanzero=lanerr*(30)
	tmp_err=trun_err*(15)
                endif

          if(restart1.eq.0)kept_max=kept_max*2/3+kept_maxf/3
		kept_min=min(3*kept_max/4, 3600)

        idx1=trun_idx
        if(restart1.le.1)then
        if(restart1.eq.0)sys_len_start=num_site/2-1
        restart1=restart2
call left_to_right(sys_len_start,Num_site-idx1-1,tmp_err,levels,.true.)
        restart1=1
        restart2=1
                endif

          if(restart1.eq.0)kept_max=kept_max*2/3+kept_maxf/3
		kept_min=min(3*kept_max/4, 3600)

                if(restart1.le.2)then
        if(restart1.le.1)sys_len_start=num_site-idx1
        restart1=restart2
	call right_to_left(sys_len_start-1,idx1,tmp_err,levels,.true.)
                        endif
        restart1=1

                kept_max=kept_maxf
		kept_min=min(3*kept_max/4, 3600)
	iter=iter+1
	do while(iter<=sweep_num)

	open(7,file="Energy_All.dat",position='append')
		write(7,*) "Sweep=",iter, ' t2=', t2
	write(7,*)"jz, jt",jz(1,1),jz(5,1),jt(1,1),jt(5,1)
	write(7,*) "kept",kept_min, kept_max, iter
	close(7)

	lanzero=lanerr
	tmp_err=trun_err

	call left_to_right(idx1,Num_site-idx1-1,tmp_err,levels,newflag)
	call right_to_left(Num_site-idx1-1,idx1,tmp_err,levels,newflag)

		call cpu_time(end_time)
		open(7,file="Energy_All.dat",position='append')
		write(7,"(A6,F12.4,2X,A3)") "Time =",end_time-bg_time,"Sec"
		write(7,*)
		close(7)

		iter=iter+1
	end do
	spectrum=.true.
	lanzero=lanerr 
	tmp_err=trun_err 
	call left_to_right(idx1,Num_site0/2+3,tmp_err,levels,newflag)
end subroutine Get_Ground_State


!Perform measurement
subroutine Measurement(trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	real(8) :: bg_time,end_time,value
	integer :: i,j,i1, j1,label,sys_len,env_len,x,y

	integer :: name_len,knum
	character(len=30) :: Filename
	type(Total_Model) :: sys,env
	type(Total_Basis) :: sys_bs,env_bs
	type(Wavefunction) :: eigvec
        character ca
        integer xc1, xc2

	real(8),dimension(Num_site) :: ave_spin,ave_spin2
	real(8),dimension(Num_site) :: ave_num_elec_up,ave_num_elec_down,ave_spin_sd2,ave_spin_sd21
	real(8),dimension(Num_site) :: ave_num_elec,ave_num_hole,ave_num_elec2,ave_num_elec_up2
	real(8),dimension(Num_site) :: ave_spin_sz,ave_spin_sz2,ave_num_elec_down2
	real(8),dimension(Num_site,Num_site) :: num_elec_cor,spin_sz_cor,spin_sd_cor,num_hole_cor
	real(8),dimension(Num_site,Num_site) :: num_elec_up_cor,num_elec_down_cor
	real(8),dimension(Num_site,Num_site) :: elec_cor,elec_down_cor,spin_cor
        integer re_meas
        re_meas=1


	!<1>: ============ Read data from disk =================
	!<1-1>: Read eigen_state from disk (Ground state only)
	call wave_from_disk(eigvec,1001,wave_name(1),7)
	sys_len=eigvec%sys_len
	env_len=eigvec%env_len

	!<1-2>: Read truncation operators from disk
	do i=trun_idx+1,sys_len
		call truns_from_disk(systruns(i),1001,i,.true.)
	end do

	do i=trun_idx+1,env_len
		call truns_from_disk(envtruns(i),1001,i,.false.)
	end do

	!<1-3>: Read basis from disk
	do i=2,sys_len
		call basis_from_disk(sys_bsm(i),1001,i,.true.)
	end do

	do i=2,env_len
		call basis_from_disk(env_bsm(i),1001,i,.false.)
	end do


        !! odd chain
   call Get_SC_Order(eigvec,trun_idx,'S0',"Sc_order.dat",1,num_site-ny)
        !! even chain
   call Get_SC_Order(eigvec,trun_idx,'S0',"Sc_order.dat",2,num_site-ny)

        if(densityc)then
	Filename="num_elec.dat"
     call measure_operator_dia(ave_num_elec,num_elec,eigvec,trun_idx,Filename)
        endif


        xc1=(Nx/6)*ny-ny+1
        xc2=xc1
        if(superc)then
   call Get_mSCOP_Operator_Cor0(eigvec,trun_idx,'S0',"Sc_cor_y.dat",xc1,xc1,1,0)
   call Get_mSCOP_Operator_Cor0(eigvec,trun_idx,'S0',"Sc_cor_x.dat",xc1,xc1,3,0)
        endif


24      continue
        xc1=(Nx/6)*ny-ny+1
        xc2=xc1

        spin_cor=0.0d0
        num_elec_cor=0.0d0
        elec_cor=0.0d0

28      format(3i8, 4f21.14)
281      format(4i8, 4f21.14)
        if(meas_line)then
        filename='density_cor.back'
        call Get_density_cor(num_elec_cor, num_elec,num_elec,eigvec,trun_idx,'B',Filename, xc1,xc1,0)
        open(251,file='density_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(num_elec_cor(i,j).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,num_elec_cor(i, j), ave_num_elec(i), ave_num_elec(j)
        endif
        enddo
        enddo

        filename='spin_cor.back'
        call Get_density_cor(spin_cor,st_sd, st_sd,eigvec,trun_idx,'B',Filename, xc1,xc2,0)

        open(251,file='spin_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(num_elec_cor(i,j).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,spin_cor(i, j) 
        endif
        enddo
        enddo

        filename='elec_cor.back'
        call Get_density_cor(elec_cor,st_elec_up, st_elec_down,eigvec,trun_idx,'F',Filename, xc1,xc1,0)

        open(251,file='elec_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(elec_cor(i,j).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,elec_cor(i, j) !!, ave_num_elec(i), ave_num_elec(j)
        endif
        enddo
        enddo
        endif



	!<6>: ========== Free space ===========
	do i=trun_idx+1,sys_len
		call deallocate_block(systruns(i))
	end do

	do i=trun_idx+1,env_len
		call deallocate_block(envtruns(i))
	end do

	do i=2,sys_len
		call deallocate_basis(sys_bsm(i))
	end do

	do i=2,env_len
		call deallocate_basis(env_bsm(i))
	end do
	call deallocate_wave(eigvec)

end subroutine Measurement

