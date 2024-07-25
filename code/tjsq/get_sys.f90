!==========================================================
!Get system hamiltonian from block and one site
!==========================================================
subroutine Get_Hamiltonian_Sys(pre,next,basis, trun)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: basis
	type(Total_Model),intent(inout) :: pre
	type(Total_Model),intent(inout) :: next
	type(Total_Block),intent(in) :: trun

	integer :: i,x,y,old_xidx,old_yidx,new_xidx,new_yidx
        integer :: i1, j, j1,ij, k, nleg1,ii
	real(8) :: coefx,coefy,coefsd,coefsz,coefsn
	type(Total_Block) :: tmp_oper, tmp_oper1
        type(Block_set) tmp_sys 


	!<1>: Update the Hamiltonian
	next%len=pre%len+1
	call update_block_dia(pre%ham,basis,next%ham)
	!<2>: Update all the operators
        do ij=1, 5,2
                j=pre%len+1
                        do k=1, neibt
                I=neib(J,k)
          if(I.le.pre%len.and.I.gt.0)then!!! coupling of this site with system
                        i1=inb1(i, pre%len) 
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, j)
	coefy = 0
	coefsd= jd(k, j)
	coefsz= jz(k, j)
	coefsn=  jn(k)
                if(ij==1)then
                if(coefx.ne.0.0)then
	call Get_hamiltonian_ndia(pre%sub(i1)%elec_up,st_elec_up, basis,coefx,next%ham,'F','S')
                endif
                endif
                if(ij==2)then
                if(coefx.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%elec_down,st_elec_down,basis,coefx,next%ham,'F','S')
                endif
                endif
                if(ij==3)then
                if(coefsd.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%spin_sd,st_sd, basis,coefsd,next%ham, 'B','S')
                endif
                endif
	endif
	endif
        enddo  

	do y=1,pre%len  
		!Update center-site operators
                I1=INB1(y, pre%len)
                if(I1.le.nleg11(pre%len))then 

                if(ij==1)then
                        call update_block_ndia_sys(pre%sub(i1)%elec_up,basis,next%sub(i1)%elec_up,'F')
                        call update_block_ndia_sys(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')

                        call deallocate_block(pre%sub(i1)%elec_up)
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==2)then
                        call update_block_ndia_sys(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==3)then
                        call update_block_ndia_sys(pre%sub(i1)%spin_sd,basis,next%sub(i1)%spin_sd,'B')
                        call deallocate_block(pre%sub(i1)%spin_sd)
                        endif
                        if(ij==4)then
                        call update_block_dia(pre%sub(i1)%spin_sz,basis,next%sub(i1)%spin_sz)
                        call deallocate_block(pre%sub(i1)%spin_sz)
                        endif
                        if(ij==5)then
                        call update_block_dia(pre%sub(i1)%num_sn,basis,next%sub(i1)%num_sn)
                        call deallocate_block(pre%sub(i1)%num_sn)
                        endif

		        endif
                        enddo

                if(ij==1)then
          call update_site_ndia_sys(st_elec_up,basis,tmp_oper,'F')
          call update_site_ndia_sys(st_elec_down,basis,tmp_oper1,'F')
                endif
                if(ij==2)then
          call update_site_ndia_sys(st_elec_down,basis,tmp_oper,'F')
                endif
        if(ij==3)then
          call update_site_ndia_sys(st_sd,basis,tmp_oper,'B')
                       ! call block_to_disk1(tmp_oper, 105)
                endif
                if(ij==4)then
          call update_site_dia(st_sz,basis,tmp_oper)
                endif
                if(ij==5)then
          call update_site_dia(num_elec,basis,tmp_oper)
                endif


                j=pre%len+1
                        do k=1, neibt
                I=neib(J,k)
          if(I.le.pre%len.and.I.gt.0)then
                        i1=inb1(i, pre%len) 
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, j)
	coefy = 0
	coefsd= jd(k, j)
	coefsz= jz(k, j)
	coefsn=  jn(k)
                if(ij==4)then
		call Get_hamiltonian_dia(next%sub(i1)%spin_sz,'N', tmp_oper,'N',coefsz,next%ham,.false.)
                endif
                if(ij==5)then
		call Get_hamiltonian_dia(next%sub(i1)%num_sn,'N', tmp_oper,'N',coefsn,next%ham,.false.)
                endif
	endif
	endif
        enddo 
!!! reorder the operators

         if(ij==5)then
                if(esite(pre%len+1).ne.0.0)then
call block_add_block_two(tmp_oper,esite(pre%len+1),next%ham)
                endif
                endif



                  nleg1=nleg11(pre%len+1)
         do I=1,pre%len   
                I1=INB1(I,pre%len) 
                J = INB1(I,pre%len+1)
                        if(J.le.nleg1.and.j.ne.i1)then

                if(ij==1)then
        call block_transfer(next%sub(i1)%elec_up,next%sub(j)%elec_up)
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==3)then
        call block_transfer(next%sub(i1)%spin_sd,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(next%sub(i1)%spin_sz,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(next%sub(i1)%num_sn,next%sub(j)%num_sn)
                endif
                                endif
                                enddo

                    j=inb1(pre%len+1,pre%len+1)
                if(ij==1)then
        call block_transfer(tmp_oper,next%sub(j)%elec_up)
        call block_transfer(tmp_oper1,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(tmp_oper,next%sub(j)%elec_down)
                        endif
                if(ij==3)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(tmp_oper,next%sub(j)%num_sn)
                endif

                        if(nleg1.lt.nleg11(pre%len))then
                        do j=nleg1+1,nleg11(pre%len)
                        if(ij==1)then
                call deallocate_block(next%sub(j)%elec_up)
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==2)then
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==3)then
                call deallocate_block(next%sub(j)%spin_sd)
                        endif
                        if(ij==4)then
                call deallocate_block(next%sub(j)%spin_sz)
                        endif
                        if(ij==5)then
                call deallocate_block(next%sub(j)%num_sn)
                        endif
                                enddo
                                endif
                
		call deallocate_block(tmp_oper)

	do ii=1, nleg11(next%len)

                        if(ij==1)then
		call update_trun_ndia(next%sub(ii)%elec_up,pre%sub(ii)%elec_up,trun) 
		call update_trun_ndia(next%sub(ii)%elec_down,pre%sub(ii)%elec_down,trun) 
                        call deallocate_block(next%sub(ii)%elec_up)
                        call deallocate_block(next%sub(ii)%elec_down)
                        endif

                        if(ij==2)then
		call update_trun_ndia(next%sub(ii)%elec_down, pre%sub(ii)%elec_down,trun)
                        call deallocate_block(next%sub(ii)%elec_down)
                                endif

                        if(ij==3)then
		call update_trun_ndia(next%sub(ii)%spin_sd,pre%sub(ii)%spin_sd,trun)
                        call deallocate_block(next%sub(ii)%spin_sd)
                        endif

                                if(ij==4)then
		call update_trun_dia(next%sub(ii)%spin_sz,pre%sub(ii)%spin_sz,trun)
                        call deallocate_block(next%sub(ii)%spin_sz)
                                endif

                        if(ij==5)then
		call update_trun_dia(next%sub(ii)%num_sn,pre%sub(ii)%num_sn,trun)
                        call deallocate_block(next%sub(ii)%num_sn)
                                endif
	end do

                        enddo
	call update_trun_dia(next%ham,tmp_oper,trun)
                        call block_transfer(tmp_oper, pre%ham)
                        call deallocate_block(tmp_oper)
                        call deallocate_block(next%ham)
                pre%len=next%len

end subroutine Get_Hamiltonian_Sys


!==========================================================
!Get system hamiltonian from block and one site
!==========================================================
subroutine Get_Hamiltonian_Env(pre,next,basis, trun)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: basis
	type(Total_Model),intent(inout) :: pre
	type(Total_Model),intent(inout) :: next
	type(Total_Block),intent(in) :: trun

	integer :: i,x,y,old_xidx,old_yidx,new_xidx,new_yidx
        integer :: i1, j, j1,ij, k, nleg1,ii
	real(8) :: coefx,coefy,coefsd,coefsz,coefsn
	type(Total_Block) :: tmp_oper, tmp_oper1
        type(Block_set) tmp_sys 


	!<1>: Update the Hamiltonian
	next%len=pre%len+1
        old_xidx=Lattice(1,nxc*nleg-pre%len+1)
        old_yidx=Lattice(2,nxc*nleg-pre%len+1)
        new_xidx=Lattice(1,nxc*nleg-pre%len+2)
        new_yidx=Lattice(2,nxc*nleg-pre%len+2)

	call update_block_dia(pre%ham,basis,next%ham)

	!<2>: Update all the operators
        do ij=1, 5, 2
                        do k=1, neibt
                j=pre%len+1+mshift
                I=neib(nposi_lat(J),k)
         II=nposi_dm(I)-mshift  
          if(II.le.pre%len.and.II.gt.0)then
                        i1=inb1(ii, pre%len) 
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, nposi_lat(j))
	coefsd= jd(k,nposi_lat(j))
	coefsz= jz(k,nposi_lat(j))
	coefsn=  jn(k)
                if(ij==1)then
                if(coefx.ne.0.0)then
	call Get_hamiltonian_ndia(pre%sub(i1)%elec_up,st_elec_up,basis,coefx,next%ham,'F','E')
                endif
                endif
                if(ij==2)then
                if(coefx.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%elec_down,st_elec_down,basis, coefx,next%ham,'F','E')
                endif
                endif
                if(ij==3)then
                if(coefsd.ne.0.0)then
		call Get_hamiltonian_ndia(pre%sub(i1)%spin_sd,st_sd, basis,coefsd,next%ham,'B','E')
                endif
                endif
	endif
	endif
        enddo

!!! reorder the operators
	do y=1,pre%len   
                I1=INB1(y, pre%len)
                if(I1.le.nleg11(pre%len))then 

                if(ij==1)then
                call update_block_ndia_env(pre%sub(i1)%elec_up,basis,next%sub(i1)%elec_up,'F')
                call update_block_ndia_env(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')
                        call deallocate_block(pre%sub(i1)%elec_up)
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==2)then
                        call update_block_ndia_env(pre%sub(i1)%elec_down,basis,next%sub(i1)%elec_down,'F')
                        call deallocate_block(pre%sub(i1)%elec_down)
                        endif
                        if(ij==3)then
                        call update_block_ndia_env(pre%sub(i1)%spin_sd,basis,next%sub(i1)%spin_sd,'B')
                        call deallocate_block(pre%sub(i1)%spin_sd)
                        endif
                        if(ij==4)then
                        call update_block_dia(pre%sub(i1)%spin_sz,basis,next%sub(i1)%spin_sz)
                        call deallocate_block(pre%sub(i1)%spin_sz)
                        endif
                        if(ij==5)then
                        call update_block_dia(pre%sub(i1)%num_sn,basis,next%sub(i1)%num_sn)
                        call deallocate_block(pre%sub(i1)%num_sn)
                        endif

		        endif
                        enddo

                if(ij==1)then
          call update_site_ndia_env(st_elec_up,basis,tmp_oper,'F')
          call update_site_ndia_env(st_elec_down,basis,tmp_oper1,'F')
                endif
                if(ij==2)then
          call update_site_ndia_env(st_elec_down,basis,tmp_oper,'F')
                endif
        if(ij==3)then
          call update_site_ndia_env(st_sd,basis,tmp_oper,'B')
                endif
                if(ij==4)then
          call update_site_dia(st_sz,basis,tmp_oper)
                endif
                if(ij==5)then
          call update_site_dia(num_elec,basis,tmp_oper)
                endif

                        do k=1, neibt
                I=neib(nposi_lat(J),k)
         II=nposi_dm(I)-mshift  
          if(II.le.pre%len.and.II.gt.0)then
                        i1=inb1(ii, pre%len) 
                if(i1.gt.0.and.i1.le.nleg11(pre%len))then
	coefx = jt(k, nposi_lat(j))
	coefsd= jd(k,nposi_lat(j))
	coefsz= jz(k,nposi_lat(j))
	coefsn=  jn(k)
                if(ij==4)then
		call Get_hamiltonian_dia(next%sub(i1)%spin_sz,'N', tmp_oper,'N',coefsz,next%ham,.false.)
                endif
                if(ij==5)then
		call Get_hamiltonian_dia(next%sub(i1)%num_sn,'N', tmp_oper,'N',coefsn,next%ham,.false.)
                endif
	endif
	endif
        enddo

                if(ij==5.and.esite(nposi_lat(j)).ne.0.0)then
            call block_add_block_two(tmp_oper,esite(nposi_lat(j)),next%ham)
                endif
11      continue


                  nleg1=nleg11(pre%len+1)
         do I=1,pre%len   
                I1=INB1(I,pre%len) 
                J = INB1(I,pre%len+1)
                        if(J.le.nleg1.and.j.ne.i1)then

                if(ij==1)then
        call block_transfer(next%sub(i1)%elec_up,next%sub(j)%elec_up)
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(next%sub(i1)%elec_down,next%sub(j)%elec_down)
                endif
                if(ij==3)then
        call block_transfer(next%sub(i1)%spin_sd,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(next%sub(i1)%spin_sz,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(next%sub(i1)%num_sn,next%sub(j)%num_sn)
                endif
                                endif
                                enddo

                    j=inb1(pre%len+1,pre%len+1)
                if(ij==1)then
        call block_transfer(tmp_oper,next%sub(j)%elec_up)
        call block_transfer(tmp_oper1,next%sub(j)%elec_down)
                endif
                if(ij==2)then
        call block_transfer(tmp_oper,next%sub(j)%elec_down)
                        endif
                if(ij==3)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sd)
                endif
                if(ij==4)then
        call block_transfer(tmp_oper,next%sub(j)%spin_sz)
                endif
                if(ij==5)then
        call block_transfer(tmp_oper,next%sub(j)%num_sn)
                endif

                        if(nleg1.lt.nleg11(pre%len))then
                        do j=nleg1+1,nleg11(pre%len)
                        if(ij==1)then
                call deallocate_block(next%sub(j)%elec_up)
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==2)then
                call deallocate_block(next%sub(j)%elec_down)
                        endif
                        if(ij==3)then
                call deallocate_block(next%sub(j)%spin_sd)
                        endif
                        if(ij==4)then
                call deallocate_block(next%sub(j)%spin_sz)
                        endif
                        if(ij==5)then
                call deallocate_block(next%sub(j)%num_sn)
                        endif
                                enddo
                                endif
		call deallocate_block(tmp_oper)

	do ii=1, nleg11(next%len)

                        if(ij==1)then
		call update_trun_ndia(next%sub(ii)%elec_up,pre%sub(ii)%elec_up,trun) 
                        call deallocate_block(next%sub(ii)%elec_up)
		call update_trun_ndia(next%sub(ii)%elec_down, pre%sub(ii)%elec_down,trun)
                        call deallocate_block(next%sub(ii)%elec_down)
                        endif

                        if(ij==2)then
		call update_trun_ndia(next%sub(ii)%elec_down, pre%sub(ii)%elec_down,trun)
                        call deallocate_block(next%sub(ii)%elec_down)
                                endif

                        if(ij==3)then
		call update_trun_ndia(next%sub(ii)%spin_sd,pre%sub(ii)%spin_sd,trun)
                        call deallocate_block(next%sub(ii)%spin_sd)
                        endif

                                if(ij==4)then
		call update_trun_dia(next%sub(ii)%spin_sz,pre%sub(ii)%spin_sz,trun)
                        call deallocate_block(next%sub(ii)%spin_sz)
                                endif

                        if(ij==5)then
		call update_trun_dia(next%sub(ii)%num_sn,pre%sub(ii)%num_sn,trun)
                        call deallocate_block(next%sub(ii)%num_sn)
                                endif
	end do

                        enddo
	call update_trun_dia(next%ham,tmp_oper,trun)
                        call block_transfer(tmp_oper, pre%ham)
                        call deallocate_block(tmp_oper)
                        call deallocate_block(next%ham)
                pre%len=next%len

end subroutine Get_Hamiltonian_Env


subroutine Get_Hamiltonian_Trun(pre,eff,trun)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: trun
	type(Total_Model),intent(in) :: pre
	type(Total_Model),intent(inout) :: eff

	integer :: ii

	call update_trun_dia(pre%ham,eff%ham,trun)
	eff%len=pre%len

	do ii=1, nleg11(pre%len)
		call update_trun_ndia(pre%sub(ii)%elec_up,eff%sub(ii)%elec_up,trun)
                        call deallocate_block(pre%sub(ii)%elec_up)
		call update_trun_ndia(pre%sub(ii)%spin_sd,eff%sub(ii)%spin_sd,trun)
                        call deallocate_block(pre%sub(ii)%spin_sd)
		call update_trun_dia(pre%sub(ii)%num_sn,eff%sub(ii)%num_sn,trun)
                        call deallocate_block(pre%sub(ii)%num_sn)
	end do

end subroutine Get_Hamiltonian_Trun


!=========================================================================
!Get hamiltonian from diagonal interaction term n_i.n_i+1
!=========================================================================
subroutine Get_hamiltonian_dia(Oper1,Flag1,Oper2,Flag2,Coef,Ham,FlagConjg)
	use pubdata
	implicit none

	real(8),intent(in) :: Coef
	logical,intent(in) :: FlagConjg
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: Oper1,Oper2
	type(Total_Block),intent(inout) :: Ham

	integer :: i,x,y
	logical :: oper1_flag,oper2_flag
	integer :: num_up,num_down,row_dim,sdim,oper1_id,oper2_id
	real(8),allocatable :: mat(:,:)
        real*8 coef1

	do i=1,ham%num
		num_up=ham%sub(i)%num_up
		num_down=ham%sub(i)%num_down
		row_dim=ham%sub(i)%row_dim
		sdim=ham%sub(i)%sdim

		oper1_flag=.false.
		do x=1,oper1%num
			if(oper1%sub(x)%num_up==num_up) then
			if(oper1%sub(x)%num_down==num_down) then
				oper1_id=x
				oper1_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		oper2_flag=.false.
		do x=1,oper2%num
			if(oper2%sub(x)%num_up==num_up) then
			if(oper2%sub(x)%num_down==num_down) then
				oper2_id=x
				oper2_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(oper1_flag.and.oper2_flag) then
			allocate(mat(sdim,sdim))

        coef1=coef/dsqrt(1.0d0+num_down)
			call DGEMM(Flag1,Flag2,sdim,sdim,sdim,coef1,oper1%sub(oper1_id)%mat&
					 &,sdim,oper2%sub(oper2_id)%mat,sdim,0.0d0,mat,sdim)

			!For the case with conjugate part
			if(FlagConjg) then
				ham%sub(i)%mat=ham%sub(i)%mat+(mat+transpose(mat))
			else
				ham%sub(i)%mat=ham%sub(i)%mat+mat
			endif
			deallocate(mat)
		endif
	end do

end subroutine Get_hamiltonian_dia


!================================================================
!Get Ham= Coef* [Oper1*(Oper2^T)+ h.c.) (i.e. C^+_1.C_2 + h.c.)
!================================================================
subroutine Get_hamiltonian_ndia(block,site,basis,coef,Ham, FlagSign, Flags)
	use pubdata
        use fact1
	implicit none

	character(len=1),intent(in) :: FlagSign, Flags
	type(Total_Block),intent(in) :: block, site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: Ham
        real(8), intent (in) :: coef

	integer :: i,x,y,new_id,block_id
        integer j1,j2,j3,j4,j5,j6
	integer :: lhs_id,lhs_num_up,lhs_num_down0, lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down, st_id

	integer :: lhs_dim,rhs_dim,up_dif,down_dif,down_dif1,sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,block_flag, st_flag
        integer st_num_down, st_num_up, new_st_num_up, new_st_num_down, st_down_dif
        real*8 coef1, coef11, Sign
	real(8),allocatable :: mat(:,:)
  real(8),external :: w6js
        

	up_dif=block%up_dif
	down_dif1=block%down_dif

	!<3>: Get new_block%sub(:)%mat
	do i=1,ham%num
		rhs_num_up=ham%sub(i)%num_up
		rhs_num_down=ham%sub(i)%num_down
		lhs_dim=ham%sub(i)%row_dim
		rhs_dim=ham%sub(i)%sdim
		
		rhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==rhs_num_up) then
			if(basis%sub(x)%new_num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue

		down_dif=ham%sub(i)%down_dif
		lhs_num_up=rhs_num_up
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		if(rhs_flag.and.lhs_flag) then
			do x=1,basis%sub(rhs_id)%num
				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
                                st_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
                                st_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                do sub_down_dif=-down_dif1, down_dif1,su
				block_flag=.false.
				do y=1,block%num
					if(block%sub(y)%num_up==rhs_sub_num_up) then
					if(block%sub(y)%num_down==rhs_sub_num_down) then
					if(block%sub(y)%down_dif==sub_down_dif) then
						block_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
				lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						block_flag=.true.
						goto 106
					endif
                                        endif
                                        endif
				end do
				106 continue

                        if(block_flag)then

        do new_st_num_down=abs(st_num_down-site%down_dif), st_num_down+site%down_dif,su

                  st_flag=.false.
                        st_down_dif=st_num_down-new_st_num_down
                        new_st_num_up=st_num_up-site%up_dif
                                             do y=1,site%num
                                       if(site%sub(y)%num_up==new_st_num_up) then
                                              if(site%sub(y)%num_down==new_st_num_down) then
                                                        if(site%sub(y)%down_dif==st_down_dif) then
                                                                st_id=y
                                                                st_flag=.true.
                                                                goto 1031
                                                        endif
                                                        endif
                                                        endif
                                                end do
                                                1031 continue


				lhs_sub_flag=.false.

				do y=1,basis%sub(lhs_id)%num
				if(basis%sub(lhs_id)%sub(y)%bl_num_up==lhs_sub_num_up) then
				if(basis%sub(lhs_id)%sub(y)%bl_num_down==lhs_sub_num_down) then
				if(basis%sub(lhs_id)%sub(y)%st_num_down==new_st_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
                                        endif
                                        endif
				end do
				105 continue

				if(lhs_sub_flag.and.block_flag.and.st_flag) then
        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}

        j1=lhs_sub_num_down 
        j2=rhs_sub_num_down  
        j3=down_dif1     
        j4=st_num_down  
        j5=new_st_num_down 
        j6=rhs_num_down     

        if(iw6j1(j1,j2,j3,j4,j5,j6)==1)then
        coef1=w6j1(j1,j2,j3,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j1(j1,j2,j3,j4,j5,j6)=1
        w6j1(j1,j2,j3,j4,j5,j6)=coef1
        endif

        coef1=coef1*dsqrt((1.0d0+j6)/(1.d0+j3))*(-1)**((j2+j5+j3+j6)/2)*site%sub(st_id)%mat(1,1)

                if(Flags=='S')then 
                                              Sign=1.0d0
                                                        if(FlagSign=='F') then
                                                         if(mod(rhs_sub_num_up, 2)==1)Sign=-1.0d0
                                                         if(mod(rhs_sub_num_down, 2)==1)Sign=-1.0d0
                                                                endif

                        else 
                        Sign=1.0d0
                                                                endif



        if(coef1.ne.0.0)then
                coef1=coef1*coef*Sign
	ham%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&=ham%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&+coef1*block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim)

        if(site%down_dif.eq.1)then
                allocate(mat(rhs_sub_dim, lhs_sub_dim))
                mat= transpose(block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim))
	ham%sub(i)%mat(rhs_pos+1:rhs_pos+rhs_sub_dim,lhs_pos+1:lhs_pos+lhs_sub_dim)&
				&=ham%sub(i)%mat(rhs_pos+1:rhs_pos+rhs_sub_dim,lhs_pos+1:lhs_pos+lhs_sub_dim)&
				&+coef1*mat

			deallocate(mat)
                        endif
                                endif 
				endif
			end do
        endif
                enddo
                enddo
		endif
        enddo
end subroutine Get_hamiltonian_ndia

subroutine Get_sys_oper_dia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,sys_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_dia(st_oper,sys_bsm(idx_one),sys_oper)
			else
				call update_site_dia(st_oper,sys_bsm(idx_one),mid_oper)
				call update_trun_dia(mid_oper,sys_oper,systruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(sys_oper,sys_bsm(idx+1),mid_oper)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_dia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(sys_oper,sys_bsm(idx+1),new_oper)
		call deallocate_block(sys_oper)
		
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_dia(st_oper,sys_bsm(idx_two),new_oper)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

end subroutine Get_sys_oper_dia


subroutine Get_env_oper_dia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,env_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_dia(st_oper,env_bsm(idx_one),env_oper)
			else
				call update_site_dia(st_oper,env_bsm(idx_one),mid_oper)
				call update_trun_dia(mid_oper,env_oper,envtruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(env_oper,env_bsm(idx+1),mid_oper)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_dia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(env_oper,env_bsm(idx+1),new_oper)
		call deallocate_block(env_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_dia(st_oper,env_bsm(idx_two),new_oper)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

end subroutine Get_env_oper_dia

subroutine Change_sys_oper_dia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,sys_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(sys_oper,sys_bsm(idx+1),mid_oper)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_dia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(sys_oper,sys_bsm(idx+1),new_oper)
		call deallocate_block(sys_oper)

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

	endif
		
	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif
end subroutine Change_sys_oper_dia


subroutine Change_env_oper_dia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns)
	use pubdata
	implicit none

	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,env_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_dia(env_oper,env_bsm(idx+1),mid_oper)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_dia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_dia(env_oper,env_bsm(idx+1),new_oper)
		call deallocate_block(env_oper)

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_dia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif
	endif
		
	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif

end subroutine Change_env_oper_dia

subroutine Get_sys_oper_ndia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,sys_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_ndia_sys(st_oper,sys_bsm(idx_one),sys_oper,Types)
			else
				call update_site_ndia_sys(st_oper,sys_bsm(idx_one),mid_oper,Types)
				call update_trun_ndia(mid_oper,sys_oper,systruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),mid_oper,Types)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_ndia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),new_oper,Types)
		call deallocate_block(sys_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_ndia_sys(st_oper,sys_bsm(idx_two),new_oper,Types)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif

end subroutine Get_sys_oper_ndia


!Get new_oper from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Get_env_oper_ndia(st_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		if(idx_one<=1) then
			call block_transfer(st_oper,env_oper)
		else
			if(idx_one<=trun_idx) then
				call update_site_ndia_env(st_oper,env_bsm(idx_one),env_oper,Types)
			else
				call update_site_ndia_env(st_oper,env_bsm(idx_one),mid_oper,Types)
				call update_trun_ndia(mid_oper,env_oper,envtruns(idx_one))
				call deallocate_block(mid_oper)
			endif
		endif

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_env(env_oper,env_bsm(idx+1),mid_oper,Types)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_ndia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_env(env_oper,env_bsm(idx+1),new_oper,Types)
		call deallocate_block(env_oper)
		
	!<2-2>: For (idx_one=idx_two)
	else if(idx_one==idx_two) then
		if(idx_one<=1) then
			call block_transfer(st_oper,new_oper)
		else
			call update_site_ndia_env(st_oper,env_bsm(idx_two),new_oper,Types)
		endif
	endif

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

end subroutine Get_env_oper_ndia


!Change operator from old_idx's to new_idx's (old_idx<=new_idx)
subroutine Change_sys_oper_ndia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,sys_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,sys_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),mid_oper,Types)
			call deallocate_block(sys_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,sys_oper)
			else
				call update_trun_ndia(mid_oper,sys_oper,systruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_sys(sys_oper,sys_bsm(idx+1),new_oper,Types)
		call deallocate_block(sys_oper)
		

	!If flag=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,sys_oper,systruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(sys_oper,new_oper)
		call deallocate_block(sys_oper)
	endif
	endif

	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif
end subroutine Change_sys_oper_ndia

subroutine Change_env_oper_ndia(old_oper,old_idx,new_oper,new_idx,trun_idx,Truns,Types)
	use pubdata
	implicit none

	character(len=1) :: Types
	logical,intent(in) :: Truns
	integer,intent(in) :: old_idx,new_idx,trun_idx
	type(Total_Block),intent(in) :: old_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: idx,idx_one,idx_two
	type(Total_Block) :: mid_oper,env_oper

	idx_one=min(old_idx,new_idx)
	idx_two=max(old_idx,new_idx)

	!<1>: For (idx_one<idx_two)
	if(idx_one<=(idx_two-1)) then
		call block_transfer(old_oper,env_oper)

		idx=idx_one
		do while((idx+1)<=(idx_two-1))
			call update_block_ndia_env(env_oper,env_bsm(idx+1),mid_oper,Types)
			call deallocate_block(env_oper)

			idx=idx+1
			if(idx<=trun_idx) then
				call block_transfer(mid_oper,env_oper)
			else
				call update_trun_ndia(mid_oper,env_oper,envtruns(idx))
			endif
			call deallocate_block(mid_oper)
		end do

		call update_block_ndia_env(env_oper,env_bsm(idx+1),new_oper,Types)
		call deallocate_block(env_oper)

	!If Truns=.true., then truncate the operator
	if(Truns) then
		call update_trun_ndia(new_oper,env_oper,envtruns(idx_two))
		call deallocate_block(new_oper)
		call block_transfer(env_oper,new_oper)
		call deallocate_block(env_oper)
	endif

	endif
		
	!<2>: For (idx_one=idx_two)
	if(idx_one==idx_two) then
		call block_transfer(old_oper,new_oper)
	endif
end subroutine Change_env_oper_ndia


!Notice : outvec=coef*sys_bl*invec
subroutine block_vec_dia(sys_bl,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none
	
	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,x,y,sys_num_up,sys_num_down
	integer :: env_dim,bl_num_up,bl_num_down,spos,sdim
        real(8):: coef1

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

        coef1=coef/dsqrt(1.0d0+bl_num_down)

					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==bl_num_up) then
						if(sys_bl%sub(k)%num_down==bl_num_down) then
							call DGEMM('N','N',sdim,env_dim,sdim,coef1,sys_bl%sub(k)%mat,sdim&
									&,invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,1.0d0&
									&,outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim)
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine block_vec_dia


!Notice : outvec=coef*env_bl*invec
subroutine vec_block_dia(env_bl,env_bs,invec,outvec,coef)
	use pubdata
	implicit none
	
	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,x,y,env_num_up,env_num_down
	integer :: sys_dim,bl_num_up,bl_num_down,spos,sdim
        real*8 coef1

	do i=1,invec%num
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					bl_num_up=env_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=env_bs%sub(x)%sub(y)%bl_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim

        coef1=coef/dsqrt(1.0d0+bl_num_down)
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==bl_num_up) then
						if(env_bl%sub(k)%num_down==bl_num_down) then
							call DGEMM('N','T',sys_dim,sdim,sdim,coef1&
									&,invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
									&,sys_dim,env_bl%sub(k)%mat,sdim,1.0d0&
									&,outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim)
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine vec_block_dia


!Notice : outvec=alpha*sys_st*invec
subroutine site_vec_dia(sys_st,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: st_flag
	integer :: i,k,x,y,st_id,spos,sdim,env_dim
	integer :: sys_num_up,sys_num_down,st_num_up,st_num_down
	real(8) :: coef_st, coef1

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==st_num_up) then
						if(sys_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 102
						endif
						endif
					end do
					102 continue

					if(st_flag) then
						coef_st=coef*sys_st%sub(st_id)%mat(1,1)
                                                coef_st=coef_st/dsqrt(1.d0+st_num_down)
						outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
							&= outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
							&+ coef_st*invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine site_vec_dia


!Notice : outvec=alpha*env_st*invec
subroutine vec_site_dia(env_st,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: st_flag
	integer :: i,k,x,y,st_id,spos,sdim,sys_dim
	integer :: env_num_up,env_num_down,st_num_up,st_num_down
	real(8) :: coef_st

	do i=1,invec%num
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					st_num_up=env_bs%sub(x)%sub(y)%st_num_up
					st_num_down=env_bs%sub(x)%sub(y)%st_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim

					st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==st_num_up) then
						if(env_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 101
						endif
						endif
					end do
					101 continue

					if(st_flag) then
						coef_st=coef*env_st%sub(st_id)%mat(1,1)
                                                coef_st=coef_st/dsqrt(1.d0+st_num_down)
						outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
							&= outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
							&+ coef_st*invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine vec_site_dia


!Notice : outvec=coef*sys_bl*sys_st*invec
subroutine block_site_vec_dia(sys_bl,sys_st,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bl_flag,st_flag
	integer :: i,k,x,y,spos,sdim
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: bl_id,bl_num_up,bl_num_down
	integer :: st_id,st_num_up,st_num_down
	real(8) :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					bl_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==bl_num_up) then
						if(sys_bl%sub(k)%num_down==bl_num_down) then
							bl_id=k
							bl_flag=.true.
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==st_num_up) then
						if(sys_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef_tmp=coef*sys_st%sub(st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))
						call DGEMM('N','N',sdim,env_dim,sdim,coef_tmp,sys_bl%sub(bl_id)%mat&
								&,sdim,invec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim&
								&,1.0d0,outvec%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine block_site_vec_dia


!Notice : outvec=coef*env_bl*env_st*invec
subroutine vec_block_site_dia(env_bl,env_st,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bl_flag,st_flag
	integer :: i,k,x,y,spos,sdim,sys_dim
	integer :: env_num_up,env_num_down
	integer :: bl_id,bl_num_up,bl_num_down
	integer :: st_id,st_num_up,st_num_down
	real(8) :: coef_tmp

	do i=1,invec%num
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					bl_num_up=env_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=env_bs%sub(x)%sub(y)%bl_num_down
					st_num_up=env_bs%sub(x)%sub(y)%st_num_up
					st_num_down=env_bs%sub(x)%sub(y)%st_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim

					bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==bl_num_up) then
						if(env_bl%sub(k)%num_down==bl_num_down) then
							bl_id=k
							bl_flag=.true.
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==st_num_up) then
						if(env_st%sub(k)%num_down==st_num_down) then
							st_id=k
							st_flag=.true.
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef_tmp=coef*env_st%sub(st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))
						call DGEMM('N','T',sys_dim,sdim,sdim,coef_tmp&
								&,invec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim)&
								&,sys_dim,env_bl%sub(bl_id)%mat,sdim,1.0d0&
								&,outvec%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine vec_block_site_dia


!Notice : outvec=coef*sys_st*env_st*invec
subroutine site_vec_site_dia(sys_st,env_st,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_st_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_st_id,env_st_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_st_num_up,sys_st_num_down,env_st_num_up,env_st_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
	real(8) :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_st_num_up=sys_bs%sub(sys_bs_id)%sub(x)%st_num_up
				sys_st_num_down=sys_bs%sub(sys_bs_id)%sub(x)%st_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_st_flag=.false.
				do k=1,sys_st%num
					if(sys_st%sub(k)%num_up==sys_st_num_up) then
					if(sys_st%sub(k)%num_down==sys_st_num_down) then
						sys_st_id=k
						sys_st_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_st_num_up=env_bs%sub(env_bs_id)%sub(y)%st_num_up
					env_st_num_down=env_bs%sub(env_bs_id)%sub(y)%st_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==env_st_num_up) then
						if(env_st%sub(k)%num_down==env_st_num_down) then
							env_st_id=k
							env_st_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_st_flag.and.env_st_flag) then
						coef_tmp=coef*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+sys_st_num_down)*(1.0d0+env_st_num_down))
						outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
							&=outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
							&+coef_tmp*invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)
					endif
				end do
			end do
		endif
	end do

end subroutine site_vec_site_dia


!Notice : outvec=coef*sys_bl*env_st*invec
subroutine block_vec_site_dia(sys_bl,env_st,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_st_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_bl_id,env_st_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,env_st_num_up,env_st_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
	real(8) :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_bl_num_up=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_up
				sys_bl_num_down=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_bl_flag=.false.
				do k=1,sys_bl%num
					if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
					if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
						sys_bl_id=k
						sys_bl_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_st_num_up=env_bs%sub(env_bs_id)%sub(y)%st_num_up
					env_st_num_down=env_bs%sub(env_bs_id)%sub(y)%st_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==env_st_num_up) then
						if(env_st%sub(k)%num_down==env_st_num_down) then
							env_st_id=k
							env_st_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_bl_flag.and.env_st_flag) then
						coef_tmp=coef*env_st%sub(env_st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+env_st_num_down)*(1.0d0+sys_bl_num_down))
						call DGEMM('N','N',sys_dim,env_dim,sys_dim,coef_tmp,sys_bl%sub(sys_bl_id)%mat,sys_dim&
								&,invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim),sys_dim&
								&,1.0d0,outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim),sys_dim)
					endif
				end do
			end do
		endif
	end do

end subroutine block_vec_site_dia


!Notice : outvec=coef*sys_st*env_bl*invec
subroutine site_vec_block_dia(sys_st,env_bl,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_bl_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_st_id,env_bl_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_st_num_up,sys_st_num_down,env_bl_num_up,env_bl_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
	real(8) :: coef_tmp

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_st_num_up=sys_bs%sub(sys_bs_id)%sub(x)%st_num_up
				sys_st_num_down=sys_bs%sub(sys_bs_id)%sub(x)%st_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_st_flag=.false.
				do k=1,sys_st%num
					if(sys_st%sub(k)%num_up==sys_st_num_up) then
					if(sys_st%sub(k)%num_down==sys_st_num_down) then
						sys_st_id=k
						sys_st_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_bl_num_up=env_bs%sub(env_bs_id)%sub(y)%bl_num_up
					env_bl_num_down=env_bs%sub(env_bs_id)%sub(y)%bl_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==env_bl_num_up) then
						if(env_bl%sub(k)%num_down==env_bl_num_down) then
							env_bl_id=k
							env_bl_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_st_flag.and.env_bl_flag) then
						coef_tmp=coef*sys_st%sub(sys_st_id)%mat(1,1)
						coef_tmp=coef_tmp/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_st_num_down))
						call DGEMM('N','T',sys_dim,env_dim,env_dim,coef_tmp&
								&,invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
								&,sys_dim,env_bl%sub(env_bl_id)%mat,env_dim,1.0d0&
								&,outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim),sys_dim)
					endif
				end do
			end do
		endif
	end do

end subroutine site_vec_block_dia


!Notice : outvec=coef*sys_bl*env_bl*invec
subroutine block_vec_block_dia(sys_bl,env_bl,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_bl_flag
	integer :: i,k,x,y,sys_bs_id,env_bs_id,sys_bl_id,env_bl_id
	integer :: sys_num_up,sys_num_down,env_num_up,env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,env_bl_num_up,env_bl_num_down
	integer :: sys_pos,sys_dim,env_pos,env_dim
	real(8),allocatable :: mat(:,:)
        real*8 coef1

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		sys_bs_flag=.false.
		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				sys_bs_id=x
				sys_bs_flag=.true.
				goto 101
			endif
			endif
		end do
		101 continue

		env_bs_flag=.false.
		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				env_bs_id=x
				env_bs_flag=.true.
				goto 102
			endif
			endif
		end do
		102 continue

		if(sys_bs_flag.and.env_bs_flag) then
			do x=1,sys_bs%sub(sys_bs_id)%num
				sys_bl_num_up=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_up
				sys_bl_num_down=sys_bs%sub(sys_bs_id)%sub(x)%bl_num_down
				sys_pos=sys_bs%sub(sys_bs_id)%sub(x)%spos
				sys_dim=sys_bs%sub(sys_bs_id)%sub(x)%sdim

				sys_bl_flag=.false.
				do k=1,sys_bl%num
					if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
					if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
						sys_bl_id=k
						sys_bl_flag=.true.
						goto 103
					endif
					endif
				end do
				103 continue

				do y=1,env_bs%sub(env_bs_id)%num
					env_bl_num_up=env_bs%sub(env_bs_id)%sub(y)%bl_num_up
					env_bl_num_down=env_bs%sub(env_bs_id)%sub(y)%bl_num_down
					env_pos=env_bs%sub(env_bs_id)%sub(y)%spos
					env_dim=env_bs%sub(env_bs_id)%sub(y)%sdim

					env_bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==env_bl_num_up) then
						if(env_bl%sub(k)%num_down==env_bl_num_down) then
							env_bl_id=k
							env_bl_flag=.true.
							goto 104
						endif
						endif
					end do
					104 continue

					if(sys_bl_flag.and.env_bl_flag) then
        if(sys_dim*env_dim==0)stop
						allocate(mat(sys_dim,env_dim))
						coef1=coef/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_bl_num_down))
						call DGEMM('N','N',sys_dim,env_dim,sys_dim,1.0d0,sys_bl%sub(sys_bl_id)%mat,sys_dim&
								&,invec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
								&,sys_dim,0.0d0,mat,sys_dim)
						call DGEMM('N','T',sys_dim,env_dim,env_dim,coef1,mat,sys_dim,env_bl%sub(env_bl_id)%mat&
								&,env_dim,1.0d0,outvec%sub(i)%vec(sys_pos+1:sys_pos+sys_dim&
								&,env_pos+1:env_pos+env_dim),sys_dim)
						deallocate(mat)
					endif
				end do
			end do
		endif
	end do

end subroutine block_vec_block_dia


!Outvec=coef*(C^+_i(sys_bl)*C_j(sys_st)+h.c.)*invec Matrix-vector multiplication
subroutine block_site_vec_ndia(sys_bl,FlagBL,sys_st,FlagST,sys_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: FlagBL,FlagST
	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bs_flag,bl_flag,st_flag
	integer :: i,k,x,y,up_dif,down_dif
	integer :: sys_num_up,sys_num_down,bl_num,st_num,tot_num
	integer :: bl_num_up,bl_num_down,new_bl_num_up,new_bl_num_down
	integer :: st_num_up,st_num_down,new_st_num_up,new_st_num_down
	integer :: bl_id,st_id,old_pos,old_dim,new_pos,new_dim,env_dim
	real(8) :: coef_tmp,Signs,SignBL,SignST
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=sys_st%up_dif).or.(down_dif1/=sys_st%down_dif)) then
		write(*,*) "Quantum number error in block_site_vec!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_dim=invec%sub(i)%env_dim

		do k=1,sys_bs%num
			if(sys_bs%sub(k)%new_num_up==sys_num_up) then
			if(sys_bs%sub(k)%new_num_down==sys_num_down) then
				do x=1,sys_bs%sub(k)%num
					bl_num_up=sys_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=sys_bs%sub(k)%sub(x)%bl_num_down
					st_num_up=sys_bs%sub(k)%sub(x)%st_num_up
					st_num_down=sys_bs%sub(k)%sub(x)%st_num_down
					old_pos=sys_bs%sub(k)%sub(x)%spos
					old_dim=sys_bs%sub(k)%sub(x)%sdim
					
        do down_dif=-down_dif1, down_dif1, su  
        do st_down_dif=-down_dif1, down_dif1, su  

						bl_flag=.false.
						do y=1,sys_bl%num
							if(sys_bl%sub(y)%num_up==bl_num_up) then
							if(sys_bl%sub(y)%num_down==bl_num_down) then
							if(sys_bl%sub(y)%down_dif==down_dif) then
								bl_id=y
								bl_flag=.true.
								goto 102
							endif
							endif
							endif
						end do
						102 continue

					new_st_num_up=st_num_up-up_dif
					new_st_num_down=st_num_down-st_down_dif

						st_flag=.false.

						do y=1,sys_st%num
							if(sys_st%sub(y)%num_up==new_st_num_up) then
							if(sys_st%sub(y)%num_down==new_st_num_down) then
							if(sys_st%sub(y)%down_dif==st_down_dif) then
								st_id=y
								st_flag=.true.
								goto 103
							endif
							endif
							endif
						end do
						103 continue

					new_bl_num_up=bl_num_up+up_dif
					new_bl_num_down=bl_num_down+down_dif
                
                if(bl_flag.and.st_flag)then

					bs_flag=.false.
					do y=1,sys_bs%sub(k)%num
						if(sys_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(sys_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
						if(sys_bs%sub(k)%sub(y)%st_num_up==new_st_num_up) then
						if(sys_bs%sub(k)%sub(y)%st_num_down==new_st_num_down) then
							bs_flag=.true.
							new_pos=sys_bs%sub(k)%sub(y)%spos
							new_dim=sys_bs%sub(k)%sub(y)%sdim
							goto 101
						endif
						endif
						endif
						endif
					end do
					101 continue

                                                if(bs_flag)then
							SignBL=1.0d0
							
							SignST=0.0d0
							if(FlagST=='B') then !For Boson
								SignST=1.0d0
							else if(FlagST=='F') then !For Fermion
								bl_num=bl_num_down !!!+bl_num_down
								tot_num=bl_num
								if(mod(tot_num,2)==0) then
									SignST=1.0d0
								else
									SignST=-1.0d0
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down 
        j6=sys_num_down
              if(iw6j1(j1,j2,j3,j4,j5,j6)==1)then
        coef1=w6j1(j1,j2,j3,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j1(j1,j2,j3,j4,j5,j6)=1
        w6j1(j1,j2,j3,j4,j5,j6)=coef1
        endif

        if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)/dsqrt(1.0d0+down_dif1)
				Signs=SignBL*SignST*coef1
			coef_tmp=Signs* coef* sys_st%sub(st_id)%mat(1,1)
	call DGEMM('N','N',new_dim,env_dim,old_dim,coef_tmp,sys_bl%sub(bl_id)%mat,new_dim&
		&,invec%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim,1.0d0&
				&,outvec%sub(i)%vec(new_pos+1:new_pos+new_dim,1:env_dim),new_dim)

if(down_dif1.eq.1)then
                call DGEMM('T','N',old_dim,env_dim,new_dim,coef_tmp,sys_bl%sub(bl_id)%mat,new_dim&
                        &,invec%sub(i)%vec(new_pos+1:new_pos+new_dim,1:env_dim),new_dim,1.0d0&
                                &,outvec%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim)

        endif



						endif
					endif

                        endif
        go to 110

					!================================================
					!For C^+_j(site).C_i(block) case
					new_bl_num_up=bl_num_up-up_dif
					new_bl_num_down=bl_num_down-down_dif
					new_st_num_up=st_num_up+up_dif
					new_st_num_down=st_num_down+down_dif

					bs_flag=.false.
					do y=1,sys_bs%sub(k)%num
						if(sys_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(sys_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
							bs_flag=.true.
							new_pos=sys_bs%sub(k)%sub(y)%spos
							new_dim=sys_bs%sub(k)%sub(y)%sdim
							goto 104
						endif
						endif
					end do
					104 continue

					if(bs_flag) then
						bl_flag=.false.
						do y=1,sys_bl%num
							if(sys_bl%sub(y)%num_up==new_bl_num_up) then
							if(sys_bl%sub(y)%num_down==new_bl_num_down) then
								bl_id=y
								bl_flag=.true.
								goto 105
							endif
							endif
						end do
						105 continue

						st_flag=.false.
						do y=1,sys_st%num
							if(sys_st%sub(y)%num_up==st_num_up) then
							if(sys_st%sub(y)%num_down==st_num_down) then
								st_id=y
								st_flag=.true.
								goto 106
							endif
							endif
						end do
						106 continue
						
						if(bl_flag.and.st_flag) then
							!<a>: For sys_block operator
							SignBL=1.0d0

							!<b>: For sys_site operator
							SignST=0.0d0
							if(FlagST=='B') then !For Boson
								SignST=1.0d0
							else if(FlagST=='F') then !For Fermion
								bl_num=new_bl_num_down !!!+new_bl_num_down
								tot_num=bl_num
								if(mod(tot_num,2)==0) then
									SignST=1.0d0
								else
									SignST=-1.0d0
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down 
        j6=sys_num_down
        coef1=1./dsqrt(1.0d0+down_dif1)*w6js(j1, j2, j3, j4,j5, j6)

        if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)
							Signs=SignBL*SignST*coef1
							coef_tmp=Signs* coef* sys_st%sub(st_id)%mat(1,1)
							call DGEMM('T','N',new_dim,env_dim,old_dim,coef_tmp,sys_bl%sub(bl_id)%mat,old_dim&
									&,invec%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim,1.0d0&
									&,outvec%sub(i)%vec(new_pos+1:new_pos+new_dim,1:env_dim),new_dim)
						endif
						endif
					endif
110     continue
                        
				end do
                        enddo
                        enddo
			endif
			endif
		end do
	end do


end subroutine block_site_vec_ndia


!Outvec=coef*(C^+_i(env_bl)*C_j(env_st)+h.c.)*invec Matrix-wavefunction multiplication
subroutine vec_block_site_ndia(env_bl,FlagBL,env_st,FlagST,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: FlagBL,FlagST
	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	logical :: bs_flag,bl_flag,st_flag
	integer :: i,k,x,y,up_dif,down_dif,env_num_up,env_num_down
	integer :: sys_num_up,sys_num_down,sys_num,bl_num,st_num,tot_num
	integer :: bl_num_up,bl_num_down,new_bl_num_up,new_bl_num_down
	integer :: st_num_up,st_num_down,new_st_num_up,new_st_num_down
	integer :: bl_id,st_id,old_pos,old_dim,new_pos,new_dim,sys_dim
	real(8) :: coef_tmp,Signs,SignBL,SignST
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=env_bl%up_dif
	down_dif1=env_bl%down_dif
	if((up_dif/=env_st%up_dif).or.(down_dif1/=env_st%down_dif)) then
		write(*,*) "Quantum number error in vec_block_site!"
		stop
	endif

	!Start multiplication mat_vec
	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
		sys_dim=invec%sub(i)%sys_dim

		do k=1,env_bs%num
			if(env_bs%sub(k)%new_num_up==env_num_up) then
			if(env_bs%sub(k)%new_num_down==env_num_down) then
				do x=1,env_bs%sub(k)%num
					bl_num_up=env_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=env_bs%sub(k)%sub(x)%bl_num_down
					st_num_up=env_bs%sub(k)%sub(x)%st_num_up
					st_num_down=env_bs%sub(k)%sub(x)%st_num_down
					old_pos=env_bs%sub(k)%sub(x)%spos
					old_dim=env_bs%sub(k)%sub(x)%sdim

					!================================================
					!For C^+_i(block).C_j(site) case

                        do bl_down_dif=-down_dif1, down_dif1, su

						bl_flag=.false.
						do y=1,env_bl%num
							if(env_bl%sub(y)%num_up==bl_num_up) then
							if(env_bl%sub(y)%num_down==bl_num_down) then
							if(env_bl%sub(y)%down_dif==bl_down_dif) then
								bl_id=y
								bl_flag=.true.
								goto 102
							endif
							endif
							endif
						end do
						102 continue

                        do st_down_dif=-down_dif1, down_dif1, su
						st_flag=.false.
					new_st_num_up=st_num_up-up_dif
					new_st_num_down=st_num_down-st_down_dif

						do y=1,env_st%num
							if(env_st%sub(y)%num_up==new_st_num_up) then
							if(env_st%sub(y)%num_down==new_st_num_down) then
							if(env_st%sub(y)%down_dif==st_down_dif) then
								st_id=y
								st_flag=.true.
								goto 103
							endif
							endif
							endif
						end do
						103 continue
					new_bl_num_up=bl_num_up+up_dif
					new_bl_num_down=bl_num_down+bl_down_dif

                if(st_flag.and.bl_flag)then
					bs_flag=.false.
					do y=1,env_bs%sub(k)%num
						if(env_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(env_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
						if(env_bs%sub(k)%sub(y)%st_num_up==new_st_num_up) then
						if(env_bs%sub(k)%sub(y)%st_num_down==new_st_num_down) then
							bs_flag=.true.
							new_pos=env_bs%sub(k)%sub(y)%spos
							new_dim=env_bs%sub(k)%sub(y)%sdim
							goto 101
						endif
						endif
						endif
						endif
					end do
					101 continue

					if(bs_flag) then
							!<a>: For env_block operator
							SignBL=0.0d0
							if(FlagBL=='B') then !For Boson
								SignBL=1.0d0
							else if(FlagBL=='F') then !For Fermion
								sys_num=sys_num_down 
								st_num=new_st_num_down 
								tot_num=sys_num+st_num
								if(mod(tot_num,2)==0) then
									SignBL=1.0d0
								else
									SignBL=-1.0d0
								endif
							endif

							!<b>: For env_site operator
							SignST=0.0d0
							if(FlagST=='B') then !For Boson
								SignST=1.0d0
							else if(FlagST=='F') then 
								sys_num=sys_num_down 
								tot_num=sys_num
								if(mod(tot_num,2)==0) then
									SignST=1.0d0
								else
									SignST=-1.0d0
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down 
        j6=env_num_down

      if(iw6j1(j1,j2,j3,j4,j5,j6)==1)then
        coef1=w6j1(j1,j2,j3,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6) 
        iw6j1(j1,j2,j3,j4,j5,j6)=1
        w6j1(j1,j2,j3,j4,j5,j6)=coef1
        endif


                if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)/dsqrt(1.0d0+down_dif1)


							Signs=SignBL*SignST*coef1
							coef_tmp=Signs* coef *env_st%sub(st_id)%mat(1,1)
				call DGEMM('N','T',sys_dim,new_dim,old_dim,coef_tmp&
								&,invec%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim),sys_dim&
									&,env_bl%sub(bl_id)%mat,new_dim,1.0d0&
									&,outvec%sub(i)%vec(1:sys_dim,new_pos+1:new_pos+new_dim),sys_dim)

if(down_dif1.eq.1)then
                        call DGEMM('N','N',sys_dim,old_dim,new_dim,coef_tmp&
                                     &,invec%sub(i)%vec(1:sys_dim,new_pos+1:new_pos+new_dim),sys_dim&
                                     &,env_bl%sub(bl_id)%mat,new_dim,1.0d0&
                                           &,outvec%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim),sys_dim)
                endif
        




						endif
					endif
                                        endif
                go to 110


					!For C_j^+(site).C_i(block) case
					new_bl_num_up=bl_num_up-up_dif
					new_bl_num_down=bl_num_down-down_dif
					new_st_num_up=st_num_up+up_dif
					new_st_num_down=st_num_down+down_dif

					bs_flag=.false.
					do y=1,env_bs%sub(k)%num
						if(env_bs%sub(k)%sub(y)%bl_num_up==new_bl_num_up) then
						if(env_bs%sub(k)%sub(y)%bl_num_down==new_bl_num_down) then
							bs_flag=.true.
							new_pos=env_bs%sub(k)%sub(y)%spos
							new_dim=env_bs%sub(k)%sub(y)%sdim
							goto 104
						endif
						endif
					end do
					104 continue

					if(bs_flag) then
						bl_flag=.false.
						do y=1,env_bl%num
							if(env_bl%sub(y)%num_up==new_bl_num_up) then
							if(env_bl%sub(y)%num_down==new_bl_num_down) then
								bl_id=y
								bl_flag=.true.
								goto 105
							endif
							endif
						end do
						105 continue

						st_flag=.false.
						do y=1,env_st%num
							if(env_st%sub(y)%num_up==st_num_up) then
							if(env_st%sub(y)%num_down==st_num_down) then
								st_id=y
								st_flag=.true.
								goto 106
							endif
							endif
						end do
						106 continue
						
						if(bl_flag.and.st_flag) then
							!<a>: For env_block operator
							SignBL=0.0d0
							if(FlagBL=='B') then !For Boson
								SignBL=1.0d0
							else if(FlagBL=='F') then !For Fermion
								sys_num=sys_num_down
								st_num=st_num_down
								tot_num=sys_num+st_num
								if(mod(tot_num,2)==0) then
									SignBL=1.0d0
								else
									SignBL=-1.0d0
								endif
							endif

							!<b>: For env_site operator
							SignST=0.0d0
							if(FlagST=='B') then !For Boson
								SignST=1.0d0
							else if(FlagST=='F') then !For Fermion
								sys_num=sys_num_down
								tot_num=sys_num
								if(mod(tot_num,2)==0) then
									SignST=1.0d0
								else
									SignST=-1.0d0
								endif
							endif

        j1=new_bl_num_down 
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down 
        j5=new_st_num_down 
        j6=env_num_down
        coef1=1./dsqrt(1.0d0+down_dif1)*w6js(j1, j2, j3, j4,j5, j6)

        if(coef1.ne.0.0)then
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)
							Signs=SignBL*SignST*coef1
							coef_tmp=Signs* coef* env_st%sub(st_id)%mat(1,1)
							call DGEMM('N','N',sys_dim,new_dim,old_dim,coef_tmp&
									&,invec%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim),sys_dim&
									&,env_bl%sub(bl_id)%mat,old_dim,1.0d0&
									&,outvec%sub(i)%vec(1:sys_dim,new_pos+1:new_pos+new_dim),sys_dim)
						endif
					endif
					endif
110             continue
				end do
				end do
				end do
			endif
			endif
		end do
	end do

end subroutine vec_block_site_ndia

!Outvec=coef*(C^+_i(sys_bl)*C_j(env_bl)+h.c.)*invec Matrix-wavefunction multiplication
subroutine block_vec_block_ndia(sys_bl,FlagSL,env_bl,FlagEL,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: FlagSL,FlagEL
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_bl_flag
	integer :: sys_bs_id,env_bs_id,sys_bl_id,env_bl_id,sys_num
	integer :: sys_st_num_down, env_st_num_up,env_st_num_down,env_st_num,sys_bl_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
        integer :: sys_st_num_up
	integer :: old_sys_pos,old_env_pos,old_sys_dim,old_env_dim
	integer :: new_sys_pos,new_env_pos,new_sys_dim,new_env_dim
	real(8) :: coef_tmp,SignSL,SignEL,Signs
	real(8),allocatable :: mat(:,:)
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j24
        integer j5, j6,j55, j66,  env_down_dif, bl_down_dif, bl_down_dif1
        integer ji1,ji2,ji3,ji4, ji5, ji6,ji7,ji8,ji9, jj, j7,j8,j9,j77,j88,j99
        integer x1, y1
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=env_bl%up_dif).or.(down_dif1/=env_bl%down_dif)) then
		write(*,*) "Quantum number error in block_vec_block!"
		stop
	endif


	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							old_sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							env_bl_num_up=env_bs%sub(y)%sub(n)%bl_num_up
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							old_env_dim=env_bs%sub(y)%sub(n)%sdim

                do down_dif=-down_dif1, down_dif1, su
							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
								if(sys_bl%sub(k)%down_dif==down_dif) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                        
				new_sys_bl_num_up=sys_bl_num_up+up_dif
				new_sys_bl_num_down=sys_bl_num_down+down_dif


                do env_down_dif=-down_dif1, down_dif1, su
				new_env_bl_num_up=env_bl_num_up-up_dif
				new_env_bl_num_down=env_bl_num_down-env_down_dif

							env_bl_flag=.false.
							do k=1,env_bl%num
					if(env_bl%sub(k)%num_up==new_env_bl_num_up) then
					if(env_bl%sub(k)%num_down==new_env_bl_num_down) then
					if(env_bl%sub(k)%down_dif==env_down_dif) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue


					if(sys_bl_flag.and.env_bl_flag) then
        do new_sys_num_down= abs(sys_num_down-down_dif1), sys_num_down+down_dif1, su
						new_sys_num_up=sys_num_up+up_dif
						new_env_num_up=env_num_up-up_dif
					new_env_num_down=new_sys_num_down
                !! total spin=0
							!<1>: For C_i^+(sys_block).C_j(env_block) case

							sys_bs_flag=.false.
					do k=1,sys_bs%num
				if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
				if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
						do l=1,sys_bs%sub(k)%num
					if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
					if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
					if(sys_bs%sub(k)%sub(l)%st_num_down==sys_st_num_down) then
					if(sys_bs%sub(k)%sub(l)%st_num_up==sys_st_num_up) then
								sys_bs_flag=.true.
					new_sys_pos=sys_bs%sub(k)%sub(l)%spos
						new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
                                        x1=k
											goto 103
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
							if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
							if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_down==env_st_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_up==env_st_num_up) then 
											env_bs_flag=.true.
									new_env_pos=env_bs%sub(k)%sub(l)%spos
									new_env_dim=env_bs%sub(k)%sub(l)%sdim
                                                        y1=k
											goto 104
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_bl_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=1.0d0

										!<b>: For env_block operator
										SignEL=0.0d0
										if(FlagEL=='B') then !For Boson
											SignEL=1.0d0
										else if(FlagEL=='F') then !For Fermion
								sys_num=sys_num_down
								env_st_num=env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=1.0d0
											else
												SignEL=-1.0d0
											endif
										endif
        coef1=0.0d0
        
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down

        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down      
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down

      if(iw6j3(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j3(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j3(j1,j2,j3-j1,j4,j5,j6)=1
        w6j3(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 105
      if(iw6j3(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j3(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  !! eq51
        iw6j3(j11,j22,j33-j11,j44,j55,j66)=1
        w6j3(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif

        if(coef11.eq.0.0)go to 105

        coef1=coef1*coef11

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))
        
        if(mod((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2, 2)==1)coef1=-coef1


										Signs=SignSL*SignEL
										coef_tmp=Signs* coef*coef1

        go to 1051

										allocate(mat(new_sys_dim,old_env_dim))
										call DGEMM('N','N',new_sys_dim,old_env_dim,old_sys_dim,1.0d0&
												&,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												&,0.0d0,mat,new_sys_dim)
										call DGEMM('N','T',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												&,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,new_env_dim,1.0d0&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)

1051    continue

                    allocate(mat(new_sys_dim,old_env_dim))
                                                       call DGEMM('N','N',new_sys_dim,old_env_dim,old_sys_dim,1.0d0&
                                                                   &,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
                                                                          &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
                                                                               &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
                                                                                    &,0.0d0,mat,new_sys_dim)
                                            call DGEMM('N','N',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
                                                    &,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,old_env_dim,1.0d0&
                                                          &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
                                                                 &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
                                                                                deallocate(mat)
if(down_dif1.eq.1)then
                allocate(mat(old_sys_dim,new_env_dim))
                        call DGEMM('T','N',old_sys_dim,new_env_dim,new_sys_dim,1.0d0&
                                        &,sys_bl%sub(sys_bl_id)%mat,new_sys_dim&
                                                &,invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
                                &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim&
                                                                &,0.0d0,mat,old_sys_dim)
                                      call DGEMM('N','T',old_sys_dim,old_env_dim,new_env_dim,coef_tmp&
                                        &,mat,old_sys_dim,env_bl%sub(env_bl_id)%mat,old_env_dim,1.0d0&
                                                &,outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
                                                &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim)
                                                                                deallocate(mat)
        endif




										goto 105
									endif
									endif
                                                                endif
								end do
							endif
							105 continue
        go to 110

							
							!<2>: For C^+_j(env_block)).C_i(sys_block) case
							new_sys_bl_num_up=sys_bl_num_up-up_dif
							new_sys_bl_num_down=sys_bl_num_down-down_dif
							new_env_bl_num_up=env_bl_num_up+up_dif
							new_env_bl_num_down=env_bl_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==new_sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==new_sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==env_bl_num_up) then
								if(env_bl%sub(k)%num_down==env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_bl_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down
        j6=new_sys_bl_num_down
        j11=new_env_num_down 
        j22=down_dif1
        j33=env_num_down    
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down
        coef1=dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2)
										SignSL=1.0d0

										!<b>: For env_block operator
										SignEL=0.0d0
										if(FlagEL=='B') then !For Boson
											SignEL=1.0d0
										else if(FlagEL=='F') then !For Fermion
								sys_num=new_sys_num_down
								env_st_num=env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=1.0d0
											else
												SignEL=-1.0d0
											endif
										endif

										Signs=SignSL*SignEL
										coef_tmp=Signs* coef*coef1
										allocate(mat(new_sys_dim,old_env_dim))
										call DGEMM('T','N',new_sys_dim,old_env_dim,old_sys_dim,1.0d0&
												 &,sys_bl%sub(sys_bl_id)%mat,old_sys_dim&
												 &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												 &,old_env_pos+1:old_env_pos+old_env_dim),old_sys_dim&
												 &,0.0d0,mat,new_sys_dim)
										call DGEMM('N','T',new_sys_dim,new_env_dim,old_env_dim,coef_tmp&
												 &,mat,new_sys_dim,env_bl%sub(env_bl_id)%mat,new_env_dim,1.0d0&
												 &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												 &,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim)
										deallocate(mat)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
                                enddo
12      continue
                endif
						end do
						end do
                                enddo
                                enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine block_vec_block_ndia

!Outvec=coef*(C^+_i(sys_bl)*C_j(env_st)+h.c.)*invec Matrix-wavefunction multiplication
subroutine block_vec_site_ndia(sys_bl,FlagSL,env_st,FlagET,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: FlagSL,FlagET
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_bl_flag,env_st_flag
	integer :: sys_bs_id,env_bs_id,sys_bl_id,env_st_id,sys_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_st_num_up,env_st_num_down,new_env_st_num_up,new_env_st_num_down
	integer :: old_sys_pos,old_env_pos,old_sys_dim,env_dim
	integer :: new_sys_pos,new_env_pos,new_sys_dim
        integer env_bl_num_down    

        integer:: sys_st_num_down
	real(8) :: coef_tmp,SignSL,SignET,Signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif
	if((up_dif/=env_st%up_dif).or.(down_dif1/=env_st%down_dif)) then
		write(*,*) "Quantum number error in block_vec_site!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
						sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							old_sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							env_dim=env_bs%sub(y)%sub(n)%sdim
							
							!<1>: For C_i^+(sys_block).C_j(env_site) case

                        do bl_down_dif=-sys_bl%down_dif, sys_bl%down_dif, su
							sys_bl_flag=.false.
				new_sys_bl_num_up=sys_bl_num_up+up_dif
				new_sys_bl_num_down=sys_bl_num_down+bl_down_dif
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
								if(sys_bl%sub(k)%down_dif==bl_down_dif) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                do st_down_dif=-env_st%down_dif, env_st%down_dif, su
							env_st_flag=.false.
				new_env_st_num_up=env_st_num_up-up_dif
				new_env_st_num_down=env_st_num_down-st_down_dif

							do k=1,env_st%num
								if(env_st%sub(k)%num_up==new_env_st_num_up) then
								if(env_st%sub(k)%num_down==new_env_st_num_down) then
								if(env_st%sub(k)%down_dif==st_down_dif) then
									env_st_id=k
									env_st_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue

							new_sys_num_up=sys_num_up+up_dif
							new_env_num_up=env_num_up-up_dif
  if(sys_bl_flag.and.env_st_flag) then
			do new_sys_num_down=abs(sys_num_down-sys_bl%down_dif), sys_num_down+sys_bl%down_dif, su
                           new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do k=1,sys_bs%num
				if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
				if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
						if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
						if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
						if(sys_bs%sub(k)%sub(l)%st_num_down==sys_st_num_down) then
									sys_bs_flag=.true.
								new_sys_pos=sys_bs%sub(k)%sub(l)%spos
								new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 103
										endif
                                                                                endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 104
										endif
                                                                                endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_bl_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=1.0d0

										!<b>: For env_site operator
										SignET=0.0d0
										if(FlagET=='B') then !For Boson
											SignET=1.0d0
										else if(FlagET=='F') then !For Fermion
								sys_num=sys_num_down
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=1.0d0
											else
												SignET=-1.0d0
											endif
										endif
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down
        j11=new_env_num_down 
        j22=down_dif1       
        j33=env_num_down    
        j44= env_st_num_down  
        j55=env_bl_num_down   
        j66=new_env_st_num_down 

      if(iw6j3(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j3(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j3(j1,j2,j3-j1,j4,j5,j6)=1
        w6j3(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j2(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j2(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j2(j11,j22,j33-j11,j44,j55,j66)=1
        w6j2(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11       

                if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.d0+j11)*(1.d0+j33)/(1.0d0+down_dif1))
        if(mod((j6+j5+j3+j2)/2+(j55++j44+j11+j22)/2,2).ne.0) coef1=-coef1
							Signs=SignSL*SignET*coef1
					coef_tmp=Signs* coef* env_st%sub(env_st_id)%mat(1,1)


										call DGEMM('N','N',new_sys_dim,env_dim,old_sys_dim,coef_tmp,sys_bl%sub(sys_bl_id)%mat&
												&,new_sys_dim,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+env_dim),old_sys_dim,1.0d0&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+env_dim),new_sys_dim)

if(down_dif1.eq.1)then
  call DGEMM('T','N',old_sys_dim,env_dim,new_sys_dim,coef_tmp,sys_bl%sub(sys_bl_id)%mat&
                                             &,new_sys_dim,invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
                                                              &,new_env_pos+1:new_env_pos+env_dim),new_sys_dim,1.0d0&
                                                                    &,outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
                                                                            &,old_env_pos+1:old_env_pos+env_dim),old_sys_dim)

                                endif
                                endif
										goto 105
									endif
									endif
								end do
							endif
							105 continue
							

                go to 110
							!<2>: For C^+_j(env_site).C_i(sys_block) case
							new_sys_bl_num_up=sys_bl_num_up-up_dif
							new_sys_bl_num_down=sys_bl_num_down-down_dif
							new_env_st_num_up=env_st_num_up+up_dif
							new_env_st_num_down=env_st_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							!new_sys_num_down=sys_num_down-down_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_bl_flag=.false.
							do k=1,sys_bl%num
								if(sys_bl%sub(k)%num_up==new_sys_bl_num_up) then
								if(sys_bl%sub(k)%num_down==new_sys_bl_num_down) then
									sys_bl_id=k
									sys_bl_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_st_flag=.false.
							do k=1,env_st%num
								if(env_st%sub(k)%num_up==env_st_num_up) then
								if(env_st%sub(k)%num_down==env_st_num_down) then
									env_st_id=k
									env_st_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(k)%sub(l)%bl_num_down==new_sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											new_sys_dim=sys_bs%sub(k)%sub(l)%sdim
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_bl_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_block operator
										SignSL=1.0d0

										!<b>: For env_site operator
										SignET=0.0d0
										if(FlagET=='B') then !For Boson
											SignET=1.0d0
										else if(FlagET=='F') then !For Fermion
							sys_num=new_sys_num_down
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=1.0d0
											else
												SignET=-1.0d0
											endif
										endif
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down 
        j6=new_sys_bl_num_down
        j11=new_env_num_down  
        j22=down_dif1        
        j33=env_num_down     
        j44= env_st_num_down 
        j55=env_bl_num_down   
        j66=new_env_st_num_down  
        coef1=dsqrt((1.d0+j11)*(1.d0+j33)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j6+j5+j3+j2)/2+(j55++j44+j11+j22)/2)
										Signs=SignSL*SignET
										coef_tmp=Signs* coef* env_st%sub(env_st_id)%mat(1,1)
										call DGEMM('T','N',new_sys_dim,env_dim,old_sys_dim,coef_tmp&
												&,sys_bl%sub(sys_bl_id)%mat,old_sys_dim&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+old_sys_dim&
												&,old_env_pos+1:old_env_pos+env_dim),old_sys_dim,1.0d0&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim&
												&,new_env_pos+1:new_env_pos+env_dim),new_sys_dim)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
						end do
                                endif
						end do
                                        enddo
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

        

end subroutine block_vec_site_ndia


!Outvec=coef*(C^+_i(sys_st)*C_j(env_bl)+h.c.)*invec Matrix-wavefunction multiplication
subroutine site_vec_block_ndia(sys_st,FlagST,env_bl,FlagEL,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: FlagST,FlagEL
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n,up_dif,down_dif
	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_bl_flag
	integer :: sys_bs_id,env_bs_id,sys_st_id,env_bl_id
	integer :: env_st_num_up,env_st_num_down,env_st_num,tot_num
	integer :: sys_bl_num_up,sys_bl_num_down,sys_bl_num,sys_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_st_num_up,sys_st_num_down,new_sys_st_num_up,new_sys_st_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
	integer :: old_sys_pos,old_env_pos,old_sys_dim,sys_dim
	integer :: new_sys_pos,new_env_pos,old_env_dim,new_env_dim
	real(8) :: coef_tmp,SignST,SignEL,Signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
	real(8) ::  coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_st%up_dif
	down_dif1=sys_st%down_dif
	if((up_dif/=env_bl%up_dif).or.(down_dif1/=env_bl%down_dif)) then
		write(*,*) "Quantum number error in site_vec_block!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							sys_dim=sys_bs%sub(x)%sub(m)%sdim

							!For envionment part
							env_bl_num_up=env_bs%sub(y)%sub(n)%bl_num_up
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							old_env_dim=env_bs%sub(y)%sub(n)%sdim
							

                do st_down_dif=-sys_st%down_dif, sys_st%down_dif, su
							sys_st_flag=.false.
							new_sys_st_num_up=sys_st_num_up+up_dif
							new_sys_st_num_down=sys_st_num_down+st_down_dif
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==sys_st_num_up) then
								if(sys_st%sub(k)%num_down==sys_st_num_down) then
								if(sys_st%sub(k)%down_dif==st_down_dif) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                do bl_down_dif=-env_bl%down_dif, env_bl%down_dif, su
							env_bl_flag=.false.
							new_env_bl_num_up=env_bl_num_up-up_dif
							new_env_bl_num_down=env_bl_num_down-bl_down_dif
							do k=1,env_bl%num
				if(env_bl%sub(k)%num_up==new_env_bl_num_up) then
					if(env_bl%sub(k)%num_down==new_env_bl_num_down) then
						if(env_bl%sub(k)%down_dif==bl_down_dif) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 102
								endif
								endif
                                                                endif
							end do
							102 continue

							!<1>: For C_i^+(sys_site).C_j(env_block) case
				new_sys_num_up=sys_num_up+up_dif
				new_env_num_up=env_num_up-up_dif
                 if(sys_st_flag.and.env_bl_flag) then
		do new_sys_num_down=abs(sys_num_down-sys_st%down_dif),sys_num_down+sys_st%down_dif, su
							new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
							if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
							if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
							if(sys_bs%sub(k)%sub(l)%bl_num_down==sys_bl_num_down) then
							if(sys_bs%sub(k)%sub(l)%bl_num_up==sys_bl_num_up) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											goto 103
										endif
										endif
										endif
                                                        endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
							if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
							if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
							if(env_bs%sub(k)%sub(l)%st_num_up==env_st_num_up) then
							if(env_bs%sub(k)%sub(l)%st_num_down==env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 104
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							if(sys_st_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
		if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
		if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down 
        j6=new_sys_st_num_down
        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down     
        j44= env_bl_num_down
        j55=env_st_num_down 
        j66=new_env_bl_num_down

      if(iw6j2(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j2(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j2(j1,j2,j3-j1,j4,j5,j6)=1
        w6j2(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j3(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j3(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66)  
        iw6j3(j11,j22,j33-j11,j44,j55,j66)=1
        w6j3(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11      

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(j33+1.d0)/(1.0d0+down_dif1))
        if(mod((j66+j55+j33+j22)/2+(j5+j4+j1+j2)/2,2).ne.0)coef1=-coef1

										SignST=0.0d0
										if(FlagST=='B') then !For Boson
											SignST=1.0d0
										else if(FlagST=='F') then !For Fermion
							sys_bl_num=sys_bl_num_down
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=1.0d0
											else
												SignST=-1.0d0
											endif
										endif

										!<b>: For env_block operator
										SignEL=0.0d0
										if(FlagEL=='B') then !For Boson
											SignEL=1.0d0
										else if(FlagEL=='F') then !For Fermion
						sys_num=sys_num_down 
						env_st_num=env_st_num_down
											tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=1.0d0
											else
												SignEL=-1.0d0
											endif
										endif

							Signs=SignST*SignEL*coef1
						coef_tmp=Signs* coef* sys_st%sub(sys_st_id)%mat(1,1)

         call DGEMM('N','N',sys_dim,new_env_dim,old_env_dim,coef_tmp&
                            &,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
                                       &,old_env_pos+1:old_env_pos+old_env_dim),sys_dim&
                                               &,env_bl%sub(env_bl_id)%mat,old_env_dim,1.0d0&
                                                    &,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
                                                          &,new_env_pos+1:new_env_pos+new_env_dim),sys_dim)
if(down_dif1.eq.1)then
                        call DGEMM('N','T',sys_dim,old_env_dim,new_env_dim,coef_tmp&
                              &,invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
                                               &,new_env_pos+1:new_env_pos+new_env_dim),sys_dim&
                                                      &,env_bl%sub(env_bl_id)%mat,old_env_dim,1.0d0&
                                                               &,outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
                                                                   &,old_env_pos+1:old_env_pos+old_env_dim),sys_dim)
                        endif
                go to 105
	call DGEMM('N','T',sys_dim,new_env_dim,old_env_dim,coef_tmp&
			&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
				&,old_env_pos+1:old_env_pos+old_env_dim),sys_dim&
				&,env_bl%sub(env_bl_id)%mat,new_env_dim,1.0d0&
					&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
						&,new_env_pos+1:new_env_pos+new_env_dim),sys_dim)
                endif
										goto 105
									endif
									endif
								end do
							endif
							105 continue
                go to 110 

							
							!<2>: For C^+_j(env_block).C_i(sys_site) case
							new_sys_st_num_up=sys_st_num_up-up_dif
							new_sys_st_num_down=sys_st_num_down-down_dif
							new_env_bl_num_up=env_bl_num_up+up_dif
							new_env_bl_num_down=env_bl_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_st_flag=.false.
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==new_sys_st_num_up) then
								if(sys_st%sub(k)%num_down==new_sys_st_num_down) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_bl_flag=.false.
							do k=1,env_bl%num
								if(env_bl%sub(k)%num_up==env_bl_num_up) then
								if(env_bl%sub(k)%num_down==env_bl_num_down) then
									env_bl_id=k
									env_bl_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
										if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(k)%sub(l)%bl_num_down==new_env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											new_env_dim=env_bs%sub(k)%sub(l)%sdim
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue


							if(sys_st_flag.and.env_bl_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_site operator
										SignST=0.0d0
										if(FlagST=='B') then !For Boson
											SignST=1.0d0
										else if(FlagST=='F') then !For Fermion
							sys_bl_num=sys_bl_num_down
								tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=1.0d0
											else
												SignST=-1.0d0
											endif
										endif

										!<b>: For env_block operator
										SignEL=0.0d0
										if(FlagEL=='B') then !For Boson
											SignEL=1.0d0
										else if(FlagEL=='F') then !For Fermion
							sys_num=new_sys_num_down
							env_st_num=env_st_num_down
									tot_num=sys_num+env_st_num
											if(mod(tot_num,2)==0) then
												SignEL=1.0d0
											else
												SignEL=-1.0d0
											endif
										endif

        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down 
        j6=new_sys_st_num_down
        j11=new_env_num_down  
        j22=down_dif1
        j33=env_num_down       
        j44= env_bl_num_down
        j55=env_st_num_down 
        j66=new_env_bl_num_down

        coef1=dsqrt((1.0d0+j11+1)*(j33+1)/(1.0d0+down_dif1))*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        coef1=coef1*(-1)**((j66+j55+j33+j22)/2+(j5+j4+j1+j2)/2)

										Signs=SignST*SignEL*coef1
										coef_tmp=Signs* coef* sys_st%sub(sys_st_id)%mat(1,1)
										call DGEMM('N','T',sys_dim,new_env_dim,old_env_dim,coef_tmp&
												&,invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim&
												&,old_env_pos+1:old_env_pos+old_env_dim),sys_dim&
												&,env_bl%sub(env_bl_id)%mat,new_env_dim,1.0d0&
												&,outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
												&,new_env_pos+1:new_env_pos+new_env_dim),sys_dim)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
						end do
                                        endif
						end do
						end do
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do
end subroutine site_vec_block_ndia


!Outvec=coef*(C^+_i(sys_st)*C_j(env_st)+h.c.)*invec Matrix-wavefunction multiplication
subroutine site_vec_site_ndia(sys_st,FlagST,env_st,FlagET,sys_bs,env_bs,invec,outvec,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: FlagST,FlagET
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,l,x,y,m,n
	logical :: sys_bs_flag,env_bs_flag,sys_st_flag,env_st_flag
	integer :: sys_bs_id,env_bs_id,sys_st_id,env_st_id,up_dif,down_dif
	integer :: sys_bl_num_up,sys_bl_num_down,sys_bl_num,sys_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_st_num_up,sys_st_num_down,new_sys_st_num_up,new_sys_st_num_down
	integer :: env_bl_num_down,env_st_num_up,env_st_num_down,new_env_st_num_up,new_env_st_num_down
	integer :: old_sys_pos,old_env_pos,sys_dim,env_dim,new_sys_pos,new_env_pos
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44
        integer j5, j6,j55, j66, st_down_dif1, st_down_dif2
          real(8),external :: w6js, w3js, w9js
	real(8) :: coef_tmp,SignST,SignET,Signs, coef1,coef11, coef12, coef13, coef14

	!Get general information
	up_dif=sys_st%up_dif
	down_dif1=sys_st%down_dif
	if((up_dif/=env_st%up_dif).or.(down_dif1/=env_st%down_dif)) then
		write(*,*) "Quantum number error in site_vec_site!"
		stop
	endif

	do i=1,invec%num
		sys_num_up=invec%sub(i)%sys_num_up
		sys_num_down=invec%sub(i)%sys_num_down
		env_num_up=invec%sub(i)%env_num_up
		env_num_down=invec%sub(i)%env_num_down
        if(sys_num_down.ne.env_num_down)write(*,*)'wrong1'

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,env_bs%num
					if(env_bs%sub(y)%new_num_up==env_num_up) then
					if(env_bs%sub(y)%new_num_down==env_num_down) then
						do m=1,sys_bs%sub(x)%num
						do n=1,env_bs%sub(y)%num
							sys_bl_num_up=sys_bs%sub(x)%sub(m)%bl_num_up
							sys_bl_num_down=sys_bs%sub(x)%sub(m)%bl_num_down
							sys_st_num_up=sys_bs%sub(x)%sub(m)%st_num_up
							sys_st_num_down=sys_bs%sub(x)%sub(m)%st_num_down
							old_sys_pos=sys_bs%sub(x)%sub(m)%spos
							sys_dim=sys_bs%sub(x)%sub(m)%sdim
							
							!For environment part
							env_bl_num_down=env_bs%sub(y)%sub(n)%bl_num_down
							env_st_num_up=env_bs%sub(y)%sub(n)%st_num_up
							env_st_num_down=env_bs%sub(y)%sub(n)%st_num_down
							old_env_pos=env_bs%sub(y)%sub(n)%spos
							env_dim=env_bs%sub(y)%sub(n)%sdim
							

        do st_down_dif1=-sys_st%down_dif, sys_st%down_dif, su
							sys_st_flag=.false.
							do k=1,sys_st%num
					if(sys_st%sub(k)%num_up==sys_st_num_up) then
					if(sys_st%sub(k)%num_down==sys_st_num_down) then
					if(sys_st%sub(k)%down_dif==st_down_dif1) then
									sys_st_id=k
									sys_st_flag=.true.
				new_sys_st_num_up=sys_st_num_up+up_dif
				new_sys_st_num_down=sys_st_num_down+st_down_dif1
									goto 101
								endif
                                                                endif
								endif
							end do
							101 continue

        do st_down_dif2=-env_st%down_dif, env_st%down_dif, su
							env_st_flag=.false.
					new_env_st_num_up=env_st_num_up-up_dif
					new_env_st_num_down=env_st_num_down-st_down_dif2

							do k=1,env_st%num
					if(env_st%sub(k)%num_up==new_env_st_num_up) then
					if(env_st%sub(k)%num_down==new_env_st_num_down) then
					if(env_st%sub(k)%down_dif==st_down_dif2) then
									env_st_id=k
									env_st_flag=.true.
									goto 102
								endif
								endif
								endif
							end do
							102 continue

				new_sys_num_up=sys_num_up+up_dif
				new_env_num_up=env_num_up-up_dif
				if(sys_st_flag.and.env_st_flag) then
                do new_sys_num_down=abs(sys_num_down-sys_st%down_dif), sys_num_down+sys_st%down_dif, su
							!<1>: For C_i^+(sys_site).C_j(env_site) case
							new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
					if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
					if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
					if(sys_bs%sub(k)%sub(l)%bl_num_down==sys_bl_num_down) then
								new_sys_pos=sys_bs%sub(k)%sub(l)%spos
									sys_bs_flag=.true.
											goto 103
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
							do l=1,env_bs%sub(k)%num
								if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
								if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
								if(env_bs%sub(k)%sub(l)%bl_num_down==env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 104
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

				if(sys_st_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
						do k=1,outvec%num
			if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
			if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_site operator
										SignST=0.0d0
										if(FlagST=='B') then !For Boson
											SignST=1.0d0
										else if(FlagST=='F') then !For Fermion
							sys_bl_num=sys_bl_num_down
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=1.0d0
											else
												SignST=-1.0d0
											endif
										endif

										!<b>: For env_block operator
										SignET=0.0d0
										if(FlagET=='B') then !For Boson
											SignET=1.0d0
										else if(FlagET=='F') then !For Fermion
								sys_num=sys_num_down 
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=1.0d0
											else
												SignET=-1.0d0
											endif
										endif

				Signs=SignST*SignET
				coef_tmp=Signs* coef

        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down 
        j6=new_sys_st_num_down
        j11=new_env_num_down  
        j22=down_dif1        
        j33=env_num_down    
        j44= env_st_num_down  
        j55=env_bl_num_down   
        j66=new_env_st_num_down  


      if(iw6j2(j1,j2,j3-j1,j4,j5,j6)==1)then
        coef1=w6j2(j1,j2,j3-j1,j4,j5,j6)
        else
        coef1=w6js(j1,j2,j3,j4, j5,j6)  
        iw6j2(j1,j2,j3-j1,j4,j5,j6)=1
        w6j2(j1,j2,j3-j1,j4,j5,j6)=coef1
        endif
        if(coef1.eq.0.0)go to 110
      if(iw6j2(j11,j22,j33-j11,j44,j55,j66)==1)then
        coef11=w6j2(j11,j22,j33-j11,j44,j55,j66)
        else
        coef11=w6js(j11,j22,j33,j44, j55,j66) 
        iw6j2(j11,j22,j33-j11,j44,j55,j66)=1
        w6j2(j11,j22,j33-j11,j44,j55,j66)=coef11
        endif
        if(coef11.eq.0.0)go to 110
        coef1=coef1*coef11      

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))
        if(mod((j55+j44+j11+j22)/2+(j5+j4+j1+j2)/2, 2).ne.0)coef1=-coef1
       coef1=coef_tmp*coef1*sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)

		outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
		&=outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
	&+coef1*invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)

if(down_dif1.eq.1)then
 outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)&
                &=outvec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)&
        &+coef1*invec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)
        endif
                endif
										goto 105
									endif
									endif
								end do
							endif
							105 continue
        go to 110
							
							!<2>: For C^+_j(env_site).C_i(sys_site) case
							new_sys_st_num_up=sys_st_num_up-up_dif
							new_sys_st_num_down=sys_st_num_down-down_dif
							new_env_st_num_up=env_st_num_up+up_dif
							new_env_st_num_down=env_st_num_down+down_dif

							new_sys_num_up=sys_num_up-up_dif
							!new_sys_num_down=sys_num_down-down_dif
							new_env_num_up=env_num_up+up_dif
							new_env_num_down=env_num_down+down_dif

							sys_st_flag=.false.
							do k=1,sys_st%num
								if(sys_st%sub(k)%num_up==new_sys_st_num_up) then
								if(sys_st%sub(k)%num_down==new_sys_st_num_down) then
									sys_st_id=k
									sys_st_flag=.true.
									goto 106
								endif
								endif
							end do
							106 continue

							env_st_flag=.false.
							do k=1,env_st%num
								if(env_st%sub(k)%num_up==env_st_num_up) then
								if(env_st%sub(k)%num_down==env_st_num_down) then
									env_st_id=k
									env_st_flag=.true.
									goto 107
								endif
								endif
							end do
							107 continue

							sys_bs_flag=.false.
							do k=1,sys_bs%num
								if(sys_bs%sub(k)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(k)%new_num_down==new_sys_num_down) then
									do l=1,sys_bs%sub(k)%num
										if(sys_bs%sub(k)%sub(l)%st_num_up==new_sys_st_num_up) then
										if(sys_bs%sub(k)%sub(l)%st_num_down==new_sys_st_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(k)%sub(l)%spos
											goto 108
										endif
										endif
									end do
								endif
								endif
							end do
							108 continue

							env_bs_flag=.false.
							do k=1,env_bs%num
								if(env_bs%sub(k)%new_num_up==new_env_num_up) then
								if(env_bs%sub(k)%new_num_down==new_env_num_down) then
									do l=1,env_bs%sub(k)%num
										if(env_bs%sub(k)%sub(l)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(k)%sub(l)%st_num_down==new_env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(k)%sub(l)%spos
											goto 109
										endif
										endif
									end do
								endif
								endif
							end do
							109 continue

							if(sys_st_flag.and.env_st_flag.and.sys_bs_flag.and.env_bs_flag) then
								do k=1,outvec%num
									if((outvec%sub(k)%sys_num_up==new_sys_num_up).and.(outvec%sub(k)%sys_num_down==new_sys_num_down)) then
									if((outvec%sub(k)%env_num_up==new_env_num_up).and.(outvec%sub(k)%env_num_down==new_env_num_down)) then
										!<a>: For sys_site operator
										SignST=0.0d0
										if(FlagST=='B') then !For Boson
											SignST=1.0d0
										else if(FlagST=='F') then !For Fermion
								sys_bl_num=sys_bl_num_down
											tot_num=sys_bl_num
											if(mod(tot_num,2)==0) then
												SignST=1.0d0
											else
												SignST=-1.0d0
											endif
										endif

										!<b>: For env_block operator
										SignET=0.0d0
										if(FlagET=='B') then !For Boson
											SignET=1.0d0
										else if(FlagET=='F') then !For Fermion
						sys_num=new_sys_num_down
											tot_num=sys_num
											if(mod(tot_num,2)==0) then
												SignET=1.0d0
											else
												SignET=-1.0d0
											endif
										endif

										Signs=SignST*SignET
										coef_tmp=Signs* coef* sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
										outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
											&=outvec%sub(k)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+env_dim)&
											&+coef_tmp*invec%sub(i)%vec(old_sys_pos+1:old_sys_pos+sys_dim,old_env_pos+1:old_env_pos+env_dim)
										goto 110
									endif
									endif
								end do
							endif
							110 continue
						end do
                        endif
                                                enddo
						end do
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine site_vec_site_ndia

!=============================================================
subroutine model_wave(sys,env,sys_bs,env_bs,invec,outvec)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Model),intent(in) :: sys,env
	type(Wavefunction),intent(in) :: invec
	type(Wavefunction),intent(inout) :: outvec

	integer :: i,k,x,y,idx,sys_len,env_len
        integer :: ii, is, js, jc, jc1, i1, j, j1,k1
	integer :: sys_xidx,sys_yidx,env_xidx,env_yidx
	real(8) :: coefx,coefy,coefsd,coefsz,coefsn
	real(8) :: maghx,maghy,maghz,bg_time,end_time
	type(Total_Block) :: sys_st,env_st, tmp_bl, tmp_bl1

	call cpu_time(bg_time)

	!<1>: Initiate outvec
	sys_len=outvec%sys_len
	env_len=outvec%env_len
	do i=1,outvec%num
		outvec%sub(i)%vec=0.0d0
	end do

	sys_xidx=Lattice(1,sys_len)
	sys_yidx=Lattice(2,sys_len)
	env_xidx=Lattice(1,Num_site-env_len+1)
	env_yidx=Lattice(2,Num_site-env_len+1)

	!<2>: For hamiltonian term
	call block_vec_dia(sys%ham,sys_bs,invec,outvec,cone)
	call vec_block_dia(env%ham,env_bs,invec,outvec,cone)
	!<3>: For all the two-body interaction term
        
        
                if(add_op.eq.1)then
                        is=nposi_lat(sys%len+1)
                        js=nposi_lat(sys%len+env%len+2)
                if(esite(is).ne.0.0)then
	call site_vec_site_dia(num_elec,st_si,sys_bs,env_bs,invec,outvec,esite(is))
                endif
                if(esite(js).ne.0.0)then
	call site_vec_site_dia(st_si,num_elec,sys_bs,env_bs,invec,outvec,esite(js))
                endif
                                ii=neiba(is,js)
                        if(ii.ne.0)then
                                jz1=jz(ii,is)
                                jd1=jd(ii, is)
                                    coefx=jt(ii, is) !! hoping in either x or y
                                      jn1=jn(ii)


                if(coefx.ne.0.0)then
	call site_vec_site_ndia(st_elec_up,'F',st_elec_up,'F',sys_bs,env_bs,invec,outvec,coefx)
                endif
                if(jd1.ne.0.0)then
	call site_vec_site_ndia(st_sd,'B',st_sd,'B',sys_bs,env_bs,invec,outvec,jd1)
                endif
	call site_vec_site_dia(num_elec,num_elec,sys_bs,env_bs,invec,outvec,jn1)
                                endif
                        if(coup_bs(1,1).eq.1)then  
	call block_site_vec_ndia(sys%snc(1)%elec_up,'F',st_elec_up,'F',sys_bs,invec,outvec,cone)
	call block_site_vec_ndia(sys%snc(1)%spin_sd,'B',st_sd,'B',sys_bs,invec,outvec,cone)
	call block_site_vec_dia(sys%snc(1)%num_sn,num_elec,sys_bs,invec,outvec,cone)
	endif
                        if(coup_bs(1,2).eq.1)then                
  call block_vec_site_ndia(sys%snc(2)%elec_up,'F',st_elec_up,'F',sys_bs,env_bs,invec,outvec,cone)
  call block_vec_site_ndia(sys%snc(2)%spin_sd,'B',st_sd,'B',sys_bs,env_bs,invec,outvec,cone)
  call block_vec_site_dia(sys%snc(2)%num_sn,num_elec,sys_bs,env_bs,invec,outvec,cone)
                                endif

                        if(coup_bs(2,1).eq.1)then
        call block_transfer_trans(env%snc(1)%elec_down, tmp_bl1)
     call site_vec_block_ndia(st_elec_up,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,cone)
                call block_transfer_trans(env%snc(1)%spin_sd, tmp_bl)
          call site_vec_block_ndia(st_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,cone)
                call site_vec_block_dia(num_elec,env%snc(1)%num_sn,sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl)
                        call deallocate_block(tmp_bl1)
                                endif

                                if(coup_bs(2,2).eq.1)then
	call vec_block_site_ndia(env%snc(2)%elec_up,'F',st_elec_up,'F',env_bs,invec,outvec,cone)
	call vec_block_site_ndia(env%snc(2)%spin_sd,'B',st_sd,'B',env_bs,invec,outvec,cone)
	call vec_block_site_dia(env%snc(2)%num_sn,num_elec,env_bs,invec,outvec,cone)
                        endif
                                jc=0
                        do i=1,sys%len+1+env%len
                        if(inb(i).ne.0)then
                                jc=jc+1
                        if(i.le.sys%len)then  !!
                                i1=inb1(i,sys%len)
        call block_transfer_trans(env%sub_s(jc)%elec_down, tmp_bl1)
        call block_vec_block_ndia(sys%sub(i1)%elec_up,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,cone)
                call block_transfer_trans(env%sub_s(jc)%spin_sd, tmp_bl)
	call block_vec_block_ndia(sys%sub(i1)%spin_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,cone)
	call block_vec_block_dia(sys%sub(i1)%num_sn,env%sub_s(jc)%num_sn,sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl)
                        call deallocate_block(tmp_bl1)
                                else
                i1=inb1(i-sys%len-1,env%len)
                        jc1=jc-coup_bb(1)
        call block_transfer_trans(env%sub(i1)%elec_down, tmp_bl1)
	call block_vec_block_ndia(sys%sub_s(jc1)%elec_up,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,cone)
                call block_transfer_trans(env%sub(i1)%spin_sd, tmp_bl)
	call block_vec_block_ndia(sys%sub_s(jc1)%spin_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,cone)
		call block_vec_block_dia(sys%sub_s(jc1)%num_sn,env%sub(i1)%num_sn,sys_bs,env_bs,invec,outvec,cone)
                        call deallocate_block(tmp_bl)
                        call deallocate_block(tmp_bl1)
                               endif
                                endif
                                enddo
101     continue

                else ! add_op=0

                        is=nposi_lat(sys%len+1)
                        js=nposi_lat(sys%len+env%len+2)
                        ii=neiba(is,js)
                if(esite(is).ne.0.0)then
	call site_vec_site_dia(num_elec,st_si,sys_bs,env_bs,invec,outvec,esite(is))
                endif
                if(esite(js).ne.0.0)then
	call site_vec_site_dia(st_si,num_elec,sys_bs,env_bs,invec,outvec,esite(js))
                endif
                        if(ii.ne.0)then

                                        coefx=jt(ii, is) !! hoping in either x or y
                                        coefsz=jz(ii, is)
                                        coefsd=jd(ii, is)
                                        coefsn=jn(ii)
                if(coefx.ne.0.0)then
	call site_vec_site_ndia(st_elec_up,'F',st_elec_up,'F',sys_bs,env_bs,invec,outvec,coefx)
                endif
                if(coefsd.ne.0.0)then
                call block_transfer_trans(st_sd, tmp_bl)
	call site_vec_site_ndia(st_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,coefsd)
                        endif
                if(coefsn.ne.0.0)then
	call site_vec_site_dia(num_elec,num_elec,sys_bs,env_bs,invec,outvec,coefsn)
                                endif
                                endif

        
                        do k1=1, neibt
                        i1=neib(sys%len+1,k1)
                                        coefx=jt(k1, sys%len+1) !! hoping in either x or y
                                        coefsz=jz(k1, sys%len+1)
                                        coefsd=jd(k1, sys%len+1)
                                        coefsn=jn(k1)
                                if(i1.le.sys%len.and.i1.gt.0)then  !!! with sys
                        j1=inb1(nposi_dm(i1), sys%len)
                        if(j1.le.nleg11(sys%len))then
        if(coefx.ne.0.0)then
		call block_site_vec_ndia(sys%sub(j1)%elec_up,'F',st_elec_up,'F',sys_bs,invec,outvec,coefx)
                endif
                if(coefsd.ne.0.0)then
		call block_site_vec_ndia(sys%sub(j1)%spin_sd,'B',st_sd,'B',sys_bs,invec,outvec,coefsd)
                        endif
                if(coefsn.ne.0.0)then
		call block_site_vec_dia(sys%sub(j1)%num_sn,num_elec,sys_bs,invec,outvec,coefsn)
                        endif
                                endif
                        	endif

                                if(i1.gt.sys%len+1)then
                        j1=inb1(nposi_dm(i1)-sys%len-1, env%len)
                        if(j1.le.nleg11(env%len).and.j1.gt.0)then
                        if(coefx.ne.0.0)then
        call block_transfer(st_elec_up, tmp_bl)
        call block_transfer_trans(env%sub(j1)%elec_down, tmp_bl1)
              call site_vec_block_ndia(tmp_bl,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,coefx)
                        endif
                if(coefsd.ne.0.0)then
                        call block_transfer_trans(env%sub(j1)%spin_sd, tmp_bl)
                call site_vec_block_ndia(st_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,coefsd)
                endif
                if(coefsn.ne.0.0)then
                call site_vec_block_dia(num_elec,env%sub(j1)%num_sn,sys_bs,env_bs,invec,outvec,coefsn)
                                endif
                                endif
                                endif
                        i1=neib(nposi_lat(sys%len+env%len+2),k1)
                                        coefx=jt(k1, nposi_lat(sys%len+env%len+2))
                                        coefsz=jz(k1, nposi_lat(sys%len+env%len+2))
                                        coefsd=jd(k1, nposi_lat(sys%len+env%len+2))
                                        coefsn=jn(k1)
                                if(i1.le.sys%len.and.i1.gt.0)then  !!! with sys
                        j1=inb1(nposi_dm(i1), sys%len)
                        if(j1.le.nleg11(sys%len))then
                        if(coefx.ne.0.0)then
		call block_vec_site_ndia(sys%sub(j1)%elec_up,'F',st_elec_up,'F',sys_bs,env_bs,invec,outvec,coefx)
                        endif
		call block_vec_site_ndia(sys%sub(j1)%spin_sd,'B',st_sd,'B',sys_bs,env_bs,invec,outvec,coefsd)
                if(coefsn.ne.0.0)then
		call block_vec_site_dia(sys%sub(j1)%num_sn,num_elec,sys_bs,env_bs,invec,outvec,coefsn)
                                endif
                                endif
                        	endif

                                if(i1.gt.sys%len+1)then

                        j1=inb1(nposi_dm(i1)-sys%len-1, env%len)
                        if(j1.le.nleg11(env%len).and.j1.gt.0)then
                        if(coefx.ne.0.0)then
		call vec_block_site_ndia(env%sub(j1)%elec_up,'F',st_elec_up,'F',env_bs,invec,outvec,coefx)
                        endif
                if(coefsd.ne.0.0)then
		call vec_block_site_ndia(env%sub(j1)%spin_sd,'B',st_sd,'B',env_bs,invec,outvec,coefsd)
                endif
                if(coefsn.ne.0.0)then
		call vec_block_site_dia(env%sub(j1)%num_sn,num_elec,env_bs,invec,outvec,coefsn)
                                endif
                                endif
                                endif

                                enddo

                                jc=0
                        do i=1,sys%len
                                i1=inb1(i,sys%len)
                                do k=1,neibt
                                j=neib(i,k)
                                        coefx=jt(k, i) !! hoping in either x or y
                                        coefsz=jz(k, i)
                                        coefsd=jd(k, i)
                                        coefsn=jn(k)
                        if(j.gt.sys%len+1)then  !!
                                j1=inb1(nposi_dm(j)-sys%len-1,env%len)
                        if(j1.gt.0.and.j1.le.nleg11(env%len))then
                if(coefx.ne.0.0)then
        call block_transfer(sys%sub(i1)%elec_up, tmp_bl)
        call block_transfer_trans(env%sub(j1)%elec_down, tmp_bl1)
	call block_vec_block_ndia(tmp_bl,'F',tmp_bl1,'F',sys_bs,env_bs,invec,outvec,coefx)
                call deallocate_block(tmp_bl)
                call deallocate_block(tmp_bl1)
                endif
                if(coefsd.ne.0.0)then
        call block_transfer_trans(env%sub(j1)%spin_sd, tmp_bl)
	call block_vec_block_ndia(sys%sub(i1)%spin_sd,'B',tmp_bl,'B',sys_bs,env_bs,invec,outvec,coefsd)
                endif
                if(coefsn.ne.0.0)then
	call block_vec_block_dia(sys%sub(i1)%num_sn,env%sub(j1)%num_sn,sys_bs,env_bs,invec,outvec,coefsn)
                                endif
                                endif
                                endif
                                enddo
                                enddo

                endif                                        
111             continue

	call cpu_time(end_time)
end subroutine model_wave

   subroutine hmerge(sys,env,sys_bs,env_bs)
	use pubdata
	implicit none
	type(Total_Basis),intent(in)::sys_bs,env_bs
	type(Total_Model),intent(inout)::sys,env
	type(Total_block)::  tmp_bl
	integer i,k,j,mbsite,ib,ib1,i1,k1,k11,imm, jc, jc1
	integer  ia, syssite, totsite, dmsite 
	real*8 h1
			totsite=(sys%len+env%len+2+nleg-1)/nleg*nleg
                        dmsite=sys%len+env%len+2
			syssite=sys%len+1
		coup_bs(1:2,1:2)=0
		do k=1,neibt
		j=neib(syssite,k) 
		if(nposi(nposi_dm(j)).eq.3)coup_bs(2,1)=1 
		if(nposi(nposi_dm(j)).eq.1)coup_bs(1,1)=1 

		j=neib(nposi_lat(dmsite),k) 
		if(nposi(nposi_dm(j)).eq.3)coup_bs(2,2)=1 
		if(nposi(nposi_dm(j)).eq.1)coup_bs(1,2)=1 
		enddo

				do ia=1,2
				if(coup_bs(1,ia).eq.1)then
        if(allocated(sys%snc(ia)%elec_up%sub))call deallocate_block(sys%snc(ia)%elec_up)
	if(allocated(sys%snc(ia)%spin_sd%sub))call deallocate_block(sys%snc(ia)%spin_sd)
	if(allocated(sys%snc(ia)%num_sn%sub))call deallocate_block(sys%snc(ia)%num_sn)
  call block_pass_info(sys%sub(1)%elec_up,sys%snc(ia)%elec_up)
  call block_pass_info(sys%sub(1)%spin_sd,sys%snc(ia)%spin_sd)
  call block_pass_info(sys%sub(1)%num_sn,  sys%snc(ia)%num_sn)
					endif

				if(coup_bs(2,ia).eq.1)then
        if(allocated(env%snc(ia)%elec_up%sub))call deallocate_block(env%snc(ia)%elec_up)
        if(allocated(env%snc(ia)%elec_down%sub))call deallocate_block(env%snc(ia)%elec_down)
	if(allocated(env%snc(ia)%spin_sd%sub))call deallocate_block(env%snc(ia)%spin_sd)
	if(allocated(env%snc(ia)%num_sn%sub))call deallocate_block(env%snc(ia)%num_sn)
  call block_pass_info(env%sub(1)%elec_up,env%snc(ia)%elec_up)
  call block_pass_info(env%sub(1)%elec_down,env%snc(ia)%elec_down)
  call block_pass_info(env%sub(1)%spin_sd,env%snc(ia)%spin_sd)
  call block_pass_info(env%sub(1)%num_sn,  env%snc(ia)%num_sn)
					endif
                                        enddo
		
		do ia=1,2 
				j=sys%len+1+(env%len+1)*(ia-1) 
			do k1=1,neibt
			Ib=neib(nposi_lat(j),k1) 
			if(nposi(nposi_dm(ib)).eq.3)then	 
			Ib1=nposi_dm(Ib)-sys%len-1 
			Ib1=inb1(Ib1,env%len)
		if(coup_bs(2,ia).ne.0)then
                jz1=jz(k1, nposi_lat(j))
                jd1=jd(k1, nposi_lat(j))   
                jt1=jt(k1,nposi_lat(j))
                jn1=jn(k1)

                if(jt1.ne.0.0)then                        
  call block_add_block_two(env%sub(ib1)%elec_up,jt1,env%snc(ia)%elec_up)
  call block_add_block_two(env%sub(ib1)%elec_down,jt1,env%snc(ia)%elec_down)
                        endif
                        if(jd1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%spin_sd,jd1,env%snc(ia)%spin_sd)
                                endif
                                if(jz1.ne.0.0)then
                                endif
                                if(jn1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%num_sn,jn1,env%snc(ia)%num_sn)
                                        endif
				endif
				endif
			if(nposi(nposi_dm(Ib)).eq.1)then	
                        ib1=nposi_dm(ib)
			Ib1=inb1(Ib1,sys%len)
	if(coup_bs(1,ia).eq.1)then
                jz1=jz(k1, nposi_lat(j))
                jd1=jd(k1, nposi_lat(j))  
                jt1=jt(k1, nposi_lat(j))
                jn1=jn(k1)
                if(jt1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%elec_up,jt1,sys%snc(ia)%elec_up)
                        endif
                        if(jd1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%spin_sd,jd1,sys%snc(ia)%spin_sd)
                                endif
                        if(jn1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%num_sn,jn1,sys%snc(ia)%num_sn)
                                endif
				endif
				endif
				enddo
				enddo

		coup_bb(1:2)=0
		k11=0
           do I=1, sys%len   
		k1=0
		inb(I)=0
		do k=1,neibt
		J=neib(I,k)
		if(nposi(nposi_dm(J)).eq.3)then 
		k1=k1+1
		endif
		enddo
		inb(I)=k1
		if(nleg11(env%len).lt.nleg11(sys%len))inb(i)=0
		if(inb(I).ne.0)k11=k11+1
		enddo
        
			coup_bb(1)=k11
				inb(sys%len+1)=0
				inb(sys%len+1+env%len+1)=0
      			do I=1, env%len  
			k1=0
			inb(i+sys%len+1)=0
			do k=1,neibt
			J=neib(nposi_lat(I+sys%len+1),k) 
			if(nposi(nposi_dm(J)).eq.1)then

			if(inb(J).eq.0)then
				k1=k1+1
			endif
				endif
			enddo
			inb(I+sys%len+1)=k1
		if(inb(I+sys%len+1).ne.0)k11=k11+1
		enddo
		coup_bb(2)=k11-coup_bb(1)

			if(coup_bb(1).ne.0)then
			do i=1,coup_bb(1)
       if(allocated(env%sub_s(i)%elec_up%sub))then
                call deallocate_block(env%sub_s(i)%elec_up)
       if(allocated(env%sub_s(i)%elec_down%sub))call deallocate_block(env%sub_s(i)%elec_down)
                call deallocate_block(env%sub_s(i)%spin_sd)
                call deallocate_block(env%sub_s(i)%num_sn)
                        endif
	call block_pass_info(env%sub(1)%elec_up,env%sub_s(i)%elec_up)
	call block_pass_info(env%sub(1)%elec_down,env%sub_s(i)%elec_down)
	call block_pass_info(env%sub(1)%spin_sd,env%sub_s(i)%spin_sd)
	call block_pass_info(env%sub(1)%num_sn,env%sub_s(i)%num_sn)
			enddo
			endif
			if(coup_bb(2).ne.0)then
			do i=1,coup_bb(2)
	if(allocated(sys%sub_s(i)%elec_up%sub))then
                call deallocate_block(sys%sub_s(i)%elec_up)
                call deallocate_block(sys%sub_s(i)%spin_sd)
                call deallocate_block(sys%sub_s(i)%num_sn)
                        endif
	call block_pass_info(sys%sub(1)%elec_up,sys%sub_s(i)%elec_up)
	call block_pass_info(sys%sub(1)%spin_sd,sys%sub_s(i)%spin_sd)
	call block_pass_info(sys%sub(1)%num_sn,sys%sub_s(i)%num_sn)
				enddo
				endif

			jc=0
		do I=1,sys%len+env%len+1
		if(inb(I).ne.0)then
		Ia=nposi_lat(i)
		if(I.le.sys%len)I1=inb1(I,sys%len)
		if(I.gt.sys%len)I1=inb1(I-sys%len-1,env%len)
		jc=jc+1
			do k1=1,neibt
			Ib=neib(Ia,k1)
		if(nposi(nposi_dm(Ib)).eq.3.and.I.le.sys%len)then	
			Ib1=inb1(nposi_dm(Ib)-sys%len-1,env%len)

                jt1=jt(k1, ia)  
                jd1=jd(k1, ia)  
                jz1=jz(k1, ia)
                jn1=jn(k1)   
        if(jt1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%elec_up,jt1,env%sub_s(jc)%elec_up)
	call block_add_block_two(env%sub(ib1)%elec_down,jt1,env%sub_s(jc)%elec_down)
                endif
        if(jd1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%spin_sd,jd1,env%sub_s(jc)%spin_sd)
                endif
        if(jn1.ne.0.0)then
	call block_add_block_two(env%sub(ib1)%num_sn,jn1,env%sub_s(jc)%num_sn)
                endif
				endif

			if(nposi(nposi_dm(Ib)).eq.1.and.I.gt.sys%len+1)then	
				if(inb(nposi_dm(ib)).eq.0)then
			Ib1=inb1(Ib,sys%len)
				jc1=jc-coup_bb(1)
                jt1=jt(k1, ia) 
                jd1=jd(k1, ia)  
                jz1=jz(k1, ia)
                jn1=jn(k1)   
                if(jt1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%elec_up,jt1,sys%sub_s(jc1)%elec_up)
                        endif

                if(jd1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%spin_sd,jd1,sys%sub_s(jc1)%spin_sd)
                        endif
                                if(jn1.ne.0.0)then
	call block_add_block_two(sys%sub(ib1)%num_sn,jn1,sys%sub_s(jc1)%num_sn)
                                endif
			endif
			endif
			enddo

	endif
	enddo

   end subroutine hmerge

