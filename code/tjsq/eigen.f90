!=============================================================================
subroutine lanczos_memo_new(sys,env,sys_bs,env_bs,super,eigvec,levels,num_idx)
	use pubdata
	implicit none

	integer,intent(in) :: levels,num_idx
	type(Super_basis),intent(in) :: super
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Model),intent(in) :: sys,env
	type(Wavefunction),intent(inout) :: eigvec(levels)

	integer,parameter :: loop_out=1,loop_in=110

	integer :: i,j,k,x,y,ii,iter,last,iter_out,tlen
	real(8) :: ritzmin,lastmin,bg_time,end_time
	real(8) :: alpha1(loop_in),beta1(loop_in)

	real(8),external :: wave_norm
	real(8),external :: wave_product
	type(Wavefunction) :: tempvec(loop_in+1)
	real(8) :: value

	integer :: ldz,info
	real(8),allocatable :: dia(:),ndia(:),dz(:,:),dwork(:)

	call cpu_time(bg_time)

	!==========================================
	!To get all level's eigenstates one by one
	do ii=1,levels
	iter_out=0
	lastmin=1024.0d0

	do while(iter_out<=loop_out)
		iter_out=iter_out+1

		last=1
		alpha1=0.0d0
		beta1=0.0d0
		ritzmin=lastmin

		!===========check===========================
		if(ii>=2) then
			if(iter_out==0) then
				call wave_initialize(eigvec(ii))
			endif
		endif
		!===========================================
		call wave_transfer(eigvec(ii),tempvec(last))
		call wave_normalize(tempvec(last))

		call allocate_wave(super,tempvec(last+1))
                if(add_op.eq.1)then
        call hmerge(sys,env, sys_bs, env_bs)
                endif
			call model_wave(sys,env,sys_bs,env_bs,tempvec(last),tempvec(last+1))

		alpha1(last)=wave_product(tempvec(last),tempvec(last+1))
		ritzmin=alpha1(last)
		
		do while(last<loop_in)
			if(last<2) then
				do i=1,tempvec(last+1)%num
					tempvec(last+1)%sub(i)%vec=tempvec(last+1)%sub(i)%vec-alpha1(last)*tempvec(last)%sub(i)%vec
				end do
			else
				do i=1,tempvec(last+1)%num
					tempvec(last+1)%sub(i)%vec=tempvec(last+1)%sub(i)%vec&
										 &-alpha1(last)*tempvec(last)%sub(i)%vec&
										 &-beta1(last-1)*tempvec(last-1)%sub(i)%vec
				end do
			endif

			beta1(last)=wave_norm(tempvec(last+1))
			if(dabs(beta1(last))<=lanzero) then
				goto 411
			endif
			call wave_normalize(tempvec(last+1))
			!Re-orthogonalization for tempvec(last+1)
			!<a>: For ground state only
			if(ii<=1) then
				call Gram_Schmidt_process(tempvec(1:last),last,tempvec(last+1))
			!<b>: For excited state
			else
				call Gram_Schmidt_Lanczos(eigvec(1:(ii-1)),ii-1,tempvec(1:last),last,tempvec(last+1))
			endif

			last=last+1
			call allocate_wave(super,tempvec(last+1))
				call model_wave(sys,env,sys_bs,env_bs,tempvec(last),tempvec(last+1))
			alpha1(last)=wave_product(tempvec(last),tempvec(last+1))

			allocate(dia(last),ndia(last),dz(last,last),dwork(2*last-2))
			dia(1:last)=alpha1(1:last)
			ndia(1:last)=beta1(1:last)
			call dstev('N',last,dia,ndia,dz,last,dwork,info)

			if((dabs(dia(1)-ritzmin)<=lanzero).or.((last==loop_in).and.(iter_out==loop_out))) then
				ritzmin=dia(1)
				lastmin=ritzmin


				deallocate(dia,ndia,dz,dwork)
				goto 411
			endif

			ritzmin=dia(1)
			lastmin=ritzmin
                        if(last.le.2)then
                write(*,*)dia(1), ritzmin, last
                                endif
			deallocate(dia,ndia,dz,dwork)
		end do

		allocate(dz(last,last),dwork(2*last-2),dia(last),ndia(last))
		dia(1:last)=alpha1(1:last)
		ndia(1:last)=beta1(1:last)
		call dstev('V',last,dia,ndia,dz,last,dwork,info)
		ritzmin=dia(1)
		lastmin=ritzmin

		do i=1,eigvec(ii)%num
			eigvec(ii)%sub(i)%vec=0.0d0
		end do
		do i=1,last
			do x=1,eigvec(ii)%num
				eigvec(ii)%sub(x)%vec=eigvec(ii)%sub(x)%vec+dz(i,1)*tempvec(i)%sub(x)%vec
			end do
		end do
		call wave_normalize(eigvec(ii))

		deallocate(dia,ndia,dz,dwork)
		do i=1,last+1
			call deallocate_wave(tempvec(i))
		end do
	end do

	411 continue
	allocate(dz(last,last),dwork(2*last-2),dia(last),ndia(last))
	dia(1:last)=alpha1(1:last)
	ndia(1:last)=beta1(1:last)
	call dstev('V',last,dia,ndia,dz,last,dwork,info)
	ritzmin=dia(1)
	Energys(ii)=ritzmin

	do i=1,eigvec(ii)%num
		eigvec(ii)%sub(i)%vec=0.0d0
	end do
	do i=1,last
		do x=1,eigvec(ii)%num
			eigvec(ii)%sub(x)%vec=eigvec(ii)%sub(x)%vec+dz(i,1)*tempvec(i)%sub(x)%vec
		end do
	end do
	call wave_normalize(eigvec(ii))

	deallocate(dz,dwork,dia,ndia)
	do i=1,last+1
		call deallocate_wave(tempvec(i))
	end do
	iter=(iter_out-1)*loop_in+last

	!For the loop of level's
	end do
	!==========================================

	!Output result
	tlen=eigvec(1)%sys_len+eigvec(1)%env_len
	call cpu_time(end_time)
	open(11,file="Energy_All.dat",position='append')
		write(11,111) "Sys_len=",eigvec(1)%sys_len,"E0=",(Energys(i),i=1,levels),"Iter=",iter,"time=",end_time-bg_time
		write(*,111) "Sys_len=",eigvec(1)%sys_len,"E0=",(Energys(i),i=1,levels),"Iter=",iter,"time=",end_time-bg_time
	close(11)

111 format(A8,I3,2X,A3,2X,F18.12,2X,A5,I4,2X,A5,2X,F12.6)
end subroutine lanczos_memo_new

