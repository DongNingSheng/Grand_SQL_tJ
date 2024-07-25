module fact1
real*8  fact(0:301)
integer nfac
end module

!=============================================================================
!Get basis between block and new added site to get the new block
!============================================================================
subroutine Get_basis(block,site,basis)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: block,site
	type(Total_Basis),intent(inout) :: basis

	integer :: i,j,ij,i1, x,y,num
	integer,allocatable :: pre_qn(:,:),res_qn(:,:)
	integer,external :: Get_Unique_QN
        integer  s1, s2, ss

	!<1>: Get general information
	basis%len=block%len+site%len
	num=block%num*3   !!! at most enlarged by 3, site%num

	allocate(pre_qn(2,num),res_qn(2,num))
	pre_qn=0
	res_qn=0

        ij=0
	do i=1,block%num
        do j=1, site%num

        s1=abs(block%sub(i)%num_down-site%sub(j)%num_down)
        s2=block%sub(i)%num_down+site%sub(j)%num_down
        do ss=s1, s2, su
        ij=ij+1
	pre_qn(1,ij)=block%sub(i)%num_up+site%sub(j)%num_up
	pre_qn(2,ij)=ss

        enddo
        enddo
        enddo

        num=ij

	!<2>: Sort pre_qn in descending order
	call QuickSort_Integer(pre_qn,num)
	basis%num=get_unique_qn(pre_qn,res_qn,num)

	allocate(basis%sub(basis%num))
	do i=1,basis%num
		basis%sub(i)%new_num_up=res_qn(1,i)
		basis%sub(i)%new_num_down=res_qn(2,i)
	end do
	deallocate(pre_qn,res_qn)

	!<3>: Get basis%sub(:)%num
	do i=1,basis%num
		basis%sub(i)%num=0
		do x=1,block%num
		do y=1,site%num
			if((block%sub(x)%num_up+site%sub(y)%num_up)==basis%sub(i)%new_num_up) then
			s1=abs(block%sub(x)%num_down-site%sub(y)%num_down)
			s2=block%sub(x)%num_down+site%sub(y)%num_down
			if(basis%sub(i)%new_num_down.ge.s1.and.basis%sub(i)%new_num_down.le.s2) then
				basis%sub(i)%num=basis%sub(i)%num+1
			endif
			endif
		end do
		101 continue
		end do
	end do

	do i=1,basis%num
		allocate(basis%sub(i)%sub(basis%sub(i)%num))
		basis%sub(i)%idx=0
		basis%sub(i)%dim=0
		do x=1,block%num
		do y=1,site%num
			if((block%sub(x)%num_up+site%sub(y)%num_up)==basis%sub(i)%new_num_up) then
			s1=abs(block%sub(x)%num_down-site%sub(y)%num_down)
			s2=block%sub(x)%num_down+site%sub(y)%num_down
			if(basis%sub(i)%new_num_down.ge.s1.and.basis%sub(i)%new_num_down.le.s2) then
				basis%sub(i)%idx=basis%sub(i)%idx+1
				basis%sub(i)%sub(basis%sub(i)%idx)%bl_num_up=block%sub(x)%num_up
				basis%sub(i)%sub(basis%sub(i)%idx)%bl_num_down=block%sub(x)%num_down
				basis%sub(i)%sub(basis%sub(i)%idx)%st_num_up=site%sub(y)%num_up
				basis%sub(i)%sub(basis%sub(i)%idx)%st_num_down=site%sub(y)%num_down
				basis%sub(i)%sub(basis%sub(i)%idx)%sdim=block%sub(x)%sdim*site%sub(y)%sdim
				basis%sub(i)%sub(basis%sub(i)%idx)%spos=basis%sub(i)%dim
				basis%sub(i)%dim=basis%sub(i)%dim+block%sub(x)%sdim*site%sub(y)%sdim
			endif
			endif
		end do
		102 continue
		end do
	end do

	!<5>: Get basis%dim
	basis%dim=0
	do i=1,basis%num
		basis%dim=basis%dim+basis%sub(i)%dim
	end do
end subroutine Get_basis


!==============================================================================
!Get super_basis from sys_bs and env_bs, according totoal quantum number
!==============================================================================
subroutine Get_super_basis(sys_bs,env_bs,super,total_up,total_down)
	use pubdata
	implicit none

	integer,intent(in) :: total_up,total_down
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Super_Basis),intent(inout) :: super

	integer :: i,x,y, s1, s2

	super%sys_len=sys_bs%len
	super%env_len=env_bs%len

	super%num=0
	do x=1,sys_bs%num
		do y=1,env_bs%num
		if(env_bs%sub(y)%new_num_up==(total_up-sys_bs%sub(x)%new_num_up)) then
			s1=abs(env_bs%sub(y)%new_num_down-sys_bs%sub(x)%new_num_down)
			s2=env_bs%sub(y)%new_num_down+sys_bs%sub(x)%new_num_down
			if(total_down.ge.s1.and.total_down.le.s2) then
				super%num=super%num+1
			endif
			endif
		end do
		101 continue
	end do

	allocate(super%sub(super%num))
	super%dim=0
	super%idx=0
	do x=1,sys_bs%num
		do y=1,env_bs%num
			if(env_bs%sub(y)%new_num_up==(total_up-sys_bs%sub(x)%new_num_up)) then

			s1=abs(env_bs%sub(y)%new_num_down-sys_bs%sub(x)%new_num_down)
			s2=env_bs%sub(y)%new_num_down+sys_bs%sub(x)%new_num_down
			if(total_down.ge.s1.and.total_down.le.s2) then
				super%idx=super%idx+1
				super%sub(super%idx)%sys_num_up=sys_bs%sub(x)%new_num_up
				super%sub(super%idx)%sys_num_down=sys_bs%sub(x)%new_num_down
				super%sub(super%idx)%env_num_up=env_bs%sub(y)%new_num_up
				super%sub(super%idx)%env_num_down=env_bs%sub(y)%new_num_down
				super%sub(super%idx)%sys_dim=sys_bs%sub(x)%dim
				super%sub(super%idx)%env_dim=env_bs%sub(y)%dim
				super%dim=super%dim+sys_bs%sub(x)%dim*env_bs%sub(y)%dim
			endif
			endif
		end do
		102 continue
	end do

end subroutine Get_super_basis


!Get unique QNs from in_qn and save to out_qn
integer function Get_Unique_QN(in_qn,out_qn,num)
	implicit none

	integer,intent(in) ::num,in_qn(2,num)
	integer,intent(inout) :: out_qn(2,num)

	logical :: flag
	integer :: i,j,start,sublen

	Get_Unique_QN=0
	flag=.false.
	start=0
	sublen=0
	out_qn=0


	do i=1,num
		do j=start+1,num
			if( (in_qn(1,j)==in_qn(1,i)).and.(in_qn(2,j)==in_qn(2,i)) ) then
				flag=.true.
				sublen=sublen+1
			endif
		end do
		if(flag) then
			start=start+sublen
			Get_Unique_QN=Get_Unique_QN+1
			out_qn(1:2,Get_Unique_QN)=in_qn(1:2,i)

			sublen=0
			flag=.false.
		endif
	end do

	return

end function Get_Unique_QN

!======================================================
!Sort interger array in descending order
!======================================================
!For integer data
subroutine QuickSort_Integer(arr,num)
	implicit none

	integer,intent(in) ::  num
	integer,intent(inout) :: arr(2,num)

	integer :: i,x,index,min,large(2),temp(2)
	integer :: tau_qn,spin_qn,tmp_tau_qn,tmp_spin_qn

	!Sort the second number in descending order
	do i=1,num
		do x=i+1,num
			if(arr(2,i)<arr(2,x)) then
				temp(1:2)=arr(1:2,i)
				arr(1:2,i)=arr(1:2,x)
				arr(1:2,x)=temp(1:2)
			endif
		end do
	end do

	!Sort the first number in descending order
	do i=1,num
		do x=i+1,num
			if(arr(2,i)==arr(2,x)) then
				if(arr(1,i)<arr(1,x)) then
					temp(1:2)=arr(1:2,i)
					arr(1:2,i)=arr(1:2,x)
					arr(1:2,x)=temp(1:2)
				endif
			endif
		end do
	end do

end subroutine QuickSort_Integer


!For real data in descending order
recursive subroutine QuickSort_Real(arr,lower,upper,len)
	implicit none

	integer,intent(in) :: len,lower,upper
	real(8),intent(inout) :: arr(len)

	integer :: i,j,index
	real(8) :: min,temp

	!start sort
	i=lower
	j=upper
	index=int((lower+upper)/2.0d0)
	min=arr(index)

	!first sequence
	do while(i<=j)
		do while(arr(i)>min)
			i=i+1
		end do
		do while(arr(j)<min)
			j=j-1
		end do
		if(i<=j) then
			temp=arr(i)
			arr(i)=arr(j)
			arr(j)=temp

			i=i+1
			j=j-1
		endif
	end do

	!recursion
	if(lower<j) then
		call QuickSort_Real(arr,lower,j,len)
	endif
	if(i<upper) then
		call QuickSort_Real(arr,i,upper,len)
	endif

end subroutine QuickSort_Real


!For real(8) data
recursive subroutine QuickSort_Double(arr,arridx,lower,upper,len)
	implicit none

	integer,intent(in) :: len,lower,upper
	integer,intent(inout) :: arridx(len)
	real(8),intent(inout) :: arr(len)

	integer :: i,j,index,tmpidx
	real(8) :: min,temp

	!start sort
	i=lower
	j=upper
	index=int((lower+upper)/2.0d0)
	min=arr(index)

	do while(i<=j)
		do while(arr(i)>min)
			i=i+1
		end do
		do while(arr(j)<min)
			j=j-1
		end do
		if(i<=j) then
			temp=arr(i)
			arr(i)=arr(j)
			arr(j)=temp

			tmpidx=arridx(i)
			arridx(i)=arridx(j)
			arridx(j)=tmpidx

			i=i+1
			j=j-1
		endif
	end do
	!recursion
	if(lower<j) then
		call QuickSort_Double(arr,arridx,lower,j,len)
	endif
	if(i<upper) then
		call QuickSort_Double(arr,arridx,i,upper,len)
	endif

end subroutine QuickSort_Double

!For complex data
recursive subroutine QuickSort_Complex(arr,lower,upper,len)
	implicit none

	integer,intent(in) :: len,lower,upper
	complex*16,intent(inout) :: arr(len)

	integer :: i,j,index
	complex*16 :: min,temp

	!start sort
	i=lower
	j=upper
	index=int((lower+upper)/2.0d0)
	min=arr(index)

	!first sequence
	do while(i<=j)
		do while(real(arr(i))>real(min))
			i=i+1
		end do
		do while(real(arr(j))<real(min))
			j=j-1
		end do
		if(i<=j) then
			temp=arr(i)
			arr(i)=arr(j)
			arr(j)=temp

			i=i+1
			j=j-1
		endif
	end do

	!recursion
	if(lower<j) then
		call QuickSort_Complex(arr,lower,j,len)
	endif
	if(i<upper) then
		call QuickSort_Complex(arr,i,upper,len)
	endif

end subroutine QuickSort_Complex


!=========================================================
!Update block for diagonal operator without truncation
!=========================================================
subroutine update_block_dia(block,basis,new_block)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: block
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_block

	integer :: i,k,x,spos,sdim
        real*8 coef1

	!<1>: Get general information of new_block from basis
	new_block%len=basis%len
	new_block%num=basis%num
	new_block%dim=basis%dim
	new_block%up_dif=0
	new_block%down_dif=0
		
	!<2>: Get basis info for new_block
	allocate(new_block%sub(new_block%num))
	do i=1,new_block%num
		new_block%sub(i)%num_up=basis%sub(i)%new_num_up
		new_block%sub(i)%num_down=basis%sub(i)%new_num_down
		new_block%sub(i)%down_dif=0
		new_block%sub(i)%row_dim=basis%sub(i)%dim
		new_block%sub(i)%sdim=basis%sub(i)%dim

		sdim=new_block%sub(i)%sdim
		allocate(new_block%sub(i)%mat(sdim,sdim))
		new_block%sub(i)%mat=0.0d0
	end do

	!<3>: Get new_block%sub(:)%mat
	do k=1,new_block%num
		do x=1,basis%sub(k)%num
			do i=1,block%num
				if(block%sub(i)%num_up==basis%sub(k)%sub(x)%bl_num_up) then
				if(block%sub(i)%num_down==basis%sub(k)%sub(x)%bl_num_down) then
					spos=basis%sub(k)%sub(x)%spos
					sdim=basis%sub(k)%sub(x)%sdim

                coef1=dsqrt((1.0d0+new_block%sub(k)%num_down)/(1.0d0+block%sub(i)%num_down))
					new_block%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
						&=new_block%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
						&+block%sub(i)%mat(1:sdim,1:sdim)*coef1
					goto 101
				endif
				endif
			end do
			101 continue
		end do
	end do

end subroutine update_block_dia


!==========================================================
!Update site for diagonal operator without truncation
!==========================================================
subroutine update_site_dia(site,basis,new_site)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_site

	integer :: i,k,x,y,spos,sdim
	real(8),allocatable :: uni_eye(:,:)
        real*8 coef1

	!<1>: Get general information of new_site from basis
	new_site%len=basis%len
	new_site%num=basis%num
	new_site%dim=basis%dim
	new_site%up_dif=0
	new_site%down_dif=0
		
	!<2>: Get basis info for new_site
	allocate(new_site%sub(new_site%num))
	do i=1,new_site%num
		new_site%sub(i)%num_up=basis%sub(i)%new_num_up
		new_site%sub(i)%num_down=basis%sub(i)%new_num_down
		new_site%sub(i)%row_dim=basis%sub(i)%dim
		new_site%sub(i)%sdim=basis%sub(i)%dim

		sdim=new_site%sub(i)%sdim
		allocate(new_site%sub(i)%mat(sdim,sdim))
		new_site%sub(i)%mat=0.0d0
	end do

	!<2>: Get new_site%sub(:)%mat
	do k=1,new_site%num
		do x=1,basis%sub(k)%num
			do i=1,site%num
				if(site%sub(i)%num_up==basis%sub(k)%sub(x)%st_num_up) then
				if(site%sub(i)%num_down==basis%sub(k)%sub(x)%st_num_down) then
					spos=basis%sub(k)%sub(x)%spos
					sdim=basis%sub(k)%sub(x)%sdim

					allocate(uni_eye(sdim,sdim))
					uni_eye=0.0d0
					do y=1,sdim
						uni_eye(y,y)=1.0d0
					end do

                coef1=dsqrt((1.0d0+new_site%sub(k)%num_down)/(1.0d0+site%sub(i)%num_down))
					new_site%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
						&=new_site%sub(k)%mat(spos+1:spos+sdim,spos+1:spos+sdim)&
						&+uni_eye(1:sdim,1:sdim)*site%sub(i)%mat(1,1)*coef1

					deallocate(uni_eye)
					goto 101
				endif
				endif
			end do
			101 continue
		end do
	end do
end subroutine update_site_dia

subroutine update_block_ndia_sys(block,basis,new_block,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: block
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_block

	integer :: i,x,y,new_id,block_id
        integer j1,j2,j3,j4,j5,j6
	integer :: lhs_id,lhs_num_up,lhs_num_down0, lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down

	integer :: lhs_dim,rhs_dim,up_dif,down_dif,down_dif1,sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,block_flag
        integer st_num_down
        real*8 coef
  real(8),external :: w6js
        

	!<1>: Get basis info for new_block
	new_block%len=basis%len
	new_block%up_dif=block%up_dif
	new_block%down_dif=block%down_dif
	up_dif=block%up_dif
	down_dif1=block%down_dif

	new_block%num=0
	do i=1,basis%num

        do down_dif=-down_dif1, down_dif1, su  
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
        if(lhs_num_down.lt.0)go to 101

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_block%num=new_block%num+1
			goto 101
			endif
			endif
		end do
		101 continue
        enddo
	end do

	!<2>: Get new_block%sub(:)%qn and new_block%sub(:)%dim
	allocate(new_block%sub(new_block%num))
	new_id=0
	new_block%dim=0
	do x=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		rhs_dim=basis%sub(x)%dim
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_block%dim=new_block%dim+rhs_dim
				
				new_block%sub(new_id)%num_up=rhs_num_up
				new_block%sub(new_id)%num_down=rhs_num_down
				new_block%sub(new_id)%down_dif=down_dif
				new_block%sub(new_id)%row_dim=lhs_dim
				new_block%sub(new_id)%sdim=rhs_dim
				allocate(new_block%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_block%sub(new_id)%mat=0.0d0
				goto 102
			endif
			end if
		end do
		102 continue
	end do
	end do

	!<3>: Get new_block%sub(:)%mat
	do i=1,new_block%num
		rhs_num_up=new_block%sub(i)%num_up
		rhs_num_down=new_block%sub(i)%num_down
		lhs_dim=new_block%sub(i)%row_dim
		rhs_dim=new_block%sub(i)%sdim
		
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

		down_dif=new_block%sub(i)%down_dif
		
		lhs_num_up=rhs_num_up+up_dif
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
                         st_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
			rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
			rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
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

				lhs_sub_flag=.false.
                        if(block_flag)then
				do y=1,basis%sub(lhs_id)%num
				if(basis%sub(lhs_id)%sub(y)%bl_num_up==lhs_sub_num_up) then
				if(basis%sub(lhs_id)%sub(y)%bl_num_down==lhs_sub_num_down) then
				if(basis%sub(lhs_id)%sub(y)%st_num_down==st_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue
                        endif

				if(lhs_sub_flag.and.block_flag) then


        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        j1=lhs_num_down 
        j2=down_dif1      
        j3=rhs_num_down  
        j4=rhs_sub_num_down  
        j5=st_num_down     
        j6=lhs_sub_num_down  
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*dsqrt((1.0d0+j1)*(1.0d0+j3))
        coef=coef*(-1)**((j6+j5+j3+j2)/2)

        if(coef.ne.0.0)then
	new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&=new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				&+coef*block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim)
                        endif
				endif
			end do
                enddo
		endif
        enddo

end subroutine update_block_ndia_sys

!================================================================
!FlagSign='F' for Fermion operator, else 'B' for Bosonic operator
!================================================================
subroutine update_block_ndia_env(block,basis,new_block,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: block
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_block

	integer :: i,x,y,new_id,block_id
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down
        integer j1,j2,j3,j4,j5,j6

	real(8) :: coef,signs
        real(8),external :: w6js
	integer :: lhs_dim,rhs_dim,up_dif,down_dif, down_dif1,sub_down_dif
	integer :: st_num_up,st_num_down,st_num
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,block_flag

	!<1>: Get basis info for new_block
	new_block%len=basis%len
	new_block%up_dif=block%up_dif
	new_block%down_dif=block%down_dif
	up_dif=block%up_dif
	down_dif1=block%down_dif

	new_block%num=0
	do i=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_block%num=new_block%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get new_block%sub(:)%qn and new_block%sub(:)%dim
	allocate(new_block%sub(new_block%num))
	new_id=0
	new_block%dim=0
	do x=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		rhs_dim=basis%sub(x)%dim
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_block%dim=new_block%dim+rhs_dim
				
				new_block%sub(new_id)%num_up=rhs_num_up
				new_block%sub(new_id)%num_down=rhs_num_down
				new_block%sub(new_id)%down_dif=down_dif
				new_block%sub(new_id)%row_dim=lhs_dim
				new_block%sub(new_id)%sdim=rhs_dim
				allocate(new_block%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_block%sub(new_id)%mat=0.0d0
				goto 102
			endif
			end if
		end do
		102 continue
	end do
	end do

	!<3>: Get new_block%sub(:)%mat
	do i=1,new_block%num
		rhs_num_up=new_block%sub(i)%num_up
		rhs_num_down=new_block%sub(i)%num_down
		lhs_dim=new_block%sub(i)%row_dim
		rhs_dim=new_block%sub(i)%sdim
		
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
		
         down_dif=new_block%sub(i)%down_dif  
		lhs_num_up=rhs_num_up+up_dif
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
				st_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				st_num_down=basis%sub(rhs_id)%sub(x)%st_num_down

				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                do sub_down_dif=-down_dif1, down_dif1,su

				block_flag=.false.
				do y=1,block%num
					if(block%sub(y)%num_up==rhs_sub_num_up) then
					if(block%sub(y)%num_down==rhs_sub_num_down) then
					if(block%sub(y)%down_dif==sub_down_dif) then
						block_flag=.true.
						block_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
		lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						goto 106
					endif
					endif
                                        endif
				end do
				106 continue

                        
			
        	lhs_sub_flag=.false.
                if(block_flag)then
				do y=1,basis%sub(lhs_id)%num
					if(basis%sub(lhs_id)%sub(y)%bl_num_up==lhs_sub_num_up) then
					if(basis%sub(lhs_id)%sub(y)%bl_num_down==lhs_sub_num_down) then
					if(basis%sub(lhs_id)%sub(y)%st_num_down==st_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue

                        endif



				if(lhs_sub_flag.and.block_flag) then

					Signs=0.0d0
					if(FlagSign=='B') then !For Boson
						Signs=1.0d0
					else if(FlagSign=='F') then !For Fermion
						st_num=st_num_down  
						if(mod(st_num,2)==0) then
							Signs=1.0d0
						else
							Signs=-1.0d0
						endif
					endif

        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        j1=lhs_num_down 
        j2=down_dif1     
        j3=rhs_num_down  
        j4=rhs_sub_num_down  
        j5=st_num_down     
        j6=lhs_sub_num_down  
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*(-1)**((j6+j5+j3+j2)/2)*dsqrt((1.0d0+j1)*(1.0d0+j3))

                if(coef.ne.0.0)Then
					coef=coef*Signs
	new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
				       &=new_block%sub(i)%mat(lhs_pos+1:lhs_pos+lhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&+coef*block%sub(block_id)%mat(1:lhs_sub_dim,1:rhs_sub_dim)
                        endif
				endif
			end do
			end do
		endif
	end do

end subroutine update_block_ndia_env




!================================================================
!FlagSign='F' for Fermion operator, else 'B' for Bosonic operator
!================================================================
subroutine update_site_ndia_sys(site,basis,new_site,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_site

	real(8) :: coef,signs
        real(8),external :: w6js
        integer j1,j2,j3,j4,j5,j6
	integer :: i,x,y,new_id,site_id
	integer :: bl_num_up,bl_num_down,bl_num
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down

	integer :: lhs_dim,rhs_dim,up_dif,down_dif, down_dif1,sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,site_flag
	real(8),allocatable :: uni_eye(:,:)

	!<1>: Get basis information for new_site
	new_site%len=basis%len
	new_site%up_dif=site%up_dif
	new_site%down_dif=site%down_dif
	up_dif=site%up_dif
	down_dif1=site%down_dif

	new_site%num=0
	do i=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_site%num=new_site%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get new_site%sub(:)%new_qn,%dim
	allocate(new_site%sub(new_site%num))
	new_id=0
	new_site%dim=0
	do x=1,basis%num
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		rhs_dim=basis%sub(x)%dim

        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_site%dim=new_site%dim+rhs_dim

				new_site%sub(new_id)%num_up=rhs_num_up
				new_site%sub(new_id)%num_down=rhs_num_down
				new_site%sub(new_id)%down_dif=down_dif
				new_site%sub(new_id)%row_dim=lhs_dim
				new_site%sub(new_id)%sdim=rhs_dim
				allocate(new_site%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_site%sub(new_id)%mat=0.0d0
				goto 102
			endif
			endif
		end do
		102 continue
	end do
	end do

	!<3>: Get new_site%sub(:)%mat
	do i=1,new_site%num
		rhs_num_up=new_site%sub(i)%num_up
		rhs_num_down=new_site%sub(i)%num_down
		rhs_dim=new_site%sub(i)%sdim

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
		
         down_dif=new_site%sub(i)%down_dif   !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
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
				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
				bl_num_up=basis%sub(rhs_id)%sub(x)%bl_num_up
				bl_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                do sub_down_dif=-down_dif1, down_dif1,su
				site_flag=.false.
				do y=1,site%num
					if(site%sub(y)%num_up==rhs_sub_num_up) then
					if(site%sub(y)%num_down==rhs_sub_num_down) then
					if(site%sub(y)%down_dif==sub_down_dif) then
						site_flag=.true.
						site_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
				lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						goto 106
					endif
                                        endif
					endif
				end do
				106 continue

				lhs_sub_flag=.false.
                if(site_flag)then
				do y=1,basis%sub(lhs_id)%num
					if(basis%sub(lhs_id)%sub(y)%st_num_up==lhs_sub_num_up) then
					if(basis%sub(lhs_id)%sub(y)%st_num_down==lhs_sub_num_down) then
					if(basis%sub(lhs_id)%sub(y)%bl_num_down==bl_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue
                        endif

				if(lhs_sub_flag.and.site_flag) then
					allocate(uni_eye(rhs_sub_dim,rhs_sub_dim))
					uni_eye=0.0d0
					do y=1,rhs_sub_dim
						uni_eye(y,y)=1.0d0
					end do

					Signs=0.0d0
					if(FlagSign=='B') then !For Boson
						Signs=1.0d0
					else if(FlagSign=='F') then !For Fermion
						bl_num=bl_num_down !!+bl_num_down
						if(mod(bl_num,2)==0) then
							Signs=1.0d0
						else
							Signs=-1.0d0
						endif
					endif

        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}

        j1=lhs_num_down 
        j2=down_dif1      
        j3=rhs_num_down  
        j4=rhs_sub_num_down  
        j5=bl_num_down     
        j6=lhs_sub_num_down  
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*(-1)**((j5+j4+j1+j2)/2)*dsqrt((1.0d0+j1)*(1.0d0+j3))

                if(coef.ne.0.0)then
					coef=coef*Signs
					new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&=new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&+coef*uni_eye(1:rhs_sub_dim,1:rhs_sub_dim)*site%sub(site_id)%mat(1,1)
                endif

					deallocate(uni_eye)
				endif
			end do
			end do
		endif
	end do
end subroutine update_site_ndia_sys

subroutine update_site_ndia_env(site,basis,new_site,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Block),intent(in) :: site
	type(Total_Basis),intent(in) :: basis
	type(Total_Block),intent(inout) :: new_site

	integer :: i,x,y,new_id,site_id
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: rhs_id,rhs_num_up,rhs_num_down
        integer bl_num_down
        integer j1,j2,j3,j4,j5,j6
        real(8),external :: w6js

	real(8) :: coef,signs
	integer :: lhs_dim,rhs_dim,up_dif,down_dif, down_dif1, sub_down_dif
	integer :: lhs_sub_num_up,lhs_sub_num_down
	integer :: rhs_sub_num_up,rhs_sub_num_down
	integer :: lhs_sub_dim,rhs_sub_dim,lhs_pos,rhs_pos
	logical :: lhs_flag,rhs_flag,lhs_sub_flag,site_flag
	real(8),allocatable :: uni_eye(:,:)

	!<1>: Get basis information for new_site
	new_site%len=basis%len
	new_site%up_dif=site%up_dif
	new_site%down_dif=site%down_dif
	up_dif=site%up_dif
	down_dif1=site%down_dif

	new_site%num=0
	do i=1,basis%num
        do down_dif=-down_dif1, down_dif1, su  
		rhs_num_up=basis%sub(i)%new_num_up
		rhs_num_down=basis%sub(i)%new_num_down
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==lhs_num_up) then
			if(basis%sub(x)%new_num_down==lhs_num_down) then
				new_site%num=new_site%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get new_site%sub(:)%new_qn,%dim
	allocate(new_site%sub(new_site%num))
	new_id=0
	new_site%dim=0
	do x=1,basis%num
		rhs_num_up=basis%sub(x)%new_num_up
		rhs_num_down=basis%sub(x)%new_num_down
		rhs_dim=basis%sub(x)%dim
		
        do down_dif=-down_dif1, down_dif1, su  !!! triangle rule
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		do y=1,basis%num
			if(basis%sub(y)%new_num_up==lhs_num_up) then
			if(basis%sub(y)%new_num_down==lhs_num_down) then
				new_id=new_id+1
				lhs_dim=basis%sub(y)%dim
				new_site%dim=new_site%dim+rhs_dim

				new_site%sub(new_id)%num_up=rhs_num_up
				new_site%sub(new_id)%num_down=rhs_num_down
				new_site%sub(new_id)%down_dif=down_dif
				new_site%sub(new_id)%row_dim=lhs_dim
				new_site%sub(new_id)%sdim=rhs_dim
				allocate(new_site%sub(new_id)%mat(lhs_dim,rhs_dim))
				new_site%sub(new_id)%mat=0.0d0
				goto 102
			endif
			endif
		end do
		102 continue
	end do
	end do

	!<3>: Get new_site%sub(:)%mat
	do i=1,new_site%num
		rhs_num_up=new_site%sub(i)%num_up
		rhs_num_down=new_site%sub(i)%num_down
		rhs_dim=new_site%sub(i)%sdim

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
		
         down_dif=new_site%sub(i)%down_dif 
		lhs_num_up=rhs_num_up+up_dif
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
				rhs_sub_num_up=basis%sub(rhs_id)%sub(x)%st_num_up
				rhs_sub_num_down=basis%sub(rhs_id)%sub(x)%st_num_down
				bl_num_down=basis%sub(rhs_id)%sub(x)%bl_num_down
				rhs_sub_dim=basis%sub(rhs_id)%sub(x)%sdim
				rhs_pos=basis%sub(rhs_id)%sub(x)%spos


                        do sub_down_dif=-down_dif1, down_dif1, su
				site_flag=.false.
				do y=1,site%num
					if(site%sub(y)%num_up==rhs_sub_num_up) then
					if(site%sub(y)%num_down==rhs_sub_num_down) then
					if(site%sub(y)%down_dif==sub_down_dif) then
						site_flag=.true.
						site_id=y
				lhs_sub_num_up=rhs_sub_num_up+up_dif
				lhs_sub_num_down=rhs_sub_num_down+sub_down_dif
						goto 106
					endif
					endif
					endif
				end do
				106 continue
				lhs_sub_flag=.false.
                        if(site_flag)then
				do y=1,basis%sub(lhs_id)%num
					if(basis%sub(lhs_id)%sub(y)%st_num_up==lhs_sub_num_up) then
					if(basis%sub(lhs_id)%sub(y)%st_num_down==lhs_sub_num_down) then
					if(basis%sub(lhs_id)%sub(y)%bl_num_down==bl_num_down) then
						lhs_sub_flag=.true.
						lhs_sub_dim=basis%sub(lhs_id)%sub(y)%sdim
						lhs_pos=basis%sub(lhs_id)%sub(y)%spos
						goto 105
					endif
					endif
					endif
				end do
				105 continue
                        endif

				if(lhs_sub_flag.and.site_flag) then
					allocate(uni_eye(rhs_sub_dim,rhs_sub_dim))
					uni_eye=0.0d0
					do y=1,rhs_sub_dim
						uni_eye(y,y)=1.0d0
					end do
        !! Wigner 6j coef {j1,j2,j3, j4,j5,j6}
        j1=lhs_num_down 
        j2=down_dif1      
        j3=rhs_num_down  
        j4=rhs_sub_num_down  
        j5=bl_num_down     
        j6=lhs_sub_num_down  
        coef=w6js(j1,j2,j3,j4, j5,j6)
        coef=coef*dsqrt((1.0d0+j1)*(1.0d0+j3))*(-1)**((j5+j4+j1+j2)/2)
                if(coef.ne.0.0)then
					new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&=new_site%sub(i)%mat(lhs_pos+1:lhs_pos+rhs_sub_dim,rhs_pos+1:rhs_pos+rhs_sub_dim)&
						&+coef*uni_eye(1:rhs_sub_dim,1:rhs_sub_dim)*site%sub(site_id)%mat(1,1)

                        endif
					deallocate(uni_eye)
				endif
			end do
			end do
		endif
	end do
end subroutine update_site_ndia_env


!================================================
!Update with truncation for number operator
!================================================
subroutine update_trun_dia(ori,eff,trun)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: ori,trun
	type(Total_Block),intent(inout) :: eff

	integer :: i,x,ori_dim,eff_dim
	integer :: eff_num_up,eff_num_down,trun_id
	real(8),allocatable :: mid(:,:)

	!<1>: Get eff%len,%num and eff%sub(:)%qn,%dim
	eff%len=ori%len
	eff%num=trun%num
	eff%dim=trun%dim
	eff%up_dif=ori%up_dif
	eff%down_dif=ori%down_dif
	allocate(eff%sub(eff%num))
	do i=1,eff%num
		eff%sub(i)%num_up=trun%sub(i)%num_up
		eff%sub(i)%num_down=trun%sub(i)%num_down
		eff%sub(i)%down_dif=0
		eff%sub(i)%row_dim=trun%sub(i)%sdim
		eff%sub(i)%sdim=trun%sub(i)%sdim
		
		eff_dim=eff%sub(i)%sdim
		allocate(eff%sub(i)%mat(eff_dim,eff_dim))
		eff%sub(i)%mat=0.0d0
	end do

	!<2>: Get eff%sub(:)%mat
	do i=1,eff%num
		eff_num_up=eff%sub(i)%num_up
		eff_num_down=eff%sub(i)%num_down
		eff_dim=eff%sub(i)%sdim
		do x=1,ori%num
			if(ori%sub(x)%num_up==eff_num_up) then
			if(ori%sub(x)%num_down==eff_num_down) then
				ori_dim=ori%sub(x)%row_dim
				allocate(mid(ori_dim,eff_dim))
				call DGEMM('N','N',ori_dim,eff_dim,ori_dim,1.0d0,ori%sub(x)%mat&
						&,ori_dim,trun%sub(i)%mat,ori_dim,0.0d0,mid,ori_dim)
				call DGEMM('T','N',eff_dim,eff_dim,ori_dim,1.0d0,trun%sub(i)%mat&
						&,ori_dim,mid,ori_dim,0.0d0,eff%sub(i)%mat,eff_dim)
				deallocate(mid)
				goto 101
			endif
			endif
		end do
		101 continue
	end do

end subroutine update_trun_dia

subroutine update_trun_ndia(ori,eff,trun)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: ori,trun
	type(Total_Block),intent(inout) :: eff

	integer :: i,j,k,x,y
	logical :: ori_flag,lhs_flag,rhs_flag
	integer :: rhs_dim,lhs_dim,up_dif,down_dif, down_dif1
	integer :: rhs_id,rhs_num_up,rhs_num_down
	integer :: lhs_id,lhs_num_up,lhs_num_down
	integer :: ori_id,eff_id,ori_row,ori_col
	real(8),allocatable :: mid(:,:)

        
	!<1>: Get eff%len and eff%num
	eff%len=ori%len
	eff%up_dif=ori%up_dif
	eff%down_dif=ori%down_dif
	up_dif=ori%up_dif
	down_dif=ori%down_dif
	down_dif1=ori%down_dif

	eff%num=0
	do x=1,trun%num
		rhs_num_up=trun%sub(x)%num_up
		rhs_num_down=trun%sub(x)%num_down
        do down_dif=-down_dif1, down_dif1, su  
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif

		do y=1,trun%num
			if(trun%sub(y)%num_up==lhs_num_up) then
			if(trun%sub(y)%num_down==lhs_num_down) then
				eff%num=eff%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do

	!<2>: Get eff%sub(:)%qn,%dim
	allocate(eff%sub(eff%num))
	eff_id=0
	eff%dim=0
	do x=1,trun%num
		rhs_num_up=trun%sub(x)%num_up
		rhs_num_down=trun%sub(x)%num_down
		rhs_dim=trun%sub(x)%sdim

        do down_dif=-down_dif1, down_dif1, su   
		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		do y=1,trun%num
			if(trun%sub(y)%num_up==lhs_num_up) then
			if(trun%sub(y)%num_down==lhs_num_down) then
				eff_id=eff_id+1
				lhs_dim=trun%sub(y)%sdim
				eff%dim=eff%dim+rhs_dim

				eff%sub(eff_id)%num_up=rhs_num_up
				eff%sub(eff_id)%num_down=rhs_num_down
				eff%sub(eff_id)%down_dif=down_dif
				eff%sub(eff_id)%row_dim=lhs_dim
				eff%sub(eff_id)%sdim=rhs_dim
				allocate(eff%sub(eff_id)%mat(lhs_dim,rhs_dim))
				eff%sub(eff_id)%mat=0.0d0
				goto 102
			endif
			endif
		end do
		102 continue
	end do
	end do

	!<3>: Get eff%sub(:)%mat
	do i=1,eff%num
		rhs_num_up=eff%sub(i)%num_up
		rhs_num_down=eff%sub(i)%num_down

		rhs_flag=.false.
		do x=1,trun%num
			if(trun%sub(x)%num_up==rhs_num_up) then
			if(trun%sub(x)%num_down==rhs_num_down) then
				rhs_flag=.true.
				rhs_id=x
				goto 103
			endif
			endif
		end do
		103 continue
         down_dif=eff%sub(i)%down_dif

		lhs_num_up=rhs_num_up+up_dif
		lhs_num_down=rhs_num_down+down_dif
		lhs_flag=.false.
		do x=1,trun%num
			if(trun%sub(x)%num_up==lhs_num_up) then
			if(trun%sub(x)%num_down==lhs_num_down) then
				lhs_flag=.true.
				lhs_id=x
				goto 104
			endif
			endif
		end do
		104 continue

		ori_flag=.false.
		do x=1,ori%num
			if(ori%sub(x)%num_up==rhs_num_up) then
			if(ori%sub(x)%num_down==rhs_num_down) then
			if(ori%sub(x)%down_dif==down_dif) then
				ori_flag=.true.
				ori_id=x
				goto 105
			endif
			endif
			endif
		end do
		105 continue

		if(rhs_flag.and.lhs_flag.and.ori_flag) then
			ori_row=ori%sub(ori_id)%row_dim
			ori_col=ori%sub(ori_id)%sdim
			rhs_dim=trun%sub(rhs_id)%sdim
			lhs_dim=trun%sub(lhs_id)%sdim
			allocate(mid(ori_row,rhs_dim))
			call DGEMM('N','N',ori_row,rhs_dim,ori_col,1.0d0,ori%sub(ori_id)%mat&
					&,ori_row,trun%sub(rhs_id)%mat,ori_col,0.0d0,mid,ori_row)
			call DGEMM('T','N',lhs_dim,rhs_dim,ori_row,1.0d0,trun%sub(lhs_id)%mat&
					&,ori_row,mid,ori_row,0.0d0,eff%sub(i)%mat,lhs_dim)
			deallocate(mid)
		endif
	end do

end subroutine update_trun_ndia


!main

SUBROUTINE FACTRL
use fact1
NFAC=142
FACT(0)=1.0
DO 10 I=1,NFAC
10 FACT(I)=FACT(I-1)*FLOAT(I)
write(*,*)fact(nfac)
RETURN
END





real*8 FUNCTION W3JS(J1,J2,J3,M1,M2,M3)
use fact1
INTEGER Z,ZMIN,ZMAX
real*8 denom, cc, cc1, cc2
W3JS=0.0d0
!w3js=1
!go to 1000
IF(M1+M2+M3.NE.0) GOTO 1000
IA=J1+J2
IF(J3.GT.IA) GOTO 1000
IB=J1-J2
IF(J3.LT.IABS(IB)) GOTO 1000
JSUM=J3+IA
IC=J1-M1
ID=J2-M2
IF(MOD(JSUM,2).NE.0) GOTO 1000
IF(MOD(IC,2).NE.0) GOTO 1000
IF(MOD(ID,2).NE.0) GOTO 1000
IF(IABS(M1).GT.J1) GOTO 1000
IF(IABS(M2).GT.J2) GOTO 1000
IF(IABS(M3).GT.J3) GOTO 1000
IE=J3-J2+M1
IF=J3-J1-M2
ZMIN=MAX0(0,-IE,-IF)
IG=IA-J3
IH=J2+M2
ZMAX=MIN0(IG,IH,IC)
CC=0.0
DO 200 Z=ZMIN,ZMAX,2
DENOM=FACT(Z/2)*FACT((IG-Z)/2)*FACT((IC-Z)/2)*FACT((IH-Z)/2)*FACT((IE+Z)/2)*FACT((IF+Z)/2)
IF(MOD(Z,4).NE.0) DENOM=-DENOM
CC=CC+1.0/DENOM
200 CONTINUE
CC1=FACT(IG/2)*FACT((J3+IB)/2)*FACT((J3-IB)/2)/FACT((JSUM+2)/2)
CC2=FACT((J1+M1)/2)*FACT(IC/2)*FACT(IH/2)*FACT(ID/2)*FACT((J3-M3)/2)*FACT((J3+M3)/2)
CC=CC*SQRT(CC1*CC2)
IF(MOD(IB-M3,4).NE.0) CC=-CC
W3JS=CC
1000 RETURN
END


real*8 FUNCTION W6JS(J1,J2,J3,L1,L2,L3)
use fact1
INTEGER W,WMIN,WMAX
INTEGER SUM1,SUM2,SUM3,SUM4
real*8 denom, omega, theta, theta1, theta2, theta3, theta4
W6JS=0.0
!w6js=1
!go to 1000

IA=J1+J2
IF(IA.LT.J3) GOTO 1000
IB=J1-J2
IF(IABS(IB).GT.J3) GOTO 1000
IC=J1+L2
IF(IC.LT.L3) GOTO 1000
ID=J1-L2
IF(IABS(ID).GT.L3) GOTO 1000
IE=L1+J2
IF(IE.LT.L3) GOTO 1000
IF=L1-J2
IF(IABS(IF).GT.L3) GOTO 1000
IG=L1+L2
IF(IG.LT.J3) GOTO 1000
IH=L1-L2
IF(IABS(IH).GT.J3) GOTO 1000
SUM1=IA+J3
SUM2=IC+L3
SUM3=IE+L3
SUM4=IG+J3
IF(MOD(SUM1,2).NE.0) GOTO 1000
IF(MOD(SUM2,2).NE.0) GOTO 1000
IF(MOD(SUM3,2).NE.0) GOTO 1000
IF(MOD(SUM4,2).NE.0) GOTO 1000
WMIN=MAX0(SUM1,SUM2,SUM3,SUM4)
II=IA+IG
IJ=J2+J3+L2+L3
IK=J3+J1+L3+L1
WMAX=MIN0(II,IJ,IK)
OMEGA=0.0
DO 200 W=WMIN,WMAX,2
DENOM=FACT((W-SUM1)/2)*FACT((W-SUM2)/2)*FACT((W-SUM3)/2)*FACT((W-SUM4)/2)*FACT((II-W)/2)*FACT((IJ-W)/2)*FACT((IK-W)/2)
IF(MOD(W,4).NE.0) DENOM=-DENOM
OMEGA=OMEGA+FACT(W/2+1)/DENOM
200 CONTINUE
THETA1=FACT((IA-J3)/2)*FACT((J3+IB)/2)*FACT((J3-IB)/2)/FACT(SUM1/2+1)
THETA2=FACT((IC-L3)/2)*FACT((L3+ID)/2)*FACT((L3-ID)/2)/FACT(SUM2/2+1)
THETA3=FACT((IE-L3)/2)*FACT((L3+IF)/2)*FACT((L3-IF)/2)/FACT(SUM3/2+1)
THETA4=FACT((IG-J3)/2)*FACT((J3+IH)/2)*FACT((J3-IH)/2)/FACT(SUM4/2+1)
THETA=THETA1*THETA2*THETA3*THETA4
W6JS=OMEGA*SQRT(THETA)
1000 RETURN
END
real*8 FUNCTION W9JS(J1,J2,J3,J4,J5,J6,J7,J8,J9)
 real(8),external :: w6js
real*8 x, x1, x2, x3, s
!w9js=1
!go to 11
KMIN=ABS(J1-J9)
KMAX=J1+J9
I=ABS(J4-J8)
IF(I.GT.KMIN) KMIN=I
I=J4+J8
IF(I.LT.KMAX) KMAX=I
I=ABS(J2-J6)
IF(I.GT.KMIN) KMIN=I
I=J2+J6
IF(I.LT.KMAX) KMAX=I
X=0.
DO 1 K=KMIN,KMAX,2
S=1.
IF(MOD(K,2).NE.0) S=-1.
X1=W6JS(J1,J9,K,J8,J4,J7)
X2=W6JS(J2,J6,K,J4,J8,J5)
X3=W6JS(J1,J9,K,J6,J2,J3)
X=X+S*X1*X2*X3*FLOAT(K+1)
1 CONTINUE
W9JS=X
11      continue
RETURN
END

subroutine block_transfer_basis(basis,outmat)
        use pubdata
        implicit none

        type(Total_Basis),intent(in) :: basis
        type(Total_Block),intent(inout) :: outmat

        integer :: i,row_dim,col_dim

        outmat%len=basis%len
        outmat%num=basis%num
        outmat%dim=basis%dim
        outmat%up_dif=0
        outmat%down_dif=0

        if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
        do i=1,outmat%num
                outmat%sub(i)%num_up=basis%sub(i)%new_num_up
                outmat%sub(i)%num_down=basis%sub(i)%new_num_down
                outmat%sub(i)%row_dim=basis%sub(i)%dim
                outmat%sub(i)%sdim=basis%sub(i)%dim

        if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
                outmat%sub(i)%mat=0.d0
        end do

end subroutine block_transfer_basis


subroutine block_transfer_trans(inmat,outmat)
        use pubdata
        implicit none

        type(Total_Block),intent(in) :: inmat
        type(Total_Block),intent(inout) :: outmat

        integer :: i,row_dim,col_dim, up_dif, down_dif

        outmat%len=inmat%len
        outmat%num=inmat%num
        outmat%dim=0
        outmat%up_dif=-inmat%up_dif
        outmat%down_dif=inmat%down_dif
        up_dif=inmat%up_dif
        down_dif=inmat%down_dif

        if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
        do i=1,outmat%num
                outmat%sub(i)%num_up=inmat%sub(i)%num_up+up_dif
                outmat%sub(i)%num_down=inmat%sub(i)%num_down+inmat%sub(i)%down_dif
                outmat%sub(i)%down_dif=-inmat%sub(i)%down_dif
                outmat%sub(i)%row_dim=inmat%sub(i)%sdim
                outmat%sub(i)%sdim=inmat%sub(i)%row_dim
                outmat%dim=outmat%dim+outmat%sub(i)%sdim
        if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
                outmat%sub(i)%mat=transpose(inmat%sub(i)%mat)
        end do
end subroutine block_transfer_trans


!Transfer the block by allocating new space
subroutine block_transfer(inmat,outmat)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: inmat
	type(Total_Block),intent(inout) :: outmat

	integer :: i,row_dim,col_dim

	outmat%len=inmat%len
	outmat%num=inmat%num
	outmat%dim=inmat%dim
	outmat%up_dif=inmat%up_dif
	outmat%down_dif=inmat%down_dif

	if(allocated(outmat%sub))deallocate(outmat%sub)
        allocate(outmat%sub(outmat%num))
	do i=1,outmat%num
		outmat%sub(i)%num_up=inmat%sub(i)%num_up
		outmat%sub(i)%num_down=inmat%sub(i)%num_down
		outmat%sub(i)%down_dif=inmat%sub(i)%down_dif
		outmat%sub(i)%row_dim=inmat%sub(i)%row_dim
		outmat%sub(i)%sdim=inmat%sub(i)%sdim

	if(allocated(outmat%sub(i)%mat))deallocate(outmat%sub(i)%mat)
        allocate(outmat%sub(i)%mat(outmat%sub(i)%row_dim,outmat%sub(i)%sdim))
		outmat%sub(i)%mat=inmat%sub(i)%mat
	end do

end subroutine block_transfer
!Transfer general information from old_block to new_block
subroutine block_pass_info(old_block,new_block)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: old_block
	type(Total_Block),intent(inout) :: new_block

	integer :: i,row_dim,col_dim

	new_block%len=old_block%len
	new_block%num=old_block%num
	new_block%dim=old_block%dim
	new_block%up_dif=old_block%up_dif
	new_block%down_dif=old_block%down_dif
	allocate(new_block%sub(new_block%num))

	do i=1,new_block%num
		new_block%sub(i)%num_up=old_block%sub(i)%num_up
		new_block%sub(i)%num_down=old_block%sub(i)%num_down
		new_block%sub(i)%down_dif=old_block%sub(i)%down_dif
		new_block%sub(i)%row_dim=old_block%sub(i)%row_dim
		new_block%sub(i)%sdim=old_block%sub(i)%sdim

		allocate(new_block%sub(i)%mat(new_block%sub(i)%row_dim,new_block%sub(i)%sdim))
		new_block%sub(i)%mat=0.0d0
	end do

end subroutine block_pass_info


!Set block data to zero
subroutine block_zero(block)
	use pubdata
	implicit none

	type(Total_Block),intent(inout) :: block

	integer :: i

	do i=1,block%num
		block%sub(i)%mat=0.0d0
	end do

end subroutine block_zero

subroutine block_add_block_two(block1,coef1,new_block)
        use pubdata
        implicit none

        real(8),intent(in) :: coef1
        type(Total_block),intent(inout) :: block1 
        type(Total_block),intent(inout) :: new_block

        logical :: flag1,flag2
        integer :: i,x,num_up,num_down,idx1,idx2, x1,down_dif

                x1=1
        do i=1,new_block%num
                num_up=new_block%sub(i)%num_up
                num_down=new_block%sub(i)%num_down
                down_dif=new_block%sub(i)%down_dif
                flag1=.false.
                do x=x1,block1%num
                        if(block1%sub(x)%num_up==num_up) then
                        if(block1%sub(x)%num_down==num_down) then
                        if(block1%sub(x)%down_dif==down_dif) then
                                flag1=.true.
                                idx1=x
                                x1=x+1
                                goto 101
                        endif
                        endif
                        endif
                end do
                101 continue

                if(flag1) then
                                new_block%sub(i)%mat=new_block%sub(i)%mat+coef1*block1%sub(idx1)%mat
                endif
        end do

end subroutine block_add_block_two

subroutine block_add_block(block1,coef1,block2,coef2,new_block,coef)
	use pubdata
	implicit none

	real(8),intent(in) :: coef1,coef2,coef
	type(Total_block),intent(inout) :: block1,block2
	type(Total_block),intent(inout) :: new_block

	logical :: flag1,flag2
	integer :: i,x,num_up,num_down,down_dif,idx1,idx2

	do i=1,new_block%num

		num_up=new_block%sub(i)%num_up
		num_down=new_block%sub(i)%num_down
		down_dif=new_block%sub(i)%down_dif

		flag1=.false.
		do x=1,block1%num
			if(block1%sub(x)%num_up==num_up) then
			if(block1%sub(x)%num_down==num_down) then
			if(block1%sub(x)%down_dif==down_dif) then
				flag1=.true.
				idx1=x
				goto 101
			endif
			endif
			endif
		end do
		101 continue


		flag2=.false.
		do x=1,block2%num
			if(block2%sub(x)%num_up==num_up) then
			if(block2%sub(x)%num_down==num_down) then
			if(block2%sub(x)%down_dif==down_dif) then
				flag2=.true.
				idx2=x
				goto 102
			endif
			endif
			endif
		end do
		102 continue

		if(flag1) then
			if(flag2) then
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat+coef1*block1%sub(idx1)%mat+coef2*block2%sub(idx2)%mat
			else
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat+coef1*block1%sub(idx1)%mat
			endif
		else
			if(flag2) then
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat+coef2*block2%sub(idx2)%mat
			else
				new_block%sub(i)%mat=coef*new_block%sub(i)%mat
			endif
		endif
	end do

end subroutine block_add_block


!Transpose the block
subroutine block_transpose(block,new_block)
	use pubdata
	implicit none

	type(Total_block),intent(in) :: block
	type(Total_Block),intent(inout) :: new_block

	integer :: i

	do i=1,block%num
		new_block%sub(i)%mat=transpose(block%sub(i)%mat)
	end do

end subroutine block_transpose


!new_block = coef*old_block
subroutine block_coef(old_block,coef,new_block)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	type(Total_Block),intent(in) :: old_block
	type(Total_Block),intent(inout) :: new_block

	integer :: i

	do i=1,new_block%num
		new_block%sub(i)%mat=coef*old_block%sub(i)%mat
	end do

end subroutine block_coef

!====================================================================
!Flag='T' means transopose; Flag='N' means no transpose
subroutine block_mul_block_dia(oper1,Flag1,oper2,Flag2,coef,new_oper)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: oper1,oper2
	type(Total_Block),intent(inout) :: new_oper
	integer :: i,sdim,num_up,num_down

	new_oper%len=oper1%len
	new_oper%num=oper1%num
	new_oper%dim=oper1%dim
	new_oper%up_dif=0
	new_oper%down_dif=0

	allocate(new_oper%sub(new_oper%num))
	do i=1,new_oper%num
		new_oper%sub(i)%num_up=oper1%sub(i)%num_up
		new_oper%sub(i)%num_down=oper1%sub(i)%num_down
		new_oper%sub(i)%down_dif=oper1%sub(i)%down_dif
		new_oper%sub(i)%row_dim=oper1%sub(i)%row_dim
		new_oper%sub(i)%sdim=oper1%sub(i)%sdim

		sdim=new_oper%sub(i)%sdim
		allocate(new_oper%sub(i)%mat(sdim,sdim))
		call DGEMM(Flag1,Flag2,sdim,sdim,sdim,coef,oper1%sub(i)%mat,sdim&
				 &,oper2%sub(i)%mat,sdim,0.0d0,new_oper%sub(i)%mat,sdim)
	end do

end subroutine block_mul_block_dia
!===========================================================================
subroutine block_mul_block_ndia(oper1,Flag1,oper2,Flag2,coef,new_oper,basis)
	use pubdata
	implicit none

	real(8),intent(in) :: coef
	character(len=1),intent(in) :: Flag1,Flag2
	type(Total_Block),intent(in) :: oper1,oper2
	type(Total_Block),intent(inout) :: new_oper
	type(Total_Basis),intent(in) :: basis

	logical :: Flags,bl_flag1,bl_flag2,bs_flag
	integer :: up_dif1,down_dif1,up_dif2,down_dif2,up_dif,down_dif
	integer :: i,j,x,y,new_dim,mid_dim,old_dim,idx1,idx2,bs_idx,new_id,LDA,LDB,LDC
	integer :: num_up,num_down,mid_num_up,mid_num_down,new_num_up,new_num_down
	integer :: mid_num_up1,mid_num_down1, down_dif11, down_dif22
        double precision :: coef1, coef11, coef12
        integer j1, j2,j3,j4,j5,j6,j11,j22,j33,j44,j55,j66, kq
        integer ji1, ji2,ji3,ji4,ji5,ji6
         real(8), external :: w3js


	!<1>: Get general information
	!<1-1>: For oper1
	up_dif1=0
	down_dif1=0
	if(Flag1=='N') then
		up_dif1=oper1%up_dif
		down_dif1=oper1%down_dif
	else if(Flag1=='T') then
		up_dif1=-oper1%up_dif
		down_dif1=oper1%down_dif
	else
		write(*,*) "Flag1 is wrong!"
		return
	endif

	!<1-2>: For oper2
	up_dif2=0
	down_dif2=0
	if(Flag2=='N') then
		up_dif2=oper2%up_dif
		down_dif2=oper2%down_dif
	else if(Flag2=='T') then
		up_dif2=-oper2%up_dif
		down_dif2=oper2%down_dif
	else
		write(*,*) "Fiag2 is wrong!"
		return
	endif

	!<1-3>: For new_oper (up_dif and down_dif)
	new_oper%len=basis%len
	new_oper%up_dif=up_dif1+up_dif2
	!!!new_oper%down_dif=down_dif1+down_dif2   we define it before calling

	up_dif=new_oper%up_dif
	down_dif=new_oper%down_dif

	!<2>: Get general information of new_oper
	new_oper%num=0
	do i=1,basis%num
		num_up=basis%sub(i)%new_num_up
		num_down=basis%sub(i)%new_num_down

		do new_num_down=abs(num_down-down_dif), num_down+down_dif, su
		new_num_up=num_up+up_dif

		do x=1,basis%num
			if(basis%sub(x)%new_num_up==new_num_up) then
			if(basis%sub(x)%new_num_down==new_num_down) then
				new_oper%num=new_oper%num+1
				goto 101
			endif
			endif
		end do
		101 continue
	end do
	end do


	!<3>: Get general information of new_oper
	allocate(new_oper%sub(new_oper%num))
	new_id=0
	new_oper%dim=0
	do i=1,basis%num
		num_up=basis%sub(i)%new_num_up
		num_down=basis%sub(i)%new_num_down
		new_num_up=num_up+up_dif

		do new_num_down=abs(num_down-down_dif), num_down+down_dif, su
		bs_flag=.false.
		bs_idx=0
		do x=1,basis%num
			if(basis%sub(x)%new_num_up==new_num_up) then
			if(basis%sub(x)%new_num_down==new_num_down) then
				new_id=new_id+1
				bs_flag=.true.
				bs_idx=x
				goto 102
			endif
			endif
		end do
		102 continue

		if(bs_flag) then
			new_oper%sub(new_id)%num_up=num_up
			new_oper%sub(new_id)%num_down=num_down
			new_oper%sub(new_id)%down_dif=new_num_down-num_down
			new_oper%sub(new_id)%row_dim=basis%sub(bs_idx)%dim
			new_oper%sub(new_id)%sdim=basis%sub(i)%dim

			old_dim=basis%sub(i)%dim
			new_dim=basis%sub(bs_idx)%dim
			allocate(new_oper%sub(new_id)%mat(new_dim,old_dim))
			new_oper%sub(new_id)%mat=0.0d0
			new_oper%dim=new_oper%dim+new_oper%sub(new_id)%sdim



        do mid_num_down= abs(num_down-down_dif2), num_down+down_dif2, su
			if(Flag2=='N') then
				mid_num_up1=num_up
				mid_num_down1=num_down
                                down_dif22= mid_num_down-num_down

			else if(Flag2=='T') then
				mid_num_up1=num_up+up_dif2
				mid_num_down1=mid_num_down
                                down_dif22= num_down-mid_num_down
			endif

			idx2=0
			bl_flag2=.false.
			do x=1,oper2%num
				if(oper2%sub(x)%num_up==mid_num_up1) then
				if(oper2%sub(x)%num_down==mid_num_down1) then
				if(oper2%sub(x)%down_dif==down_dif22) then
					bl_flag2=.true.
					idx2=x
					goto 103
				endif
                                endif
				endif
			end do
			103 continue

			!<b>: ============ For oper1 ===================
			if(Flag1=='N') then
				mid_num_up1=num_up+up_dif2
				mid_num_down1=mid_num_down
                                down_dif11= new_num_down-mid_num_down
			else if(Flag1=='T') then
				mid_num_up1=new_num_up
				mid_num_down1=new_num_down
                                down_dif11= mid_num_down-new_num_down
			endif

			idx1=0
			bl_flag1=.false.
			do x=1,oper1%num
				if(oper1%sub(x)%num_up==mid_num_up1) then
				if(oper1%sub(x)%num_down==mid_num_down1) then
				if(oper1%sub(x)%down_dif==down_dif11) then
					bl_flag1=.true.
					idx1=x
					goto 104
				endif
				endif
				endif
			end do
			104 continue

			if(bl_flag1.and.bl_flag2) then
				!<a>: For oper1
				if(Flag1=='N') then
					mid_dim=oper1%sub(idx1)%sdim
					LDA=new_dim
				else if(Flag1=='T') then
					mid_dim=oper1%sub(idx1)%row_dim
					LDA=mid_dim
				endif

				!<b>: For oper2
				if(Flag2=='N') then
					LDB=mid_dim
				else if(Flag2=='T') then
					LDB=old_dim
				endif

        j1=new_num_down    
        j4=-j1
        j3=num_down
        j6=j3
        j2=new_oper%down_dif
        j5=-j4-j6  
        kq=j5    
        coef11=w3js(j1,j2,j3,j4, j5, j6)
        if(coef11.eq.0.0)write(*,*)'w3js wrong', j1, j2, j3,j4,j5,j6
        if(coef11.eq.0.0)stop

        coef1=0.d0
        do j6=-mid_num_down, mid_num_down, su  
        j1=new_num_down    
        j4=-j1

        j2=oper1%down_dif
        j3=mid_num_down
        j5=-j4-j6
        coef12=w3js(j1,j2,j3,j4,j5,j6)

        j11=mid_num_down   
        j44=-j6
        j22=oper2%down_dif
        j33=num_down
        j66=j33
        j55=-j44-j66
        coef12=coef12*w3js(j11,j22,j33,j44,j55,j66)*(-1)**((j11+j44)/2)

        ji1=oper1%down_dif     
        ji2=oper2%down_dif    
        ji3=new_oper%down_dif   
        ji4=j5  
        ji5=j55
        ji6=-kq
        coef12=coef12*w3js(ji1,ji2,ji3,ji4, ji5, ji6)*(-1)**((-ji1+ji2+kq)/2)*dsqrt(ji3+1.0d0)
        coef1=coef1+coef12
        enddo

        if(coef1.ne.0.0)then
        coef1=coef*coef1/coef11

		call DGEMM(Flag1,Flag2,new_dim,old_dim,mid_dim,coef1,oper1%sub(idx1)%mat,LDA,&
			  &oper2%sub(idx2)%mat,LDB,1.0d0,new_oper%sub(new_id)%mat,new_dim)
			endif
		endif
	end do
        endif

	end do
	end do

end subroutine block_mul_block_ndia


!Deallocate block
subroutine deallocate_block(block)
	use pubdata
	implicit none

	type(Total_block),intent(inout) :: block

	integer :: i

	do i=1,block%num
		deallocate(block%sub(i)%mat)
	end do
	deallocate(block%sub)
end subroutine deallocate_block

!Transfer model data by allocating new space
subroutine model_transfer(inmod,outmod)
	use pubdata
	implicit none

	type(Total_Model),intent(in) :: inmod
	type(Total_Model),intent(inout) :: outmod

	integer :: ii

	!Transfer index and Hamiltonian
	outmod%len=inmod%len
	call block_transfer(inmod%ham,outmod%ham)

	!For operators in the middle of system
	do ii=1,nleg11(inmod%len) !!!Ny
		call block_transfer(inmod%sub(ii)%elec_up,outmod%sub(ii)%elec_up)
		call block_transfer(inmod%sub(ii)%elec_down,outmod%sub(ii)%elec_down)
		call block_transfer(inmod%sub(ii)%spin_sd,outmod%sub(ii)%spin_sd)
		call block_transfer(inmod%sub(ii)%num_sn,outmod%sub(ii)%num_sn)
	end do


end subroutine model_transfer


!deallocate model data
subroutine deallocate_model(model)
	use pubdata
	implicit none

	type(Total_Model),intent(inout) :: model
	integer :: ii

	!Deallocate Hamiltonian
        if(allocated(model%ham%sub))then
        call deallocate_block(model%ham)
                endif


	!For operators in the middle of system
	do ii=1,nleg11(model%len)
                if(allocated(model%sub(ii)%elec_up%sub))then
		call deallocate_block(model%sub(ii)%elec_up)
		call deallocate_block(model%sub(ii)%elec_down)
		call deallocate_block(model%sub(ii)%spin_sd)
		call deallocate_block(model%sub(ii)%num_sn)
                        endif
	end do

end subroutine deallocate_model


!=========================================================================
!Allocate new wavefunction according to the superbasis information
!=========================================================================
subroutine allocate_wave(super,wave)
	use pubdata
	implicit none

	type(Super_Basis),intent(in) :: super
	type(Wavefunction),intent(inout) :: wave

	integer :: i,sys_dim,env_dim

	wave%sys_len=super%sys_len
	wave%env_len=super%env_len
	wave%num=super%num
	wave%dim=super%dim

	allocate(wave%sub(wave%num))
	do i=1,wave%num
		wave%sub(i)%sys_num_up=super%sub(i)%sys_num_up
		wave%sub(i)%sys_num_down=super%sub(i)%sys_num_down
		wave%sub(i)%env_num_up=super%sub(i)%env_num_up
		wave%sub(i)%env_num_down=super%sub(i)%env_num_down
		wave%sub(i)%sys_dim=super%sub(i)%sys_dim
		wave%sub(i)%env_dim=super%sub(i)%env_dim

		allocate(wave%sub(i)%vec(wave%sub(i)%sys_dim,wave%sub(i)%env_dim))
		wave%sub(i)%vec=0.0d0
	end do

end subroutine allocate_wave


!Transfer wavefunction by allocating new space
subroutine wave_transfer(inwave,outwave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: inwave
	type(Wavefunction),intent(inout) :: outwave

	integer :: i,sys_dim,env_dim

	outwave%sys_len=inwave%sys_len
	outwave%env_len=inwave%env_len
	outwave%num=inwave%num
	outwave%dim=inwave%dim

	allocate(outwave%sub(outwave%num))
	do i=1,outwave%num
		outwave%sub(i)%sys_num_up=inwave%sub(i)%sys_num_up
		outwave%sub(i)%sys_num_down=inwave%sub(i)%sys_num_down
		outwave%sub(i)%env_num_up=inwave%sub(i)%env_num_up
		outwave%sub(i)%env_num_down=inwave%sub(i)%env_num_down
		outwave%sub(i)%sys_dim=inwave%sub(i)%sys_dim
		outwave%sub(i)%env_dim=inwave%sub(i)%env_dim

		allocate(outwave%sub(i)%vec(outwave%sub(i)%sys_dim,outwave%sub(i)%env_dim))
		outwave%sub(i)%vec=inwave%sub(i)%vec
	end do

end subroutine wave_transfer


!Notice that inwave and outwave have the same configuration
subroutine wave_pass(inwave,outwave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: inwave
	type(Wavefunction),intent(inout) :: outwave

	integer :: i

	outwave%sys_len=inwave%sys_len
	outwave%env_len=inwave%env_len
	outwave%num=inwave%num
	outwave%dim=inwave%dim

	do i=1,outwave%num
		outwave%sub(i)%sys_num_up=inwave%sub(i)%sys_num_up
		outwave%sub(i)%sys_num_down=inwave%sub(i)%sys_num_down
		outwave%sub(i)%env_num_up=inwave%sub(i)%env_num_up
		outwave%sub(i)%env_num_down=inwave%sub(i)%env_num_down
		outwave%sub(i)%sys_dim=inwave%sub(i)%sys_dim
		outwave%sub(i)%env_dim=inwave%sub(i)%env_dim
		outwave%sub(i)%vec=inwave%sub(i)%vec
	end do

end subroutine wave_pass

!Deallocate wavefunction
subroutine deallocate_wave(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(inout) :: wave

	integer :: i
	
	do i=1,wave%num
		deallocate(wave%sub(i)%vec)
	end do
	deallocate(wave%sub)

end subroutine deallocate_wave


!Note: left and right have the same configuration
real(8) function wave_product(left,right)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: left,right

	integer :: i,x

	wave_product=0.0d0
	do i=1,left%num
		do x=1,left%sub(i)%env_dim
			wave_product=wave_product+dot_product(left%sub(i)%vec(:,x),right%sub(i)%vec(:,x))
		end do
	end do
	return

end function wave_product


!Get the norm of the wavefunction
real(8) function wave_norm(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: wave

	integer :: i,x

	wave_norm=0.0d0
	do i=1,wave%num
		do x=1,wave%sub(i)%env_dim
			wave_norm=wave_norm+dot_product(wave%sub(i)%vec(:,x),wave%sub(i)%vec(:,x))
		end do
	end do
	wave_norm=dsqrt(wave_norm)
	return

end function wave_norm


!Nomalize the wavefunction
subroutine wave_normalize(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(inout) :: wave

	integer :: i,j
	real(8) :: vec_norm
	real(8),external :: wave_norm

	vec_norm=wave_norm(wave)
	do i=1,wave%num
		wave%sub(i)%vec=wave%sub(i)%vec/vec_norm
	end do

end subroutine wave_normalize


!Initialize the wavefunction
subroutine wave_initialize(wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(inout) :: wave

	integer :: i, j

	do i=1,wave%num
		call random_number(wave%sub(i)%vec)
	end do
	call wave_normalize(wave)

end subroutine wave_initialize



subroutine wave_initialize_zero(wave)
        use pubdata
        implicit none

       double precision,allocatable :: mid_r(:,:),mid_r1(:,:)
        type(Wavefunction),intent(inout) :: wave
        integer :: i, sys_dim, env_dim
      do i=1,wave%num
     wave%sub(i)%vec=0.0d0  !!!!mid_r) !,(thetay-thetax)*mid_r)
        end do
end subroutine wave_initialize_zero

double precision function wave_overlap(left, right)
        use pubdata
        implicit none

        type(Wavefunction),intent(in) :: left, right
        integer :: i,j, x
        wave_overlap=0.0d0
        do i=1,left%num
                do j=1, right%num
                   if(left%sub(i)%sys_num_up==right%sub(j)%sys_num_up)then
                   if(left%sub(i)%sys_num_down==right%sub(j)%sys_num_down)then
                do x=1,   left%sub(i)%env_dim
                if(left%sub(i)%env_dim.ne.right%sub(j)%env_dim)write(*,*)'overlap wrong'
        wave_overlap=wave_overlap+dot_product(left%sub(i)%vec(:,x),right%sub(i)%vec(:,x))
                end do
                endif
                endif
        end do
        end do
        return
end function wave_overlap

!Get the Gram¨CSchmidt projection operator
!proj_u(v) = <u,v>*u/<u,u> = wave_new
subroutine Gram_Schmidt_WF(wave_u,wave_v,wave_new)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: wave_u,wave_v
	type(Wavefunction) :: wave_new

	integer :: i,j,k
	real(8) :: prod_uv,prod_uu
	real(8),external :: wave_product

	!<1>: Get <u,v> and <u,u>
	prod_uv=wave_product(wave_u,wave_v)
	prod_uu=wave_product(wave_u,wave_u)

	!<2>: Get proj_u(v)=<u,v>*u/<u,u> = wave_new
	wave_new%sys_len=wave_u%sys_len
	wave_new%env_len=wave_u%env_len
	do i=1,wave_new%num
		wave_new%sub(i)%vec = wave_u%sub(i)%vec*(prod_uv/prod_uu)
	end do

end subroutine Gram_Schmidt_WF


subroutine Gram_Schmidt_process(wave_prev,levels,wave_new)
	use pubdata
	implicit none

	integer,intent(in) :: levels
	type(Wavefunction) :: wave_prev(levels)
	type(Wavefunction) :: wave_new

	integer :: i,j,x,y
	real(8) :: prod_uv(levels),prod_uu(levels)
	real(8),external :: wave_product
	type(Wavefunction) :: tmp_vec

	do j=1,levels
		prod_uv(j) = wave_product(wave_prev(j),wave_new)
		prod_uu(j) = wave_product(wave_prev(j),wave_prev(j))
	end do

	call wave_transfer(wave_new,tmp_vec)
	do i=1,tmp_vec%num
		tmp_vec%sub(i)%vec=0.0d0
	end do

	do j=1,levels
	do i=1,tmp_vec%num
		tmp_vec%sub(i)%vec = tmp_vec%sub(i)%vec + (prod_uv(j)/prod_uu(j))*wave_prev(j)%sub(i)%vec
	end do
	end do

	!<3>: Normalize wave_new
	do i=1,wave_new%num
		wave_new%sub(i)%vec = wave_new%sub(i)%vec - tmp_vec%sub(i)%vec
	end do
	call wave_normalize(wave_new)
	call deallocate_wave(tmp_vec)

end subroutine Gram_Schmidt_process


subroutine Gram_Schmidt_Lanczos(wave_prev,levels,wave_lanc,lanc_num,wave_new)
	use pubdata
	implicit none

	integer,intent(in) :: levels,lanc_num
	type(Wavefunction) :: wave_prev(levels),wave_lanc(lanc_num)
	type(Wavefunction) :: wave_new

	integer :: i,j,x,y
	real(8) :: prod_uv(levels+lanc_num),prod_uu(levels+lanc_num)
	real(8),external :: wave_product
	type(Wavefunction) :: prev_vec,lanc_vec

	do j=1,levels
		prod_uv(j) = wave_product(wave_prev(j),wave_new)
		prod_uu(j) = wave_product(wave_prev(j),wave_prev(j))
	end do

	do j=1,lanc_num
		prod_uv(levels+j) = wave_product(wave_lanc(j),wave_new)
		prod_uu(levels+j) = wave_product(wave_lanc(j),wave_lanc(j))
	end do

	call wave_transfer(wave_new,prev_vec)
	do i=1,prev_vec%num
		prev_vec%sub(i)%vec=0.0d0
	end do

	do j=1,levels
	do i=1,prev_vec%num
		prev_vec%sub(i)%vec = prev_vec%sub(i)%vec + wave_prev(j)%sub(i)%vec*(prod_uv(j)/prod_uu(j))
	end do
	end do

	call wave_transfer(wave_new,lanc_vec)
	do i=1,lanc_vec%num
		lanc_vec%sub(i)%vec=0.0d0
	end do

	do j=1,lanc_num
	do i=1,lanc_vec%num
		lanc_vec%sub(i)%vec = lanc_vec%sub(i)%vec + wave_lanc(j)%sub(i)%vec*(prod_uv(levels+j)/prod_uu(levels+j))
	end do
	end do

	!<3>: Normalize wave_new
	do i=1,wave_new%num
		wave_new%sub(i)%vec = wave_new%sub(i)%vec - prev_vec%sub(i)%vec - lanc_vec%sub(i)%vec
	end do
	call wave_normalize(wave_new)
	call deallocate_wave(prev_vec)
	call deallocate_wave(lanc_vec)

end subroutine Gram_Schmidt_Lanczos


!==================================================================
!Transfer basis information from the old one to the new one
!==================================================================
subroutine basis_transfer(old_basis,new_basis)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: old_basis
	type(Total_Basis),intent(inout) :: new_basis

	integer :: i,x

         if(allocated(new_basis%sub))call deallocate_basis(new_basis)

	new_basis%len=old_basis%len
	new_basis%num=old_basis%num
	new_basis%dim=old_basis%dim

	allocate(new_basis%sub(new_basis%num))
	do i=1,new_basis%num
		new_basis%sub(i)%new_num_up=old_basis%sub(i)%new_num_up
		new_basis%sub(i)%new_num_down=old_basis%sub(i)%new_num_down
		new_basis%sub(i)%idx=old_basis%sub(i)%idx
		new_basis%sub(i)%num=old_basis%sub(i)%num
		new_basis%sub(i)%dim=old_basis%sub(i)%dim
		allocate(new_basis%sub(i)%sub(new_basis%sub(i)%num))

		do x=1,new_basis%sub(i)%num
			new_basis%sub(i)%sub(x)%spos=old_basis%sub(i)%sub(x)%spos
			new_basis%sub(i)%sub(x)%sdim=old_basis%sub(i)%sub(x)%sdim
			new_basis%sub(i)%sub(x)%bl_num_up=old_basis%sub(i)%sub(x)%bl_num_up
			new_basis%sub(i)%sub(x)%bl_num_down=old_basis%sub(i)%sub(x)%bl_num_down
			new_basis%sub(i)%sub(x)%st_num_up=old_basis%sub(i)%sub(x)%st_num_up
			new_basis%sub(i)%sub(x)%st_num_down=old_basis%sub(i)%sub(x)%st_num_down
		end do
	end do

end subroutine basis_transfer


!Output the basis information
subroutine output_basis(basis)
	use pubdata
	implicit none

	type(Total_basis),intent(in) :: basis

	integer :: i,x

	write(*,"(A35)") "******** Basis Information ********"
	write(*,111) "Len=",basis%len,"Num=",basis%num,"Dim=",basis%dim
	do i=1,basis%num
		write(*,112) "Num=",basis%sub(i)%num,"Dim=",basis%sub(i)%dim,&
					&"new_num_up=",basis%sub(i)%new_num_up,"new_num_down=",basis%sub(i)%new_num_down

	end do

111 format(A4,I3,2X,A4,I3,4x,A4,I6)
112 format(A4,I4,2X,A4,I4,2X,A11,I4,2X,A13,I4)
113 format(A10,I4,2X,A12,I4)
114 format(A4,I4,2X,A4,I4)

end subroutine output_basis


!Deallocate the basis information
subroutine deallocate_basis(basis)
	use pubdata
	implicit none

	type(Total_basis),intent(inout) :: basis

	integer :: i

	do i=1,basis%num
		deallocate(basis%sub(i)%sub)
	end do
	deallocate(basis%sub)

end subroutine deallocate_basis


!Deallocate the superbasis information
subroutine deallocate_super_basis(super)
	use pubdata
	implicit none

	type(Super_basis),intent(inout) :: super

	deallocate(super%sub)

end subroutine deallocate_super_basis


!==================================
!For filename operation
!==================================
!Translate number into string
subroutine num_to_str(num,num_str)
	implicit none

	integer,intent(in) :: num
	character(len=6) :: num_str

	integer :: i,idx,ra,res
	character(len=6) :: str1

	!<1>: Initiate
	num_str(1:6)=" "

	!<2>: Get num_str
	idx=1
	res=num
	do while(res>0)
		if(abs(res)<10) then
			str1(idx:idx)=char(res+48)
			goto 111
		endif
		ra=res-res/10*10
		str1(idx:idx)=char(ra+48)
		res=res/10
		idx=idx+1
	end do
	111 continue

	do i=1,idx
		num_str(i:i)=str1(idx-i+1:idx-i+1)
	end do

end subroutine num_to_str


!Link string with string
subroutine str_link_str(str_name,str_len,num)
	implicit none

	integer,intent(in) :: str_len,num
	character(len=20) :: str_name

	integer :: num_len
	character(len=6) :: num_str

	call num_to_str(num,num_str)

	num_len=len_trim(num_str)
	str_name(str_len+1:str_len+num_len)=num_str(1:num_len)

end subroutine str_link_str

!Save wavefunction to disk
subroutine wave_to_disk(wave,fileid,filename,name_len)
	use pubdata
	implicit none

	integer,intent(in) :: fileid,name_len
	character(len=name_len) :: filename
	type(Wavefunction),intent(in) :: wave

	integer :: i,x,y
	
	open(fileid,file=filename,form='unformatted')

	write(fileid) wave%sys_len
	write(fileid) wave%env_len
	write(fileid) wave%num
	write(fileid) wave%dim

	do i=1,wave%num
		write(fileid)
		write(fileid) wave%sub(i)%sys_num_up,wave%sub(i)%sys_num_down
		write(fileid) wave%sub(i)%env_num_up,wave%sub(i)%env_num_down
		write(fileid) wave%sub(i)%sys_dim,wave%sub(i)%env_dim
				write(fileid) wave%sub(i)%vec  
	end do

	close(fileid)

end subroutine wave_to_disk


!Read wavefunction from disk
subroutine wave_from_disk(wave,fileid,filename,name_len)
	use pubdata
	implicit none

	integer,intent(in) :: fileid,name_len
	character(len=name_len) :: filename
	type(Wavefunction),intent(inout) :: wave

	integer :: i,x,y

	open(fileid,file=filename,form='unformatted')

	read(fileid) wave%sys_len
	read(fileid) wave%env_len
	read(fileid) wave%num
	read(fileid) wave%dim

	allocate(wave%sub(wave%num))
	do i=1,wave%num
		read(fileid)
		read(fileid) wave%sub(i)%sys_num_up,wave%sub(i)%sys_num_down
		read(fileid) wave%sub(i)%env_num_up,wave%sub(i)%env_num_down
		read(fileid) wave%sub(i)%sys_dim,wave%sub(i)%env_dim

		allocate(wave%sub(i)%vec(wave%sub(i)%sys_dim,wave%sub(i)%env_dim))
				read(fileid) wave%sub(i)%vec  !!!!(x,y)
	end do

	close(fileid)

end subroutine wave_from_disk


!Read wavefunction from disk
!Save block to disk
subroutine block_to_disk(block,fileid)
	use pubdata
	implicit none

	integer,intent(in) :: fileid
	type(Total_Block),intent(in) :: block

	integer :: i,x,y,row_dim,col_dim

	write(fileid) block%len
	write(fileid) block%num
	write(fileid) block%dim
	write(fileid) block%up_dif
	write(fileid) block%down_dif

	do i=1,block%num
		write(fileid)
		write(fileid) block%sub(i)%num_up,block%sub(i)%num_down, block%sub(i)%down_dif
		write(fileid) block%sub(i)%row_dim,block%sub(i)%sdim
				write(fileid) block%sub(i)%mat  
	end do
end subroutine block_to_disk

!Read block from disk
subroutine block_from_disk(block,fileid)
	use pubdata
	implicit none
	integer,intent(in) :: fileid
	type(Total_Block),intent(inout) :: block
	integer :: i,x,y
	read(fileid) block%len
	read(fileid) block%num
	read(fileid) block%dim
	read(fileid) block%up_dif
	read(fileid) block%down_dif

	allocate(block%sub(block%num))
	do i=1,block%num
		read(fileid)
		read(fileid) block%sub(i)%num_up,block%sub(i)%num_down, block%sub(i)%down_dif
		read(fileid) block%sub(i)%row_dim,block%sub(i)%sdim
		allocate(block%sub(i)%mat(block%sub(i)%row_dim,block%sub(i)%sdim))
				read(fileid) block%sub(i)%mat  !!(x,y)
	end do
end subroutine block_from_disk

!Note: flag=.true. for system, flag=.false. for environment
subroutine model_to_disk(model,fileid,mod_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,mod_idx
	type(Total_Model),intent(in) :: model
	
	integer :: ii

	!<1>: Save index
	if(flag) then
		open(fileid,file=sysname(mod_idx,1),form='unformatted')
	else
		open(fileid,file=envname(mod_idx,1),form='unformatted')
	endif

	write(fileid) model%len

	!<2>: Save general data
	call block_to_disk(model%ham,fileid)

	!For operators in the middle of the system
	do ii=1,nleg11(model%len)
		call block_to_disk(model%sub(ii)%elec_up,fileid)
		call block_to_disk(model%sub(ii)%elec_down,fileid)
		call block_to_disk(model%sub(ii)%spin_sd,fileid)
		call block_to_disk(model%sub(ii)%num_sn,fileid)
	end do


	close(fileid)

end subroutine model_to_disk

!Read model from disk
!Note: flag=.true. for system, flag=.false. for environment
subroutine model_from_disk(model,fileid,mod_idx,Flag)
	use pubdata
	implicit none

	logical,intent(in) :: Flag
	integer,intent(in) :: fileid,mod_idx
	type(Total_Model),intent(inout) :: model
	
	integer :: ii

	!<1>: Read index
	if(flag) then
		open(fileid,file=sysname(mod_idx,1),form='unformatted')
	else
		open(fileid,file=envname(mod_idx,1),form='unformatted')
	endif

	read(fileid) model%len

	!<2>: Save general data
	call block_from_disk(model%ham,fileid)

	!For operators in the middle of system
	do ii=1, nleg11(model%len)
		call block_from_disk(model%sub(ii)%elec_up,fileid)
		call block_from_disk(model%sub(ii)%elec_down,fileid)
		call block_from_disk(model%sub(ii)%spin_sd,fileid)
		call block_from_disk(model%sub(ii)%num_sn,fileid)
	end do


	!==================== For space saving ======================

        if((infi_del.ne.1).or.(mod_idx.gt.num_site0/2+4))then
	if(Flag) then
		open(fileid,file=sysname(mod_idx+1,1),form='unformatted')
		close(fileid,status='delete')
	else
		open(fileid,file=envname(mod_idx+1,1),form='unformatted')
		close(fileid,status='delete')
	endif
        endif
end subroutine model_from_disk

!Save truncation matrix to disk:
!Note: flag=.true. for system, flag=.false. for environment
subroutine truns_to_disk(truns,fileid,tr_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,tr_idx
	type(Total_Block),intent(in) :: truns

	!Open file
	if(flag) then
		open(fileid,file=sysname(tr_idx,2),form='unformatted')
	else
		open(fileid,file=envname(tr_idx,2),form='unformatted')
	endif
	call block_to_disk(truns,fileid)
	close(fileid)
end subroutine truns_to_disk

!Read truncation matrix to disk:
!Note: flag=.true. for system
!		 flag=.false. for environment
subroutine truns_from_disk(truns,fileid,tr_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,tr_idx
	type(Total_Block),intent(inout) :: truns

	!Open file
	if(flag) then
		open(fileid,file=sysname(tr_idx,2),form='unformatted')
	else
		open(fileid,file=envname(tr_idx,2),form='unformatted')
	endif

	call block_from_disk(truns,fileid)

	close(fileid)

end subroutine truns_from_disk


!Save basis to disk
subroutine basis_to_disk(basis,fileid,bs_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,bs_idx
	type(Total_Basis),intent(in) :: basis

	integer :: i,x,y

	!Open file
	if(flag) then
		open(fileid,file=sysname(bs_idx,3))
	else
		open(fileid,file=envname(bs_idx,3))
	endif

	write(fileid,101) "len=",basis%len
	write(fileid,101) "num=",basis%num
	write(fileid,101) "dim=",basis%dim

	do i=1,basis%num
		write(fileid,*)
		write(fileid,101) "idx=",basis%sub(i)%idx
		write(fileid,101) "num=",basis%sub(i)%num
		write(fileid,101) "dim=",basis%sub(i)%dim
		write(fileid,102) "new_num_up=",basis%sub(i)%new_num_up
		write(fileid,103) "new_num_down=",basis%sub(i)%new_num_down

		do x=1,basis%sub(i)%num
			write(fileid,101) "pos=",basis%sub(i)%sub(x)%spos
			write(fileid,101) "dim=",basis%sub(i)%sub(x)%sdim
			write(fileid,104) "bl_num_up=",basis%sub(i)%sub(x)%bl_num_up,"bl_num_down=",basis%sub(i)%sub(x)%bl_num_down
			write(fileid,104) "st_num_up=",basis%sub(i)%sub(x)%st_num_up,"st_num_down=",basis%sub(i)%sub(x)%st_num_down
		end do
	end do
	close(fileid)

101 format(A4,I6)
102 format(A11,I4)
103 format(A13,I4)
104 format(A10,I4,2X,A12,I4)

end subroutine basis_to_disk


!Read basis from disk
subroutine basis_from_disk(basis,fileid,bs_idx,flag)
	use pubdata
	implicit none

	logical,intent(in) :: flag
	integer,intent(in) :: fileid,bs_idx
	type(Total_Basis),intent(inout) :: basis

	integer :: i,x,y
	character(len=4) :: str4
	character(len=10) :: str10
	character(len=11) :: str11
	character(len=12) :: str12
	character(len=13) :: str13

	!Open file
	if(flag) then
		open(fileid,file=sysname(bs_idx,3))
	else
		open(fileid,file=envname(bs_idx,3))
	endif

	read(fileid,101) str4,basis%len
	read(fileid,101) str4,basis%num
	read(fileid,101) str4,basis%dim

	allocate(basis%sub(basis%num))
	do i=1,basis%num
		read(fileid,*)
		read(fileid,101) str4,basis%sub(i)%idx
		read(fileid,101) str4,basis%sub(i)%num
		read(fileid,101) str4,basis%sub(i)%dim
		read(fileid,102) str11,basis%sub(i)%new_num_up
		read(fileid,103) str13,basis%sub(i)%new_num_down

		allocate(basis%sub(i)%sub(basis%sub(i)%num))
		do x=1,basis%sub(i)%num
			read(fileid,101) str4,basis%sub(i)%sub(x)%spos
			read(fileid,101) str4,basis%sub(i)%sub(x)%sdim
			read(fileid,104) str10,basis%sub(i)%sub(x)%bl_num_up,str12,basis%sub(i)%sub(x)%bl_num_down
			read(fileid,104) str10,basis%sub(i)%sub(x)%st_num_up,str12,basis%sub(i)%sub(x)%st_num_down
		end do
	end do
	close(fileid)

101 format(A4,I6)
102 format(A11,I4)
103 format(A13,I4)
104 format(A10,I4,2X,A12,I4)

end subroutine basis_from_disk
