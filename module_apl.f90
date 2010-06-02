module apl
    !use kinds
    use readtraj
    !use trajop
    use util
    implicit none
    integer(kind=ik),allocatable ::&
    apl_atoms(:),apl_atoms_invatom(:),apl_atoms_invmol(:),grid(:,:,:)
    contains

    subroutine apl_atomlist(atoms)
        character(kind=1,len=1) :: atoms(:,:)
        integer(kind=ik) :: i,j,k,imol
        j=0
        do i=1,size(atoms,2)
            j=j+molt(moltypeofuatom(atomindex(trim(stringconv(atoms(:,i))))))%nmol
        end do
        if(allocated(apl_atoms))deallocate(apl_atoms)
        if(allocated(apl_atoms_invatom))deallocate(apl_atoms_invatom,apl_atoms_invmol)
        allocate(apl_atoms(j),apl_atoms_invatom(j),apl_atoms_invmol(j))
        k=0
        do i=1,size(atoms,2)
            do imol=1,molt(moltypeofuatom(atomindex(trim(stringconv(atoms(:,i))))))%nmol
                k=k+1
                apl_atoms(k)=cind(atomindex(trim(stringconv(atoms(:,i)))),imol)
                apl_atoms_invatom(k)=atomindex(trim(stringconv(atoms(:,i))))
                apl_atoms_invmol(k)=imol
            end do
        end do
    end subroutine apl_atomlist
    
    subroutine apl_grid(instr)
        type(instruct) :: instr
        real(kind=rk) :: rmin
        integer(kind=ik),allocatable ::&
        apl_lower(:),apl_lower_inv(:),apl_upper(:),apl_upper_inv(:)!,grid(:,:,:)
        integer(kind=ik) ::i,j,kl,ku,dmpc,mgdg,g
        i=size(apl_atoms)
        allocate(apl_lower(i),apl_lower_inv(i),apl_upper(i),apl_upper_inv(i))
        kl=0;ku=0;dmpc=0;mgdg=0
        do i=1,size(apl_atoms)
            if(dot_product(center_of_molecule(moltypeofuatom(apl_atoms_invatom(i))&
            ,apl_atoms_invmol(i))-centerofmembrane,director)<0)then
                kl=kl+1
                apl_lower(kl)=apl_atoms(i)
                apl_lower_inv(kl)=i
            else
                ku=ku+1
                apl_upper(ku)=apl_atoms(i)
                apl_upper_inv(ku)=i
            end if
        end do
        call reallocate(apl_lower,kl)
        call reallocate(apl_lower_inv,kl)
        call reallocate(apl_upper,ku)
        call reallocate(apl_upper_inv,ku)

       if(.NOT.allocated(grid))allocate(grid(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2),1:2))
       if(size(grid,1)/=instr%set%aplgrid(1) .OR.&
       size(grid,2)/=instr%set%aplgrid(2))then
       deallocate(grid);allocate(grid(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2),1:2))
       end if
        do i=1,size(grid,2)
            do j=1,size(grid,1)
                grid(j,i,1)=nneighbour(j,i,apl_upper,apl_upper_inv)
                grid(j,i,2)=nneighbour(j,i,apl_lower,apl_lower_inv)
            end do
        end do
       contains

     !  function nneighbour(a,b,apl_side,apl_side_inv) result(k)
     !      integer(kind=ik) :: a,b,i,k,apl_side(:),apl_side_inv(:)
     !      real(kind=4) :: gridcoor(2),r2,rmin
     !      gridcoor=[(box(1)/real(size(grid,1),rk))*real(a,rk),(box(2)/real(size(grid,2),rk))*real(b,rk)]
     !      rmin=huge(rmin)
     !      k=-1
     !      do i=1,size(apl_side)
     !         r2=sum(mymodulo(real(coor(1:2,apl_side(i)),4) -&
     !          gridcoor,real(box(1:2),4))**2)
     !          if(r2<rmin)then
     !              rmin=r2
     !              k=apl_side_inv(i)!moltypeofuatom(apl_atoms_invatom(apl_upper_inv(i)))
     !          end if
     !      end do
     !  end function nneighbour

       function nneighbour(a,b,apl_side,apl_side_inv) result(k)
           integer(kind=ik) :: a,b,i,k,apl_side(:),apl_side_inv(:)
           real(kind=rk) :: gridcoor(2),r2,rmin
           gridcoor=[(box(1)/real(size(grid,1),rk))*real(a,rk),(box(2)/real(size(grid,2),rk))*real(b,rk)]
           rmin=huge(rmin)
           k=-1
           do i=1,size(apl_side)
              r2=sum(mymodulo(coor(1:2,apl_side(i)) -&
               gridcoor,box(1:2))**2)
               if(r2<rmin)then
                   rmin=r2
                   k=apl_side_inv(i)
               end if
           end do
       end function nneighbour


!        elemental function mymodulo(a,b) result(c)
!            real(kind=4),intent(in) :: a,b
!            real(kind=4) :: c
!            c=a-nint(a/b)*b !ok
!        end function mymodulo

!        elemental function mymodulo(a,b) result(c)
!            real(kind=rk),intent(in) :: a,b
!            real(kind=rk) :: c
!            c=a-nint(a/b)*b !ok
!        end function mymodulo

    end subroutine apl_grid

    subroutine apl_calc(instr,frame)
        type(instruct) :: instr
        integer(kind=ik),intent(in) :: frame
        integer(kind=ik) :: i,j,k(instr%nmolop),side,imol!,atom,side
        if(.NOT.allocated(instr%apl_side))allocate(instr%apl_side(instr%nmolop))
        k=0;side=1
        do side=1,2
            do j=1,size(grid,2)
                do i=1,size(grid,1)
                    if(ANY(instr%atoms==apl_atoms_invatom(grid(i,j,side))))then
                        imol=apl_atoms_invmol(grid(i,j,side))
                        instr%apl_side(imol)=side
                        k(imol)=k(imol)+1
                    end if
                end do
            end do
        end do
        instr%datam(:,frame)=real(k,rk)*box(1)*box(2)/(real(size(grid,1),rk)*real(size(grid,2),rk))
    end subroutine apl_calc

    subroutine apl_matrix_out(frame)
        integer(kind=ik) :: i,j,side,frame
        integer(kind=ik),allocatable :: v1(:),v2(:)
        allocate(v1(size(grid,1)),v2(size(grid,1)))
        do side=1,2
            if(side==1)open(43,file='APL_polygons_upper_grid'&
            //trim(adjustl(intstr(size(grid,1))))//'.'//trim(adjustl(intstr(size(grid,2))))&
            //'_frame'//trim(adjustl(intstr(frame))))
            if(side==2)open(43,file='APL_polygons_lower_grid'&
            //trim(adjustl(intstr(size(grid,1))))//'.'//trim(adjustl(intstr(size(grid,2))))&
            //'_frame'//trim(adjustl(intstr(frame))))
            if(side==1)open(44,file='APL_polygons_upper_grid'&
            //trim(adjustl(intstr(size(grid,1))))//'.'//trim(adjustl(intstr(size(grid,2))))&
            //'_frame'//trim(adjustl(intstr(frame)))//'_mtrx')
            if(side==2)open(44,file='APL_polygons_lower_grid'&
            //trim(adjustl(intstr(size(grid,1))))//'.'//trim(adjustl(intstr(size(grid,2))))&
            //'_frame'//trim(adjustl(intstr(frame)))//'_mtrx')

            do j=1,size(grid,2)
                do i=1,size(v1)
                    v1(i)=apl_atoms_invmol(grid(i,j,side))
                    v2(i)=moltypeofuatom(apl_atoms_invatom(grid(i,j,side)))
                    write(43,*)i,j,v1(i),v2(i)
                    v1(i)=v1(i)+sum(molt(1:v2(i)-1)%nmol)
                    end do
                    write(43,*)
                    write(44,*)v1(:)
            end do
            close(43);close(44)
        end do
        deallocate(v1,v2)
    end subroutine apl_matrix_out

end module apl
