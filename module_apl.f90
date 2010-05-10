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
            !write(*,*)"MOL",molt(moltypeofuatom(atomindex(trim(stringconv(atoms(:,i))))))%nmol
                !write(*,*)i,trim(stringconv(atoms(:,i))),cind(atomindex(trim(stringconv(atoms(:,i)))),imol)
                k=k+1
                apl_atoms(k)=cind(atomindex(trim(stringconv(atoms(:,i)))),imol)
                apl_atoms_invatom(k)=atomindex(trim(stringconv(atoms(:,i))))
                apl_atoms_invmol(k)=imol
            end do
        end do
    end subroutine apl_atomlist
    
    subroutine apl_grid
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
        !call mindist(apl_upper)
        !allocate(grid(1:ceiling((box(1)/rmin)*(sqrt(12._rk)+2)),1:ceiling((box(2)/rmin)*(sqrt(12._rk)+2))))
        g=50
       ! open(44,file="areatest_dmpc");open(45,file="areatest_mgdg")
       ! do
       ! g=g+50
       ! if(g==1050)exit
       ! dmpc=0;mgdg=0
       if(.NOT.allocated(grid))allocate(grid(1:g,1:g,1:2))
        do i=1,size(grid,2)
            do j=1,size(grid,1)
                grid(j,i,1)=nneighbour(j,i,apl_upper,apl_upper_inv)
                grid(j,i,2)=nneighbour(j,i,apl_lower,apl_lower_inv)
               ! if(grid(j,i)==1)then
               !     dmpc=dmpc+1
               ! else
               !     mgdg=mgdg+1
               ! end if
            end do
        end do
        !write(44,*)g,(box(1)*box(2)*(real(dmpc,rk)/(real(size(grid,1),rk)*real(size(grid,2),rk))))/48._rk !DMPC 72 48
        !write(45,*)g,(box(1)*box(2)*(real(mgdg,rk)/(real(size(grid,1),rk)*real(size(grid,2),rk))))/40._rk !MGDG 18 40

       ! if(allocated(grid))deallocate(grid)
       ! end do
       ! close(44);close(45)
        !stop 'END LOOP'
       contains

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
                   k=apl_side_inv(i)!moltypeofuatom(apl_atoms_invatom(apl_upper_inv(i)))
               end if
           end do
       end function nneighbour

        elemental function mymodulo(a,b) result(c)
            real(kind=rk),intent(in) :: a,b
            real(kind=rk) :: c
            c=a-nint(a/b)*b !ok
        end function mymodulo

      ! subroutine mindist(atoms)
      ! integer(kind=ik) :: atoms(:),i,j
      ! real(kind=rk) :: r2
      ! rmin=huge(rmin)
      ! do i=2,size(atoms)
      !     do j=1,i-1
      !         r2=sum(modulo(coor(1:2,atoms(i)) -&
      !         coor(1:2,atoms(j)),box(1:2)/2)**2)
      !         rmin=min(rmin,r2)
      !     end do
      ! end do
      ! rmin=sqrt(rmin)
      ! end subroutine mindist

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
                    !if(apl_atoms(grid(i,j,side))==cind(atom,imol))then
                        imol=apl_atoms_invmol(grid(i,j,side))
                        instr%apl_side(imol)=side
                        k(imol)=k(imol)+1
                    end if
                end do
            end do
        end do
        instr%datam(:,frame)=real(k,rk)*box(1)*box(2)/(real(size(grid,1),rk)*real(size(grid,2),rk))
    end subroutine apl_calc

    !subroutine apl_matrix_out(k)
    !    integer(kind=ik) :: i,j,k
    !    do i=1,size(grid,2)
    !        write(42,*)real(grid(:,i,k),rk)
    !        do j=1,size(grid,1) !Gnuplot splot matrix, set dm3d ...
    !            write(43,*)i,j,real(grid(j,i,k),rk)
    !        end do
    !    end do

    !end subroutine apl_matrix_out

end module apl
