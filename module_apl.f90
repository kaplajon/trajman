module apl
    !use kinds
    use readtraj
    use trajop
    use util
    implicit none
    integer(kind=ik),allocatable ::&
    apl_atoms(:),apl_atoms_invatom(:),apl_atoms_invmol(:)
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
        apl_lower(:),apl_lower_inv(:),apl_upper(:),apl_upper_inv(:),grid(:,:)
        integer(kind=ik) ::i,j,kl,ku,dmpc,mgdg
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
        call mindist(apl_upper)
        !allocate(grid(1:ceiling((box(1)/rmin)*(sqrt(12._rk)+2)),1:ceiling((box(2)/rmin)*(sqrt(12._rk)+2))))
        allocate(grid(1:1000,1:1000))
        write(*,*)
        write(*,*)'Hej Grid!',shape(grid),size(apl_upper),size(apl_lower),size(apl_atoms)
        do i=1,size(grid,2)
            do j=1,size(grid,1)
                grid(j,i)=nneighbour(j,i,apl_upper)
                if(grid(j,i)==1)then
                    dmpc=dmpc+1
                else
                    mgdg=mgdg+1
                end if
            end do
        end do
        write(*,*)'AREA DMPC: ',(box(1)*box(2)*(real(dmpc,rk)/(real(size(grid,1),rk)*real(size(grid,2),rk))))/72._rk
        write(*,*)'AREA MGDG: ',(box(1)*box(2)*(real(mgdg,rk)/(real(size(grid,1),rk)*real(size(grid,2),rk))))/18._rk

        do i=1,size(grid,2)
            write(42,*)real(grid(:,i),rk)
            do j=1,size(grid,1)
            write(43,*)i,j,real(grid(j,i),rk)
            end do
        end do
        !stop 'END LOOP'
       contains

       function nneighbour(a,b,apl_side) result(k)
           integer(kind=ik) :: a,b,i,k,apl_side(:)
           real(kind=rk) :: gridcoor(2),r2,rmin
           gridcoor=[(box(1)/real(size(grid,1),rk))*real(a,rk),(box(2)/real(size(grid,2),rk))*real(b,rk)]
           rmin=huge(rmin)
           k=-1
           do i=1,size(apl_side)
              r2=sum(mymodulo(coor(1:2,apl_side(i)) -&
               gridcoor,box(1:2))**2)
               if(r2<rmin)then
                   rmin=r2
                   k=moltypeofuatom(apl_atoms_invatom(apl_upper_inv(i)))
               end if
           end do
       end function nneighbour

        elemental function mymodulo(a,b) result(c)
            real(kind=rk),intent(in) :: a,b
            real(kind=rk) :: c
            c=a-nint(a/b)*b !ok
        end function mymodulo

       subroutine mindist(atoms)
       integer(kind=ik) :: atoms(:),i,j
       real(kind=rk) :: r2
       rmin=huge(rmin)
       do i=2,size(atoms)
           do j=1,i-1
               r2=sum(modulo(coor(1:2,atoms(i)) -&
               coor(1:2,atoms(j)),box(1:2)/2)**2)
               rmin=min(rmin,r2)
           end do
       end do
       rmin=sqrt(rmin)
       end subroutine mindist
    end subroutine apl_grid


end module apl
