module apl
    use kinds
    use readtraj
    implicit none
    integer(kind=ik),allocatable :: apl_atoms(:),apl_atoms_inv(:)
    contains

    subroutine apl_atomlist(atoms)
        character(kind=1,len=1) :: atoms(:,:)
        integer(kind=ik) :: i,j,k,imol
        j=0
        do i=1,size(atoms,2)
            j=j+molt(moltypeofuatom(atomindex(trim(stringconv(atoms(:,i))))))%nmol
        end do
        if(allocated(apl_atoms))deallocate(apl_atoms)
        if(allocated(apl_atoms_inv))deallocate(apl_atoms_inv)
        allocate(apl_atoms(j),apl_atoms_inv(j))
        k=0
        do i=1,size(atoms,2)
            do imol=1,molt(moltypeofuatom(atomindex(trim(stringconv(atoms(:,i))))))%nmol
                !write(*,*)i,trim(stringconv(atoms(:,i))),cind(atomindex(trim(stringconv(atoms(:,i)))),imol)
                k=k+1
                apl_atoms(k)=cind(atomindex(trim(stringconv(atoms(:,i)))),imol)
                apl_atoms_inv(k)=atomindex(trim(stringconv(atoms(:,i))))
            end do
        end do
    end subroutine apl_atomlist
    
    subroutine apl_grid
        real(kind=rk) :: rmin,r2
        integer(kind=ik) ::i,j
        rmin=huge(rmin)
        do i=2,size(apl_atoms)
            do j=1,i-1
                r2=sum(modulo(coor(1:2,apl_atoms(i)) -&
                coor(1:2,apl_atoms(j)),box(1:2)/2)**2)
                rmin=min(rmin,r2)
            end do
        end do
        write(*,*)sqrt(rmin)
!! Glöm inte att separera bilagren!!
    end subroutine apl_grid


end module apl
