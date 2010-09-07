!---LICENSE-------------------------------------------------------!{{{
! This file is part of
!
!  Trajman: A MD Trajectory Analysis Tool
!
!  Copyright (C) 2010, Jon Kapla, jon.kapla@mmk.su.se
!  Department of Material and Environmental Chemistry,
!  Stockholm University
!
! Trajman is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Trajman is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Trajman.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------!}}}
module apl
    use readtraj
    use util
    implicit none
    integer(kind=ik),allocatable ::&
    apl_atoms(:),apl_atoms_invatom(:),apl_atoms_invmol(:),grid(:,:)&
    ,apl_side(:),apl_side_inv(:),apl_side_invmol(:)
    real(kind=rk) :: meanbox(4)=0
    contains
    subroutine apl_atomlist(atoms)!{{{
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
    end subroutine apl_atomlist!}}}
    
    subroutine apl_grid(instr)!{{{
        type(instruct) :: instr
        real(kind=rk) :: rmin
        integer(kind=ik) ::i,j,kl,ku,dmpc,mgdg,g
        i=size(apl_atoms)
        if(allocated(apl_side))deallocate(apl_side,apl_side_inv,apl_side_invmol)
        allocate(apl_side(i),apl_side_inv(i),apl_side_invmol(i))
        kl=0;ku=0;dmpc=0;mgdg=0
        select case(instr%set%leaflet)
            case(0) !BOTH
                apl_side=apl_atoms
                apl_side_invmol=apl_atoms_invmol
                apl_side_inv=[(i,i=1,size(apl_atoms))]
                kl=size(apl_atoms)
            case(1) !LOWER
                do i=1,size(apl_atoms)
                    if(dot_product(center_of_molecule(moltypeofuatom(apl_atoms_invatom(i))&
                    ,apl_atoms_invmol(i))-centerofmembrane,director)<0._rk)then
                        kl=kl+1
                        apl_side(kl)=apl_atoms(i)
                        apl_side_invmol(kl)=apl_atoms_invmol(i)
                        apl_side_inv(kl)=i
                    end if
                end do
            case(2) !UPPER
                do i=1,size(apl_atoms)
                    if(dot_product(center_of_molecule(moltypeofuatom(apl_atoms_invatom(i))&
                    ,apl_atoms_invmol(i))-centerofmembrane,director)>0._rk)then
                        kl=kl+1
                        apl_side(kl)=apl_atoms(i)
                        apl_side_invmol(kl)=apl_atoms_invmol(i)
                        apl_side_inv(kl)=i
                    end if
                end do
            end select
        call reallocate(apl_side,kl)
        call reallocate(apl_side_inv,kl)
        call reallocate(apl_side_invmol,kl)
        if(allocated(grid))deallocate(grid)
        allocate(grid(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
        do i=1,size(grid,2)
            do j=1,size(grid,1)
                grid(j,i)=nneighbour(j,i,apl_side,apl_side_inv)
            end do
        end do
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
                   k=apl_side_inv(i)
               end if
           end do
       end function nneighbour
    end subroutine apl_grid!}}}

    subroutine apl_calc(instr,frame)!{{{
        type(instruct) :: instr
        integer(kind=ik),intent(in) :: frame
        integer(kind=ik) :: i,j,l,imol,k(instr%nmolop)
        k=0
            do j=1,size(grid,2)
                do i=1,size(grid,1)
                    if(ANY(instr%atoms==apl_atoms_invatom(grid(i,j))))then
                        imol=apl_atoms_invmol(grid(i,j))
                        do l=1,size(instr%molind)
                            if(imol==instr%molind(l))then
                                k(l)=k(l)+1
                            end if
                        end do
                    end if
                end do
            end do
        !instr%datam(:,frame)=real(k,rk)*box(1)*box(2)/(real(size(grid,1),rk)*real(size(grid,2),rk))
        instr%datam(:,frame)=real(k,rk)/(real(size(grid,1),rk)*real(size(grid,2),rk))
    end subroutine apl_calc!}}}

    subroutine apl_matrix_out(frame,instr)!{{{
        integer(kind=ik) :: i,j,frame
        type(instruct) :: instr 
        integer(kind=ik),allocatable :: v1(:),v2(:)
        character(len=5) :: side
        allocate(v1(size(grid,1)),v2(size(grid,1)))
        side=''
        select case(instr%set%leaflet)
            case(1)
                side="lower"
            case(2)
                side="upper"
            case(0)
                side="both"
        end select
        open(43,file='apl_polygons_'//trim(side)//'_grid'&
        //trim(adjustl(intstr(size(grid,1))))//'.'//trim(adjustl(intstr(size(grid,2))))&
        //'_frame'//trim(adjustl(intstr(frame))))
        open(44,file='apl_polygons_'//trim(side)//'_grid'&
        //trim(adjustl(intstr(size(grid,1))))//'.'//trim(adjustl(intstr(size(grid,2))))&
        //'_frame'//trim(adjustl(intstr(frame)))//'_mtrx')
        do j=1,size(grid,2)
            do i=1,size(v1)
                v1(i)=apl_atoms_invmol(grid(i,j))
                v2(i)=moltypeofuatom(apl_atoms_invatom(grid(i,j)))
                write(43,*)i,j,v1(i),v2(i)
                v1(i)=v1(i)+sum(molt(1:v2(i)-1)%nmol)
            end do
                write(43,*)
                write(44,*)v1(:)
        end do
        close(43);close(44)
        deallocate(v1,v2)
    end subroutine apl_matrix_out!}}}
end module apl
