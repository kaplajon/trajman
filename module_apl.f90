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
    apl_atoms(:),apl_atoms_invatom(:),apl_atoms_invmol(:)&
    ,apl_both(:),apl_both_inv(:),apl_both_invmol(:)&
    ,apl_lower(:),apl_lower_inv(:),apl_lower_invmol(:)&
    ,apl_upper(:),apl_upper_inv(:),apl_upper_invmol(:)
    integer(kind=ik),target,allocatable ::grid_both(:,:),grid_lower(:,:),grid_upper(:,:)
    real(kind=rk) :: meanbox(4)=0
    contains
    subroutine apl_atomlist(atoms)!{{{
        character(kind=1,len=1) :: atoms(:,:)
        integer(kind=ik) :: i,j,k,imol,kl,ku
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
        integer(kind=ik),pointer :: grid(:,:)
        integer(kind=ik) ::i,j,kl,ku
        
        i=size(apl_atoms)
        if(allocated(apl_both))deallocate(apl_both,apl_both_inv,apl_both_invmol)
        if(allocated(apl_lower))deallocate(apl_lower,apl_lower_inv,apl_lower_invmol)
        if(allocated(apl_upper))deallocate(apl_upper,apl_upper_inv,apl_upper_invmol)
        allocate(apl_both(i),apl_both_inv(i),apl_both_invmol(i))
        allocate(apl_lower(i),apl_lower_inv(i),apl_lower_invmol(i))
        allocate(apl_upper(i),apl_upper_inv(i),apl_upper_invmol(i))
        kl=0;ku=0
                apl_both=apl_atoms
                apl_both_invmol=apl_atoms_invmol
                apl_both_inv=[(i,i=1,size(apl_atoms))]
                do i=1,size(apl_atoms)
                    if(dot_product(center_of_molecule(moltypeofuatom(apl_atoms_invatom(i))&
                    ,apl_atoms_invmol(i))-centerofmembrane,director)<0._rk)then
                        kl=kl+1
                        apl_lower(kl)=apl_atoms(i)
                        apl_lower_invmol(kl)=apl_atoms_invmol(i)
                        apl_lower_inv(kl)=i
                    else
                        ku=ku+1
                        apl_upper(ku)=apl_atoms(i)
                        apl_upper_invmol(ku)=apl_atoms_invmol(i)
                        apl_upper_inv(ku)=i
                    end if
                end do
        call reallocate(apl_lower,kl)
        call reallocate(apl_lower_inv,kl)
        call reallocate(apl_lower_invmol,kl)
        call reallocate(apl_upper,ku)
        call reallocate(apl_upper_inv,ku)
        call reallocate(apl_upper_invmol,ku)
        
        select case(instr%set%apl)
        case(.TRUE.)
            if(allocated(grid_both))deallocate(grid_both)
            allocate(grid_both(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
            do i=1,size(grid_both,2)
                do j=1,size(grid_both,1)
                    grid_both(j,i)=nneighbour(j,i,apl_both_invmol,apl_both_inv,grid_both)
                end do
            end do
            if(allocated(grid_lower))deallocate(grid_lower)
            allocate(grid_lower(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
            do i=1,size(grid_lower,2)
                do j=1,size(grid_lower,1)
                    grid_lower(j,i)=nneighbour(j,i,apl_lower_invmol,apl_lower_inv,grid_lower)
                end do
            end do
            if(allocated(grid_upper))deallocate(grid_upper)
            allocate(grid_upper(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
            do i=1,size(grid_upper,2)
                do j=1,size(grid_upper,1)
                    grid_upper(j,i)=nneighbour(j,i,apl_upper_invmol,apl_upper_inv,grid_upper)
                end do
            end do
        case(.FALSE.)
            if(allocated(grid_both))deallocate(grid_both)
            allocate(grid_both(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
            grid=>grid_both
            call countmol(apl_both_inv,apl_both_invmol,grid)
            if(allocated(grid_lower))deallocate(grid_lower)
            allocate(grid_lower(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
            grid=>grid_lower
            call countmol(apl_lower_inv,apl_lower_invmol,grid)
            if(allocated(grid_upper))deallocate(grid_upper)
            allocate(grid_upper(1:instr%set%aplgrid(1),1:instr%set%aplgrid(2)))
            grid=>grid_upper
            call countmol(apl_upper_inv,apl_upper_invmol,grid)
        end select

       contains
       function nneighbour(a,b,apl_side_invmol,apl_side_inv,grid) result(k)
           integer(kind=ik) :: a,b,i,k,apl_side_invmol(:),apl_side_inv(:)
           integer(kind=ik),allocatable :: grid(:,:)
           real(kind=rk) :: gridcoor(2),r2,rmin
           gridcoor=[(box(1)/real(size(grid,1),rk))*real(a,rk),(box(2)/real(size(grid,2),rk))*real(b,rk)]
           rmin=huge(rmin)
           k=-1
           do i=1,size(apl_side_inv)
              r2=sum(mymodulo(atom(apl_side_inv(i))%coor(1:2,apl_side_invmol(i)) -&
               gridcoor,box(1:2))**2)
               if(r2<rmin)then
                   rmin=r2
                   k=apl_side_inv(i)
               end if
           end do
       end function nneighbour
    end subroutine apl_grid!}}}

    subroutine countmol(apl_side_inv,apl_side_invmol,g)!{{{
           integer(kind=ik) :: apl_side_inv(:),apl_side_invmol(:),i&
           ,ma1,mi1,ma2,mi2,bi1,bi2
           integer(kind=ik) :: g(:,:)
           real(kind=rk) :: bin1,bin2
           g=0
           ma1=box(1);mi1=0
           ma2=box(2);mi2=0
           bin1=(ma1-mi1)/real(size(g,1),rk)
           bin2=(ma2-mi2)/real(size(g,2),rk)
           do i=1,size(apl_side_inv)
               bi1=modulo(int(((atom(apl_side_inv(i))%coor(1,apl_side_invmol(i))-mi1)/bin1+1._rk))-1,size(g,1))+1
               bi2=modulo(int(((atom(apl_side_inv(i))%coor(2,apl_side_invmol(i))-mi2)/bin2+1._rk))-1,size(g,2))+1
               g(bi1,bi2)=g(bi1,bi2)+1
           end do
    end subroutine countmol!}}}

    subroutine apl_calc(instr,frame)!{{{
        type(instruct) :: instr
        integer(kind=ik),pointer :: grid(:,:)
        integer(kind=ik),intent(in) :: frame
        integer(kind=ik) :: i,j,l,imol,k(instr%nmolop)
        select case(instr%set%leaflet)
        case(0)
            grid=>grid_both
        case(1)
            grid=>grid_lower
        case(2)
            grid=>grid_upper
        case default
            stop 'APL'
        end select
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
        nullify(grid)
    end subroutine apl_calc!}}}

    function griddiff(grid1,grid2) result(res)!{{{
        real(kind=rk) :: res
        integer(kind=ik) :: grid1(:,:),grid2(:,:)
        real(kind=rk) :: rg1(size(grid1,1),size(grid1,2)),rg2(size(grid2,1),size(grid2,2))
        rg1=real(grid1,rk)/real(sum(grid1),rk)
        rg2=real(grid2,rk)/real(sum(grid2),rk)
        !res=(sqrt(real(sum((grid1-grid2)**2),rk))*sum(grid1)/real(size(grid1),rk))
        res=sum(abs(rg1-rg2))/2._rk
        !res=sum(abs(real(((real(grid1,rk)/real(sum(grid1),rk))-(real(grid2,rk)/real(sum(grid2),rk))),rk)))/2._rk
    end function griddiff!}}}

    subroutine apl_matrix_out(frame,instr)!{{{
        integer(kind=ik) :: i,j,frame
        type(instruct) :: instr 
        integer(kind=ik),pointer :: grid(:,:)
        integer(kind=ik),allocatable :: v1(:),v2(:)
        character(len=5) :: side
        allocate(v1(size(grid_lower,1)),v2(size(grid_lower,1)))
        side=''
        nullify(grid)
        select case(instr%set%leaflet)
            case(1)
                side="lower"
                grid=>grid_lower
            case(2)
                side="upper"
                grid=>grid_upper
            case(0)
                side="both"
                grid=>grid_both
            case default
                stop 'Matrix'
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
                v2(i)=atom(apl_atoms_invatom(grid(i,j)))%moltype
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
