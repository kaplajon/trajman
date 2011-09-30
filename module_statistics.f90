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
module statistics
    use readtraj
    use util
    implicit none
    integer(kind=ikr),parameter :: kk=100,ll=37,mm=2_ikr**(8_ikr*ikr-2_ikr)
    integer(kind=ikr) :: ranx(kk)
    contains
    subroutine rnarry(aa)!{{{
!     FORTRAN 77 version of "ran_array"
!     from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
!     including the MODIFICATIONS made in the 9th printing (2002)
!     ********* see the book for explanations and caveats! *********
!     Translated to fortran90 2010-09-15 by Jon Kapla.
    integer(kind=ikr) :: j,aa(:),n
    n=size(aa)
    do j=1,kk
        aa(j)=ranx(j)
    end do !bY jk 
    do j=kk+1,n
       aa(j)=aa(j-kk)-aa(j-ll)
       if (aa(j) .lt. 0) aa(j)=aa(j)+mm
    end do !bY jk
    do j=1,ll
       ranx(j)=aa(n+j-kk)-aa(n+j-ll)
       if (ranx(j) .lt. 0) ranx(j)=ranx(j)+mm
    end do !bY jk
    do j=ll+1,kk
       ranx(j)=aa(n+j-kk)-ranx(j-ll)
       if (ranx(j) .lt. 0) ranx(j)=ranx(j)+mm
    end do !bY jk
    end subroutine rnarry!}}}
    subroutine rnstrt(seed)!{{{
    integer(kind=ikr),parameter :: tt=70,kkk=kk+kk-1
    integer(kind=ikr) :: x(kkk),seed,j,ss,sseed,t
    if (seed .lt. 0) then
       sseed=mm-1-mod(-1-seed,mm)
    else
       sseed=mod(seed,mm)
    end if
    ss=sseed-mod(sseed,2)+2
    do j=1,kk
       x(j)=ss
       ss=ss+ss
       if (ss .ge. mm) ss=ss-mm+2
    end do
    x(2)=x(2)+1
    ss=sseed
    t=tt-1
    do while(t .gt. 0)
    do j=kk,2,-1
       x(j+j-1)=x(j)
       x(j+j-2)=0
    end do
    do j=kkk,kk+1,-1
       x(j-(kk-ll))=x(j-(kk-ll))-x(j)
       if (x(j-(kk-ll)) .lt. 0) x(j-(kk-ll))=x(j-(kk-ll))+mm
       x(j-kk)=x(j-kk)-x(j)
       if (x(j-kk) .lt. 0) x(j-kk)=x(j-kk)+mm
    end do
    if (mod(ss,2) .eq. 1) then
       do j=kk,1,-1
       x(j+1)=x(j)
       end do
       x(1)=x(kk+1)
       x(ll+1)=x(ll+1)-x(kk+1)
       if (x(ll+1) .lt. 0) x(ll+1)=x(ll+1)+mm
    end if
    if (ss .ne. 0) then
        ss=ss/2
    else
        t=t-1
    end if
    end do !while
    do j=1,ll
      ranx(j+kk-ll)=x(j)
    end do
    do j=ll+1,kk
      ranx(j-ll)=x(j)
    end do
    do j=1,10
      call rnarry(x)
    end do
    end subroutine rnstrt!}}}
    subroutine randint(V)!{{{
        integer(kind=ik) :: V(:)
        integer(kind=ikr) :: VR(max(1009,size(V)))
        do
        call RNARRY(VR)
        if(.not.any(VR>=size(V)*int(MM/size(V))))exit
        end do
        V=modulo(VR(1:size(V)),int(size(V),ikr))+1
    end subroutine randint!}}}
    subroutine randint_s(V)!{{{
        integer(kind=ik) :: V(:),i
        integer(kind=ikr) :: VR(max(1009,size(V)))
        a: do
        call RNARRY(VR)
        do i=1,size(V)
        if(.not.VR(i)>=(size(V)+1-i)*int(MM/(size(V)+1-i)))exit a
        end do
        end do a
        do i=1,size(V)
        V(i)=modulo(VR(i),int(size(V)+1-i,ikr))+1
        end do
        !V=modulo(VR(1:size(V)),int(size(V),ikr))+1
    end subroutine randint_s!}}}
    subroutine shuffle(V,VSHUF)!{{{
        integer(kind=ik),intent(in) :: V(:)
        integer(kind=ik),intent(out) :: VSHUF(size(V))
        integer(kind=ik) :: RAND(size(V)),i,j,k
        logical :: VLOGIC(size(V)) 
        VLOGIC=.false.
        VSHUF=-1
        call randint_s(RAND)
       ! write(*,*)RAND
        do i=1,size(V)
            k=0
            do j=1,size(V)
                if(.not.VLOGIC(j))k=k+1
                if(k==RAND(i))then
                    VLOGIC(j)=.true.
                    VSHUF(i)=V(j)
                    exit
                end if
            end do
        end do
    end subroutine shuffle!}}}

    subroutine atomshuffle(shuffle_atoms,atoms,side,shuffled_atom_cindexes)!{{{
        integer(kind=ik) :: shuffle_atoms(:),atoms(:),i,j,k,l,m,imol,side
        integer(kind=ik),allocatable :: cindexes(:),cindexes_invatom(:),shuffled_cindexes(:),shuffled_atom_cindexes(:)
        if(allocated(cindexes))deallocate(cindexes,cindexes_invatom,shuffled_cindexes)
        if(allocated(shuffled_atom_cindexes))deallocate(shuffled_atom_cindexes)
        j=0
        do i=1,size(shuffle_atoms)
            if(side==0)j=j+molt(moltypeofuatom(shuffle_atoms(i)))%nmol
            if(side==1)j=j+size(molt(moltypeofuatom(shuffle_atoms(i)))%lower)
            if(side==2)j=j+size(molt(moltypeofuatom(shuffle_atoms(i)))%upper)
        end do
        m=0
        do l=1,size(atoms)
            if(side==0)m=m+molt(moltypeofuatom(atoms(l)))%nmol
            if(side==1)m=m+size(molt(moltypeofuatom(atoms(l)))%lower)
            if(side==2)m=m+size(molt(moltypeofuatom(atoms(l)))%upper)
        end do
        allocate(cindexes(j),shuffled_cindexes(j),cindexes_invatom(j))
        allocate(shuffled_atom_cindexes(m))
        m=1
        !Storing molecular indexes in cindexes temporarily
        select case(side)
        case(0)
            do i=1,size(shuffle_atoms)
                cindexes(m:)=[molt(moltypeofuatom(shuffle_atoms(i)))%lower,molt(moltypeofuatom(shuffle_atoms(i)))%upper]
                cindexes_invatom(m:)=shuffle_atoms(i)
                m=m+(size(molt(moltypeofuatom(shuffle_atoms(i)))%lower)+size(molt(moltypeofuatom(shuffle_atoms(i)))%upper))
            end do
        case(1)
            do i=1,size(shuffle_atoms)
                cindexes(m:)=molt(moltypeofuatom(shuffle_atoms(i)))%lower
                cindexes_invatom(m:)=shuffle_atoms(i)
                m=m+size(molt(moltypeofuatom(shuffle_atoms(i)))%lower)
            end do
            !write(*,*)cindexes
            !write(*,*)            
            !call shuffle(cindexes,shuffled_cindexes)
            !do i=1,size(shuffled_cindexes)
            !write(33,*)shuffled_cindexes(i)
            !end do
        case(2)
            do i=1,size(shuffle_atoms)
                cindexes(m:)=molt(moltypeofuatom(shuffle_atoms(i)))%upper
                cindexes_invatom(m:)=shuffle_atoms(i)
                m=m+size(molt(moltypeofuatom(shuffle_atoms(i)))%upper)
            end do
           ! write(*,*)'-------------------------------------------'
           ! write(*,*)cindexes
           ! write(*,*)            
           ! call shuffle(cindexes,shuffled_cindexes)
           ! do i=1,size(shuffled_cindexes)
           ! write(34,*)shuffled_cindexes(i)
           ! end do
           ! stop
        end select
        !Transfer the molecular indexes to atomic indexes in
        do i=1,size(cindexes)
            cindexes(i)=cind(cindexes_invatom(i),cindexes(i))
            !write(*,*)cindexes(i),cindexes_invatom(i),i
        end do
        call shuffle(cindexes,shuffled_cindexes) !Shuffle the cindexes
        k=0
        do i=1,size(atoms) !Pick randomized atom indexes
            do j=1,size(shuffled_cindexes)
            if(cindexes_invatom(j)==atoms(i))then
                k=k+1
                shuffled_atom_cindexes(k)=shuffled_cindexes(j)
            end if
            end do
        end do
    end subroutine atomshuffle!}}}

end module statistics
