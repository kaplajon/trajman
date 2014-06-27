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
module readtraj
    use kinds
    use util
    implicit none
    !external :: f77_molfile_read_next,f77_molfile_init,f77_molfile_open_read,&
    !f77_molfile_close_read,f77_molfile_finish
    real (kind=rk),allocatable :: coor(:,:),box(:),masses(:),mgratios(:)
    character(kind=1, len=11),allocatable :: temp(:)
!    character(kind=1, len=11),allocatable :: atomnames(:)
   ! character(kind=1, len=5),allocatable :: moltypenames(:)
   ! character(kind=1, len=1),allocatable :: trajfile(:)
    character(kind=1,len=3) :: trajtype
    character(kind=1,len=*), parameter :: dlpoly3histfile='HISTORY'
    character(kind=1,len=*), parameter :: dlpoly3histtag='dlpoly3hist'
    integer (kind=ik) :: mols=0,nunit=13,atot,rowsperframe!,tunit=12
    integer(kind=4),allocatable :: tunit(:)!=12!,runit
    integer (kind=ik),allocatable ::&
    natoms(:),nmolsatoms(:),shift(:),moltypeofuatom(:)
    real(rk) :: director(1:3)=[0._rk,0._rk,1._rk],centerofmembrane(1:3)=0
    private :: atomindex_a,atomindex_b
    interface atomindex
        module procedure atomindex_a,atomindex_b
    end interface
!    type trjf
!        character(kind=1,len=1),allocatable :: filename(:)
!    end type trjf
!    type(trjf),allocatable :: trajfile(:)
!    type atomdata
!        character(kind=1,len=len(atomnames)) :: aname
!        real(kind=rk) :: mass,mgratio
!    end type atomdata
!    type(atomdata),allocatable :: atomd(:)
    interface operator(.str.)
        module procedure stringconv
    end interface operator(.str.)
contains

    subroutine wf_gro(filename,wf,funit)!{{{
    character(kind=1,len=*) :: filename,wf
    integer(kind=ik) :: funit,i,j,k,l
        open(funit,file=trim(filename))
        write(funit,*)'frame=',trim(adjustl(wf))
        write(funit,*)size(coor,2)
        l=1
        do i=1,size(molt)
            do j=1,molt(i)%nmol
                do k=molt(i)%firstatom,molt(i)%lastatom
                    !write(*,'(i5,a5,a10,i5,3F10.7)')j,molt(i)%molname,atomnames(k),l,getatom(k,j)
                    write(funit,'(i5,2a5,i5,3F8.3)')j,molt(i)%molname,atom(k)%aname,l,atom(k)%coor(:,j)
                    l=l+1
                end do
            end do
        end do
        write(funit,'(3F8.3)')box
        close(funit)
    end subroutine wf_gro!}}}

    subroutine wf_xyz(filename,funit)!{{{
    character(kind=1,len=*) :: filename
    integer(kind=ik) :: funit,i,j,k
        open(funit,file=trim(filename))
        write(funit,*)size(coor,2)
        write(funit,*)
        do i=1,size(molt)
            do j=1,molt(i)%nmol
                do k=molt(i)%firstatom,molt(i)%lastatom
                    write(funit,*)atom(k)%aname,10*atom(k)%coor(:,j)
                end do
            end do
        end do
        close(funit)
    end subroutine wf_xyz!}}}

subroutine reallocinstruct(v,i)!{{{
    type(instruct),intent(inout),allocatable :: v(:)
    type(instruct),allocatable ::copy(:)
    integer(kind=ik) :: i,j
    if (allocated(v))then
        j=min(i,size(v))
        allocate(copy(i))
        copy(1:j)=v(1:j)
        call move_alloc(copy,v)
        if(allocated(copy))deallocate(copy)
    else
        allocate(v(i))
    end if
    !do j=1,size(v)
    !write(*,*)j,allocated(v(j)%datam),allocated(v(j)%atoms)
    !end do
end subroutine reallocinstruct!}}}

!subroutine reallocinstruct(v,i)!{{{
!    class(instruct),allocatable :: v(:),copy(:)
!    !type(instruct),allocatable ::copy(:)
!    integer(kind=ik) :: i,j
!    if (allocated(v))then
!        j=min(i,size(v))
!        allocate(copy(i))
!        copy(1:j)=v(1:j)
!        call move_alloc(copy,v)
!    else
!        allocate(v(i))
!        
!    end if
!end subroutine reallocinstruct!}}}

    subroutine readline(onerow,ios,runit)!{{{
        character(kind=1,len=1),allocatable :: onerow(:) 
        integer(kind=ik) :: n,ios,readunit
        integer(kind=ik),optional,intent(in) :: runit
        readunit=stdin
        if (present(runit))readunit=runit
            comment:do
                n=0
                do
                    if(n+1>=size(onerow))call reallocate(onerow,int(2*size(onerow)+1,ik))
                    read(readunit,fmt="(1A1)",advance='NO',iostat=ios)onerow(n+1)!oneletter
                    if(ios==endf)exit comment
                    if(onerow(n+1)=='#' .or. ios==endr )then
                        if(n==0)then
                            if(ios/=endr)read(readunit,fmt='()',advance='YES',iostat=ios)
                            cycle comment
                        else
                            exit comment
                        endif
                    end if
                n=n+1
                enddo
            end do comment
        if(ios==endf)return       
        n=max(n,1)
       call reallocate(onerow,n)
    end subroutine readline !}}}

subroutine globals!{{{
    character(kind=1, len=30),allocatable :: defmass(:)
    integer(kind=ik) :: i,ios
    if(.NOT.allocated(atomd))then
        global_setflags%autofilename=.TRUE.
        global_setflags%folding=.FALSE.
        global_setflags%filename=''
        global_setflags%distbin=100
        global_setflags%fileprefix='auto_'
        global_setflags%filesuffix='.out'
    !    global_setflags%writeframe%framenumber=0
        global_setflags%wftot=0
        global_setflags%apl=.FALSE.
        global_setflags%gd=.FALSE.
        global_setflags%aplgrid=[250,250]
        global_setflags%leaflet=0
        global_setflags%whole=.FALSE.
        global_setflags%centerofmembrane=.FALSE.
        global_setflags%foldcenterofmembrane=.FALSE.
        global_setflags%leaflets_defined=.FALSE.
        global_setflags%molaverage=.FALSE.
        global_setflags%instructionsum=.FALSE.
        global_setflags%tshift=0
        common_setflags%traj_cscale=1
        global_setflags%rdf_binsize=0.2
        global_setflags%xyrdf=.FALSE.
        global_setflags%zrdf=.FALSE.
        global_setflags%const%value=0
        global_setflags%const%switch=.FALSE.
        global_setflags%Vsnorm=.TRUE.
        global_setflags%distminmax%switch=.FALSE.
        global_setflags%scaling%switch=.FALSE.
        global_setflags%ch_bondlength=1.09_rk
        global_setflags%slice%switch=.FALSE.
        global_setflags%slice%switch2=.FALSE.
        global_setflags%karplus=.FALSE.
        global_setflags%karplus_fnc=0
        global_setflags%coorsys=.FALSE.
    end if
    !maxframes=0;minframe=0
    !if(.NOT.allocated(global_setflags%calc))then
    !    allocate(global_setflags%calc(1))
    !    global_Setflags%calc=''
    !end if
    !global_setflags%silent=.FALSE.
    ! Default atom masses
    allocate(defmass(5))
    if(.not.allocated(atomd))allocate(atomd(size(defmass)))
    ! Masses and magnetogyric ratios from webelements.com 
    defmass=''
    defmass=[&
    !   Atom mass      magnetogyric
        'O   15.9994   -3.62808e7     ',&
        'C   12.0107   6.728284e7     ',&
        'H   1.00794   26.7522128e7   ',&
        'N   14.0067   -2.71261804e7  ',&
        'P   30.973762 10.8394e7      ' &
        ]
    do i=1,size(defmass)    
        read(defmass(i),*,iostat=ios)atomd(i)%aname,atomd(i)%mass,atomd(i)%mgratio
    end do
    deallocate(defmass)

   ! global_setflags%calc=''
end subroutine globals!}}}

    subroutine initgro(fname)!{{{
        ! Bestämmer antalet molekyltyper, atomer och indexering av dessa utifrån
        ! en frame.
        character(kind=1,len=1),intent(in) :: fname(:)
        character(kind=1,len=1),allocatable :: charvec(:)
   !     character(kind=1,len=11),allocatable :: temp(:)
        character(kind=1,len=10),allocatable :: moltype_atom(:,:)
        character(kind=1,len=10) :: sdr
        integer(kind=ik) :: ios,ia,i
       !stop 'HEJ'
        allocate(charvec(1))
        allocate(tunit(1))
        tunit=12
        open(unit=tunit(1),file=trim(stringconv(fname)),position='rewind',status='old',iostat=ios)

        if(ios/=0)then
                write(*,*)"Error, cannot open ",&
                "file ",":",trim(stringconv(fname)),":"
                stop
        endif

        read(unit=tunit(1),fmt=*,iostat=ios)
        if(ios==endf)return
        read(unit=tunit(1),fmt=*,iostat=ios)atot ! Totala antalet atomer i en frame
        allocate(moltype_atom(2,atot),coor(1:3,atot),box(1:3))
        moltype_atom=''
        do ia=1,atot
            read(tunit(1),"(5x,A10,5x)",advance='no')sdr
            sdr=adjustl(sdr)
            i=scan(sdr,' ')
            if(i>10.or.i<1)i=6
            moltype_atom(1,ia)=trim(sdr(1:i-1))
            moltype_atom(2,ia)=trim(adjustl(sdr(i:)))
            read(tunit(1),*)coor(1:3,ia)
            !write(*,*)coor(1:3,ia),'COOR'
        end do
            read(tunit(1),*)box
        !rewind(tunit)
        coor=coor*10._rk ! Scale to Å (default in .trr files)
        box=box*10._rk ! Scale to Å (default in .trr files)
        close(tunit(1))
        rowsperframe=atot+3
        call trajindex(moltype_atom)
        deallocate(charvec)
        deallocate(tunit)
    end subroutine initgro!}}}

    subroutine inithist(fname,frames)!{{{
        ! Bestämmer antalet molekyltyper, atomer och indexering av dessa utifrån
        ! en frame definderad av filen HISTORY från DL_POLY
        character(kind=1,len=1),intent(in) :: fname(:)
        integer(kind=ik), optional :: frames
        character(kind=1,len=10),allocatable :: moltype_atom(:,:)
        character(kind=1,len=10) :: sdr
        integer(kind=ik) :: ios,ia,i,iframe
        real(kind=rk) :: boxvec(3,3)
        character(kind=1,len=72) :: header
        integer(kind=ik) :: keytrj,imcon,megatm,frame,records,iatm,nstep
        character(kind=1,len=8) :: atmnam,timestep
        real(kind=rk) :: weight,charge,rsd,velocity(1:3),forces(1:3),tstep,time
        if(present(frames))then
           open(unit=11,file=trim(stringconv(fname)),position='rewind',status='old',iostat=ios)
           read(unit=11,fmt=*,iostat=ios)header
           if(ios==endf)return
           read(unit=11,fmt=*,iostat=ios)keytrj,imcon,megatm,frame,records
           frames=frame
           close(unit=11)
           return
        end if
        allocate(tunit(1))
        tunit(1)=12
        open(unit=tunit(1),file=trim(stringconv(fname)),position='rewind',status='old',iostat=ios)
        if(ios/=0)then
            write(*,*)"Error, cannot open ",&
                 "file ",":",trim(stringconv(fname)),":"
            stop
        endif
        read(unit=tunit(1),fmt=*,iostat=ios)header
        if(ios==endf)return
        read(unit=tunit(1),fmt=*,iostat=ios)keytrj,imcon,megatm,frame,records
        if(imcon/=1)write(*,*)"Warning, trajectory has not an orthorombic symmetry"
        atot=megatm
        allocate(moltype_atom(2,atot),coor(1:3,atot),box(1:3))
        read(tunit(1),fmt=*)timestep,nstep,megatm,keytrj,imcon,tstep,time
        read(unit=tunit(1),fmt=*,iostat=ios)boxvec(1,1),boxvec(2,1),boxvec(3,1)
        read(unit=tunit(1),fmt=*,iostat=ios)boxvec(1,2),boxvec(2,2),boxvec(3,2)
        read(unit=tunit(1),fmt=*,iostat=ios)boxvec(1,3),boxvec(2,3),boxvec(3,3)
        box(1)=boxvec(1,1)
        box(2)=boxvec(2,2)
        box(3)=boxvec(3,3)
        moltype_atom=''
        do ia=1,atot
            read(tunit(1),*)atmnam,iatm,weight,charge,rsd
            sdr=''
            sdr=adjustl(atmnam)
            moltype_atom(1,ia)=trim(sdr)//'M'! M means molecule
            moltype_atom(2,ia)=trim(sdr)
            read(tunit(1),*)coor(1:3,ia)! unit is already in Ångström
            if(keytrj>=1)read(tunit(1),*)velocity(1:3)
            if(keytrj>=2)read(tunit(1),*)forces(1:3)
        end do
! the rest of the HISTORY file has the following format:
!        read(tunit(1),*)timestep,nstep,megatm,keytrj,imcon,tstep,time
!        read(unit=tunit(1),fmt=*,iostat=ios)boxvec(1,1),boxvec(2,1),boxvec(3,1)
!        read(unit=tunit(1),fmt=*,iostat=ios)boxvec(1,2),boxvec(2,2),boxvec(3,2)
!        read(unit=tunit(1),fmt=*,iostat=ios)boxvec(1,3),boxvec(2,3),boxvec(3,3)
!        do iframe=1,frame
!           do ia=1,megatm
!              read(tunit(1),*)coor(1:3,ia)! unit is already in Ångström
!            if(keytrj>=1)read(tunit(1),*)velocity(1:3)
!            if(keytrj>=2)read(tunit(1),*)forces(1:3)
!           end do
!        end do
        close(tunit(1))
        deallocate(tunit)
        call trajindex(moltype_atom)
    end subroutine inithist!}}}

    function readgro(fhandle) result(i)!{{{
        integer(kind=4) :: fhandle
        integer(kind=ik) :: i,ios,ia
        read(unit=fhandle,fmt=*,iostat=ios)
        select case(ios)
        case(0)
            i=ios
        case(endf)
            if(fhandle==tunit(size(tunit)))then
                i=1 !for the frameloop to work as expected.
            else
                i=2
            end if
        end select
        if(i==1)return
        read(unit=fhandle,fmt=*,iostat=ios)atot
        if(allocated(coor))stop 'READGRO: Coor is already allocated, why?'
        allocate(coor(1:3,atot))
        do ia=1,atot
            read(unit=fhandle,fmt='(20x)',iostat=ios,advance='no')!coor(:,ia)
            read(unit=fhandle,fmt=*,iostat=ios)coor(:,ia)
        end do
        read(unit=fhandle,fmt=*,iostat=ios)box(:)
        coor=coor*10._rk ! Scale to Å (default in .trr files)
        box=box*10._rk ! Scale to Å (default in .trr files)

    end function readgro!}}}

    function readtrr(fhandle) result(ios) !MOLFILEPLUGIN FROM VMD!{{{
        integer(kind=4) :: fhandle,statf,natm
        integer(kind=ik) :: ios
        real(kind=4) :: bx(6)
        real(kind=4),allocatable :: coorv(:)
        allocate(coorv(atot*3))
        natm=int(atot,4)
        statf=1
        call f77_molfile_read_next(fhandle,natm,coorv,bx,statf)
        if(allocated(coor))stop 'READTRR: Coor is already allocated, why?'
        allocate(coor(1:3,atot))
        coor=reshape(real(coorv,kind(coor)),[3,size(coor,2)],PAD=0._rk*coor)
        box=bx(1:3)
        select case(statf)
        case(0)
            if(fhandle==tunit(size(tunit)))then
                ios=1 !for the frameloop to work as expected.
            else
                ios=2
            end if
        case default
            ios=0 !molfile stat and iostat conversion 
        end select
    end function readtrr!}}}

    subroutine closetraj(funit)!{{{
        integer(kind=4) :: stat,funit
        select case(trajtype)
        case('gro')
            close(funit)
        case('trr','dcd','pdb')
            call f77_molfile_close_read(funit,stat)
            if(funit==tunit(size(tunit)))call f77_molfile_finish
        case default
            call f77_molfile_close_read(funit,stat)
            if(funit==tunit(size(tunit)))call f77_molfile_finish
        end select
    end subroutine closetraj!}}}

    function readframe(tunit) result(ios)!{{{
        integer(kind=ik) :: ios,i,j
        integer(kind=4) :: tunit
        select case(trajtype)
        case('gro')
            ios=readgro(tunit)
        case('trr','dcd','pdb')
            ios=readtrr(tunit)
        case default
            ios=readtrr(tunit)
        end select
        do i=1,size(atom)-atomsdefined
            do j=1,molt(atom(i)%moltype)%nmol
                atom(i)%coor(:,j)=getatom(i,j)
            end do
        end do
        ! Make sure everything uses the new coordinate handling
        deallocate(coor) 
    end function readframe        !}}}

    pure function stringconv(vector) result(string)!{{{
    ! Convert a vector to a string

        character(kind=1,len=1),intent(in) :: vector(:)
        character(kind=1,len=size(vector)) :: string
        integer(kind=ik) :: i
        do i=1,len(string)
            string(i:i)=vector(i)(1:1)        
        enddo
    end function stringconv!}}}

    function atomindex_a(a,b,n) result(ind)!{{{
        character(kind=1,len=*),intent(in) :: a,b(:)
        integer (kind=ik) :: i,n,ind
        ind=0
        do i=1,n
            if(b(i)==a)then
                if(len_trim(b(i))==len_trim(a))then
                    ind=i
                    exit
                end if
            endif
        end do
        
    end function atomindex_a!}}}

    function atomindex_b(a) result(ind)!{{{
        character(kind=1,len=*),intent(in) :: a
        integer (kind=ik) :: i,ind
        ind=0
        do i=1,size(atomnames)
            if(trim(atomnames(i))==trim(a))then
                if(len_trim(atomnames(i))==len_trim(a))then
                    ind=i
                    exit
                end if
            endif
        end do
        
    end function atomindex_b!}}}

    elemental function cind(iuatom,imol) result(ind)!{{{

        integer(kind=ik),intent(in) :: imol,iuatom
        integer(kind=ik) :: ind
        ! Index i coor för atom iuatom i molekyl imol
        ind=shift(iuatom)+natoms(iuatom)*imol
    end function cind!}}}

    pure function getatom(aindex,imol) result(acoor)!{{{

        integer (kind=ik),intent(in) :: aindex,imol
        real (kind=rk) :: acoor(1:3)
        !if(.not.allocated(coor))stop 'getatom: coor is not allocated!'

        acoor(1:3)=coor(1:3,cind(aindex,imol))
        
    end function  getatom!}}}

    subroutine trajindex(moltype_atom)!{{{
        character(kind=1,len=11) :: uatom
   !     character(kind=1,len=11),allocatable :: temp(:)
        character(kind=1,len=10),allocatable ::&
        temp2(:),temp3(:),moltype_atom(:,:)
        integer(kind=ik),allocatable :: natomsoftype(:),nmols(:)
        character(kind=1, len=5),allocatable :: moltypenames(:),test(:)
        integer(kind=ik) :: ia,i,j,atoms=0       
        allocate(atomnames(atot),moltypenames(atot),temp(atot),temp2(atot),temp3(atot),test(3))
        atomnames="";moltypenames="";temp="";temp2="";temp3=""
!test=['OME','0LB','3YB']
!write(*,*)test
!        do i=1,atot
!        write(*,*)moltype_atom(1,i),moltype_atom(2,i)
!        enddo
!        stop
        do ia=1,atot

            uatom=""
            uatom=trim(moltype_atom(2,ia))//"_"//trim(adjustl(moltype_atom(1,ia))) ! Sätter unikt atomnamn
            if(allocated(groupres))then
                if(ANY(groupres(1:size(groupres)-1)==trim(moltype_atom(1,ia))))then
                    moltype_atom(1,ia)=trim(groupres(size(groupres)))
                end if
            end if
             ! Alla atomnamn i traj (atot)
            temp(ia)=uatom
             ! Alla förekomster av molekyltyp i traj (atot) 
             temp2(ia)=trim(moltype_atom(1,ia))
            
            if(atomindex(uatom,atomnames,atoms)==0)then
                atoms=atoms+1 ! Räknar unika atomnamn
                atomnames(atoms)=uatom ! Alla unika atomnamn (atoms)
                temp3(atoms)=trim(moltype_atom(1,ia))
            endif

            if(atomindex(trim(moltype_atom(1,ia)),moltypenames,mols)==0)then
                mols=mols+1 ! Räknar molekyltyper
                !Alla unika molekylnamn (mols)
                moltypenames(mols)=trim(moltype_atom(1,ia))
            endif
        end do
        allocate(nmols(atoms),natoms(atoms),shift(atoms),nmolsatoms(mols),&
        natomsoftype(mols),moltypeofuatom(atoms),molt(mols))
        shift=0
        nmols=0;natoms=0;nmolsatoms=0
        do ia=1,atot
         ! Vilket index har atom ia i den unika atomnamnslistan?
            i=atomindex(temp(ia),atomnames,atoms)
            nmols(i)=nmols(i)+1 ! Räknar molekyler med atom ia

            ! Vilket index har motsvarande molekyltyp ia? 
            j=atomindex(temp2(ia),moltypenames,mols)
            nmolsatoms(j)=nmolsatoms(j)+1 ! Räknar totala antalet atomer i molekyltyp

        end do
        do ia=1,atoms
            !i=atomindex(atomnames(ia),atomnames,atoms) 
            j=atomindex(temp3(ia),moltypenames,mols)
            natoms(ia)=nmolsatoms(j)/nmols(ia) ! Atomer ia per molekyl ia
            natomsoftype(j)=natoms(ia) ! Atomer per molekyltyp

            ! Skiftet som krävs i coor för att hitta samma atom i alla
            ! molekyler. Första summan skiftar tidigare molekyltyper. natoms
            ! skiftar en molekyl. Sista summan skiftar ia med tidigare
            ! molekyltyper.
            shift(ia)=sum(nmolsatoms(1:j-1))-natoms(ia)+ia-sum(natomsoftype(1:j-1))
            moltypeofuatom(ia)=j
        end do

        do j=1,mols
        do i=1,atoms
            if(moltypeofuatom(i)==j)then
                molt(j)%firstatom=i
                exit
            end if
        end do
        molt(j)%lastatom=molt(j)%firstatom+natomsoftype(j)-1
        molt(j)%nmol=nmols(molt(j)%lastatom)
        molt(j)%molname=moltypenames(j)
        molt(j)%natoms=natomsoftype(j)
        !write(*,*)molt(j)%molname
        end do

        temp=atomnames;temp2=moltypenames
        deallocate(atomnames,moltypenames)
        allocate(atomnames(atoms))
        atomnames=temp(1:atoms)
        deallocate(temp,temp2,temp3,natomsoftype,moltype_atom)
!        allocate(masses(size(atomnames)),mgratios(size(atomnames)))
!        masses=0;mgratios=0
        !Get away from this harsh indexing by storing coordinates in a new type
        allocate(atom(size(atomnames)))
        do i=1,size(atom)
            j=moltypeofuatom(i)
            allocate(atom(i)%coor(3,molt(j)%nmol))
            atom(i)%aname=atomnames(i)
            atom(i)%mname=molt(j)%molname
            atom(i)%moltype=j
            do j=1,size(atom(i)%coor,2)
                atom(i)%coor(:,j)=getatom(i,j)
            end do
        end do
        deallocate(coor) ! To make sure everything uses the new coordinate handling
    end subroutine trajindex!}}}

    pure function center_of_molecule(umol,imol) result(centerofmolecule)!{{{
        integer(kind=ik),intent(in) :: umol,imol
        integer(kind=ik) :: i,j
        real(kind=rk) :: centerofmolecule(3)
        centerofmolecule=0
        j=umol
        do i=molt(j)%firstatom,molt(j)%lastatom
            centerofmolecule=centerofmolecule+atom(i)%coor(:,imol)*masses(i)
        end do
        centerofmolecule=centerofmolecule/sum(masses(molt(j)%firstatom:molt(j)%lastatom))
       ! write(*,*)centerofmolecule,sum(masses(molt(j)%firstatom:molt(j)%lastatom))
    end function center_of_molecule!}}}
end module readtraj
