module readtraj
    use kinds
    use util
    implicit none
    real (kind=rk),allocatable :: coor(:,:),box(:),masses(:),mgratios(:)
    character(kind=1, len=11),allocatable :: temp(:)
    character(kind=1, len=11),allocatable :: atomnames(:)
   ! character(kind=1, len=5),allocatable :: moltypenames(:)
    character(kind=1, len=1),allocatable :: trajfile(:)
    integer (kind=ik) :: mols=0,tunit=12,nunit=13,atot,rowsperframe
    integer (kind=ik),allocatable ::&
    natoms(:),nmolsatoms(:),shift(:),moltypeofuatom(:)
    real(rk) :: director(1:3)=[0._rk,0._rk,1._rk],centerofmembrane(1:3)=0
    private :: atomindex_a,atomindex_b
    
    interface atomindex
        module procedure atomindex_a,atomindex_b
    end interface

!    type write_frame
!        integer(kind=ik) :: framenumber
!        character(kind=1,len=3) :: outformat
!    end type write_frame
!
!    type setflags
!        logical :: autofilename,cbl_switch,folding,apl
!        integer(kind=ik) :: distbin,ounit,wftot,aplgrid(2),leaflet !,writeframe
!        character(kind=1,len=255) :: filename,fileprefix,filesuffix
!        type(write_frame),allocatable :: writeframe(:)
!        character(kind=1,len=100),allocatable :: calc(:)
!        real(kind=rk) :: constant_bl
!    end type setflags
!    
!    type natom
!        character(kind=1,len=100) :: atomname,from_mol_prop,molecule
!    end type natom
!
!    type setcommon
!        integer(kind=ik),allocatable :: membrane_moltypes(:)
!    end type
!    type(setcommon) :: common_setflags
!
!    type(setflags) :: global_setflags
!
!    type calcval
!        real(kind=rk) :: mean,meandev,entropy,entropymutual,pearsoncoeff
!    end type calcval

!    type instruct
!        integer(kind=ik) :: atoms_bak(20),findex,nmolop,average_count
!        integer(kind=ik),allocatable :: atoms(:),apl_side(:),molind(:)
!        logical :: setapl
!        character(kind=1, len=50) :: instructionstring
!        real(kind=rk),allocatable :: datam(:,:)
!        type(setflags) :: set
!        type(calcval) :: cv
!        type(natom) :: newatom
!    end type instruct

    type atomdata
        character(kind=1,len=len(atomnames)) :: aname
        real(kind=rk) :: mass,mgratio
    end type atomdata

    type(atomdata),allocatable :: atomd(:)

!    type moltype
!        integer(kind=ik) :: firstatom,lastatom,nmol,natoms
!        integer(kind=ik),allocatable :: upper(:),lower(:)
!        character(kind=1,len=255) :: molname
!    end type moltype
!
!    type(moltype),allocatable :: molt(:)

    interface operator(.str.)
        module procedure stringconv
    end interface operator(.str.)

contains

    subroutine wf_gro(filename,wf,funit)!{{{
    character(kind=1,len=*) :: filename,wf
    integer(kind=ik) :: funit,i,j,k,l
        open(funit,file=trim(filename))
        write(funit,*)'frame=',trim(adjustl(wf))
        write(funit,*)atot
        l=1
        do i=1,size(molt)
            do j=1,molt(i)%nmol
                do k=molt(i)%firstatom,molt(i)%lastatom
                    !write(*,'(i5,a5,a10,i5,3F10.7)')j,molt(i)%molname,atomnames(k),l,getatom(k,j)
                    write(funit,'(i5,2a5,i5,3F8.3)')j,molt(i)%molname,atomnames(k),l,getatom(k,j)
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
        write(funit,*)atot
        write(funit,*)
        do i=1,size(molt)
            do j=1,molt(i)%nmol
                do k=molt(i)%firstatom,molt(i)%lastatom
                    write(funit,*)atomnames(k),10*getatom(k,j)
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
    global_setflags%autofilename=.TRUE.
    global_setflags%folding=.FALSE.
    global_setflags%filename=''
    global_setflags%distbin=100
    global_setflags%fileprefix='auto_'
    global_setflags%filesuffix='.out'
!    global_setflags%writeframe%framenumber=0
    global_setflags%wftot=0
    global_setflags%apl=.FALSE.
    global_setflags%aplgrid=[250,250]
    global_setflags%leaflet=0
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
        character(kind=1,len=5),allocatable :: moltype_atom(:,:)
        character(kind=1,len=200) :: sdr
        integer(kind=ik) :: ios,ia
       !stop 'HEJ'
        allocate(charvec(1))
        open(unit=tunit,file=trim(stringconv(fname)),iostat=ios)

        if(ios/=0)then
                write(*,*)"Error, cannot open ",&
                "file ",":",trim(stringconv(fname)),":"
                stop
        endif

        read(unit=tunit,fmt=*,iostat=ios)
        if(ios==endf)return
        read(unit=tunit,fmt=*,iostat=ios)atot ! Totala antalet atomer i en frame
        allocate(moltype_atom(2,atot),coor(1:3,atot),box(1:3))
        do ia=1,atot
            !read(tunit,"(5x,2A5)",iostat=ios)moltype_atom(1:2,ia)
            !write(*,*)size(charvec),'charvec'
            call readline(charvec,ios,tunit)!,*,iostat=ios)sdr
            !write(*,*)size(charvec),'charvec'
            !stop
            sdr=stringconv(charvec)
            read(sdr(6:16),"(2A5)")moltype_atom(1:2,ia)
            read(sdr(17:),*)coor(1:3,ia)
        end do
            read(tunit,*)box
        rewind(tunit)
        rowsperframe=atot+3
        call trajindex(moltype_atom)
        deallocate(charvec)   
    end subroutine initgro!}}}

    function readgro(tunit) result(ios)!{{{

        integer(kind=ik) :: tunit
        integer(kind=ik) :: ios,ia
        read(unit=tunit,fmt=*,iostat=ios)
        if(ios==endf)return
        read(unit=tunit,fmt=*,iostat=ios)atot
        do ia=1,atot
            read(unit=tunit,fmt='(20x,3f12.7)',iostat=ios)coor(:,ia)
        end do
        read(unit=tunit,fmt=*,iostat=ios)box(:)

    end function readgro!}}}

    function stringconv(vector) result(string)!{{{
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

    function cind(iuatom,imol) result(ind)!{{{

        integer(kind=ik) :: ind,imol,iuatom
        ! Index i coor för atom iuatom i molekyl imol
        ind=shift(iuatom)+natoms(iuatom)*imol
    end function cind!}}}

    function getatom(aindex,imol) result(acoor)!{{{

        integer (kind=ik) :: aindex,imol
        real (kind=rk) :: acoor(1:3)

        acoor(1:3)=coor(1:3,cind(aindex,imol))
        
    end function  getatom!}}}

    subroutine trajindex(moltype_atom)!{{{
        character(kind=1,len=11) :: uatom
   !     character(kind=1,len=11),allocatable :: temp(:)
        character(kind=1,len=5),allocatable ::&
        temp2(:),temp3(:),moltype_atom(:,:)
        integer(kind=ik),allocatable :: natomsoftype(:),nmols(:)
        character(kind=1, len=5),allocatable :: moltypenames(:)
        integer(kind=ik) :: ia,i,j,atoms=0       
        allocate(atomnames(atot),moltypenames(atot),temp(atot),temp2(atot),temp3(atot))
        atomnames="";moltypenames="";temp="";temp2="";temp3=""
        do ia=1,atot

            uatom=""
            uatom=trim(moltype_atom(1,ia))//"_"//trim(adjustl(moltype_atom(2,ia))) ! Sätter unikt atomnamn
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

    end subroutine trajindex!}}}

    function center_of_molecule(umol,imol) result(centerofmolecule)!{{{
        integer(kind=ik) :: umol,i,j,imol
        real(kind=rk) :: centerofmolecule(3)
        centerofmolecule=0
        j=umol
        do i=molt(j)%firstatom,molt(j)%lastatom
            centerofmolecule=centerofmolecule+getatom(i,imol)*masses(i)
        end do
        centerofmolecule=centerofmolecule/sum(masses(molt(j)%firstatom:molt(j)%lastatom))
       ! write(*,*)centerofmolecule,sum(masses(molt(j)%firstatom:molt(j)%lastatom))
    end function center_of_molecule!}}}
end module readtraj
