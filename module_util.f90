module util
    use kinds
!    use readtraj
!    use input
    implicit none
    integer(kind=ik) :: maxframes
    type moltype
        integer(kind=ik) :: firstatom,lastatom,nmol,natoms
        integer(kind=ik),allocatable :: upper(:),lower(:)
        character(kind=1,len=255) :: molname
    end type moltype

    type(moltype),allocatable :: molt(:)

    type write_frame
        integer(kind=ik) :: framenumber
        character(kind=1,len=3) :: outformat
    end type write_frame

    type setflags
        logical ::&
        autofilename,cbl_switch,folding,apl,whole,leaflets_defined,centerofmembrane,&
        molaverage
        integer(kind=ik) :: distbin,ounit,wftot,aplgrid(2),leaflet !,writeframe
        character(kind=1,len=255) :: filename,fileprefix,filesuffix,corrindex(2)
        type(write_frame),allocatable :: writeframe(:)
        character(kind=1,len=100),allocatable :: calc(:)
        real(kind=rk) :: constant_bl
    end type setflags
    
    type natom
        character(kind=1,len=100) :: atomname,from_mol_prop,molecule
    end type natom

    type setcommon
        logical :: silent,centerofmembrane
        !integer(kind=ik),allocatable :: membrane_moltypes(:)
    end type
    type(setcommon) :: common_setflags

    type(setflags) :: global_setflags

    type calcval
        real(kind=rk) :: mean,meandev,entropy,entropymutual,pearsoncoeff
    end type calcval

    type instruct
        integer(kind=ik) :: findex,nmolop,average_count,define!atoms_bak(20),
        integer(kind=ik),allocatable ::&
        atoms(:),apl_side(:),molind(:),membrane_moltypes(:)
        logical :: setapl
        character(kind=1, len=50) :: instructionstring,ref
        real(kind=rk),allocatable :: datam(:,:)
        type(setflags) :: set
        type(calcval) :: cv
        type(natom) :: newatom
    end type instruct

    interface reallocate
        module procedure &
        reallocatepointerchar,reallocateint,reallocatemoltype,reallocatewriteframe,reallocatereal
    end interface

    contains

    subroutine subtract(a,b,c)!{{{
        real(kind=rk) :: a(:,:),b(:,:),c(:,:)

        c(:,:)=a(:,:)-b(:,:)

    end subroutine subtract!}}}

    function normalize(vec) result(nvec)!{{{
    real(kind=rk),intent(in) :: vec(:)
    real(kind=rk) :: nvec(1:3)

    nvec=vec/sqrt(sum(vec**2))

    end function normalize!}}}

    elemental function mymodulo(a,b) result(c)!{{{
        real(kind=rk),intent(in) :: a,b
        real(kind=rk) :: c
        c=a-nint(a/b)*b !ok
    end function mymodulo!}}}

    function cross_product(v1,v2) result(v3)!{{{
    real(kind=rk) :: v1(:),v2(:),v3(3)

           v3 = [ v1(2)*v2(3)-v1(3)*v2(2),&
                  v1(3)*v2(1)-v1(1)*v2(3),&
                  v1(1)*v2(2)-v1(2)*v2(1) &
                ]

    end function cross_product!}}}

    subroutine argparse(i)!{{{

        character(kind=1,len=100) :: BUFFER,H
        integer(kind=4) :: i,n
        !n=iargc() ! Räknar antalet argument
        !write(*,*)n
        n=2
        call getarg(i,buffer)
        if (n<i)then
            write(*,*)"too many arguments!"
            stop
        else
            write(*,*)buffer
        end if

    end subroutine argparse!}}}

    subroutine reallocatepointerchar(vector,n)!{{{
        character(kind=1,len=1),allocatable,intent(inout) :: vector(:)
        character(kind=1,len=1),allocatable :: copy(:)
        integer(kind=ik),intent(in) :: n
        if (allocated(vector))then
            allocate(copy(1:min(size(vector),n)))
            copy=vector(1:size(copy))
            deallocate(vector)
            allocate(vector(1:n))
            vector(1:size(copy))=copy
            deallocate(copy)
        else
            allocate(vector(1:n))
        endif
    end subroutine reallocatepointerchar !}}}

    subroutine reallocatechar(vector,n)!{{{
        character(kind=1,len=*),allocatable,intent(inout) :: vector(:)
        character(kind=1,len=len(vector)),allocatable :: copy(:)
        integer(kind=ik),intent(in) :: n
        if (allocated(vector))then
            allocate(copy(1:n))
            copy=''
            copy(1:min(n,size(vector)))=vector(1:min(n,size(vector)))
            call move_alloc(copy,vector)
        else
            allocate(vector(1:n))
        endif
    end subroutine reallocatechar !}}}

    subroutine reallocatecharle(vector,n,le)!{{{
    integer(kind=ik) :: le    
    character(kind=1,len=le),allocatable,intent(inout) :: vector(:)
        character(kind=1,len=le),allocatable :: copy(:)
        integer(kind=ik),intent(in) :: n
        if (allocated(vector))then
            allocate(copy(1:n))
            copy=''
            copy(1:min(n,size(vector)))=vector(1:min(n,size(vector)))
            call move_alloc(copy,vector)
        else
            allocate(vector(1:n))
        endif
    end subroutine reallocatecharle !}}}

    elemental function intstr(i) result(str)!{{{
    integer(kind=ik),intent(in) :: i
    character(len=40) :: str
    write(str,*)i
    end function intstr!}}}

    subroutine reallocateint(vector,n)!{{{
        integer(kind=ik),allocatable,intent(inout) :: vector(:)
        integer(kind=ik),allocatable :: copy(:)
        integer(kind=ik),intent(in) :: n
        if (allocated(vector))then
            allocate(copy(1:n))
            copy=0
            copy(1:min(n,size(vector)))=vector(1:min(n,size(vector)))
            call move_alloc(copy,vector)
        else
            allocate(vector(1:n))
        endif
    end subroutine reallocateint !}}}

    subroutine reallocatereal(vector,n)!{{{
        real(kind=rk),allocatable,intent(inout) :: vector(:)
        real(kind=rk),allocatable :: copy(:)
        integer(kind=ik),intent(in) :: n
        if (allocated(vector))then
            allocate(copy(1:n))
            copy=0
            copy(1:min(n,size(vector)))=vector(1:min(n,size(vector)))
            call move_alloc(copy,vector)
        else
            allocate(vector(1:n))
        endif
    end subroutine reallocatereal !}}}

    subroutine reallocatemoltype(v,i)!{{{
        type(moltype),intent(inout),allocatable :: v(:)
        type(moltype),allocatable ::copy(:)
        integer(kind=ik) :: i,j
        if (allocated(v))then
            j=min(i,size(v))
            allocate(copy(i))
            copy(1:j)=v(1:j)
            call move_alloc(copy,v)
        else
            allocate(v(i))
        
        end if
    end subroutine reallocatemoltype!}}}

    subroutine reallocatewriteframe(v,i)!{{{
        type(write_frame),intent(inout),allocatable :: v(:)
        type(write_frame),allocatable ::copy(:)
        integer(kind=ik) :: i,j
        if (allocated(v))then
            j=min(i,size(v))
            allocate(copy(i))
            copy(1:j)=v(1:j)
            call move_alloc(copy,v)
        else
            allocate(v(i))
        
        end if
    end subroutine reallocatewriteframe!}}}

    subroutine reallocinstratoms(v,i)!{{{
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
    end subroutine reallocinstratoms!}}}

function strvecindex(refvec,teststr) result(j)
character (len=*) :: refvec(:),teststr
integer (kind=ik) :: i,j
!if(ANY(refvec==trim(teststr)))then
    j=0
    do i=1,size(refvec)
        if(trim(refvec(i))==trim(teststr))j=i
    end do
!else
!    i=0
!end if
end function strvecindex

    elemental function readint(str) result(i)
        integer :: i
        character(len=*),intent(in) :: str
        read(str,*)i
    end function readint

end module util
