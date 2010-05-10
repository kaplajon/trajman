module util
    use kinds
    use readtraj
!    use input
    implicit none
    integer(kind=ik) :: maxframes

    interface reallocate
        module procedure &
        reallocatepointerchar,reallocateint,reallocatemoltype,reallocatewriteframe
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

!    subroutine APL(grid,apl_atoms)
!        integer(kind=ik) :: grid(200,200),apl_atoms(:)
!
!        do i=1,size(grid,1)
!            do j=1,size(grid,2)
!                do k=1,size(apl_atoms)
!                    
!                    rsq= (modulo( getatomx - (box(1)/(size(grid,1)))*i,box(1)/2 ))**2&
!                    + (modulo( getatomy - (box(2)/size(grid,2))*j,box(2)/2 ))**2
!                    if(rsq<rmin)then
!                        minind=k
!                        rmin=rsq
!                    end if
!                end do
!                grid(i,j)=minind
!            end do
!        end do
!
!    end subroutine APL

end module util
