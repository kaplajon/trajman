module util
    use kinds
    use readtraj
!    use input
    implicit none
    integer(kind=ik) :: maxframes
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
