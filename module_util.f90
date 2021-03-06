!---LICENSE-------------------------------------------------------!{{{
! This file is part of
!
!  Trajman: A MD Trajectory Analysis Tool
!
!  Copyright (C) 2010-2014, Jon Kapla, jon.kapla@mmk.su.se
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
module util
    use kinds
!    use readtraj
!    use input
    implicit none
    integer(kind=ik) :: maxframes,skipframes,atomsdefined
    character(kind=1, len=11),allocatable :: atomnames(:)
    character(kind=1, len=3),allocatable :: groupres(:)
    type moltype!{{{
        integer(kind=ik) :: firstatom,lastatom,nmol,natoms
        integer(kind=ik),allocatable :: upper(:),lower(:)
        character(kind=1,len=255) :: molname
    end type moltype!}}}
    type(moltype),allocatable :: molt(:)
    type trjf
        character(kind=1,len=1),allocatable :: filename(:)
    end type trjf
    type(trjf),allocatable :: trajfile(:)
    type uatoms!{{{
        character(kind=1,len=255) :: aname
        character(kind=1,len=255) :: mname
        integer(kind=ik) :: moltype
        real (kind=rk),allocatable :: coor(:,:)
    end type!}}}
    type(uatoms),allocatable :: atom(:)
    type write_frame!{{{
        integer(kind=ik) :: framenumber
        character(kind=1,len=3) :: outformat
    end type write_frame!}}}
    type constant!{{{
        real(kind=rk) :: value
        logical :: switch
    end type!}}}
    type distmima!{{{
        real(kind=rk) :: mi,ma
        logical :: switch
    end type!}}}
    type scaletype!{{{
        !real(kind=rk) ::
        character(kind=1,len=3) :: typ
        real(kind=rk) :: val
        logical :: switch
    end type!}}}
    type slicetype!{{{
        !real(kind=rk) ::
        character(kind=1,len=4) :: typ
        real(kind=rk) :: upper,lower,upper2,lower2
        logical :: switch,switch2
        logical,allocatable :: bintest(:,:)
    end type!}}}
    type setflags!{{{
        logical ::&
        autofilename,cbl_switch,folding,apl,gd,whole,leaflets_defined,centerofmembrane,&
        molaverage,instructionsum,xyrdf,zrdf,VSnorm,foldcenterofmembrane,karplus,coorsys
        integer(kind=ik) :: &
        distbin,ounit,wftot,aplgrid(2),leaflet,tshift,karplus_fnc,coorsys_type,coorsys_helpers(3)
        character(kind=1,len=255) :: filename,fileprefix,filesuffix,corrindex(6)
        character(kind=1,len=255),allocatable :: averagetags(:),dealloctags(:)
        type(write_frame),allocatable :: writeframe(:)
        type(constant) :: const
        type(distmima) :: distminmax
        type(scaletype) :: scaling
        type(slicetype) :: slice
        character(kind=1,len=100),allocatable :: calc(:)
        integer(kind=ik),allocatable :: &
        sph1(:),sph2(:),spt1(:),spt2(:)
        real(kind=rk) :: constant_bl,ch_bondlength,rdf_binsize,karplus_params(5)
    end type setflags!}}}
    type natom!{{{
        character(kind=1,len=100) :: atomname,atype,molecule
        real(kind=rk) :: angle,bl
        integer(kind=ik),allocatable :: helpers(:)
    end type natom!}}}
    type setcommon!{{{
        logical :: silent,centerofmembrane
        real(kind=rk) :: traj_cscale
        integer(kind=ik),allocatable :: membrane_moltypes(:),shuffle_atoms(:)
        !integer(kind=ik),allocatable :: membrane_moltypes(:)
    end type!}}}
    
    type(setcommon) :: common_setflags
    type(setflags) :: global_setflags

    type atomdata
        character(kind=1,len=len(atomnames)) :: aname
        real(kind=rk) :: mass,mgratio
    end type atomdata
    type(atomdata),allocatable :: atomd(:)
    type calcval!{{{
        real(kind=rk) :: mean,meandev,entropy,entropymutual,pearsoncoeff
        integer(kind=ik) :: n
    end type calcval!}}}
    type instruct!{{{
        integer(kind=ik) :: findex,nmolop,average_count,define,dpctype!atoms_bak(20),
        integer(kind=ik),allocatable ::&
        atoms(:),apl_side(:),molind(:)!,membrane_moltypes(:)
        logical :: setapl
        character(kind=1, len=50) :: instructionstring,ref
        real(kind=rk),allocatable :: datam(:,:),rdf_dist(:),rdf_pairs(:)
        real(kind=rk) :: rdf_bin
        type(setflags) :: set
        type(calcval) :: cv
        type(natom) :: newatom
    end type instruct!}}}
    interface reallocate
        module procedure &
        reallocatepointerchar,reallocateint,reallocatemoltype,&
        reallocatewriteframe,reallocatereal,reallocatetrajfile,&
        reallocateatomd
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

    elemental function realstr(i) result(str)!{{{
    real(kind=rk),intent(in) :: i
    character(len=40) :: str
    write(str,*)i
    end function realstr!}}}

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

    subroutine reallocatecoor(mtrx,n)!{{{
        real(kind=rk),allocatable,intent(inout) :: mtrx(:,:)
        real(kind=rk),allocatable :: copy(:,:)
        integer(kind=ik),intent(in) :: n
        if (allocated(mtrx))then
            allocate(copy(1:3,1:n))
            copy=0
            copy(1:3,1:min(n,size(mtrx,2)))=mtrx(1:3,1:min(n,size(mtrx,2)))
            call move_alloc(copy,mtrx)
        else
            allocate(mtrx(1:3,1:n))
        endif
    end subroutine reallocatecoor !}}}

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


    subroutine reallocateatomd(v,i)!{{{
        type(atomdata),intent(inout),allocatable :: v(:)
        type(atomdata),allocatable ::copy(:)
        integer(kind=ik) :: i,j
        if (allocated(v))then
            j=min(i,size(v))
            allocate(copy(i))
            copy(1:j)=v(1:j)
            call move_alloc(copy,v)
        else
            allocate(v(i))
        
        end if
    end subroutine reallocateatomd!}}}

    subroutine reallocateatom(v,i)!{{{
        type(uatoms),intent(inout),allocatable :: v(:)
        type(uatoms),allocatable ::copy(:)
        integer(kind=ik) :: i,j
        if (allocated(v))then
            j=min(i,size(v))
            allocate(copy(i))
            copy(1:j)=v(1:j)
            call move_alloc(copy,v)
        else
            allocate(v(i))
        
        end if
    end subroutine reallocateatom!}}}

    subroutine reallocatetrajfile(v,i)!{{{
        type(trjf),intent(inout),allocatable :: v(:)
        type(trjf),allocatable ::copy(:)
        integer(kind=ik) :: i,j
        if (allocated(v))then
            j=min(i,size(v))
            allocate(copy(i))
            copy(1:j)=v(1:j)
            call move_alloc(copy,v)
        else
            allocate(v(i))
        
        end if
    end subroutine reallocatetrajfile!}}}

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

function strvecindex(refvec,teststr) result(j)!{{{
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
end function strvecindex!}}}

    elemental function readint(str) result(i)!{{{
        integer :: i
        character(len=*),intent(in) :: str
        read(str,*)i
    end function readint!}}}

    elemental function readreal(str) result(r)!{{{
        real :: r
        character(len=*),intent(in) :: str
        read(str,*)r
    end function readreal!}}}

function r8vec_norm ( n, a )!{{{

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer ( kind = ik ) n

  real    ( kind = rk ) a(n)
  real    ( kind = rk ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end function r8vec_norm!}}}

subroutine vector_rotate_3d ( v1, axis, angle, v2 )!{{{

!*****************************************************************************80
!
!! VECTOR_ROTATE_3D rotates a vector around an axis vector in 3D.
!
!  Discussion:
!
!    Thanks to Cody Farnell for correcting some errors in a previous
!    version of this routine!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), the vector to be rotated.
!
!    Input, real ( kind = 8 ) AXIS(3), the vector about which the
!    rotation is to be carried out.
!
!    Input, real ( kind = 8 ) ANGLE, the angle, in radians, of the rotation
!    to be carried out.
!
!    Output, real ( kind = 8 ) V2(3), the rotated vector.
!
  implicit none

  real    ( kind = rk ) angle
  real    ( kind = rk ) axis(3)
  real    ( kind = rk ) dot
  real    ( kind = rk ) norm
  real    ( kind = rk ) norm_vn
  real    ( kind = rk ) normal2(3)
!  real    ( kind = 8 ) r8vec_norm
  real    ( kind = rk ) v1(3)
  real    ( kind = rk ) v2(3)
  real    ( kind = rk ) vn(3)
  real    ( kind = rk ) vp(3)
  real    ( kind = rk ) vr(3)
!
!  Compute the length of the rotation axis.
!
  norm = r8vec_norm ( 3, axis )

  if ( norm == 0.0_rk ) then
    v2(1:3) = v1(1:3)
    return
  end if
!
!  Compute the dot product of the vector and the (unit) rotation axis.
!
  dot = dot_product ( v1(1:3), axis(1:3) ) / norm
!
!  Compute the parallel component of the vector.
!
  vp(1:3) = dot * axis(1:3) / norm
!
!  Compute the normal component of the vector.
!
  vn(1:3) = v1(1:3) - vp(1:3)

  norm_vn = r8vec_norm ( 3, vn )

  if ( norm == 0.0_rk ) then
    v2(1:3) = vp(1:3)
    return
  end if

  vn(1:3) = vn(1:3) / norm_vn
!
!  Compute a second vector, lying in the plane, perpendicular
!  to V1 and VN, and forming a right-handed system.
!
  normal2(1) = axis(2) * vn(3) - axis(3) * vn(2)
  normal2(2) = axis(3) * vn(1) - axis(1) * vn(3)
  normal2(3) = axis(1) * vn(2) - axis(2) * vn(1)

  norm = r8vec_norm ( 3, normal2 )
  normal2(1:3) = normal2(1:3) / norm
!
!  Rotate the normal component by the angle.
!
!  v2(1:3) = ( cos ( angle ) * v1(1:3) + sin ( angle ) * cross_product(v1,axis)!<Up>(normal2(1:3) )
!
!  The rotated vector is the parallel component plus the rotated component.
!
 

  return
end subroutine vector_rotate_3d !}}}

function v2q(v,theta) result(q)
!Create a quaternion from the vector v and angle theta
real(kind=rk),intent(in) :: v(:),theta ! Rad
real(kind=rk) :: q(4)
q=[cos(theta/2),sin(theta/2)*normalize(v)]
end function v2q

function RV(q) result(m)
!Rotational matrix: B.Stevensson et. al. 2011
real(kind=rk),intent(in) :: q(:) !quat
real(kind=rk) :: m(3,3) !quat
!First col
m(1,1)=q(1)**2+q(2)**2-q(3)**2-q(4)**2
m(2,1)=2*(q(2)*q(3)+q(1)*q(4))
m(3,1)=2*(q(2)*q(4)-q(1)*q(3))
!Second col
m(1,2)=2*(q(2)*q(3)-q(1)*q(4))
m(2,2)=q(1)**2-q(2)**2+q(3)**2-q(4)**2
m(3,2)=2*(q(3)*q(4)+q(1)*q(2))
!Third col
m(1,3)=2*(q(2)*q(4)+q(1)*q(3))
m(2,3)=2*(q(3)*q(4)-q(1)*q(2))
m(3,3)=q(1)**2-q(2)**2-q(3)**2+q(4)**2
end function RV

end module util
