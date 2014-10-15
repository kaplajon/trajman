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
module kinds
   ! use ifport
    implicit none
    integer(kind=4) ,parameter::ik=4
    integer(kind=ik),parameter :: rk=8,ikr=ik,endf=-1,endr=-2,stdin=5,stdout=6
    real(kind=rk),parameter :: planck=6.62606957e-34 !(29) CODATA2010
 !6.62606896e-34 !Js (33)
    real(kind=rk),protected :: pi,magnetic,hbar!=log(real(cmplx(-1.,0.,rk)**cmplx(0.,-1.,rk)))
!complex(kind=rk) :: pi
    contains
    subroutine constants
        pi=log(real(cmplx(-1.,0.,rk)**cmplx(0.,-1.,rk)))
        magnetic=4*pi*1.e-7 !NA^-2
        hbar=planck/(2*pi)
        !pi=cmplx(-1._rk,0._rk)**cmplx(0.,-1.,rk)
    end subroutine constants
end module kinds
