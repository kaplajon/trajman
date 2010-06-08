module kinds
   ! use ifport
    implicit none
    integer(kind=4) ,parameter::ik=4
    integer(kind=ik),parameter :: rk=8,endf=-1,endr=-2,stdin=5,stdout=6
    real(kind=rk),parameter :: planck=6.62606896e-34 !Js (33)
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
