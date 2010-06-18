!-----------------------------------------------------------------
! This file is part of
!
!  Trajman: A MD Trajectory Analysis Tool
!
!  Copyright (C) 2010, Jon Kapla, jon.kapla@mmk.su.se
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
!-----------------------------------------------------------------

module trajop
    use kinds
    use readtraj
    use apl
!    use input
    use util
    implicit none
    
    contains

    function dirangle(a,b,imol) result(teta)!{{{
        ! Angle: a-b  against director in degrees
        real(kind=rk) :: teta,dir(3)
        integer(kind=ik) :: a,b,imol
        dir=director*sign(1._rk,dot_product(center_of_molecule(moltypeofuatom(a),imol)-centerofmembrane,director))
        teta=acos(dot_product(normalize(getatom(b,imol)-getatom(a,imol)),dir))*180._rk/pi
    end function dirangle!}}}

    function valenceangle(a,b,c,imol) result(teta)!{{{
        ! Angle a-b against a-c in degrees
        real(kind=rk) :: teta
        integer(kind=ik) :: a,b,c,imol
        teta=acos(dot_product(normalize(getatom(a,imol)-getatom(b,imol))&
                ,normalize(getatom(a,imol)-getatom(c,imol))))*180._rk/pi
    end function valenceangle!}}}

    function torsionangle(a,b,c,d,imol) result(teta)!{{{
        ! Torsion angle: a-b-c-d in degrees 
        real(kind=rk) :: Vba(1:3),Vcd(1:3),Vcb(1:3),teta,teta2
        integer(kind=ik) :: a,b,c,d,imol
        Vcb=normalize(getatom(b,imol) - getatom(c,imol))
        Vba = normalize(cross_product( &
                        cross_product(getatom(a,imol) - getatom(b,imol) ,Vcb) &
                        , Vcb))
        Vcd = normalize(cross_product( &
                        cross_product(getatom(d,imol) - getatom(c,imol) ,Vcb) &
                        ,Vcb))
        !teta = modulo(atan2(dot_product(Vcb,cross_product(Vcd,Vba)),dot_product(Vba,Vcd)) * 180._rk / pi,360._rk)
        teta = atan2(dot_product(Vcb,cross_product(Vcd,Vba)),dot_product(Vba,Vcd)) * 180._rk / pi
    end function torsionangle!}}}

    function bond_length(a,b,imol) result(c)!{{{
        integer(kind=ik) :: a,b,imol
        real(kind=rk) :: c,vector(1:3)
        vector=getatom(b,imol)-getatom(a,imol)
        c = sqrt(sum(vector**2)) 
    end function bond_length!}}}

    function distance_com(a,imol) result(distance)!{{{
        integer(kind=ik) :: a,imol
        real(kind=rk) :: dir(3),distance
        dir=director*sign(1._rk,dot_product(center_of_molecule(moltypeofuatom(a),imol)-centerofmembrane,director))
        distance=dot_product(getatom(a,imol)-centerofmembrane,dir)
    end function distance_com!}}}

    function order_parameter(a,b,imol) result(c)!{{{
        integer(kind=ik) :: a,b,imol
        real(kind=rk) :: c
        c=1.5_rk*cos(pi/180._rk*&
        dirangle(a,b,imol))**2-0.5
    end function order_parameter!}}}

    subroutine procop(instr,frame)!{{{
    ! Hanterar operationerna hämtade från input
        integer(kind=ik) :: imol,jmol,i,j,frame
        real(kind=rk) :: teta,bl
        type(instruct) :: instr(:)
          do i=1,size(instr)
            select case(instr(i)%findex)
                case(0)
                    !if(instr(i)%setapl)call apl_grid(instr(i))
                case(1) !DIRECTOR ANGLE
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=&
                        cos(pi/180._rk*dirangle(instr(i)%atoms(1),instr(i)%atoms(2),imol))
                    end do
                case(2) !VALENCE ANGLE
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=valenceangle(instr(i)%atoms(1),instr(i)%atoms(2)&
                        ,instr(i)%atoms(3),imol)
                    end do
                case(3) !TORSION ANGLE
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=torsionangle(instr(i)%atoms(1),instr(i)%atoms(2)&
                        ,instr(i)%atoms(3),instr(i)%atoms(4),imol)
                        
                    end do
                case(4) !BOND LENGTH
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol)
                    end do
                case(5) !ORDER PARAMETER
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        !instr(i)%datam(imol,frame)=1.5_rk*cos(pi/180._rk*&
                        !dirangle(instr( Ti)%atoms(1),instr(i)%atoms(2),imol))**2-0.5
                        instr(i)%datam(jmol,frame)=order_parameter(instr(i)%atoms(1),instr(i)%atoms(2),imol)
                    end do
                case(6) ! DISTANCE CENTER OF MEMBRANE
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=distance_com(instr(i)%atoms(1),imol)
                    end do
                case(7) ! CORRELATE
                        ! Allt sköts i postproc
                case(8) ! DIPOLE COUPLING
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        if(instr(i)%set%cbl_switch)then
                            bl=instr(i)%set%constant_bl
                        else
                            bl=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol)
                        end if
                        instr(i)%datam(jmol,frame) =-(magnetic/(4*pi))&
                        *(mgratios(instr(i)%atoms(1))*mgratios(instr(i)%atoms(2))*hbar/(2*pi))&
                        *order_parameter(instr(i)%atoms(1),instr(i)%atoms(2),imol)&
                        *(bl*1.e-9)**(-3)/1000.!(1e-9 for meters and 1000 for kHz) 
                    end do
                 case(9) ! AVERAGE
                         ! Allt sköts i postproc
                 case(10) !DEFINE
                     select case(instr(i)%define)
                        case(1) !CENTEROFMEMBRANE
                            call center_of_membrane(instr(i)%membrane_moltypes)
                        case(2) !ATOM
                            select case(instr(i)%newatom%from_mol_prop)
                                case('com','center_of_mass')
                                ! Define atom from center of mass of a defined submolecule.
                                    do imol=1,molt(moltypeofuatom(instr(i)%atoms(1)))%nmol
                                        coor(:,cind(instr(i)%atoms(1),imol))=&
                                        center_of_molecule(moltypeofuatom(instr(i)%atoms(1)),imol)
                                    end do
                            end select
                         case(3) !LEAFLET
                                    if(global_setflags%apl)call apl_grid(instr(i))
                       end select
                 case(11) ! AREA PER LIPID
                     call apl_calc(instr(i),frame)
                     if(allocated(global_setflags%writeframe))then
                     if(ANY(instr(i)%set%writeframe(:)%framenumber==frame)&
                     .AND.instr(i)%set%writeframe(strvecindex(instr(i)%set%writeframe(:)%outformat,'apl'))%framenumber&
                     ==frame)then
                        call apl_matrix_out(frame)
                    end if
                    end if
            end select
        end do 
    end subroutine procop!}}}

    subroutine mean_var(vec,mean,var)!{{{
        real(kind=rk) :: vec(:),mean,var
        mean=sum(vec)/real(size(vec),rk)
        var=sum((vec(:)-mean)**2)/real(size(vec)-1,rk)
    end subroutine mean_var!}}}

    function std(meant) result(meandev)!{{{
        real(kind=rk) :: meant(:),mean,var,meandev
        call mean_var(meant,mean,var)
        !mean=sum(meant)/real(size(meant),rk)
        !var=sum((meant(:)-mean)**2)/real(size(meant)-1,rk)
        meandev=sqrt(var/real(size(meant),rk))
    end function std!}}}

    function average(datam) result(string2)!{{{
        implicit none
        character(len=24) :: string
        character(len=300) :: string2
        integer(kind=ik) :: imol,i
        real(kind=rk) :: datam(:,:),meant(1:size(datam,1)),mean,meandev,var
        do imol=1,size(datam,1)
            meant(imol)=sum(datam(imol,:))/real(size(datam,2),rk)
        end do
        mean=sum(meant)/real(size(meant),rk)
        var=sum((meant(:)-mean)**2)/real(size(meant)-1,rk)
        meandev=sqrt(var/real(size(meant),rk))
        string=getmeanwithdev(mean,sqrt(var))
        write(string2,*)mean,meandev,sqrt(var)
        !write(*,*)string2
        string2=trim(adjustl(string))//" "//trim(string2)
    end function average!}}}

    function test_var(vector) result(var)!{{{
        real(kind=rk) :: vector(:),vector2(size(vector)),var,mean
        integer(kind=ik) :: imol
        mean=sum(vector)
        do imol=1,size(vector)
            vector2(imol)=(mean-vector(imol))/(real(size(vector),rk)-1._rk)
        end do
        mean=sum(vector(:))/real(size(vector),rk)
        var=sum((vector2(:)-mean)**2)*real(size(vector)-1,rk)/real(size(vector),rk)
    end function test_var!}}}

    subroutine instruction_average(instr,ind,lastind)!{{{
        type(instruct),intent(inout) :: instr(:)
        integer(kind=ik) :: ounit,ind,j,k,lastind,mols
        j=ind;k=0;mols=0
        if(instr(ind)%set%molaverage)then
        do 
            j=j+1
            if(instr(j)%findex/=0 .AND. instr(j)%findex/=10)then
                k=k+1
                mols=mols+instr(j)%nmolop
                if(k>=instr(ind)%average_count)exit
            endif
        end do
        allocate(instr(ind)%datam(mols,size(instr(j)%datam,2)))
        else
        allocate(instr(ind)%datam(size(instr(lastind)%datam,1),size(instr(lastind)%datam,2)))!*instr(ind)%average_count))
        end if
        mols=1;j=ind;k=0
        instr(ind)%datam=0
        do !j=i+1,i+instr(i)%average_count
            j=j+1
            if(instr(j)%findex/=0 .AND. instr(j)%findex/=10)then
                k=k+1
                if(instr(ind)%set%molaverage)then
                    mols=mols+instr(j)%nmolop
                    instr(ind)%datam(mols-instr(j)%nmolop:mols,:)=instr(j)%datam
                    write(*,*)mols,instr(j)%nmolop,j
                else
                    instr(ind)%datam=(instr(ind)%datam+instr(j)%datam)!/real(instr(ind)%average_count,rk)!(:,size(instr(j)%datam,2)*(k-1)+1:size(instr(j)%datam,2)*k)=instr(j)%datam
                end if
                if(k>=instr(ind)%average_count)exit
            endif
        end do
        if(.NOT.instr(ind)%set%molaverage)instr(ind)%datam=instr(ind)%datam/instr(ind)%average_count
       ! instr(ind)%datam=instr(ind)%datam/real(instr(ind)%average_count,rk)
    end subroutine instruction_average!}}}

    function getmeanwithdev(mean,meandev) result(str)!{{{
        implicit none ! write mean with meandeviation in a compact form.
        real(kind=rk), intent(in) :: mean,meandev
        integer(kind=ik), parameter :: lenlength=20
        character(kind=1,len=lenlength+4) :: str
        character(kind=1,len=lenlength) :: str1
        character(kind=1,len=8) :: ffmtstr,esfmtstr
        integer(kind=ik) :: i,n,p,decdiff
        do i=1,len(str)
           str(i:i)=' ' !blank_char(1:1)
        end do
        n=len(str1)
        do i=1,n
           str1(i:i)=' '!blank_char(1:1)
        end do
        esfmtstr='(ES00.1)'
        write(unit=esfmtstr(4:5),fmt='(I2)')lenlength
        write(unit=str1(1:n),fmt=esfmtstr)meandev
        p=0
        do i=1,n
           if(str1(i:i)=='.')then
              p=i
              exit
           end if
        end do
        decdiff=ceiling(-log(meandev)/log(10._rk),ik)+1_ik
        select case(decdiff)
        case(:-1_ik)
           str(1:11)='unit is NaN'
           return
        case(0_ik:9_ik)
           ffmtstr='(F00.0) '
           write(unit=ffmtstr(6:6),fmt='(I1)')decdiff
        case(10_ik:99_ik)
           ffmtstr='(F00.00)'
           write(unit=ffmtstr(6:7),fmt='(I2)')decdiff
        end select
        write(unit=ffmtstr(3:4),fmt='(I2)')lenlength
        write(unit=str(1:n),fmt=ffmtstr)mean
        str(n+1:n+1)='('
        str(n+2:n+2)=str1(p-1:p-1)
        str(n+3:n+3)=str1(p+1:p+1)
        str(n+4:n+4)=')'
      end function getmeanwithdev!}}}

    subroutine distrib(instr,ind,calcind)!{{{
        type(instruct),intent(inout) :: instr
        integer(kind=ik) :: bi,i,j,ind,calcind,indold=0
        !integer(kind=ik),optional :: ounit
        real(kind=rk) :: vec(size(instr%datam))
        real(kind=rk) :: mi,ma,x,bin,dp,dpm,isoentropy
        real(kind=rk),allocatable :: dist(:),distmol(:,:),molentropy(:)
        save indold,dist,distmol,mi,bin

        if(ind>indold)then
            !if(allocated(dist))deallocate(dist);
            allocate(dist(instr%set%distbin))
            !if(allocated(distmol)deallocate(distmol);
            allocate(distmol(size(instr%datam,1),size(dist)))
            !if(allocated(molentropy)deallocate(molentropy);
            allocate(molentropy(size(distmol,1)))
            vec=reshape(instr%datam,[size(instr%datam)])
            dist=0
            distmol=0
            select case(instr%findex)
            case(1)
                mi=-1._rk
                ma=1._rk
            case(2)
                mi=0._rk
                ma=180._rk
            case(3)
               ! mi=0._rk
               ! ma=360._rk
                mi=-180._rk
                ma=180._rk
            case default
                mi=minval(vec)
                ma=maxval(vec)
            end select
            bin=(ma-mi)/real(size(dist),rk)
            dp=1._rk/(real(size(vec),rk)*bin)
            dpm=1._rk/(real(size(instr%datam,2),rk)*bin)
            do i=1,size(vec)
                bi=int((vec(i)-mi)/bin+1._rk)
                if(bi==size(dist)+1)bi=size(dist)
                if(bi<=0.or.bi>size(dist))then
                    write(*,*)bi,"bi",vec(i),"vec",i,"i",ma,"ma",mi,"mi",dp,"dp",bin,"bin"
                    stop "subroutine distrib"   
                else
                    dist(bi)=dist(bi)+dp
                    j=mod(i-1,size(distmol,1))+1
                    distmol(j,bi)=distmol(j,bi)+dpm
                end if
            end do
           !Räkna ut entropin för distributionen
           isoentropy=log(ma-mi)
           instr%cv%entropy=-log(real(size(dist),rk)*bin)-bin*sum(dist(:)*log(dist(:)),MASK=dist(:)/=0)&
                         -isoentropy
           !Standardavvikelse för entropin
           do i=1,size(distmol,1)
           molentropy(i)=-log(real(size(distmol,2),rk)*bin)&
           -bin*sum(distmol(i,:)*log(distmol(i,:)),MASK=distmol(i,:)/=0)-isoentropy
           end do
           instr%cv%entropymutual=instr%cv%entropy-sum(molentropy)/real(size(distmol,1),rk)
           deallocate(molentropy)
        end if !ind>indold
        select case(instr%set%calc(calcind))
        case('distrib')
            !Skriv ut distribution till disk
            do bi=1,size(dist)
                x=mi+(real(bi,rk)-1._rk)*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,dist(bi),std(distmol(:,bi))
                x=mi+(real(bi,rk))*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,dist(bi),std(distmol(:,bi))
            end do
        case('distribm')
            do i=1,size(distmol,1)
             do bi=1,size(dist)
                x=mi+(real(bi,rk)-1._rk)*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,i,distmol(i,bi)!std(distmol(:,bi))
                x=mi+(real(bi,rk))*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,i,distmol(i,bi)!std(distmol(:,bi))
            end do
                write(instr%set%ounit,*)
            end do

        end select
        if(COUNT(instr%set%calc(:)(1:7)=='distrib')==1.OR.ind==indold)deallocate(dist,distmol)
        indold=ind
    end subroutine distrib!}}}

    subroutine corr_distrib(instr0,instr1,instr2)!{{{
        type(instruct),intent(inout) :: instr0,instr1,instr2
        integer(kind=ik) :: bi1,bi2,i,j,imol,ios
        character(len=len(instr0%set%filename)) :: filename
        !integer(kind=ik),optional :: ounit
        real(kind=rk) :: vec1(size(instr1%datam)),vec2(size(instr2%datam))
        real(kind=rk) :: mi1,mi2,ma1,ma2,x1,x2,bin1,bin2,dp,dpm,isoentropy&
        ,m1,m2,v1,v2
        real(kind=rk),allocatable :: dist(:,:),distmol(:,:,:),molentropy(:)
        allocate(dist(instr1%set%distbin,instr2%set%distbin))
        allocate(distmol(size(instr1%datam,1),size(dist,1),size(dist,2)))
        allocate(molentropy(size(distmol,1)))
        vec1=reshape(instr1%datam,[size(instr1%datam)])
        vec2=reshape(instr2%datam,[size(instr2%datam)])

        dist=0
        distmol=0
        select case(instr1%findex)
            case(1)
                mi1=-1._rk
                ma1=1._rk
            case(2)
                mi1=0._rk
                ma1=180._rk
            case(3)
               ! mi1=0._rk
               ! ma1=360._rk
                mi1=-180._rk
                ma1=180._rk
            case default
                mi1=minval(vec1)
                ma1=maxval(vec1)
        end select

        select case(instr2%findex)
            case(1)
                mi2=-1._rk
                ma2=1._rk
            case(2)
                mi2=0._rk
                ma2=180._rk
            case(3)
               ! mi2=0._rk
               ! ma2=360._rk
                mi2=-180._rk
                ma2=180._rk
            case default
                mi2=minval(vec2)
                ma2=maxval(vec2)
        end select

        bin1=max((ma1-mi1)/real(size(dist,1),rk),epsilon(bin1))
        bin2=(ma2-mi2)/real(size(dist,2),rk)

        dp=1._rk/(real(size(vec1),rk)*bin1*bin2)
        dpm=1._rk/(real(size(instr1%datam,2),rk)*bin1*bin2)
        do i=1,size(vec1)
            bi1=int((vec1(i)-mi1)/bin1+1._rk)
            bi2=int((vec2(i)-mi2)/bin2+1._rk)

            if(bi1==size(dist,1)+1)bi1=size(dist,1)
            if(bi2==size(dist,2)+1)bi2=size(dist,2)
            if(bi1<=0.or.bi1>size(dist,1).OR.bi2<=0.or.bi2>size(dist,2))then
                write(*,*)bi1,"bi1",vec1(i),"vec1",i,"i",ma1,"ma1",mi1,"mi1",dp,"dp",bin1,"bin1"
                write(*,*)bi2,"bi2",vec2(i),"vec2",i,"i",ma2,"ma2",mi2,"mi2",dp,"dp",bin2,"bin2"
            stop "subroutine corr_distrib"   
            else
                    dist(bi1,bi2)=dist(bi1,bi2)+dp

                    j=mod(i-1,size(distmol,1))+1
                    distmol(j,bi1,bi2)=distmol(j,bi1,bi2)+dpm
                end if
        end do
        if(COUNT(instr0%set%calc(:)(1:7)=='distrib')>0)then
        do i=1,size(instr0%set%calc)
        select case(instr0%set%calc(i))
        case('distrib')
            filename=trim(instr0%set%filename)//'_'//trim(instr0%set%calc(i))//trim(instr0%set%filesuffix)
        if(instr0%set%ounit/=stdout)open(unit=instr0%set%ounit,file=trim(filename),status="unknown",iostat=ios)
            !Skriv ut distribution till disk
            do bi1=1,size(dist,1)
                do bi2=1,size(dist,2)
                    x1=mi1+(real(bi1,rk)-.5_rk)*bin1
                    x2=mi2+(real(bi2,rk)-.5_rk)*bin2
                    if(instr1%findex==1)x1=acos(x1)*180._rk/pi
                    if(instr2%findex==1)x2=acos(x2)*180._rk/pi
                    write(instr0%set%ounit,*)x1,x2,dist(bi1,bi2),std(distmol(:,bi1,bi2))
                    !x=mi+(real(bi,rk))*bin
                    !if(instr%findex==1)x=acos(x)*180._rk/pi
                    !write(ounit,*)x,dist(bi),std(distmol(:,bi))
                end do
                write(instr0%set%ounit,*)
            end do
            if(instr0%set%ounit/=stdout)close(instr0%set%ounit)
        case('distribm')
            filename=trim(instr0%set%filename)//'_'//trim(instr0%set%calc(i))//trim(instr0%set%filesuffix)
            if(instr0%set%ounit/=stdout)open(unit=instr0%set%ounit,file=trim(filename),status="unknown",iostat=ios)
            do imol=1,size(distmol,1)
            do bi1=1,size(dist,1)
                do bi2=1,size(dist,2)
                    x1=mi1+(real(bi1,rk)-.5_rk)*bin1
                    x2=mi2+(real(bi2,rk)-.5_rk)*bin2
                    if(instr1%findex==1)x1=acos(x1)*180._rk/pi
                    if(instr2%findex==1)x2=acos(x2)*180._rk/pi
                    write(instr0%set%ounit,*)x1,x2,imol,distmol(imol,bi1,bi2)!,std(distmol(:,bi1,bi2))
                    !x=mi+(real(bi,rk))*bin
                    !if(instr%findex==1)x=acos(x)*180._rk/pi
                    !write(ounit,*)x,dist(bi),std(distmol(:,bi))
                end do
                write(instr0%set%ounit,*)
            end do
            !write(ounit,*)
            end do
            if(instr0%set%ounit/=stdout)close(instr0%set%ounit)
        end select
        end do
        end if

      !   Pearson's correlation coefficiant

        call mean_var(vec1,m1,v1)
        call mean_var(vec2,m2,v2)
        instr0%cv%pearsoncoeff= sum((vec1(:)-m1)*(vec2(:)-m2))&
                /(real(size(vec1)-1,rk)*sqrt(v1*v2))

        !Räkna ut entropin för distributionen
        isoentropy=log((ma1-mi1)*(ma2-mi2))
        instr0%cv%entropy=-log(real(size(dist),rk)*bin1*bin2)&
        -bin1*bin2*sum(dist(:,:)*log(dist(:,:)),MASK=dist(:,:)/=0)-isoentropy

        !Standardavvikelse för entropin
        do i=1,size(distmol,1)
            molentropy(i)=-log(real(size(distmol,2)*size(distmol,3),rk)*bin1*bin2)&
            -bin1*bin2*sum(distmol(i,:,:)*log(distmol(i,:,:)),MASK=distmol(i,:,:)/=0)&
            -isoentropy
        end do
        instr0%cv%entropymutual=instr0%cv%entropy-sum(molentropy)/real(size(distmol,1),rk)
       ! instr%cv%entropydev=instr%cv%entropy-&
        !(-log(real(size(distmol),rk))-sum(distmol(:,:)*log(distmol(:,:)),MASK=distmol(:,:)/=0))

        deallocate(dist,distmol,molentropy)
    end subroutine corr_distrib!}}}

    subroutine center_of_membrane(molecules)!{{{
        integer(kind=ik) :: molecules(:),i,imol,j
        real(kind=rk) :: n
        ! Calculate center of membrane based on mass density of
        ! all atoms between atom a and b
        centerofmembrane=0;n=0
        do j=1,size(molecules)
           i=molecules(j)
            do imol=1,molt(i)%nmol
                centerofmembrane=centerofmembrane+center_of_molecule(i,imol)&
                *sum(masses(molt(i)%firstatom:molt(i)%lastatom))
            end do
            n=n+sum(masses(molt(i)%firstatom:molt(i)%lastatom))*molt(i)%nmol
        end do
        centerofmembrane=centerofmembrane/n

    end subroutine center_of_membrane!}}}

    subroutine whole!{{{
        integer(kind=ik) :: i,imol,j,k
        real(kind=rk) :: com(3),shift
        ! Fold back molecules that is outside the box
        do i=1,size(molt)
            do imol=1,molt(i)%nmol
                com=getatom(molt(i)%firstatom,imol) !center_of_molecule(i,imol) !Ändrade molt(i)%firstatom till i
                do k=molt(i)%firstatom,molt(i)%lastatom
                    coor(:,cind(k,imol))=com+mymodulo(coor(:,cind(k,imol))-com(:),box(:))
                end do
            end do
        end do
    end subroutine whole!}}}

   subroutine foldmol!{{{
        integer(kind=ik) :: i,imol,j,k
        real(kind=rk) :: com(3),shift
        ! Fold back molecules that is outside the box
        do i=1,size(molt)
            do imol=1,molt(i)%nmol
                com=center_of_molecule(i,imol) !Ändrade molt(i)%firstatom till i
                do j=1,3
                    shift=com(j)-modulo(com(j),box(j))
                    if(abs(shift)>box(j)/2._rk)then
                        do k=molt(i)%firstatom,molt(i)%lastatom
                            coor(j,cind(k,imol))=coor(j,cind(k,imol))-shift
                        end do
                    end if
                end do
                com=center_of_molecule(i,imol)
                if(com(3)<0 .OR. com(3)>box(3))then
                write(*,*)com,'COM'
                stop
                end if
            end do
        end do
    end subroutine foldmol!}}}

    function totalmoi(umol,imol) result(m)!{{{
        integer(kind=ik) :: umol,imol,i
        real(kind=rk) :: com(3),m
        m=0
        com=center_of_molecule(umol,imol)
        do i=molt(umol)%firstatom,molt(umol)%lastatom
            m=m+sum((getatom(i,imol)-com)**2)*(masses(i))
        end do
        
    end function totalmoi!}}}

!    subroutine test(a)!{{{
!    use f95_lapack
!!    use solver
!    integer(kind=ik) :: a,imol
!    integer(kind=4) ::info=0,lwork=102
!    real(kind=rk) :: tens(3,3),eigenvalues(3),work(1:102)
!    imol=1
!    tens=moi(a,imol)
!    write(*,*)tens
!    eigenvalues=0;work=0
!   ! tens=reshape([1,0,0,0,2,0,0,0,3],[3,3])
!    call&
!!    DSYEV('V','U',size(tens,1),tens,size(tens,2),eigenvalues,work,lwork,info)
!    LA_SYEVD(tens,eigenvalues,'V','U',info)
!    write(*,*)eigenvalues,info
!    write(*,*)tens
!    write(*,*)'**************************************************'
!    eigenvalues=0;tens=0
!!    call automatic_evd2_noorder(moi(a,imol),eigenvalues,tens)
!    write(*,*)eigenvalues
!    write(*,*)tens
!
!    stop
!    end subroutine test!}}}
   
!    function center_of_molecule(umol,imol) result(centerofmolecule)!{{{
!        integer(kind=ik) :: umol,i,j,imol
!        real(kind=rk) :: centerofmolecule(3)
!        centerofmolecule=0
!        j=umol
!        do i=molt(j)%firstatom,molt(j)%lastatom
!            centerofmolecule=centerofmolecule+getatom(i,imol)*masses(i)
!        end do
!        centerofmolecule=centerofmolecule/sum(masses(molt(j)%firstatom:molt(j)%lastatom))
!    end function center_of_molecule!}}}

    function moi(umol,imol) result(tensor)!{{{

    integer(kind=ik) :: umol,imol,i,j
    real(kind=rk) :: com(3),ri(3), tensor(3,3)
    real(kind=rk) :: Ixx,Iyy,Izz,Ixy,Ixz,Iyz

    j=umol
    com=center_of_molecule(j,imol)

    Ixx=0;Iyy=0;Izz=0
    Ixy=0;Ixz=0;Iyz=0

    do i=molt(j)%firstatom,molt(j)%lastatom
    
        ri=getatom(i,imol)-com
        Ixx=Ixx+masses(i)*(ri(2)**2+ri(3)**2)
        Iyy=Iyy+masses(i)*(ri(1)**2+ri(3)**2)
        Izz=Izz+masses(i)*(ri(1)**2+ri(2)**2)
        Ixy=Ixy-masses(i)*ri(1)*ri(2)
        Ixz=Ixz-masses(i)*ri(1)*ri(3)
        Iyz=Iyz-masses(i)*ri(2)*ri(3)

    end do
        tensor=reshape([&
                [Ixx,Ixy,Ixz],&
                [Ixy,Iyy,Iyz],&
                [Ixz,Iyz,Izz] &
                ],[3,3])

    end function moi!}}}

    subroutine traj(instr,j)!{{{
        type(instruct) :: instr
        real(kind=rk) :: eaverage(size(instr%datam,2))
        integer(kind=ik) :: frame,i,j
        select case(instr%set%calc(j))
        case('traj')
            do frame=1,size(instr%datam,2)
                eaverage(frame)=sum(instr%datam(:,frame))/size(instr%datam,1)
            end do
        
            do frame=1,size(instr%datam,2)
                write(instr%set%ounit,*)frame,eaverage(frame),std(instr%datam(:,frame))
            end do
        case('trajm')
            do i=1,size(instr%datam,1)
                do frame=1,size(instr%datam,2)
                    write(instr%set%ounit,*)frame,i,instr%datam(i,frame)
                end do
                write(instr%set%ounit,*)
            end do
        end select
    end subroutine traj!}}}

    subroutine postproc(instr)!{{{
        type(instruct) :: instr(:)
        integer(kind=ik) :: ounit,ios,i,j,k,n,corr1,corr2
        real(kind=rk) :: a
        character(kind=1,len=len(global_setflags%filename)) :: filename
        logical :: exists
        character(kind=1,len=100) :: pos,avfilename
        a=0._rk
        instr%cv%entropy=a/a
        n=0
        do i=1,size(instr)
            if(instr(i)%findex/=0)n=max(len(trim(instr(i)%instructionstring)),n)
        end do
        !open(78,file=trim(global_setflags%fileprefix)//'averages'//trim(global_setflags%filesuffix))
        do i=1,size(instr) ! Loop över instruktionsrader
        select case(instr(i)%findex)
            case(0,10)
            case(7) ! CORRELATE
                corr1=strvecindex(instr(:)%ref,instr(i)%set%corrindex(1))
                corr2=strvecindex(instr(:)%ref,instr(i)%set%corrindex(2))
                if(instr(i)%set%filename/='')then
                    instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%set%filename)!//trim(instr(i)%set%filesuffix)
                    filename=instr(i)%set%filename
                    instr(i)%set%ounit=nunit
                else
                    if(instr(i)%set%autofilename)then
                        instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%instructionstring)&
                                //trim(instr(corr1)%instructionstring)//'_'//&
                                trim(instr(corr2)%instructionstring)!//trim(instr(i)%set%filesuffix)
                        filename=instr(i)%set%filename
                        instr(i)%set%ounit=nunit
                    else
                        instr(i)%set%ounit=stdout
                        filename="stdout"
                    end if
                end if
                call corr_distrib(instr(i),instr(corr1),instr(corr2))
                avfilename=trim(instr(i)%set%fileprefix)//'averages'//trim(instr(i)%set%filesuffix)
                inquire(file=avfilename,exist=exists)
                if(.NOT. exists)then
                    pos='ASIS'
                else
                    pos='APPEND'
                end if
                open(78,file=trim(avfilename)&
                ,position=trim(pos),status='unknown')
                write(78,*)int(i,2),filename(len_trim(instr(i)%set%fileprefix)+1:len_trim(filename))," = ",&
                    instr(i)%cv%pearsoncoeff,instr(i)%cv%entropy,instr(i)%cv%entropymutual
                close(78)

            case(9) ! AVERAGE
                if(instr(i)%set%filename/='')then
                    instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%set%filename)
                    j=i;k=0
                    do 
                        j=j+1
                        if(instr(j)%findex/=0 .AND. instr(j)%findex/=10)then
                            k=k+1
                            if(k>=instr(i)%average_count)exit
                        endif
                    end do
                    instr(i)%set%ounit=nunit
                else
                    if(instr(i)%set%autofilename)then
                    instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%instructionstring(1:2))&
                                //"_"
                    j=i;k=0
                    do 
                        j=j+1
                        if(instr(j)%findex/=0 .AND. instr(j)%findex/=10)then
                            k=k+1
                            if(k==1)instr(i)%set%filename=trim(instr(i)%set%filename)//trim(instr(j)%instructionstring(1:2))
                            instr(i)%set%filename=trim(instr(i)%set%filename)&
                            //'_'//trim(instr(j)%instructionstring(4:))
                            if(k>=instr(i)%average_count)exit
                        endif
                    end do
                    instr(i)%set%ounit=nunit
                else
                    instr(i)%set%ounit=stdout
                    filename="stdout"
                    j=i;k=0
                    do 
                        j=j+1
                        if(instr(j)%findex/=0 .AND. instr(j)%findex/=10)then
                            k=k+1
                            if(k>=instr(i)%average_count)exit
                        endif
                    end do
                end if
                end if
               call instruction_average(instr,i,j)
               instr(i)%findex=instr(j)%findex

            case default ! DEFAULT

                    if(instr(i)%set%filename/='')then
                        instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%set%filename)
                        instr(i)%set%ounit=nunit
                        
                    else
                        if(instr(i)%set%autofilename)then
                        instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%instructionstring)
                        instr(i)%set%ounit=nunit
                    else
                        instr(i)%set%ounit=stdout
                        filename="stdout"
                    end if
                end if
            end select
        end do
        do i=1,size(instr) ! Loop över instruktionsrader
            avfilename=trim(instr(i)%set%fileprefix)//'averages'//trim(instr(i)%set%filesuffix)
            select case(instr(i)%findex)
                case(0,7,10)
                case default
                    !write(*,*)size(instr(i)%set%calc),allocated(instr(i)%set%calc)
                    !if(size(instr(i)%set%calc)/=0)then ! Make sure even Intel knows what to do
                    if(allocated(instr(i)%set%calc))then
                    do j=1,size(instr(i)%set%calc) ! Loop över postberäkningar
                        if(instr(i)%set%autofilename)then
                            filename=trim(instr(i)%set%filename)&
                            //'_'//trim(instr(i)%set%calc(j))//trim(instr(i)%set%filesuffix)
                        else 
                            filename=trim(instr(i)%set%filename)//trim(instr(i)%set%filesuffix)
                        end if
                        if(instr(i)%set%ounit/=stdout)open(unit=instr(i)%set%ounit,file=trim(filename),status="unknown",iostat=ios)
                        select case(trim(instr(i)%set%calc(j)))
                            case('distrib','distribm')
                                !write(*,*)'I',i,instr(i)%datam
                                call distrib(instr(i),i,j)
                            case('traj','trajm')
                                call traj(instr(i),j)
                            case('acorr','acorrm')
                                call autocorr(instr(i),i,j)
                        end select
                        if(instr(i)%set%ounit/=stdout)close(instr(i)%set%ounit,iostat=ios)
                    end do
                    end if
                    !if(size(instr(i)%set%calc)/=0)then ! Make sure even Intel knows what to do
                    select case(instr(i)%findex)
                        case(0,7,10)
                        case(1)
                            inquire(file=avfilename,exist=exists)
                            if(.NOT. exists)then
                                pos='ASIS'
                            else
                                pos='APPEND'
                            end if
                            open(78,file=trim(avfilename)&
                            ,position=trim(pos),status='unknown')
                            write(78,*)int(i,2),instr(i)%instructionstring(1:n)," = ",&
                            trim(adjustl(average(acos(instr(i)%datam)*180._rk/pi)))!,&
                            close(78)
                            !instr(i)%cv%entropy,instr(i)%cv%entropymutual
                        ! case(9)

                        case default
                            inquire(file=avfilename,exist=exists)
                            if(.NOT. exists)then
                                pos='ASIS'
                            else
                                pos='APPEND'
                            end if
                            open(78,file=trim(avfilename)&
                            ,position=trim(pos),status='unknown')
                            write(78,*)int(i,2),instr(i)%instructionstring(1:n)," = ",&
                            trim(adjustl(average(instr(i)%datam)))!,&
                            close(78)
!                            instr(i)%cv%entropy,instr(i)%cv%entropymutual
                    end select
                    !end if
            end select
            if(allocated(instr(i)%datam))deallocate(instr(i)%datam)
        end do
        if(.NOT.allocated(global_setflags%calc))then
            write(*,*)'     Nothing except averages is calculated.'
            write(*,*)'     To enable postcalcualtions (e.g. distributions), '
            write(*,*)'     use the "set calc" function.'
        end if
    end subroutine postproc!}}}
    
    subroutine autocorr(instr,ind,calcind)!{{{
        type(instruct) :: instr
        integer(kind=ik) :: i,j,imol,ind,calcind,indold=0
        real(kind=rk) :: norm
        real(kind=rk),allocatable :: acorr(:,:)
        save acorr,indold
        if(ind>indold)then
            allocate(acorr(size(instr%datam,1),0:size(instr%datam,2)-1))
        acorr=0
        do imol=1,size(instr%datam,1)
            call directtcf(instr%datam(imol,:),acorr(imol,:))
           !call directzero(instr%datam(imol,:),acorr(imol,:))
             norm=1._rk/acorr(imol,0)
            acorr(imol,:)=acorr(imol,:)*norm
        end do
        end if !ind>indold
        select case(instr%set%calc(calcind))
        case('acorr')
            do i=0,ubound(acorr,2)
                write(instr%set%ounit,*)i,sum(acorr(:,i))/real(size(acorr,1),rk),std(acorr(:,i))
            end do
        case('acorrm')
            do imol=1,size(instr%datam,1)
                do i=0,ubound(acorr,2)
                    write(instr%set%ounit,*)i,imol,acorr(imol,i)!,std(acorr(:,i))
                end do
                write(instr%set%ounit,*)
            end do
         end select
        if(COUNT(instr%set%calc(:)(1:4)=='traj')==1.OR.ind==indold)deallocate(acorr)
        indold=ind
        contains
        subroutine directtcf(datav,acorrv)
            real(kind=rk) :: datav(:),acorrv(0:size(datav)-1)
            integer(kind=ik) :: i,j
            acorrv=0
            do j=0,size(datav)-1
                !do i=1,size(datav)-j
                !    acorrv(j)=acorrv(j)+(datav(i)*datav(i+j))
                !end do
                !n=size(datav)-j
                acorrv(j)=sum(datav(1:size(datav)-j)*datav(j+1:size(datav)))/real(size(datav)-j,rk)
            end do
        end subroutine directtcf

        subroutine directtcfrs(datav,acorrv)
            real(kind=rk) :: datav(:),acorrv(0:size(datav)-1)
            integer(kind=ik) :: i,j
            acorrv=0
            do j=0,size(datav)-1
                !do i=1,size(datav)-j
                    acorrv(j)=acorrv(j)+(datav(1)*datav(1+j))/&
                    real(size(datav)-j,rk)
                !end do
            end do
        end subroutine directtcfrs


        subroutine directzero(datav,acorrv)
            real(kind=rk) ::&
            datav(:),acorrv(0:size(datav)-1),dataz(1:2**ceiling(log(real(2*size(datav)-1,rk))/log(2._rk)))
            integer(kind=ik) :: i,j
            acorrv=0;dataz=0
            dataz(1:size(datav))=datav(:)
            do j=0,size(datav)-1
                do i=1,size(datav)
                    acorrv(j)=acorrv(j)+(dataz(i)*dataz(i+j))
                end do
            end do
            do j=0,size(datav)-1
                acorrv(j)=acorrv(j)/real(size(acorrv)-j,rk)
            end do
        end subroutine directzero

        subroutine myft(datav,acorrv)
            real(kind=rk) ::&
            datav(:),acorrv(0:size(datav)-1),dataz(1:2**ceiling(log(real(2*size(datav)-1,rk))/log(2._rk))),&
            ftdataz(-size(dataz)/2+1:size(dataz)/2)
            integer(kind=ik) :: i,j
            acorrv=0;dataz=0
            dataz(1:size(datav))=datav(:)
            do j=0,size(datav)-1
                do i=1,size(datav)
                    acorrv(j)=acorrv(j)+(dataz(i)*dataz(i+j))
                end do
            end do
            do j=0,size(datav)-1
                acorrv(j)=acorrv(j)/real(size(acorrv)-j,rk)
            end do
        end subroutine myft
    end subroutine autocorr!}}}
end module trajop
