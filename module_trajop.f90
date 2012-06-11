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

module trajop
    use kinds
    use readtraj
    use apl
    use statistics
!    use input
    use util
    implicit none
   ! real(kind=rk),allocatable :: centers(:,:)
    
    contains

    function dirangle(a,b,imol) result(teta)!{{{
        ! Angle: a-b  against director in degrees
        real(kind=rk) :: teta,dir(3),xtz
        integer(kind=ik) :: a,b,imol
        !dir=director*sign(1._rk,dot_product(center_of_molecule(moltypeofuatom(a),imol)-centerofmembrane,director))
            dir=director
        if(global_setflags%centerofmembrane)then
            if(global_setflags%leaflets_defined)then
                if(ANY(imol==molt(moltypeofuatom(a))%lower))dir=-1._rk*director
            end if
        end if
        xtz=dot_product(atom(b)%coor(:,imol)-atom(a)%coor(:,imol),dir)
        teta=atan2(sqrt(abs(sum((atom(b)%coor(:,imol)-atom(a)%coor(:,imol))**2)*sum(dir**2)-xtz**2)),xtz)*180._rk/pi
        ! teta=(atom(b)%coor(3,imol)-atom(a)%coor(3,imol))/sqrt(sum((atom(b)%coor(:,imol)-atom(a)%coor(:,imol))**2))
        !write(300,*)teta,imol,atom(b)%coor(3,imol),atom(a)%coor(3,imol),sqrt(sum((atom(b)%coor(:,imol)-atom(a)%coor(:,imol))**2))
       ! teta=acos(dot_product(normalize(atom(b)%coor(:,imol)-atom(a)%coor(:,imol)),dir))*180._rk/pi
    end function dirangle!}}}

    function valenceangle(a,b,c,imol) result(teta)!{{{
        ! Angle a-b against a-c in degrees
        real(kind=rk) :: teta
        integer(kind=ik) :: a,b,c,imol
        teta=acos(dot_product(normalize(atom(a)%coor(:,imol)-atom(b)%coor(:,imol))&
                ,normalize(atom(a)%coor(:,imol)-atom(c)%coor(:,imol))))*180._rk/pi
    end function valenceangle!}}}

    function torsionangle(a,b,c,d,imol) result(teta)!{{{
        ! Torsion angle: a-b-c-d in degrees 
        real(kind=rk) :: Vba(1:3),Vcd(1:3),Vcb(1:3),teta,teta2
        integer(kind=ik) :: a,b,c,d,imol
        Vcb=normalize(atom(b)%coor(:,imol) - atom(c)%coor(:,imol))
        Vba = normalize(cross_product( &
                        cross_product(atom(a)%coor(:,imol) - atom(b)%coor(:,imol) ,Vcb) &
                        , Vcb))
        Vcd = normalize(cross_product( &
                        cross_product(atom(d)%coor(:,imol) - atom(c)%coor(:,imol) ,Vcb) &
                        ,Vcb))
        !teta = modulo(atan2(dot_product(Vcb,cross_product(Vcd,Vba)),dot_product(Vba,Vcd)) * 180._rk / pi,360._rk)
        teta = atan2(dot_product(Vcb,cross_product(Vcd,Vba)),dot_product(Vba,Vcd)) * 180._rk / pi
    end function torsionangle!}}}

    function bond_length(a,b,imol,jmol,instr) result(c)!{{{
        integer(kind=ik) :: a,b,imol
        integer(kind=ik),optional :: jmol
        type(instruct),optional :: instr
        real(kind=rk) :: c,vector(1:3),ac(3),bc(3)!,com(3)
        if(present(jmol))then
        if(present(instr).AND.instr%findex==13)then
            vector=mymodulo(atom(b)%coor(:,jmol)-atom(a)%coor(:,imol),box)
        else
           ! write(*,*)getatom(b,jmol)
           ! write(*,*)getatom(a,imol)
           ! stop
           
           
         !   if(any(imol==molt(moltypeofuatom(a))%upper).neqv.any(jmol==molt(moltypeofuatom(b))%upper)&
         !   .or.any(imol==molt(moltypeofuatom(a))%lower).neqv.any(jmol==molt(moltypeofuatom(b))%lower))then
            ac=atom(a)%coor(:,imol)
            bc=atom(b)%coor(:,jmol)
          !  if(ac(3)<0._rk)then
          !      if(any(imol==molt(moltypeofuatom(a))%upper))ac(3)=box(3)-ac(3)
          !      if(any(imol==molt(moltypeofuatom(a))%lower))ac(3)=ac(3)
          !  else if(ac(3)>box(3))then
          !      if(any(imol==molt(moltypeofuatom(a))%upper))ac(3)=ac(3)
          !      if(any(imol==molt(moltypeofuatom(a))%lower))ac(3)=box(3)-ac(3)
          !  end if
          !  if(bc(3)<0._rk)then
          !      if(any(jmol==molt(moltypeofuatom(b))%upper))bc(3)=box(3)-bc(3)
          !      if(any(jmol==molt(moltypeofuatom(b))%lower))bc(3)=bc(3)
          !  else if(bc(3)>box(3))then
          !      if(any(jmol==molt(moltypeofuatom(b))%upper))bc(3)=bc(3)
          !      if(any(jmol==molt(moltypeofuatom(b))%lower))bc(3)=box(3)-bc(3)
          !  end if
        !end if
            vector=bc-ac
        end if
        if(present(instr).AND.instr%set%xyrdf)vector(3)=0._rk
        if(present(instr).AND.instr%set%zrdf)vector(1:2)=0._rk
        c=sqrt(sum(vector**2))
        else
        vector=atom(b)%coor(:,imol)-atom(a)%coor(:,imol)
        c = sqrt(sum(vector**2))
        end if
    end function bond_length!}}}

    function distance_com(a,imol) result(distance)!{{{
        integer(kind=ik) :: a,imol
        real(kind=rk) :: dir(3),distance
        dir=director*sign(1._rk,dot_product(center_of_molecule(moltypeofuatom(a),imol)-centerofmembrane,director))
        distance=dot_product(atom(a)%coor(:,imol)-centerofmembrane,dir)
    end function distance_com!}}}

    function order_parameter(a,b,imol) result(c)!{{{
        integer(kind=ik) :: a,b,imol
        real(kind=rk) :: c
        c=1.5_rk*cos(pi/180._rk*&
        dirangle(a,b,imol))**2-0.5

    end function order_parameter!}}}

    function imoi(a,imol) result(c)!{{{
        integer(kind=ik) :: a,i,imol
        real(kind=rk) :: c,com(3)
        com=center_of_molecule(a,imol)
        c=0
        do i=molt(a)%firstatom,molt(a)%lastatom
            c=c+sum((atom(i)%coor(:,imol)-com)**2)*masses(i) 
        end do
            c=c*(common_setflags%traj_cscale/1.e-9_rk)**2
    end function imoi!}}}

    function center_of_mass(atoms,imol,a) result(center)
        integer(kind=ik) :: atoms(:),i,imol
        real(kind=rk) :: center(3),massweight,mass
        logical,optional :: a
        mass=1;massweight=0;center=0
        do i=1,size(atoms)
            if(.not.present(a))mass=masses(i)
            center=center+(atom(atoms(i))%coor(:,imol)*mass)
            massweight=massweight+mass
        end do
        center=center/massweight
    end function center_of_mass 

    subroutine add_hydrogen(helpers,hcoor,imol,func,bondlength,theta)
        integer(kind=ik),intent(in) :: helpers(:),func,imol
        integer(kind=ik) :: i
        real(kind=rk) ::&
        v1(3),v2(3),v3(3),v4(3),v5(3),u(3),bondlength,thetal
        real(kind=rk),intent(inout) :: hcoor(:)
        real(kind=rk),optional :: theta
        if(.not.present(theta) .and. func==3)stop 'CH2x needs rotation angle'
        select case(func)
        case(1) ! CH
            v1=atom(helpers(1))%coor(:,imol)
            v2=0._rk
            do i=2,size(helpers)
            ! Center of the nomalized vectors
            v2=v2+normalize(atom(helpers(i))%coor(:,imol)-v1)
            end do
            ! Normalize by number and shift it back to v1
            v2=v2/(size(helpers)-1)+v1
            hcoor=bondlength*normalize(v1-v2)+v1
        case(2) !CH double bond
            v1=atom(helpers(1))%coor(:,imol)
            v2=normalize(atom(helpers(2))%coor(:,imol)-v1)
            v3=normalize(atom(helpers(3))%coor(:,imol)-v1)
            ! thetal is the angle 2pi - C-C-C devided by 2
            ! to ensure equal (C-C-H) angles from both directions
            thetal=pi*(2-valenceangle(helpers(1),&
                    helpers(2),&
                    helpers(3),imol)/180._rk)/2
            u=normalize(cross_product(v2,v3))
            ! RV returns a rotational 3x3 matrix from a quaternion
            ! v2q turns a vector into a quaternion with angle theta
            ! Rotate by matrix multiplication RV*V3 where V3 is rotated
            hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v3))+v1
        case(3) !CH2
            thetal=theta*pi/180._rk
            v1=atom(helpers(1))%coor(:,imol)
            v2=normalize(atom(helpers(2))%coor(:,imol)-v1)
            v3=normalize(atom(helpers(3))%coor(:,imol)-v1)
            ! Perpendicular to the C-C-C plane
            v4=normalize(cross_product(v3,v2))
            ! Vector to rotate around. In the plane.
            u=normalize(v2-v3) 
            ! The vector to be rotated theta/2, perpendicular to u and v4
            v4=normalize(cross_product(v4,u))
            ! Rotate the new v4 around u by thetal 
            hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v4))+v1
        case(4) !CH3e
            thetal=acos(-1._rk/3) ! 109.47 deg in radians
            v1=atom(helpers(1))%coor(:,imol)
            v2=normalize(atom(helpers(2))%coor(:,imol)-v1)
            v3=normalize(atom(helpers(3))%coor(:,imol)-v1)
            u=normalize(cross_product(v3,v2))
            hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v2))+v1
        case(5) !CH3r
            thetal=120._rk*pi/180._rk
            v1=atom(helpers(1))%coor(:,imol)
            v2=atom(helpers(2))%coor(:,imol)
            v3=atom(helpers(3))%coor(:,imol) !First hydrogen added by CH3e
            u=normalize(v2-v1) ! C-C bond to rotate around
            v4=normalize(v3-v1) 
            hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v4))+v1
        case(6) !CH3s
            thetal=-120._rk*pi/180._rk
            v1=atom(helpers(1))%coor(:,imol)
            v2=atom(helpers(2))%coor(:,imol)
            v3=atom(helpers(3))%coor(:,imol) !First hydrogen added by CH3e
            u=normalize(v2-v1) ! C-C bond to rotate around
            v4=normalize(v3-v1)
            hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v4))+v1
        end select
    end subroutine add_hydrogen
    subroutine procop(instr,frame)!{{{
    ! Processing of operations from input
        integer(kind=ik) :: imol,jmol,i,j,k,l,m,frame,steps
        integer(kind=ik),allocatable ::&
        cindexes_l(:),cindexes_u(:),grid1(:,:),grid2(:,:)
        logical :: moltest
        real(kind=rk) :: teta,bl,v(3),x,scaling
        type(instruct) :: instr(:)
          do i=1,size(instr)
            select case(instr(i)%findex)
                case(0)
                    !if(instr(i)%setapl)call apl_grid(instr(i))
                case(1) !DIRECTOR ANGLE
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        !instr(i)%datam(jmol,frame)=&
                        !cos((pi/180._rk)*dirangle(instr(i)%atoms(1),instr(i)%atoms(2),imol))
                        instr(i)%datam(jmol,frame)=&
                        cos((pi/180._rk)*dirangle(instr(i)%atoms(1),instr(i)%atoms(2),imol))
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
                        instr(i)%datam(jmol,frame)=modulo(torsionangle(instr(i)%atoms(1),instr(i)%atoms(2)&
                        ,instr(i)%atoms(3),instr(i)%atoms(4),imol)-real(instr(i)%set%tshift,rk),360._rk)&
                        +real(instr(i)%set%tshift,rk)
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
                        ! Everything is handled in subroutine postproc
                case(8) ! DIPOLAR COUPLING
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
                        *(bl*common_setflags%traj_cscale)**(-3)/1000.!(1e-9 for meters and 1000 for kHz) 
                    end do
                 case(9) ! AVERAGE
                        ! Everything is handled in subroutine postproc
                 case(10) !DEFINE
                     select case(instr(i)%define)!{{{
                     case(1) !CENTEROFMEMBRANE
                        call center_of_membrane(common_setflags%membrane_moltypes)
                     case(2) !ATOM
                        do imol=1,molt(moltypeofuatom(instr(i)%atoms(1)))%nmol
                            select case(instr(i)%newatom%atype)
                            case('com','center_of_mass')
                            ! Define atom from center of mass of a defined submolecule.
                               atom(instr(i)%atoms(1))%coor(:,imol)=&
                                center_of_molecule(atom(instr(i)%atoms(1))%moltype,imol)
                            case('CH1')
                                call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,1,instr(i)%set%ch_bondlength)
                            case('CH1db')
                                call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,2,instr(i)%set%ch_bondlength)
                            case('CH2r','CH2s')
                                ! The following papers explain the choice of
                                ! theta:
                                !(1) MASTRYUKOV, V.; OSINA, E. Journal of Molecular Structure. 1977, 36, 127-132.
                                !(2) Dakkouri, M. Journal of Molecular Structure. 1997, 413-414, 133-152.
                                teta=126.1_rk-0.175_rk*valenceangle(&
                                instr(i)%newatom%helpers(1),&
                                instr(i)%newatom%helpers(2),&
                                instr(i)%newatom%helpers(3),imol)
                                select case(instr(i)%newatom%atype)
                                case('CH2r')
                                    call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,3,instr(i)%set%ch_bondlength,teta/2._rk)
                                case('CH2s')
                                    call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,3,instr(i)%set%ch_bondlength,-teta/2._rk)
                                end select
                            case('CH3e')
                                call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,4,instr(i)%set%ch_bondlength)
                            case('CH3r')
                                call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,5,instr(i)%set%ch_bondlength)
                            case('CH3s')
                                call add_hydrogen(instr(i)%newatom%helpers,&
                                atom(instr(i)%atoms(1))%coor(:,imol),imol,6,instr(i)%set%ch_bondlength)
                            end select
                        end do
                     case(3) !LEAFLET
                            if(global_setflags%apl .OR. global_setflags%gd)call apl_grid(instr(i))
                     case(4)
                         call foldmol('com')
                     end select!}}}
                 case(11) ! AREA PER LIPID
                     call apl_calc(instr(i),frame)
                     if(allocated(global_setflags%writeframe))then
                        do j=1,size(instr(i)%set%writeframe)
                        if(trim(instr(i)%set%writeframe(j)%outformat)=='apl'&
                        .AND.instr(i)%set%writeframe(j)%framenumber==frame)&
                        call apl_matrix_out(frame,instr(i))
                        end do
                    end if
                case(12) ! ISOTROPIC MOMENT OF INERTIA
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=imoi(atom(instr(i)%atoms(1))%moltype,imol)
                    end do
                case(13) ! RADIAL DISTRIBUTION FUNCTION (RDF, g(r))
                    l=0
                    instr(i)%rdf_pairs=0
                    moltest=(atom(instr(i)%atoms(1))%moltype==atom(instr(i)%atoms(2))%moltype)
                    select case(instr(i)%set%leaflet)!{{{
                    case(0)!Both
                        do j=1,molt(atom(instr(i)%atoms(1))%moltype)%nmol
                            imol=j
                            do k=merge(j+1,1,moltest),molt(atom(instr(i)%atoms(2))%moltype)%nmol
                                jmol=k
                                l=l+1
                                instr(i)%rdf_pairs(l)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol,jmol,instr(i)) !%set%xyrdf)
                            end do
                        end do
                    case(1)!Lower
                        do j=1,size(molt(atom(instr(i)%atoms(1))%moltype)%lower)
                            imol=molt(atom(instr(i)%atoms(1))%moltype)%lower(j)
                            do k=merge(j+1,1,moltest),size(molt(atom(instr(i)%atoms(2))%moltype)%lower)
                                jmol=molt(atom(instr(i)%atoms(2))%moltype)%lower(k)
                                l=l+1
                                instr(i)%rdf_pairs(l)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol,jmol,instr(i)) !%set%xyrdf)
                            end do
                        end do
                    case(2)!Upper

                        do j=1,size(molt(atom(instr(i)%atoms(1))%moltype)%upper)
                            imol=molt(atom(instr(i)%atoms(1))%moltype)%upper(j)
                            do k=merge(j+1,1,moltest),size(molt(atom(instr(i)%atoms(2))%moltype)%upper)
                                jmol=molt(moltypeofuatom(instr(i)%atoms(2)))%upper(k)
                                l=l+1
                                instr(i)%rdf_pairs(l)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol,jmol,instr(i)) !%set%xyrdf)
                            end do
                        end do
                    end select!}}}
                    call rdf_distrib(instr(i),frame)
                    instr(i)%datam(1,frame)=sum((1._rk/instr(i)%rdf_pairs)**3)/size(instr(i)%rdf_pairs)
                case(14) !XYZ
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        v=[atom(instr(i)%atoms(1))%coor(:,imol)]
                        select case(instr(i)%instructionstring(1:1))
                        case('X')
                            instr(i)%datam(jmol,frame)=v(1)
                        case('Y')
                            instr(i)%datam(jmol,frame)=v(2)
                        case('Z')
                            instr(i)%datam(jmol,frame)=v(3)
                        end select
                    end do
                case(15) !Griddiff
                    instr(i)%datam(1,frame)=griddiff(grid_lower,grid_upper)
                    if(.not.allocated(instr(i)%rdf_dist))then
                        allocate(instr(i)%rdf_dist(-sum(grid_upper):sum(grid_lower)))
                        instr(i)%rdf_dist=0
                    end if
                    do j=1,size(grid_lower,1)
                    do k=1,size(grid_lower,2)
                        l=grid_lower(j,k)-grid_upper(j,k)
                        instr(i)%rdf_dist(l)=instr(i)%rdf_dist(l)+1/real((maxframes-skipframes),rk)
                    end do
                    end do

                case(16) ! shuffle_GD
                    stop 'SHUFFLE (16): atomshuffle is still using coor. FIX!!'
                    allocate(grid1(1:instr(i)%set%aplgrid(1),1:instr(i)%set%aplgrid(2)))!{{{
                    allocate(grid2(1:instr(i)%set%aplgrid(1),1:instr(i)%set%aplgrid(2)))
!                    if(.not.allocated(instr(i)%rdf_dist))then
!                        allocate(instr(i)%rdf_dist(-sum(grid_upper):sum(grid_lower)))
!                        instr(i)%rdf_dist=0
!                    end if
                    steps=1
                    do j=1,steps
                        call atomshuffle(common_setflags%shuffle_atoms,instr(i)%atoms,1,cindexes_l)
                        call atomshuffle(common_setflags%shuffle_atoms,instr(i)%atoms,2,cindexes_u)
                    !    call countmol(cindexes_l,grid1) ! DISABLED DUE TO NEW HANDLING OF COORDINATES!!!
                    !    call countmol(cindexes_u,grid2) ! DISABLED DUE TO NEW HANDLING OF COORDINATES!!!
                    if(.not.allocated(instr(i)%rdf_dist))then
                        allocate(instr(i)%rdf_dist(-sum(grid1):sum(grid1)))
                        instr(i)%rdf_dist=0
                    end if

                     !   !x=griddiff(grid1,grid2)
                     !   !instr(i)%cv%mean=instr(i)%cv%mean+x
                     !   !instr(i)%cv%n=instr(i)%cv%n+1
                     !   !instr(i)%cv%meandev=instr(i)%cv%meandev+x**2
                     !   !k=int(x*size(instr(i)%rdf_dist)+1)
                        !write(*,*)k,griddiff(grid1,grid2)
                     !   !if(k==(size(instr(i)%rdf_dist)+1))k=size(instr(i)%rdf_dist)
                        
                       ! write(*,*)k,'k',griddiff(grid1,grid2),allocated(instr(i)%rdf_dist)
                    !    !instr(i)%rdf_dist(k)=instr(i)%rdf_dist(k)+real(size(cindexes_l),rk)/real(((maxframes-skipframes)*steps),rk)
                    do k=1,size(grid1,1)
                    do l=1,size(grid2,2)
                        m=grid1(k,l)-grid2(k,l)
                        instr(i)%rdf_dist(m)=instr(i)%rdf_dist(m)+1/real(((maxframes-skipframes)*steps),rk)
                    end do
                    end do
                    end do
                    
                   ! if(frame==maxframes)then
                   ! do j=1,size(instr(i)%rdf_dist)
                   !     x=(real(j,rk)-0.5_rk)*(1._rk/size(instr(i)%rdf_dist))
                   !     write(instr(i)%set%ounit,*)x,instr(i)%rdf_dist(j)
                   ! end do
                   ! write(instr(i)%set%ounit,*)
                   ! end if

                    deallocate(grid1,grid2)!}}}
                case(17)
                    select case(instr(i)%instructionstring(1:2))
                        case('BX')
                            instr(i)%datam(1,frame)=box(1)
                        case('BY')
                            instr(i)%datam(1,frame)=box(2)
                        case('BZ')
                            instr(i)%datam(1,frame)=box(3)
                        case('XY')
                            instr(i)%datam(1,frame)=box(1)*box(2)
                        case('XZ')
                            instr(i)%datam(1,frame)=box(1)*box(3)
                        case('YZ')
                            instr(i)%datam(1,frame)=box(2)*box(3)
                        case('BV')
                            instr(i)%datam(1,frame)=box(1)*box(2)*box(3)
                    end select
                case(18)
                    instr(i)%datam(1,frame)=(box(1)*box(2))/(molt(instr(i)%atoms(1))%nmol/2._rk)
                case(21) !Constant datam
                    do jmol=1,instr(i)%nmolop
                        imol=instr(i)%molind(jmol)
                        instr(i)%datam(jmol,frame)=&
                        instr(i)%set%const%value
                    end do
                case(22) !MBL distance between atoms in different molecules
                    l=0!{{{
                    moltest=(atom(instr(i)%atoms(1))%moltype==atom(instr(i)%atoms(2))%moltype)
                    select case(instr(i)%set%leaflet)
                    case(0)!Both
                        do j=1,molt(atom(instr(i)%atoms(1))%moltype)%nmol
                            imol=j
                            do k=merge(j+1,1,moltest),molt(atom(instr(i)%atoms(2))%moltype)%nmol
                                jmol=k
                                l=l+1
                                !write(*,*)l,size(instr(i)%datam,1)
                                instr(i)%datam(l,frame)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol,jmol,instr(i)) !%set%xyrdf)
                            end do
                        end do
                    case(1)!Lower
                        do j=1,size(molt(atom(instr(i)%atoms(1))%moltype)%lower)
                            imol=molt(atom(instr(i)%atoms(1))%moltype)%lower(j)
                            do k=merge(j+1,1,moltest),size(molt(atom(instr(i)%atoms(2))%moltype)%lower)
                                jmol=molt(atom(instr(i)%atoms(2))%moltype)%lower(k)
                                l=l+1
                                instr(i)%datam(l,frame)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol,jmol,instr(i)) !%set%xyrdf)
                            end do
                        end do
                    case(2)!Upper

                        do j=1,size(molt(atom(instr(i)%atoms(1))%moltype)%upper)
                            imol=molt(atom(instr(i)%atoms(1))%moltype)%upper(j)
                            do k=merge(j+1,1,moltest),size(molt(atom(instr(i)%atoms(2))%moltype)%upper)
                                jmol=molt(atom(instr(i)%atoms(2))%moltype)%upper(k)
                                l=l+1
                                instr(i)%datam(l,frame)=bond_length(instr(i)%atoms(1),instr(i)%atoms(2),imol,jmol,instr(i)) !%set%xyrdf)
                            end do
                        end do
                    end select!}}}
            end select

            if(instr(i)%set%slice%switch .and. allocated(instr(i)%datam))then
            ! SLICE THE BOX AND CHECK MOLECULES IN THE SLICE
                if(.not.allocated(instr(i)%set%slice%bintest))&
                allocate(instr(i)%set%slice%bintest(instr(i)%nmolop,skipframes+1:maxframes))
                do jmol=1,instr(i)%nmolop
                    imol=instr(i)%molind(jmol)
                    v=center_of_molecule(atom(instr(i)%atoms(1))%moltype,imol)
                    select case(instr(i)%set%slice%typ)
                    case('Z_ou')
                    if(v(3)<=instr(i)%set%slice%lower.or.v(3)>=instr(i)%set%slice%upper)then
                        instr(i)%set%slice%bintest(jmol,frame)=.TRUE.
                    else
                        instr(i)%set%slice%bintest(jmol,frame)=.FALSE.
                    end if
                    case('Z_in')
                    if(v(3)>=instr(i)%set%slice%lower.and.v(3)<=instr(i)%set%slice%upper)then
                        instr(i)%set%slice%bintest(jmol,frame)=.TRUE.
                    else
                        instr(i)%set%slice%bintest(jmol,frame)=.FALSE.
                    end if
                    case default
                        stop 'SLICE: something is wrong (subroutine procop)'
                    end select
                end do
            end if

            if(instr(i)%set%scaling%switch .and. allocated(instr(i)%datam))then
                scaling=0
                select case(trim(instr(i)%set%scaling%typ))
                ! Scale the function results with the box in different ways
                case('xy','yx')
                    scaling=1._rk/(box(1)*box(2))
                case('xz','zx')
                    scaling=1._rk/(box(1)*box(3))
                case('yz','zy')
                    scaling=1._rk/(box(2)*box(3))
                case('xyz')
                    scaling=1._rk/(box(1)*box(2)*box(3))
                case('x');scaling=1._rk/box(1)
                case('y');scaling=1._rk/box(2)
                case('z');scaling=1._rk/box(3)
                end select
                !write(*,*)i,frame,scaling,trim(instr(i)%set%scaling%typ),instr(i)%instructionstring,instr(i)%findex
                instr(i)%datam(:,frame)=instr(i)%datam(:,frame)*scaling
                !write(*,*)i,frame,scaling
            end if
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

    subroutine rdf_distrib(instr,frame)!{{{
        type(instruct),intent(inout) :: instr
        integer(kind=ik) :: bi,i,j
        integer(kind=ik),optional :: frame
        real(kind=rk) :: mi,ma,x,dp,bvol,numberdensity,modbox(3)
            mi=0._rk
        if(present(frame))then
            modbox=box
            if(trajtype=='gro')then !Rescale distances to ångström
                instr%rdf_pairs=instr%rdf_pairs*10._rk
                modbox=modbox*10._rk
            end if
            if(instr%set%xyrdf)modbox(3)=1._rk !Normalize with area if xyrdf
            !if(instr%set%zrdf)modbox(1:2)=1._rk !ZRDF
            ma=instr%set%rdf_binsize*instr%set%distbin
            !if(instr%set%zrdf)ma=modbox(3) !!ZRDF
            instr%rdf_bin=(ma-mi)/real(size(instr%rdf_dist,1),rk)
            numberdensity=(size(instr%rdf_pairs)/product(modbox(1:3))) 
            do i=1,size(instr%rdf_pairs)
                bi=int((instr%rdf_pairs(i)-mi)/instr%rdf_bin+1.5_rk)
                if(instr%set%xyrdf)then
                    bvol=2*pi*((mi+(real(bi,rk)-1._rk)*instr%rdf_bin)*instr%rdf_bin)
                    if(bi==1)bvol=(instr%rdf_bin)**2*pi/4._rk
                    ! special integration domain for the first point
                else
                    bvol=(4*pi*((mi+(real(bi,rk)-1._rk)*instr%rdf_bin)**2+instr%rdf_bin**2/12._rk))*instr%rdf_bin
                    if(bi==1)bvol=(instr%rdf_bin)**3*pi/6._rk
                    ! special integration domain for the first point
                    !if(instr%set%zrdf)bvol=instr%rdf_bin**2 !ZRDF
                end if

                if(bi<=0.or.bi>size(instr%rdf_dist,1))then
                else
                    instr%rdf_dist(bi)=instr%rdf_dist(bi)+1._rk/(bvol*numberdensity*(maxframes-skipframes))
                    !&
                    !((4*pi*((mi+(real(bi,rk)-0.5_rk)*instr%rdf_bin)**2+instr%rdf_bin**2/12._rk))&
                    !*instr%rdf_bin*(maxframes-skipframes)*(size(instr%rdf_pairs)/product(modbox(1:3))))
                end if
            end do
        else !.not.present(frame)
            !Skriv ut distribution till disk
!                x=mi+(real(1,rk)-1._rk)*instr%rdf_bin
!                write(instr%set%ounit,*)x,instr%rdf_dist(1)
            do bi=1,size(instr%rdf_dist,1)
                x=mi+(real(bi,rk)-1._rk)*instr%rdf_bin
                write(instr%set%ounit,*)x,instr%rdf_dist(bi)
            end do
!                x=mi+(real(size(instr%rdf_dist,1),rk)-1._rk)*instr%rdf_bin
!                write(instr%set%ounit,*)x,instr%rdf_dist(size(instr%rdf_dist,1))
                write(instr%set%ounit,*)
                if(allocated(instr%rdf_dist))deallocate(instr%rdf_dist)
                if(allocated(instr%rdf_pairs))deallocate(instr%rdf_pairs)
        end if !present(frame)
    end subroutine rdf_distrib!}}}

    subroutine write_distrib(instr)!{{{
        type(instruct) :: instr
        integer(kind=ik) :: bi
        real(kind=rk) :: x,mean,var,meandev
        select case(instr%findex)
        case(0)
            mean=instr%cv%mean/instr%cv%n
            var=(1/(real(instr%cv%n,rk)-1))*(instr%cv%meandev-&
            (1/real(instr%cv%n,rk))*(instr%cv%mean)**2)
            meandev=sqrt(var/real(instr%cv%n,rk))
            !Skriv ut distribution till disk
            do bi=1,size(instr%rdf_dist,1)
                x=0._rk+(real(bi,rk)-0.5_rk)*instr%rdf_bin
                write(instr%set%ounit,*)x,instr%rdf_dist(bi),(1/sqrt(2*pi*var))*exp(-((x-mean)**2)/(2*var))
            end do
            write(instr%set%ounit,*)
        case(15,16)
            do bi=lbound(instr%rdf_dist,1),ubound(instr%rdf_dist,1)
                write(*,*)lbound(instr%rdf_dist,1),size(instr%rdf_dist),bi
                write(instr%set%ounit,*)bi,instr%rdf_dist(bi)!,(1/sqrt(2*pi*var))*exp(-((x-mean)**2)/(2*var))
            end do
            write(instr%set%ounit,*)
        end select
            
            
    end subroutine write_distrib!}}}

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
        if(size(datam,1)==1)then !Added 2011-08-31: Corrected 120329
            var=sum((datam(1,:)-mean)**2)/real(size(datam,2)-1,rk)
        else
            var=sum((meant(:)-mean)**2)/real(size(meant)-1,rk)
        end if
        meandev=sqrt(var/real(size(meant),rk))
        string=getmeanwithdev(mean,sqrt(var))
        write(string2,*)mean,meandev,sqrt(var)
        !write(*,*)string2
        string2=trim(adjustl(string))//" "//trim(string2)
    end function average!}}}

    function slice_average(instr) result(string2)!{{{
        implicit none
        type(instruct) :: instr
        character(len=24) :: string
        character(len=300) :: string2
        integer(kind=ik) :: frame,i
        real(kind=rk) :: meant(1:size(instr%datam,2)),mean,meandev,var
        do frame=1,size(instr%datam,2)
            meant(frame)=&
            sum(instr%datam(:,skipframes+frame),MASK=instr%set%slice%bintest(:,skipframes+frame))&
            /real(count(instr%set%slice%bintest(:,skipframes+frame)),rk)
        end do
        mean=sum(meant)/real(size(meant),rk)
        !if(size(datam,1)==1)then !Added 2011-08-31
        !    var=sum((datam(:,2)-mean)**2)/real(size(datam,2)-1,rk)
        !else
            var=sum((meant(:)-mean)**2)/real(size(meant)-1,rk)
        !end if
        meandev=sqrt(var/real(size(meant),rk))
        string=getmeanwithdev(mean,sqrt(var))
        write(string2,*)mean,meandev,sqrt(var)
        !write(*,*)string2
        string2=trim(adjustl(string))//" framedev. "//trim(string2)
    end function slice_average!}}}

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
        if(instr(ind)%set%molaverage)then !COMBINE
        do 
            j=j+1
            if(instr(j)%findex/=0 .AND. instr(j)%findex/=10)then
                k=k+1
                mols=mols+instr(j)%nmolop
                if(k>=instr(ind)%average_count)exit
            endif
        end do
        allocate(instr(ind)%datam(mols,size(instr(j)%datam,2)))
        else ! AVERAGE
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
                    !write(*,*)mols,instr(j)%nmolop,j
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
        integer(kind=ik) :: bi,i,j,k,l,ind,calcind,indold=0
        !integer(kind=ik),optional :: ounit
        real(kind=rk) :: vec(size(instr%datam))
        real(kind=rk) :: mi,ma,x,bin,dp,dpm,isoentropy
        real(kind=rk),allocatable :: dist(:),distmol(:,:),molentropy(:)
        save indold,dist,distmol,mi,bin

        if(ind/=indold)then
            !if(allocated(dist))deallocate(dist);
            allocate(dist(instr%set%distbin))
            !if(allocated(distmol)deallocate(distmol);
            allocate(distmol(size(instr%datam,1),size(dist)))
            !if(allocated(molentropy)deallocate(molentropy);
            allocate(molentropy(size(distmol,1)))
            vec=reshape(instr%datam,[size(instr%datam)])
            dist=0
            distmol=0
            if(instr%set%distminmax%switch)then
                mi=instr%set%distminmax%mi
                ma=instr%set%distminmax%ma
            else
                select case(instr%findex)
                case(1)
                    mi=-1._rk
                    ma=1._rk
                case(2)
                    mi=0._rk
                    ma=180._rk
                !case(3)
                ! mi=0._rk
                ! ma=360._rk
                !   mi=-180._rk
                !  ma=180._rk
                case default
                    mi=minval(vec)
                    ma=maxval(vec)
                end select
            end if
            bin=(ma-mi)/real(size(dist),rk)
            dp=1._rk/(real(size(vec),rk)*bin)
            dpm=1._rk/(real(size(instr%datam,2),rk)*bin)
            k=0;l=1
            do i=1,size(vec)
                bi=int((vec(i)-mi)/bin+1._rk)
                if(bi==size(dist)+1)bi=size(dist)
                if(bi<=0.or.bi>size(dist))then
                        if(.not.instr%set%distminmax%switch)then
                        write(*,*)bi,"bi",vec(i),"vec",i,"i",ma,"ma",mi,"mi",dp,"dp",bin,"bin"
                        stop "subroutine distrib"
                    end if
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
            x=mi
            if(instr%findex==1)x=acos(x)*180._rk/pi
            write(instr%set%ounit,*)x,dist(1),std(distmol(:,1))
            do bi=1,size(dist)
                !x=mi+(real(bi,rk)-1._rk)*bin
                !if(instr%findex==1)x=acos(x)*180._rk/pi
                !write(instr%set%ounit,*)x,dist(bi),std(distmol(:,bi))
                x=mi+(real(bi,rk)-.5_rk)*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,dist(bi),std(distmol(:,bi))
            end do
            x=ma
            if(instr%findex==1)x=acos(x)*180._rk/pi
            write(instr%set%ounit,*)x,dist(size(dist)),std(distmol(:,size(dist)))
       ! case('distribm')
       !     do i=1,size(distmol,1)
       !      do bi=1,size(dist)
       !         x=mi+(real(bi,rk)-1._rk)*bin
       !         if(instr%findex==1)x=acos(x)*180._rk/pi
       !         write(instr%set%ounit,*)x,i,distmol(i,bi)!std(distmol(:,bi))
       !         x=mi+(real(bi,rk))*bin
       !         if(instr%findex==1)x=acos(x)*180._rk/pi
       !         write(instr%set%ounit,*)x,i,distmol(i,bi)!std(distmol(:,bi))
       !     end do
       !         write(instr%set%ounit,*)
       !     end do

        case('distribm')
            do i=1,size(distmol,1)
                !x=mi
                !if(instr%findex==1)x=acos(x)*180._rk/pi
                !write(instr%set%ounit,*)x,i,distmol(i,1)!std(distmol(:,bi))
             do bi=1,size(dist)
                x=mi+(real(bi,rk)-.5_rk)*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,i,distmol(i,bi)!std(distmol(:,bi))
                !x=mi+(real(bi,rk))*bin
                !if(instr%findex==1)x=acos(x)*180._rk/pi
                !write(instr%set%ounit,*)x,i,distmol(i,bi)!std(distmol(:,bi))
            end do
                !x=ma
                !if(instr%findex==1)x=acos(x)*180._rk/pi
                !write(instr%set%ounit,*)x,i,distmol(i,size(distmol))!std(distmol(:,bi))
                write(instr%set%ounit,*)
                write(instr%set%ounit,*)
            end do

        end select
        if(COUNT(instr%set%calc(:)(1:7)=='distrib')==1.OR.ind==indold)deallocate(dist,distmol)
        indold=ind
    end subroutine distrib!}}}

    subroutine distrib_test(instr,ind,calcind)!{{{
        type(instruct),intent(inout) :: instr
        integer(kind=ik) :: bi,i,j,ind,calcind,indold=0
        !integer(kind=ik),optional :: ounit
        real(kind=rk) :: vec(size(instr%datam))
        real(kind=rk) :: mi,ma,x,bin,dp,dpm,isoentropy
        real(kind=rk),allocatable :: norm(:),dist(:),distmol(:,:),molentropy(:)
        save indold,dist,distmol,mi,bin

        if(ind>indold)then
            !if(allocated(dist))deallocate(dist);
            allocate(dist(instr%set%distbin))
            allocate(norm(size(dist)))
            !if(allocated(distmol)deallocate(distmol);
            allocate(distmol(size(instr%datam,1),size(dist)))
            !if(allocated(molentropy)deallocate(molentropy);
            allocate(molentropy(size(distmol,1)))
            vec=reshape(instr%datam,[size(instr%datam)])
            dist=0;norm=0
            distmol=0
            select case(instr%findex)
            case(1)
                mi=-1._rk
                ma=1._rk
            case(2)
                mi=0._rk
                ma=180._rk
            case(13)
                mi=0._rk
                ma=30._rk
            !case(3)
               ! mi=0._rk
               ! ma=360._rk
             !   mi=-180._rk
              !  ma=180._rk
            case default
                mi=minval(vec)
                ma=maxval(vec)
            end select
            bin=(ma-mi)/real(size(dist),rk)
            dp=1._rk/(real(size(vec),rk)*bin)
            dpm=1._rk/(real(size(instr%datam,2),rk)*bin)
            do i=1,size(vec)
                bi=int((vec(i)-mi)/bin+1._rk)
           ! write(*,*)bi,'BI'
                !if(bi==size(dist)+1)bi=size(dist)
                if(bi<=0.or.bi>size(dist))then
                   ! write(*,*)bi,"bi",vec(i),"vec",i,"i",ma,"ma",mi,"mi",dp,"dp",bin,"bin"
                   ! stop "subroutine distrib"   
                else
                    if(instr%findex==13)then
                    dist(bi)=dist(bi)+1/((4*pi*(mi+(real(bi,rk)-0.5_rk)*bin)**2)*bin)!vec(i)**2*bin)
                    norm(bi)=norm(bi)+1
                   ! dist(bi)=dist(bi)+1/(2*pi*vec(i))
                    j=mod(i-1,size(distmol,1))+1
                    distmol(j,bi)=distmol(j,bi)+1/(4*pi*vec(i)**2)
                    else
                    dist(bi)=dist(bi)+dp
                    j=mod(i-1,size(distmol,1))+1
                    distmol(j,bi)=distmol(j,bi)+dpm
                    end if
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
!            dist=dist/sum(dist(1:300)*4*pi*[((mi+(real(bi,rk)-0.5_rk)*bin)**2,bi=1,300)])
!dist(301:)=0._rk
!isoentropy=sum(norm(1:300))/((4*pi*(mi+300._rk*bin)**3)/3)
isoentropy=(size(molt(atom(instr%atoms(1))%moltype)%upper)&
+size(molt(atom(instr%atoms(2))%moltype)%upper))/(box(1)*box(2)*box(3))*1001
dist(1:300)=dist(1:300)/(isoentropy)!*bin*4*pi*[((mi+(real(bi,rk)-0.5_rk)*bin)**2,bi=1,300)])

            !Skriv ut distribution till disk
            do bi=1,size(dist)
                x=mi+(real(bi,rk)-0.5_rk)*bin
                if(instr%findex==1)x=acos(x)*180._rk/pi
                write(instr%set%ounit,*)x,dist(bi),std(distmol(:,bi))
                !x=mi+(real(bi,rk))*bin
                !if(instr%findex==1)x=acos(x)*180._rk/pi
                !write(instr%set%ounit,*)x,dist(bi),std(distmol(:,bi))
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
    end subroutine distrib_test!}}}

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
        if(instr1%set%distminmax%switch)then
            mi1=instr1%set%distminmax%mi
            ma1=instr1%set%distminmax%ma
        else
            select case(instr1%findex)
                case(1)
                    mi1=-1._rk
                    ma1=1._rk
                case(2)
                    mi1=0._rk
                    ma1=180._rk
                !case(3)
                ! mi1=0._rk
                ! ma1=360._rk
                !   mi1=-180._rk
                !  ma1=180._rk
                case default
                    mi1=minval(vec1)
                    ma1=maxval(vec1)
            end select
        end if
        if(instr2%set%distminmax%switch)then
            mi2=instr2%set%distminmax%mi
            ma2=instr2%set%distminmax%ma
        else
            select case(instr2%findex)
                case(1)
                    mi2=-1._rk
                    ma2=1._rk
                case(2)
                    mi2=0._rk
                    ma2=180._rk
                !case(3)
                ! mi2=0._rk
                ! ma2=360._rk
                !   mi2=-180._rk
                !  ma2=180._rk
                case default
                    mi2=minval(vec2)
                    ma2=maxval(vec2)
            end select
        end if
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
                if(.not.instr1%set%distminmax%switch .and. bi1<=0 .or. bi1>size(dist,1))then
                        write(*,*)bi1,"bi1",vec1(i),"vec1",i,"i",ma1,"ma1",mi1,"mi1",dp,"dp",bin1,"bin1"
                        stop 'subroutine corr_distrib'
                else if(.not.instr2%set%distminmax%switch)then
                    write(*,*)bi2,"bi2",vec2(i),"vec2",i,"i",ma2,"ma2",mi2,"mi2",dp,"dp",bin2,"bin2"
                    stop "subroutine corr_distrib"
                end if
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

    subroutine vs_distrib(instr0,instr1,instr2)!{{{
        type(instruct),intent(inout) :: instr0,instr1,instr2
        integer(kind=ik) :: bi,i,j,l,imol,ios
        character(len=len(instr0%set%filename)) :: filename
        real(kind=rk) :: mi,ma,x,bin
        real(kind=rk),allocatable :: dist(:,:),distmol(:,:,:)
        !write(*,*)allocated(instr1%datam),"allocated"
        !write(*,*)instr0%set%filename
        allocate(dist(instr2%set%distbin,2))
        allocate(distmol(size(instr1%datam,1),size(dist,1),2))
        dist=0
        distmol=0

        if(instr2%set%distminmax%switch)then
            mi=instr2%set%distminmax%mi
            ma=instr2%set%distminmax%ma
        else
            select case(instr2%findex)
                case(1)
                    mi=-1._rk
                    ma=1._rk
                case(2)
                    mi=0._rk
                    ma=180._rk
                case default
                    mi=minval(instr2%datam)
                    ma=maxval(instr2%datam)
            end select
        end if
        bin=(ma-mi)/real(size(dist,1),rk)
        do l=lbound(instr2%datam,1),ubound(instr2%datam,1)
            do i=lbound(instr2%datam,2),ubound(instr2%datam,2)
                bi=int(((instr2%datam(l,i)-mi)/bin)+1._rk)
                if(bi==size(dist,1)+1)bi=size(dist,1)
                if(bi<=0.or.bi>size(dist,1))then
                    if(.not.instr2%set%distminmax%switch)then
                    write(*,*)bi,"bi",instr2%datam(l,i),"instr2",i,"i",ma,"ma",mi,"mi",bin,"bin"
                    stop "subroutine corr_distrib"
                    end if
                else
                    dist(bi,1)=dist(bi,1)+instr1%datam(l,i)
                    dist(bi,2)=dist(bi,2)+1
                    distmol(l,bi,1)=distmol(l,bi,1)+instr1%datam(l,i)
                    distmol(l,bi,2)=distmol(l,bi,2)+1
                end if
            end do
        end do
        if(instr0%set%VSnorm)dist(:,1)=dist(:,1)/dist(:,2)
        if(.not.instr0%set%VSnorm)dist(:,1)=dist(:,1)/(bin*size(instr1%datam,2))
        !dist(:,1)=merge(dist(:,1)/dist(:,2),0._rk,MASK=dist(:,2)/=0)
        if(instr0%set%VSnorm)distmol(:,:,1)=distmol(:,:,1)/distmol(:,:,2)
        if(.not.instr0%set%VSnorm)distmol(:,:,1)=distmol(:,:,1)/(bin*size(instr1%datam,2))
        !distmol(:,:,1)=merge(distmol(:,:,1)/distmol(:,:,2),0._rk,MASK=distmol(:,:,2)/=0)
        !s=sum((distmol(:,bi,1)-dist(bi,1))**2/(dist(bi,2)-1._rk),MASK=distmol(:,bi,2)/=0)

        if(COUNT(instr0%set%calc(:)(1:7)=='distrib')>0)then
        do i=1,size(instr0%set%calc)
        select case(instr0%set%calc(i))
        case('distrib')
            filename=trim(instr0%set%filename)//'_'//trim(instr0%set%calc(i))//trim(instr0%set%filesuffix)
        if(instr0%set%ounit/=stdout)open(unit=instr0%set%ounit,file=trim(filename),status="unknown",iostat=ios)
            !Skriv ut distribution till disk
            do bi=1,size(dist,1)
                    x=mi+(real(bi,rk)-.5_rk)*bin
                    if(instr2%findex==1)x=acos(x)*180._rk/pi
                    write(instr0%set%ounit,*)x,dist(bi,1),merge(sqrt(sum((distmol(:,bi,1)-dist(bi,1))**2&
                    ,MASK=distmol(:,bi,2)/=0)/(real(COUNT(distmol(:,bi,2)/=0._rk),rk)-1._rk)),&
                    dist(bi,1),MASK=.not.all(distmol(:,bi,2)==0._rk))
            end do
            if(instr0%set%ounit/=stdout)close(instr0%set%ounit)
        case('distribm')
            filename=trim(instr0%set%filename)//'_'//trim(instr0%set%calc(i))//trim(instr0%set%filesuffix)
            if(instr0%set%ounit/=stdout)open(unit=instr0%set%ounit,file=trim(filename),status="unknown",iostat=ios)
            do imol=1,size(distmol,1)
            do bi=1,size(dist,1)
                    x=mi+(real(bi,rk)-.5_rk)*bin
                    if(instr2%findex==1)x=acos(x)*180._rk/pi
                    write(instr0%set%ounit,*)x,imol,distmol(imol,bi,1)
            end do
                write(instr0%set%ounit,*)
                write(instr0%set%ounit,*)
            end do
            if(instr0%set%ounit/=stdout)close(instr0%set%ounit)
        end select
        end do
        end if
        deallocate(dist,distmol)
    end subroutine vs_distrib!}}}

    subroutine center_of_membrane(molecules)!{{{
        integer(kind=ik) :: molecules(:),i,imol,j,k,l
        real(kind=rk),allocatable,save :: centers(:,:)
        real(kind=rk) :: n,shift(3),com(3)
        ! Calculate center of membrane based on mass density of
        ! all atoms between atom a and b
        if(.not.allocated(centers))then
            allocate(centers(3,sum(molt(molecules(:))%nmol)))
            k=0
            do j=1,size(molecules)
                i=molecules(j)
                do imol=1,molt(i)%nmol
                    k=k+1
                    centers(:,k)=center_of_molecule(i,imol)!com-shift
                end do
            end do
            !centers=centers*10._rk !SCALE INITFILE INSTEAD!
        end if
        centerofmembrane=0;n=0;k=0;shift=0
        do j=1,size(molecules)
           i=molecules(j)
            do imol=1,molt(i)%nmol
                k=k+1
                com=center_of_molecule(i,imol)
                !if(allocated(molt(i)%upper))shift(3)=(com(3)-centers(3,k))-mymodulo(com(3)-centers(3,k),box(3))
                shift(:)=(com(:)-centers(:,k))-mymodulo(com(:)-centers(:,k),box(:))
                do l=molt(i)%firstatom,molt(i)%lastatom
                    atom(l)%coor(:,imol)=atom(l)%coor(:,imol)-shift(:)
                end do
                !if(allocated(molt(i)%upper))centers(:,k)=center_of_molecule(i,imol)!com-shift
                centers(:,k)=center_of_molecule(i,imol)!com-shift
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
        ! Fold back atoms that are outside the box
        do i=1,size(molt)
            do imol=1,molt(i)%nmol
                com=atom(molt(i)%firstatom)%coor(:,imol) !center_of_molecule(i,imol) !Ändrade molt(i)%firstatom till i
                do k=molt(i)%firstatom,molt(i)%lastatom
                   ! coor(:,cind(k,imol))=com+mymodulo(coor(:,cind(k,imol))-com(:),box(:))
                   atom(k)%coor(:,imol)=com+mymodulo(atom(k)%coor(:,imol)-com(:),box(:))
                end do
            end do
        end do
    end subroutine whole!}}}

   subroutine foldmol(a)!{{{
        integer(kind=ik) :: i,imol,j,k
        real(kind=rk) :: com(3),shift
        character(len=*),optional,intent(in) :: a

        do i=1,size(molt)
            do imol=1,molt(i)%nmol
                if(present(a))then ! Fold and shift the molecules with respect to centerofmembrane
                    com=center_of_molecule(i,imol)-(centerofmembrane) !box/2 till box
                    do j=1,3
                        shift=0
                        if(ALL(i/=common_setflags%membrane_moltypes))shift=com(j)-mymodulo(com(j),box(j))
                        shift=shift+centerofmembrane(j)
                        do k=molt(i)%firstatom,molt(i)%lastatom
                          !  coor(j,cind(k,imol))=coor(j,cind(k,imol))-shift
                          atom(k)%coor(j,imol)=atom(k)%coor(j,imol)-shift
                        end do
                    end do
                else ! Fold back molecules that are outside the box
                    com=center_of_molecule(i,imol) !Ändrade molt(i)%firstatom till i
                    do j=1,3
                        shift=com(j)-modulo(com(j),box(j))
                        if(abs(shift)>box(j)/2._rk)then
                            do k=molt(i)%firstatom,molt(i)%lastatom
                          !      coor(j,cind(k,imol))=coor(j,cind(k,imol))-shift
                                 atom(k)%coor(j,imol)=atom(k)%coor(j,imol)-shift
                            end do
                        end if
                    end do
                    if(.not.present(a))then
                        com=center_of_molecule(i,imol)
                        if(com(3)<0 .OR. com(3)>box(3))then
                            write(*,*)com,'COM'
                            stop
                        end if
                    end if
                end if
            end do
        end do
        if(present(a))centerofmembrane=0._rk
    end subroutine foldmol!}}}

    function totalmoi(umol,imol) result(m)!{{{
        integer(kind=ik) :: umol,imol,i
        real(kind=rk) :: com(3),m
        m=0
        com=center_of_molecule(umol,imol)
        do i=molt(umol)%firstatom,molt(umol)%lastatom
            m=m+sum((atom(i)%coor(:,imol)-com)**2)*(masses(i))
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
    
        ri=atom(i)%coor(:,imol)-com
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
        real(kind=rk) :: eaverage(lbound(instr%datam,2):ubound(instr%datam,2))
        integer(kind=ik) :: frame,i,j,XYZ
        select case(instr%set%calc(j))
        case('traj')
            do frame=lbound(instr%datam,2),ubound(instr%datam,2)
           ! write(*,*)instr%datam(1,frame+minframe),"frame+minframe"
           ! write(*,*)instr%datam(1,frame),"frame"
            select case(instr%findex)
            case(14)! A try to average coordinates in the periodic box 
                select case(instr%instructionstring(1:1))
                case('X')
                    XYZ=1
                case('Y')
                    XYZ=2
                case('Z')
                    XYZ=3
                end select
                eaverage(frame)=&
                modulo(aimag(log(sum(exp(cmplx(0,instr%datam(:,frame)&
                *(2*pi/box(XYZ)))))/size(instr%datam,1)))*(box(XYZ)/(2*pi)),box(XYZ))
                !eaverage(frame)=&
                !abs(sum(exp(cmplx(0,instr%datam(:,frame)&
                !*(2*pi/box(XYZ)))))/size(instr%datam,1))!*(box(XYZ)/(2*pi))
            case default
                eaverage(frame)=sum(instr%datam(:,frame))/size(instr%datam,1)
            end select
            end do
        
            do frame=lbound(instr%datam,2),ubound(instr%datam,2)
                write(instr%set%ounit,*)frame,eaverage(frame),std(instr%datam(:,frame))
            end do
        case('trajm')
            do i=1,size(instr%datam,1)
                do frame=lbound(instr%datam,2),ubound(instr%datam,2)
                    write(instr%set%ounit,*)frame,i,instr%datam(i,frame)
                end do
                write(instr%set%ounit,*)
            end do
        end select
    end subroutine traj!}}}

    subroutine postproc(instr)!{{{
        type(instruct) :: instr(:)
        integer(kind=ik) :: ounit,ios,i,j,k,m,n,corr(6)!1,corr2
        real(kind=rk) :: a,mean,meandev,var
        character(kind=1,len=len(global_setflags%filename)) :: filename
        character(len=24) :: string
        character(len=300) :: string2
        logical :: exists
        character(kind=1,len=100) :: pos,avfilename
        meanbox=meanbox/(maxframes-skipframes)
        a=0._rk
        instr%cv%entropy=a/a
        n=0
        do i=1,size(instr)
            if(instr(i)%findex/=0)n=max(len(trim(instr(i)%instructionstring)),n)
            if(instr(i)%findex==11)instr(i)%datam=instr(i)%datam*(meanbox(1)*meanbox(2)) !APL
        end do
        do i=1,size(instr) ! Loop över instruktionsrader
        select case(instr(i)%findex)
            case(0,10)
            case(7,19,20) ! CORRELATE,VERSUS, ROTCORR
                corr(1)=strvecindex(instr(:)%ref,instr(i)%set%corrindex(1))
                corr(2)=strvecindex(instr(:)%ref,instr(i)%set%corrindex(2))
                if(instr(i)%findex==20)then !ROTCORR
                    corr(3)=strvecindex(instr(:)%ref,instr(i)%set%corrindex(3))
                    corr(4)=strvecindex(instr(:)%ref,instr(i)%set%corrindex(4))
                    corr(5)=strvecindex(instr(:)%ref,instr(i)%set%corrindex(5))
                    corr(6)=strvecindex(instr(:)%ref,instr(i)%set%corrindex(6))
                end if
                if(instr(i)%set%filename/='')then
                    instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%set%filename)!//trim(instr(i)%set%filesuffix)
                    filename=instr(i)%set%filename
                    instr(i)%set%ounit=nunit
                else
                    if(instr(i)%set%autofilename)then
                        select case(instr(i)%findex)
                        case(7,19)
                        instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%instructionstring)&
                                //trim(instr(corr(1))%instructionstring)//'_'//&
                                trim(instr(corr(2))%instructionstring)!//trim(instr(i)%set%filesuffix)
                        case(20)
                        instr(i)%set%filename=trim(instr(i)%set%fileprefix)//trim(instr(i)%instructionstring)&
                                //trim(instr(corr(1))%instructionstring)//'_'//&
                                trim(instr(corr(4))%instructionstring)!//trim(instr(i)%set%filesuffix)
                        end select
                        filename=instr(i)%set%filename
                        instr(i)%set%ounit=nunit
                    else
                        instr(i)%set%ounit=stdout
                        filename="stdout"
                    end if
                end if
                if(instr(i)%findex==7)then !CORRELATE
                    call corr_distrib(instr(i),instr(corr(1)),instr(corr(2)))
                    avfilename=trim(instr(i)%set%fileprefix)//'averages'//trim(instr(i)%set%filesuffix)
                    inquire(file=avfilename,exist=exists)
                    if(.NOT. exists)then
                        pos='ASIS'
                    else
                        pos='APPEND'
                    end if
                    open(78,file=trim(avfilename)&
                    ,position=trim(pos),status='unknown')
                    write(78,*)int(i,2),filename(len_trim(instr(i)%set%fileprefix)+1:len_trim(filename))," ",&
                    instr(i)%cv%pearsoncoeff,instr(i)%cv%entropy,instr(i)%cv%entropymutual
                    close(78)
                end if
                !!!!
                if(instr(i)%findex==19)call vs_distrib(instr(i),instr(corr(1)),instr(corr(2))) !VERSUS
                !!!!
                if(instr(i)%findex==20)then !ROTCORR
                    call rotcorr(instr(i),instr(corr(1)),instr(corr(2)),instr(corr(3)),&
                    instr(corr(4)),instr(corr(5)),instr(corr(6)))
                end if

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
                case(0,7,10,19,20,21)
                case default
                    if(allocated(instr(i)%set%calc))then
                    do j=1,size(instr(i)%set%calc) ! Loop över postberäkningar
                        write(0,'(5X,A)',advance="no")'Instruction '//trim(adjustl(intstr(i)))//&
                        ' / '//trim(adjustl(intstr(size(instr))))
                        if(instr(i)%set%autofilename)then
                            filename=trim(instr(i)%set%filename)&
                            //'_'//trim(instr(i)%set%calc(j))//trim(instr(i)%set%filesuffix)
                        else 
                            filename=trim(instr(i)%set%filename)//trim(instr(i)%set%filesuffix)
                        end if
                        if(instr(i)%set%ounit/=stdout)open(unit=instr(i)%set%ounit,file=trim(filename),status="unknown",iostat=ios)
                        select case(trim(instr(i)%set%calc(j)))
                            case('distrib','distribm')
                                write(0,'(A10)',advance="no")'distrib'
                                write(0,'(A)',advance="no")' '//trim(instr(i)%instructionstring)
                                !write(*,*)'I',i,instr(i)%datam
                                if(instr(i)%findex/=13.AND.instr(i)%findex/=16.AND.instr(i)%findex/=15)call distrib(instr(i),i,j)
                                if(instr(i)%findex==13.AND.trim(instr(i)%set%calc(j))=='distrib')&
                                call rdf_distrib(instr(i))
                                if(instr(i)%findex==16.OR.instr(i)%findex==15.AND.trim(instr(i)%set%calc(j))=='distrib')&
                                call write_distrib(instr(i))
                            case('traj','trajm')
                                write(0,'(A10)',advance="no")'traj'
                                write(0,'(A)',advance="no")' '//trim(instr(i)%instructionstring)
                                if(instr(i)%findex/=16)call traj(instr(i),j)
                            case('acorr','acorrm')
                                write(0,'(A10)',advance="no")'acorr'
                                write(0,'(A)',advance="no")' '//trim(instr(i)%instructionstring)
                                if(instr(i)%findex/=13.AND.instr(i)%findex/=16)call autocorr(instr(i),i,j)
                        end select
                        m=len_trim(instr(i)%instructionstring)&
                        +len_trim('Instruction '//trim(adjustl(intstr(i)))//' / '//trim(adjustl(intstr(size(instr)))))+16

                        write(0,'(A)',advance="no")char(13)
                        write(0,'(A'//trim(adjustl(intstr(m)))//')',advance="no")' '!('A',k=1,m)
                        write(0,'(A)',advance="no")char(13)
                        if(instr(i)%set%ounit/=stdout)close(instr(i)%set%ounit,iostat=ios)
                    end do
                    end if
                    select case(instr(i)%findex)
                    case(0,7,10,13,19,20)
                    case default
                        inquire(file=avfilename,exist=exists)
                        if(.NOT. exists)then
                            pos='ASIS'
                        else
                            pos='APPEND'
                        end if
                        open(78,file=trim(avfilename)&
                        ,position=trim(pos),status='unknown')
                    end select

                    select case(instr(i)%findex)
                        case(0,7,10,13,19,20)
                        case(1)
                            write(78,*)int(i,2),instr(i)%instructionstring(1:n)," ",&
                            trim(adjustl(average(acos(instr(i)%datam)*180._rk/pi)))!,&
                            close(78)
                        case(16)
                            mean=instr(i)%cv%mean/instr(i)%cv%n
                            var=(1/(real(instr(i)%cv%n,rk)-1))*(instr(i)%cv%meandev-&
                            (1/real(instr(i)%cv%n,rk))*(instr(i)%cv%mean)**2)
                            meandev=sqrt(var/real(instr(i)%cv%n,rk))
                            string=getmeanwithdev(mean,meandev)
                            write(string2,*)mean,meandev,sqrt(var)
                            string2=trim(adjustl(string))//" "//trim(string2)
                            write(78,*)int(i,2),instr(i)%instructionstring(1:n)," ",&
                            trim(string2)
                            close(78)
                        case default
                            write(78,*)int(i,2),instr(i)%instructionstring(1:n)," ",&
                            trim(adjustl(average(instr(i)%datam)))!,&
                            !!SLICE
                            if(instr(i)%set%slice%switch)then
                            write(78,*)int(i,2),instr(i)%instructionstring(1:n)//"_SLICE_"//&
                            realstr(instr(i)%set%slice%lower)//"-"//realstr(instr(i)%set%slice%upper)," ",&
                            trim(adjustl(slice_average(instr(i)))) !,&
                            end if
                            !if(instr(i)%molaverage)then
                            !    do k=1,size(instr(i)%datam,1)
                            !        write(42,*)k,instr(i)%instructionstring(1:n),sum(instr(i)%datam(k,:))/size(instr(i)%datam,2)
                            !    end do
                            !end if
                            close(78)
                    end select
            end select
            if(allocated(instr(i)%datam))deallocate(instr(i)%datam)
        end do
        write(*,*)
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

subroutine rotcorr(instr,xa,ya,za,xb,yb,zb)!{{{
    type(instruct),intent(inout) :: instr,xa,ya,za,xb,yb,zb
    real(kind=rk),allocatable :: dataa(:,:,:),datab(:,:,:),acorr(:,:)
    real(kind=rk) :: norm
    integer(kind=ik) :: imol,calcind,i,ios
    character(kind=1,len=len(global_setflags%filename)) :: filename
    allocate(dataa(3,lbound(xa%datam,2):ubound(xa%datam,2),size(xa%datam,1)))
    allocate(datab(3,lbound(xb%datam,2):ubound(xb%datam,2),size(xb%datam,1)))
    allocate(acorr(size(xa%datam,1),0:size(xa%datam,2)-1))
    dataa(1,:,:)=transpose(xa%datam(:,:))
    dataa(2,:,:)=transpose(ya%datam(:,:))
    dataa(3,:,:)=transpose(za%datam(:,:))
    datab(1,:,:)=transpose(xb%datam(:,:))
    datab(2,:,:)=transpose(yb%datam(:,:))
    datab(3,:,:)=transpose(zb%datam(:,:))
    !write(*,*)lbound(xa%datam,2),lbound(dataa,2)
    if(size(xa%datam,1)/=size(xb%datam,1))then
        write(*,*)'WARNING RC: NMOLa /= NMOLb'
        write(*,*)instr%set%filename,' <-- Doing nothing!!'
        return
    else if(size(xa%datam,2)/=size(xb%datam,2))then
        write(*,*)'WARNING RC: NFRAMESa /= NFRAMESb'
        write(*,*)instr%set%filename,' <-- Doing nothing!!'
        return
    end if
        
    do imol=1,size(xa%datam,1)
        call rotacf(dataa(:,:,imol),datab(:,:,imol),acorr(imol,0:))
        !norm=1._rk/acorr(imol,0)
        !acorr(imol,:)=acorr(imol,:)*norm
    end do
    do calcind=1,size(instr%set%calc)
        select case(instr%set%calc(calcind))
        case('acorr')
            filename=trim(instr%set%filename)//'_'//trim(instr%set%calc(calcind))//trim(instr%set%filesuffix)
            if(instr%set%ounit/=stdout)open(unit=instr%set%ounit,file=trim(filename),status="unknown",iostat=ios)

            do i=0,ubound(acorr,2)
                write(instr%set%ounit,*)i,sum(acorr(:,i))/real(size(acorr,1),rk),std(acorr(:,i))
            end do
            if(instr%set%ounit/=stdout)close(instr%set%ounit)
        case('acorrm')
            filename=trim(instr%set%filename)//'_'//trim(instr%set%calc(calcind))//trim(instr%set%filesuffix)
            if(instr%set%ounit/=stdout)open(unit=instr%set%ounit,file=trim(filename),status="unknown",iostat=ios)
            do imol=1,size(dataa,3)
                do i=0,ubound(acorr,2)
                    write(instr%set%ounit,*)i,imol,acorr(imol,i)!,std(acorr(:,i))
                end do
                write(instr%set%ounit,*)
            end do
            if(instr%set%ounit/=stdout)close(instr%set%ounit)
         end select
     end do
     contains
     subroutine rotacf(dataam,databm,acorrv)
         real(kind=rk) :: dataam(:,:),databm(:,:),datav(3,size(dataam,2)),acorrv(0:)
         integer(kind=ik) :: i,j
         acorrv=0
         datav=0
         do i=lbound(dataam,2),ubound(dataam,2)
             datav(:,i)=normalize(databm(:,i)-dataam(:,i))
         end do
             do j=0,size(datav,2)-1
                 do i=1,size(datav,2)-j
                     !write(*,*)i,j,j+i,lbound(acorrv)
                     acorrv(j)=acorrv(j)+P2(datav(:,i),datav(:,j+i))
                 end do
                 acorrv(j)=acorrv(j)/real(size(datav,2)-j,rk)
             end do
     end subroutine rotacf
     function P2(a,b) result(c)
         real(kind=rk) :: a(:),b(:),c
         c=1.5_rk*dot_product(a,b)**2-.5_rk
         !write(*,*)c
     end function P2
     end subroutine rotcorr!}}}
end module trajop
