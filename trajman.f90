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
program trajman
    use util
    use kinds
    use input
    use readtraj
    use trajop
    use statistics
    implicit none
    character(kind=1,len=1),allocatable :: inputrad(:),words(:,:)
    character(kind=1,len=255) :: sysstr,pid,wf,filename
    integer(kind=ik) :: i,j,k,l,m,n,ios,p,frame
    integer(kind=4) :: runit
    integer(kind=8) :: trs,frs
    integer(kind=4) :: natm,statf
    real(kind=4),allocatable :: coorv(:),bx(:) !MOLFILEPLUGIN FROM VMD
    type(instruct),allocatable :: troptype(:)
    integer(kind=ik) :: v(1:4),vs(1:4)
    call RNSTRT(310952_ikr) !init. random seed
    skipframes=0;maxframes=0
    atomsdefined=0
    call arguments(runit)
    call constants
    call globals
    i=0
    allocate(inputrad(1))
    do !Input processing 
        call readline(inputrad,ios,runit)
        if(ios==endf)exit
        i=i+1
        call reallocinstruct(troptype,i)
        call procinp(inputrad,troptype(i))
    end do
    troptype(:)%nmolop=0
    do i=1,size(troptype)
        select case(troptype(i)%findex)
        case(0,7,9,10,15,16,17,18,19,20,23,24)
        case(13,22) !RDF
            j=atom(troptype(i)%atoms(1))%moltype
            k=atom(troptype(i)%atoms(2))%moltype
            select case(troptype(i)%set%leaflet)
            case(0)!Both
                l=molt(j)%nmol*molt(k)%nmol
            case(1)!Lower
                l=size(molt(j)%lower)*size(molt(k)%lower)
            case(2)!Upper
                l=size(molt(j)%upper)*size(molt(k)%upper)
            end select
            !if(j==k)l=nint(sqrt(real(l,rk))*(sqrt(real(l,rk))-1))/2
            if(j==k)l=nint(sqrt(real(l,rk))*(sqrt(real(l,rk))-1)/2)
            if(troptype(i)%findex==13)allocate(troptype(i)%rdf_pairs(l))
            if(troptype(i)%findex==22)troptype(i)%nmolop=l
           ! write(*,*)troptype(i)%nmolop,troptype(i)%findex,j,k,molt(1)%nmol
           ! stop

                
        case default
            j=atom(troptype(i)%atoms(1))%moltype
            select case(troptype(i)%set%leaflet)
                case(0)
                    call reallocate(troptype(i)%molind,molt(j)%nmol)
                    do k=1,molt(j)%nmol
                        troptype(i)%molind(k)=k
                    end do
                case(1)
                    call reallocate(troptype(i)%molind,size(molt(j)%lower))
                    troptype(i)%molind=molt(j)%lower
                case(2)
                    call reallocate(troptype(i)%molind,size(molt(j)%upper))
                    troptype(i)%molind=molt(j)%upper
            end select
            troptype(i)%nmolop=size(troptype(i)%molind)
        end select
    end do
    ! Allocate data matrix to a size based on input
    if(maxframes==0)then
        select case(trajtype)
        case('gro')
            !Sum the number of frames of all trajfiles
            j=0
            do i=1,size(trajfile)
            p=getpid() ! Returns pid of the running instance of trajman
            write(pid,*)p
            !Make a system call to stat to get filesize of trajectory file
            write(sysstr,*)'stat -c%s '//stringconv(trajfile(i)%filename)//' >',trim(adjustl(pid)),'trsize'
            ios=system(trim(sysstr))
            !Make a system call to extract one frame from trajectory file (based on
            !rowsperframe)
            write(sysstr,*)'tail ',-rowsperframe,' ',stringconv(trajfile(i)%filename),' >',trim(adjustl(pid)),'oneframe'
            ios=system(trim(sysstr))
            ! Make a system call to get the size of one frame
            sysstr='stat -c%s '//trim(adjustl(pid))//'oneframe >>'//trim(adjustl(pid))//'trsize'
            ios=system(trim(sysstr))
            !Get sizes and calculate number of frames
            open(42,file=trim(adjustl(pid))//'trsize')
            read(42,*)trs
            read(42,*)frs
            j=((trs-mod(trs,frs))/frs+nint(real(mod(trs,frs),8)/real(frs,8)))
            maxframes=maxframes+j
            if(trs==0)maxframes=0
            close(42)
            ! Remove temporary files
            ios=system('rm '//trim(adjustl(pid))//'trsize '//trim(adjustl(pid))//'oneframe')
            if(.NOT. common_setflags%silent)write(0,"(A,5x,A15,I6)",advance="no")&
           char(13),"Counting frames: ",j!,char(13)     
            write(*,*)stringconv(trajfile(i)%filename)
            end do
            
        case('trr','xtc','dcd','pdb')!VMD MOLFILEPLUGIN
            statf=1
            allocate(coorv(atot*3),bx(6))
            open(42,file="framecount.dat")
            do i=1,size(tunit)
                j=0
                do 
                    call f77_molfile_read_next(tunit(i),int(atot,4),coorv,bx,statf)
                    if(.NOT. common_setflags%silent)write(0,"(A,5x,A15,I6)",advance="no")&
                    char(13),"Counting frames: ",j!,char(13)     
                    write(42,"(A,5x,A15,I6)",advance="no")&
                    char(13),"Counting frames: ",j!,char(13)     
                    if(statf==0)exit
                    maxframes=maxframes+1
                    j=j+1
                end do
                if(.NOT. common_setflags%silent)write(*,*)stringconv(trajfile(i)%filename)
                call f77_molfile_close_read(tunit(i),statf)
                call f77_molfile_open_read(tunit(i),natm,stringconv(trajfile(i)%filename),trajtype)
            end do
            write(42,*)"FRAMES: ",trim(intstr(maxframes))
            close(42)
            deallocate(coorv,bx)
         case(dlpoly3histfile(1:min(len(dlpoly3histfile),len(trajtype))))
            statf=1
            maxframes=0
            do i=1,size(tunit)
               call inithist(trajfile(i)%filename,j)
               maxframes=maxframes+j
               if(.NOT. common_setflags%silent)then
                  write(0,"(A,5x,A15,I6)",advance="no")&
                       char(13),"Counting frames: ",j
                  write(0,*)stringconv(trajfile(i)%filename)
               end if
            end do
        case default
            statf=1
            allocate(coorv(atot*3),bx(6))
            do i=1,size(tunit)
                j=0
                do 
                    call f77_molfile_read_next(tunit(i),int(atot,4),coorv,bx,statf)
                    if(.NOT. common_setflags%silent)write(0,"(A,5x,A15,I6)",advance="no")&
                    char(13),"Counting frames: ",j!,char(13)     
                    if(statf==0)exit
                    maxframes=maxframes+1
                    j=j+1
                end do
                if(.NOT. common_setflags%silent)write(*,*)stringconv(trajfile(i)%filename)
                call f77_molfile_close_read(tunit(i),statf)
                call f77_molfile_open_read(tunit(i),natm,stringconv(trajfile(i)%filename),'auto')
            end do
            deallocate(coorv,bx)
        end select
    end if
    do i=1,size(troptype)
        if(troptype(i)%nmolop*maxframes/=0)then
            allocate(troptype(i)%datam(troptype(i)%nmolop,skipframes+1:maxframes))
            troptype(i)%datam=0
        else
            select case(troptype(i)%findex)
            case(13)
                allocate(troptype(i)%rdf_dist(troptype(i)%set%distbin))
                allocate(troptype(i)%datam(1,skipframes+1:maxframes))
                troptype(i)%rdf_dist=0
            case(15)
                allocate(troptype(i)%datam(1,skipframes+1:maxframes))
                troptype(i)%datam=0
            case(16)
                !allocate(troptype(i)%rdf_dist(troptype(i)%set%distbin))
                !troptype(i)%rdf_dist=0
                !troptype(i)%rdf_bin=1._rk/troptype(i)%set%distbin
            case(17,18)
                allocate(troptype(i)%datam(1,skipframes+1:maxframes))
                troptype(i)%datam=0
            end select
        end if
    end do
    !write(*,*)size(troptype(18)%datam,2),'DATAM';stop
    if(.NOT. common_setflags%silent)then
        write(*,*)"input processing done."
        call summary
        write(*,'(5X,A19,I6,A2,I6)')'Frames to process: '&
        ,maxframes-skipframes," /",maxframes
    end if
!###########################FRAME LOOP##################################################    
j=1;i=0
do !i=1,skipframes
    if(i==skipframes)exit
    i=i+1
    if(.NOT. common_setflags%silent)write(0,"(A,5x,A15,I6)",advance="no")&
    char(13),"Skipping frame: ",i!,char(13) 
    ios=readframe(tunit(j))
    if(ios==2)then
        j=j+1 ! Change file!
        i=i-1 ! Rewind the counter so that we get the last frame
    end if
end do
if(skipframes/=0.and..NOT.common_setflags%silent)write(*,*)&
"    Starting calculation on frame: ",trim(adjustl(intstr(skipframes+1)))
frame=skipframes
    do k=j,size(tunit) !Loop over trajectory files
    do while (readframe(tunit(k))/=2)! Trajectory processing, frameloop
   ! write(*,*)tunit(k),'TUNIT'
        frame=frame+1
        if(modulo(frame-skipframes,max((maxframes-skipframes)/100,1))==0.OR.&
        frame==skipframes+1)then
            if(.NOT. common_setflags%silent)then
                write(0,'(A,5X,A10,I3,A2,I6)',advance='no')&
                char(13),'Progress: ',nint(real(100*(frame-skipframes),rk)/real(maxframes-skipframes,rk)),'% ',frame!,char(13)
            end if
        end if
        if(global_setflags%whole)call whole
        if(global_setflags%folding)call foldmol
        meanbox(1:3)=meanbox(1:3)+box
        meanbox(4)=meanbox(4)+1/(box(1)*box(2)*box(3))
        call procop(troptype,frame) ! CALCULATIONS: Perform instructions on frame
        if(frame==skipframes+1)then !Write atomnames and coordinates for the first molecules
                        !to a .xyz file
            open(37,file='atoms.xyz')
            write(37,*)size(atomnames)
            write(37,*)
            do i=1,size(atomnames)
                select case(trajtype)
                case('gro')
                    write(37,*)trim(atom(i)%aname),atom(i)%coor(:,1)
                case('trr','dcd','pdb')
                    write(37,*)trim(atom(i)%aname),atom(i)%coor(:,1)
                case default
                    write(37,*)trim(atom(i)%aname),atom(i)%coor(:,1)
                end select
            end do
            close(37)
        end if

        do m=1,global_setflags%wftot !BEGIN WRITEFRAME: Write specified frame to file
            if(allocated(global_setflags%writeframe) .AND. frame==global_setflags%writeframe(m)%framenumber)then
                write(wf,*)global_setflags%writeframe(m)%framenumber
                
                select case(trim(global_setflags%writeframe(m)%outformat))
                case('gro')
                    filename='frame'//trim(adjustl(wf))//'.gro'
                    call wf_gro(trim(filename),wf,37)
                case('xyz') ! Defaults to writing a .xyz file with long atomnames,
                     ! viewable by gOpenmol.
                    filename='frame'//trim(adjustl(wf))//'.xyz'
                    call wf_xyz(trim(filename),37)
                end select
            end if
        end do !END WRITEFRAME       
        if (frame==maxframes)exit
        end do ! Frameloop
        if(frame==maxframes)exit
        end do ! Trajectory file loop
!#######################################################################################
    do i=1,size(tunit)
        call closetraj(tunit(i))
    end do
    if(.NOT. common_setflags%silent)then
        write(*,*)
        write(*,'(5x,A)')'Postprocessing...'
    end if
    call postproc(troptype) ! POSTPROCESSING
    if(.NOT. common_setflags%silent)write(*,'(5X,A)')'Done!'
end program trajman
