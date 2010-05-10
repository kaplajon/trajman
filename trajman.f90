program trajman
    use util
    use kinds
    use input
    use readtraj
    use trajop
    implicit none
    character(kind=1,len=1),allocatable :: inputrad(:),words(:,:)
    character(kind=1,len=255) :: sysstr,pid,wf,filename
    integer(kind=ik) :: i,j,k,l,m,n,ios,p,frame,runit
    integer(kind=8) :: trs,frs
    logical,allocatable :: logicmolt(:)
    type(instruct),allocatable :: troptype(:)
    call arguments(runit)
    call constants
    call globals
    i=0
    do !Input processing 
        call readline(inputrad,ios,runit)
        !write(*,*)inputrad
        if(ios==endf)exit
        i=i+1
        call reallocinstruct(troptype,i)
        call procinp(inputrad,troptype(i))
    end do
    troptype(:)%nmolop=0
    do i=1,size(troptype)
        if(troptype(i)%findex/=0)then!troptype(i)%nmolop=molt(moltypeofuatom(troptype(i)%atoms(1)))%nmol
            allocate(logicmolt(size(molt)));logicmolt=.FALSE.
            do j=1,size(troptype(i)%atoms)
                logicmolt(moltypeofuatom(troptype(i)%atoms(j)))=.TRUE.
            end do
            troptype(i)%nmolop=sum(molt(:)%nmol,MASK=logicmolt)
            deallocate(logicmolt)
        end if
    end do

    ! Allocate data matrix to a size based on input
    if(maxframes==0)then !maxframes=1001
        p=getpid() ! Returns pid of the running instance of trajman
        write(pid,*)p
        !Make a system call to stat to get filesize of trajectory file
        write(sysstr,*)'stat -c%s '//stringconv(trajfile)//' >',trim(adjustl(pid)),'trsize'
        ios=system(trim(sysstr))
        !Make a system call to extract one frame from trajectory file (based on
        !rowsperframe)
        write(sysstr,*)'tail ',-rowsperframe,' ',stringconv(trajfile),' >',trim(adjustl(pid)),'oneframe'
        ios=system(trim(sysstr))
        ! Make a system call to get the size of one frame
        sysstr='stat -c%s '//trim(adjustl(pid))//'oneframe >>'//trim(adjustl(pid))//'trsize'
        ios=system(trim(sysstr))
        !Get sizes and calculate number of frames
        open(42,file=trim(adjustl(pid))//'trsize')
        read(42,*)trs
        read(42,*)frs
        maxframes=(trs-mod(trs,frs))/frs+nint(real(mod(trs,frs),8)/real(frs,8))
        if(trs==0)maxframes=0
        close(42)
        ! Remove temporary files
        ios=system('rm '//trim(adjustl(pid))//'trsize '//trim(adjustl(pid))//'oneframe')
    end if

    do i=1,size(troptype)
        if(troptype(i)%nmolop*maxframes/=0)then
            allocate(troptype(i)%datam(troptype(i)%nmolop,maxframes))
            troptype(i)%datam=0
        end if
    end do

    write(*,*)"input done"
    write(*,'(5X,A19,I4)')'Frames to process: ',maxframes
!#######################################################################################    
    frame=0
    do while (readgro(tunit)==0)! Trajectory processing
        frame=frame+1
        if(modulo(frame-1,max(maxframes/100,1))==0)then
            write(*,'(5X,A10,I3,2A)',advance='no')'Progress: ',nint(real(100*frame,rk)/real(maxframes,rk)),'%',char(13)
        end if
        if(global_setflags%folding)call foldmol
        if(allocated(common_setflags%membrane_moltypes))&
        call center_of_membrane(common_setflags%membrane_moltypes)
        if(global_setflags%apl)call apl_grid
        call procop(troptype,frame) ! Perform instructions on frame

        if(frame==1)then !Write atomnames and coordinates for the first molecules
                        !to a .xyz file
            open(37,file='atoms.xyz')
            write(37,*)size(atomnames)
            write(37,*)
            do i=1,size(atomnames)
                    write(37,*)atomnames(i),10*coor(1:3,cind(atomindex(atomnames(i),atomnames,size(atomnames)),1_ik))
            end do
            close(37)
        end if

        do m=1,global_setflags%wftot !BEGIN WRITEFRAME: Write specified frame to file
            if(allocated(global_setflags%writeframe) .AND. frame==global_setflags%writeframe(m)%framenumber)then
                write(wf,*)global_setflags%writeframe(m)%framenumber
                
                if(trim(global_setflags%writeframe(m)%outformat)=='gro')then
                    filename='frame'//trim(adjustl(wf))//'.gro'
                    call wf_gro(trim(filename),wf,37)
                else ! Defaults to writing a .xyz file with long atomnames,
                     ! viewable by gOpenmol.
                    filename='frame'//trim(adjustl(wf))//'.xyz'
                    call wf_xyz(trim(filename),37)
                end if

            end if
        end do !END WRITEFRAME       
        if (frame==maxframes)exit
        end do
!#######################################################################################    
    write(*,*)
    write(*,'(5x,A)')'Postprocessing...'
    call postproc(troptype)
    write(*,'(5X,A5)')'Done!'
end program trajman
