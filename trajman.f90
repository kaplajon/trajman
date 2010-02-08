program trajman
!TESTTESTTEST
    use util
    use kinds
    use input
    use readtraj
    use trajop
    implicit none
    character(kind=1,len=1),pointer :: inputrad(:),words(:,:)
    character(kind=1,len=255) :: sysstr,pid,wf
    integer(kind=ik) :: i,j,k,m,n,ios,p,frame,runit
    integer(kind=8) :: trs,frs
    type(instruct),allocatable :: troptype(:)
    call arguments(runit)
    call constants
    call globals
    i=0
    do !Input processing 
        call readline(inputrad,ios,runit)   
        if(ios==endf)exit
        i=i+1
        call reallocinstruct(troptype,i)
        call procinp(inputrad,troptype(i))
    end do
    ! Maximala antalet molekyler som processas per frame är maxmols
    !allocate(nmolop(size(trop(1,:))))
    troptype(:)%nmolop=0
    do i=1,size(troptype)
        if(troptype(i)%findex/=0)troptype(i)%nmolop=molt(moltypeofuatom(troptype(i)%atoms(1)))%nmol
    end do
!    maxmols=maxval(troptype(:)%nmolop)

!write(*,fmt='(1A11)')(atomnames(i),i=1,size(atomnames));stop

    ! Allokera datamatris till storlek baserad på input
    if(maxframes==0)then !maxframes=1001
        p=getpid() ! Returns pid of the running instance of trajman
        write(pid,*)p
        !Make a system call to stat to get filesize of trajectory file
        write(sysstr,*)'stat -c%s '//stringconv(trajfile)//' >',trim(adjustl(pid)),'trsize'
        !write(*,*)sysstr
        !call fstat(tunit,buff)
        !write(*,*)buff(8),'buffert(8)'
        call system(trim(sysstr))
        !Make a system call to extract one frame from trajectory file (based on
        !rowsperframe)
        write(sysstr,*)'tail ',-rowsperframe,' ',stringconv(trajfile),' >',trim(adjustl(pid)),'oneframe'
        call system(trim(sysstr))
        ! Make a system call to get the size of one frame
        sysstr='stat -c%s '//trim(adjustl(pid))//'oneframe >>'//trim(adjustl(pid))//'trsize'
        call system(trim(sysstr))
        !Get sizes and calculate number of frames
        open(42,file=trim(adjustl(pid))//'trsize')
        read(42,*)trs
        read(42,*)frs
        !write(*,*)trs,frs
        !stopa
        maxframes=(trs-mod(trs,frs))/frs+nint(real(mod(trs,frs),8)/real(frs,8))
        if(trs==0)maxframes=0
        !maxframes=nint(real(trs,8)/real(frs,8))
        close(42)
        ! Remove temporary files
        call system('rm '//trim(adjustl(pid))//'trsize '//trim(adjustl(pid))//'oneframe')
    end if

    do i=1,size(troptype)
        if(troptype(i)%nmolop*maxframes/=0)then
            allocate(troptype(i)%datam(troptype(i)%nmolop,maxframes))
            troptype(i)%datam=0
        end if
    end do

    write(*,*)"input done"
    write(*,'(5X,A19,I4)')'Frames to process: ',maxframes
    frame=0
    do while (readgro(tunit)==0)! Trajectory processing
       frame=frame+1
       if(modulo(frame-1,max(maxframes/100,1))==0)then
           write(*,'(5X,A10,I3,2A)',advance='no')'Progress: ',nint(real(100*frame,rk)/real(maxframes,rk)),'%',char(13)
       end if

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

       if(frame==global_setflags%writeframe)then
           write(wf,*)global_setflags%writeframe
           open(37,file='frame'//trim(adjustl(wf))//'.xyz')
           write(37,*)atot
           write(37,*)
            do i=1,size(molt)
                do j=1,molt(i)%nmol
                    do k=molt(i)%firstatom,molt(i)%lastatom
                        write(37,*)atomnames(k),10*getatom(k,j)
                    end do
                end do
            end do
            close(37)
        end if
                
      
      ! flush(6)
       call procop(troptype,frame) ! Perform instructions on frame
       if (frame==maxframes)exit
    end do
    write(*,*)
    write(*,'(5x,A)')'Postprocessing...'
!            do i=1,size(atomnames)
!        write(37,*)atomnames(i),10*coor(1:3,cind(atomindex(atomnames(i),atomnames,atoms),1_ik))
!        end do
    call postproc(troptype)
    write(*,'(5X,A5)')'Done!'
end program trajman
