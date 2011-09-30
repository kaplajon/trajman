program trr_test
implicit none
 integer(kind=4) :: natom,handle,stat,i,j
 real(kind=4) :: box(6)
 real(kind=8),allocatable :: coor(:,:)
 real(kind=4),allocatable :: coorv(:)
 character(len=200) :: infile
 character(len=16) :: intype
handle=-1
natom=0
infile='tedd.st2.trr'
intype=trim(infile(scan(infile,'.',BACK=.TRUE.)+1:))
write(*,*)intype,scan(infile,'.'),scan(infile,'.',BACK=.TRUE.)
!stop
intype='trr'
 call f77_molfile_init
 call f77_molfile_open_read(handle,natom,trim(infile),trim(intype))
if(handle/=-1)then
write(*,*)'ATOMS: ',natom,' Handle: ',handle
else
stop 'CANNOT READ FILE'
end if
allocate(coorv(natom*3),coor(1:3,natom))
stat=1
j=0
do
 call f77_molfile_read_next(handle,natom,coorv,box,stat)
 call f77_molfile_read_next(handle,natom,coorv,box,stat)
if(stat==0)exit
!if(ioconv(stat)/=0)exit
j=j+1
coor=reshape(real(coorv,kind(coor)),[3,natom])
do i=1,size(coor,2)
write(20,*)coor(1:3,i)
end do
exit
!write(*,*)'Handle: ',handle,' STAT: ',stat
end do
write(*,*)'Handle: ',handle,' STAT: ',stat,' FRAMES: ',j,box
!coor=reshape(real(coorv,rk),[3,natom])
!do i=1,200
!write(*,*)coor(1:3,i)
!end do
 call f77_molfile_close_read(handle,stat)
write(*,*)'CLOSE STAT: ',stat
 call f77_molfile_finish
end program trr_test
