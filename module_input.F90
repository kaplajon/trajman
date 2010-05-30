module input
    use kinds
    use version
    use util
    use readtraj
    use apl
    use trajop
    implicit none ! Routines for input file treatment

    contains

    subroutine arguments(runit)!{{{
        integer(kind=ik) :: i,ios,runit
        character(kind=1,len=255) :: carg,carg2,ctime

        !If no arguments, print help info
        if(command_argument_count()==0)then
            call print_help
        end if

        i=1
        do !Loop over command line arguments 
           call get_command_argument(i,carg)
           if (len_trim(carg)==0)exit
           select case(trim(carg))
            case('-i','--input')
                call get_command_argument(i+1,carg2)
                if(len_trim(carg2)==0)stop 'Input flag needs a filename!'
                i=i+1
                runit=14
                open(runit,file=trim(carg2),status='old',iostat=ios)
                if(ios/=0)then
                    write(*,*)'Cannot open inputfile :',trim(carg2),':'
                    stop
                end if
            case('-v','--version')
                !ctime is set by the compiler using c type macro and cpp. Make
                !sure the file is suffixed with capital F or that the compiler
                !is flagged to use the preprocessor.
                ctime=CINFO
                write(*,*)'     Trajman experimental :: ',trim(ctime)
                write(*,*)'                   Branch :: ',branch
                write(*,*)'                 Revision :: r',revision,', committed ',revdate
                write(*,*)'     (c) Jon Kapla, 2009  :: Contact: jon.kapla@mmk.su.se'
                stop
            case('-h','--help')
               call print_help 
            case('-')
                runit=stdin

            case default
                write(*,*)':',trim(carg),':'
                write(*,*)'I cannot do that! Try "trajman -h" for more information.'
                stop
           end select
           i=i+1

        end do
    end subroutine arguments!}}}

    subroutine print_help!{{{
        character(len=100) :: pid
        integer(kind=ik) :: p
        write(*,*)
        write(*,*)'Usage: trajman [FLAGS] [FILE]'
        write(*,*)
        write(*,*)' FLAGS:'
        write(*,*)
        write(*,*)'  -v, --version'
        write(*,*)'          Print program version information.'
        write(*,*)'  -h, --help'       
        write(*,*)'          Print this help and quit.'
        write(*,*)'  -i, --input FILE'
        write(*,*)'          Use FILE as inputfile'
        write(*,*)'  -'         
        write(*,*)'          Read from standard input, STDIN'
        write(*,*)
        p=getpid()
        write(pid,*)p
        p=system("cat $(readlink /proc/"//trim(adjustl(pid))//"/exe)_documentation.txt")
        stop
    end subroutine print_help!}}}

!    subroutine readline(onerow,ios,runit)!{{{
!        character(kind=1,len=1),allocatable :: onerow(:) 
!        integer(kind=ik) :: n,ios,readunit
!        integer(kind=ik),optional,intent(in) :: runit
!        readunit=stdin
!        if (present(runit))readunit=runit
!            comment:do
!                n=0
!                do
!                    if(n+1>=size(onerow))call reallocate(onerow,int(2*size(onerow)+1,ik))
!                    read(readunit,fmt="(1A1)",advance='NO',iostat=ios)onerow(n+1)!oneletter
!                    if(ios==endf)exit comment
!                    if(onerow(n+1)=='#' .or. ios==endr )then
!                        if(n==0)then
!                            if(ios/=endr)read(readunit,fmt='()',advance='YES',iostat=ios)
!                            cycle comment
!                        else
!                            exit comment
!                        endif
!                    end if
!                n=n+1
!                enddo
!            end do comment
!        if(ios==endf)return       
!        n=max(n,1)
!       call reallocate(onerow,n)
!    end subroutine readline !}}}

    function wordcount(onerow) result(words)!{{{
        character(kind=1,len=1) :: onerow(:)
        integer :: i,words
        logical :: whitespace
        whitespace=.TRUE.
        words=0
        do i=1,size(onerow)
            select case(iachar(onerow(i)))
                case(iachar(' '),iachar(','),9,0) ! char(9) är tab
                    whitespace=.TRUE.
                    onerow(i)=" "
                case default
                    if(whitespace)words=words+1
                    whitespace=.FALSE.
            end select
        end do
    end function wordcount!}}}

    subroutine getwords(vector,words)!{{{
        character(len=1),allocatable :: vector(:)
        character(len=1),allocatable :: words(:,:)
        integer :: i,nwords,k
        logical :: whitespace
        nwords=wordcount(vector)
        
        if (allocated(words))deallocate(words)
        allocate(words(1:size(vector),1:nwords))
        words(:,:)=" "
        whitespace=.TRUE.
        nwords=0
        do i=1,size(vector)
            select case(iachar(vector(i)))
                case(iachar(' '),iachar(','),9) ! char(9) är tab
                    whitespace=.TRUE.
                    vector(i)=" "
                case default
                    if(whitespace)then
                        nwords=nwords+1
                        k=1
                    else
                        k=k+1
                    endif
                    words(k,nwords)=vector(i)
                    whitespace=.FALSE.
            end select
        end do
 
    end subroutine getwords!}}}

    function concatargs(charmat) result(st)!{{{
        integer(kind=ik) :: i
        character(kind=1,len=1) :: charmat(:,:)
        character(kind=1,len=size(charmat,2)*(size(charmat,1)+1)) :: st
        st=""
        do i=1,size(charmat,2)
            st=trim(st)//trim(.str.charmat(:,i))//"-"
        end do
        i=len(trim(st))
        st(i:)=""
    end function concatargs!}}}

    subroutine procinp(charvector,trajop)!{{{
        implicit none    
        character(kind=1,len=1),allocatable ::charvector(:),arguments(:,:)
        character(kind=1,len=3) :: funcstr
        character(kind=1,len=20) :: arg2
        integer(kind=ik) ::&
        ios,i,j,aind1,aind2,aind3,aind4,findex,p!,trajop(:,:)
        type(instruct) :: trajop
        call getwords(charvector,arguments)
        trajop%findex=0
        trajop%set=global_setflags
        trajop%setapl=.FALSE.
        p=0
        select case(trim(stringconv(arguments(:,1)))) ! Arg 1
            case('traj')
                call initgro(arguments(:,2))
                allocate(trajfile(len_trim(stringconv(arguments(:,2)))))
                trajfile=arguments(1:size(trajfile),2)

            case('set','SET')
                if(trim(stringconv(arguments(:,2)))=='apl')trajop%setapl=.TRUE.
                call set(arguments)

            case('dirangle','DIRANGLE','da','DA')

                !findex=1 ! Funktionsindex, för att kunna kalla på rätt funktion
                trajop%findex=1
                p=3 ! Antalet argument (operation, atomnamn etc.).
                ! Argument p+1 är filnamn för output
                funcstr='DA_'

            case('valenceangle','VALENCEANGLE','va','VA')

                trajop%findex=2
                p=4
                funcstr='VA_'

            case('torsionangle','ta','TA')

                trajop%findex=3
                p=5
                funcstr='TA_'

            case('bondlength','BL','bl')

                trajop%findex=4
                p=3 ! Antalet argument (operation, atomnamn etc.).
                ! Argument p+1 är filnamn för output
                funcstr='BL_'

            case('orderparameter','Sv','SV')
                trajop%findex=5
                p=3
                funcstr='SV_'

            case('membraneposition','MP')
                trajop%findex=6
                p=2
                funcstr='MP_'
            case('corr','co','CO')
                trajop%findex=7
                p=1
                funcstr='CO_'
            case('dc','DC','dipolecoupling')
                trajop%findex=8
                p=3
                funcstr='DC_'
            case('average')
                trajop%findex=9
                p=2 !faktiskt 2, men andra arg ej molekylnamn
                funcstr='AV_'
            case('define')
                 select case(trim(stringconv(arguments(:,2))))
                    case('atom')
                        trajop%findex=10
                        p=4
                 end select
            case('al','AL','apl')
                trajop%findex=11
                p=size(arguments,2)
                funcstr='AL_'
                
            case('exit')
                stop
                    
                        
                !end select

            case default 
                write(*,*)"Not a valid input, ",":",&
                trim(stringconv(arguments(:,1))),":"
                write(*,*)len(trim(stringconv(arguments(:,1)))),size(arguments(:,1)),arguments(:,1)
                stop
            end select
            if(p>=2)allocate(trajop%atoms(p-1))
            trajop%instructionstring=''
            if(trajop%findex/=0)then
                if(trajop%findex/=9 .AND. trajop%findex/=10)then
                do i=1,p-1
                    trajop%atoms(i)=atomindex(trim(stringconv(arguments(:,i+1))),molt(:)%molname,size(molt))&
                    +atomindex(trim(stringconv(arguments(:,i+1))))
                    if(trajop%atoms(i)==0)then
                        write(*,*)'Input is not an atom or a molecule '&
                        ,trim(stringconv(arguments(:,i+1)))
                        stop
                    endif
                end do
                end if
                if(p>=2)then
                    if(trajop%findex==9)then ! AVERAGE
                    trajop%instructionstring=funcstr//&
                    trim(concatargs(arguments(:,2:size(arguments,2))))
                        arg2=trim(stringconv(arguments(:,p)))
                        read(arg2,*,iostat=ios)trajop%average_count
                        if(ios/=0)then
                            write(*,*)'Input is not an integer',&
                            trim(stringconv(arguments(:,p-1)))
                            stop
                        endif
                    else if(trajop%findex==10)then ! DEFINE ATOM
                        trajop%newatom%atomname=trim(stringconv(arguments(:,3)))
                        trajop%newatom%from_mol_prop=trim(stringconv(arguments(:,4)))
                        trajop%newatom%molecule=trim(stringconv(arguments(:,5)))
                       ! This indexing works for one and only one atom per
                       ! submolecule!
                       ! Reallocate indexing vectors to add one atom:
                        call reallocatechar(atomnames,size(atomnames)+1)
                        call reallocate(shift,size(shift)+1)
                        call reallocate(natoms,size(natoms)+1)
                        call reallocate(moltypeofuatom,size(moltypeofuatom)+1)
                        call reallocate(masses,size(masses)+1)
                        ! Handle the indexing for the new atom:
                        masses(size(masses))=1._rk
                        i=size(atomnames)
                        j=atomindex(trim(trajop%newatom%molecule),molt(:)%molname,size(molt))
                        moltypeofuatom(i)=j
                        molt(j)%natoms=molt(j)%natoms+1
                        natoms(size(atomnames))=molt(j)%natoms
                        atomnames(size(atomnames))=trajop%newatom%atomname
                        trajop%atoms(1)=size(atomnames)
                        atot=atot+molt(j)%nmol
                        shift(i)=sum(molt(1:j-1)%nmol*molt(1:j-1)%natoms)+i-sum(molt(1:j)%natoms)
                        deallocate(coor)
                        allocate(coor(1:3,atot))


                    else
                        trajop%instructionstring=funcstr//trim(concatargs(arguments(:,2:p)))
                    endif

                else
                    trajop%instructionstring=trim(funcstr)
                endif
            end if

    end subroutine procinp!}}}

    subroutine set(args)!arg2,arg3)!{{{
        character(kind=1,len=1) :: args(:,:)
        character(kind=1,len=size(args,1)) :: arg2,arg3,arg4,arg5
        integer(kind=ik) :: i,j,p,ios,imol,kl,ku
        arg2='';arg3='';arg4='';arg5=''
        if(size(args,2)>=2)arg2=trim(stringconv(args(:,2)))
        if(size(args,2)>=3)arg3=trim(stringconv(args(:,3)))
        if(size(args,2)>=4)arg4=trim(stringconv(args(:,4)))
        if(size(args,2)>=5)arg5=trim(stringconv(args(:,5)))

        select case(arg2)
            case ('autofilename','AUTOFILENAME','Autofilename')
                select case(arg3)
                    case('on','On','ON')
                        global_setflags%autofilename=.TRUE.
                    case('off','Off','OFF')
                        global_setflags%autofilename=.FALSE.
                    case default
                        global_setflags%autofilename=.TRUE.
                end select

            case('maxframes','mf')
                read(arg3,*,iostat=ios)maxframes
                if (ios/=0)stop 'SET: maxframes should be of type integer'

            case('atomnames')
                open(nunit,file=arg3,action='read',status='old',iostat=ios)
                if(ios/=0)then
                    write(*,*)"Error, cannot open ",&
                    "file ",":",arg3,":"
                    stop
                endif
                do i=1,size(atomnames)
                read(nunit,*,iostat=ios)atomnames(i)
                if(ios/=0)then
                    write(*,*)"ERROR: (",i,")",arg3,ios
                    stop
                endif
                end do
                close(nunit)

            case('filename')
                if(arg3=='off')then
                    global_setflags%filename=''
                else
                    global_setflags%filename=''
                    global_setflags%filename=arg3
                endif

            case('calc')

                if (allocated(global_setflags%calc))deallocate(global_setflags%calc)
                allocate(global_setflags%calc(size(args,2)-2))

                do i=1,size(global_setflags%calc)
                    global_setflags%calc(i)=trim(stringconv(args(:,i+2)))
                end do

            case('distbin')
                read(arg3,*,iostat=ios)global_setflags%distbin
                if (ios/=0 .OR. global_setflags%distbin<=0)&
                stop 'SET: distbin should be of type integer >0'

            case('atommasses')
                allocate(masses(size(atomnames)),mgratios(size(atomnames)))
                masses=0;mgratios=0

                if(len_trim(arg3)/=0)then
                open(nunit,file=arg3,action='read', status='old',iostat=ios)
                if(ios==0)then
                    !Läs in atommassor från fil:
                    j=size(atomd)
                    i=j
                    !Räkna rader i filen
                    do while(ios/=endf)
                      read(nunit,*,iostat=ios)
                      i=i+1
                    end do
                    allocate(atomd(i))
                    call globals
                    i=j+1
                    
                    do while(ios/=endf)
                      read(nunit,*,iostat=ios)atomd(i)%aname,atomd(i)%mass,atomd(i)%mgratio
                      if(ios>0)stop 'Cannot read atom masses'
                      i=i+1
                    end do
                end if
                close(nunit,iostat=ios)
                end if
                !Ge varje unikt atomnamn en massa: 
                do i=1,size(atomd)
                    p=len_trim(atomd(i)%aname)
                    do j=1,size(atomnames)
                        if(atomd(i)%aname(1:p)==atomnames(j)(1:p))then
                            masses(j)=atomd(i)%mass
                            mgratios(j)=atomd(i)%mgratio
                        endif
                     end do
                end do
            case('centerofmembrane')
               allocate(common_setflags%membrane_moltypes(size(args,2)-2))
               do i=1,size(common_setflags%membrane_moltypes)
                   common_setflags%membrane_moltypes(i)=atomindex(trim(stringconv(args(:,2+i))),molt(:)%molname,size(molt))
               end do
            case('constant_bondlength','cbl')
                select case(arg3)
                    case('on','ON')
                       global_setflags%cbl_switch=.TRUE.
                       read(arg4,*,iostat=ios)global_setflags%constant_bl
                       if(ios/=0)stop 'SET: cannot read constant bond length'
                    case('off','OFF')
                        global_setflags%cbl_switch=.FALSE.
                        global_setflags%constant_bl=0
                end select
            case('fileprefix')
                if(size(args,2)>=3)then
                    global_setflags%fileprefix=arg3
                else
                    global_setflags%fileprefix=''
                end if
            case('filesuffix')
                if(size(args,2)>=3)then
                    global_setflags%filesuffix=arg3
                else
                    global_setflags%filesuffix=''
                end if
            case('writeframe')
                if(size(args,2)>=3)then
                    global_setflags%wftot=global_setflags%wftot+1
                    call reallocate(global_setflags%writeframe,global_setflags%wftot)
                    read(arg3,*)global_setflags%writeframe(global_setflags%wftot)%framenumber
                    if(size(args,2)==4)then
                        read(arg4,*)global_setflags%writeframe(global_setflags%wftot)%outformat
                    else
                        global_setflags%writeframe(global_setflags%wftot)%outformat='xyz'
                    end if
                else
                    write(*,*)'ERROR:SET:writeframe: Needs an integer'
                    stop
                end if
            case('area_per_lipid','apl')
                if(size(args,2)>=3)then
                    call apl_atomlist(args(:,3:))
                    global_setflags%apl=.TRUE.
                else
                    write(*,*)'SET: Area per lipid: Requires atomnames'
                end if
                
                
               !allocate(common_setflags%apl_moltypes(size(args,2)-2))
            case('submoltype')
                call reallocate(molt,size(molt)+1)
                mols=mols+1
                molt(size(molt))%molname=arg3
                if(arg4=='')stop 'ERROR:SET:NEWMOLTYPE: Needs start atom!'
                if(arg5=='')stop 'ERROR:SET:NEWMOLTYPE: Needs end atom!'
                if(moltypeofuatom(atomindex(arg4))/=moltypeofuatom(atomindex(arg5)))stop &
                'ERROR:SET:NEWMOLTYPE: You cannot mix atoms from different &
                molecules'
                molt(size(molt))%firstatom=atomindex(arg4)
                molt(size(molt))%lastatom=atomindex(arg5)
                molt(size(molt))%nmol=&
                molt(moltypeofuatom(molt(size(molt))%firstatom))%nmol
                molt(size(molt))%natoms=0
                if(molt(size(molt))%firstatom > molt(size(molt))%lastatom)then
                    write(*,*)'ERROR:SET:NEWMOLTYPE: Illegal atom order!'
                    stop
                endif
            case('folding')
                global_setflags%folding=.TRUE.

            case('aplgrid')
                read(arg3,*,iostat=ios)global_setflags%aplgrid(1)
                if(ios/=0)stop 'SET:aplgrid: Arg1'
                read(arg4,*,iostat=ios)global_setflags%aplgrid(2)
                if(ios/=0)stop 'SET:aplgrid: Arg2'

            case('leaflets')
                
                select case(arg3)
                    case('both')
                        global_setflags%leaflet=0
                    case('lower')
                        global_setflags%leaflet=1
                    case('upper')
                        global_setflags%leaflet=2
                    case default
                        !write(*,*)coor(:,10)
                        call whole
                        if(allocated(masses))then
                            call foldmol
                        else
                            stop 'LEAFLETS:FOLDMOL: Atommasses not initiated'
                        end if
                        if(allocated(common_setflags%membrane_moltypes))then
                            call&
                            center_of_membrane(common_setflags%membrane_moltypes)
                        else
                            stop 'LEAFLETS: Needs center of membrane'
                        end if
                        !write(*,*)masses
                        do i=1,size(molt)
                        kl=0;ku=0
                            do imol=1,molt(i)%nmol
                         !   write(*,*)center_of_molecule(i,imol),centerofmembrane,director
                            if(dot_product(center_of_molecule(i,imol)-centerofmembrane,director) < 0._rk)then
                                    kl=kl+1
                                    call reallocate(molt(i)%lower,kl)
                                    molt(i)%lower(kl)=imol
                                else
                                    ku=ku+1
                                    call reallocate(molt(i)%upper,ku)
                                    molt(i)%upper(ku)=imol
                                end if
                            end do
                           ! write(*,*)size(molt(i)%upper)
                        end do
                        
                end select


            case default
                if(size(args,2)>=2)then
                    write(*,*)'SET: >',trim(arg2),'<  is not a valid argument'
                else
                    write(*,*)'SET: No valid arguments'
                end if
                stop
        end select
    end subroutine set!}}}
end module input
