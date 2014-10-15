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
        integer(kind=ik) :: i,ios!,runit
        integer(kind=4) :: runit
        character(kind=1,len=255) :: carg,carg2,ctime
        common_setflags%silent=.FALSE.
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
            case('-s','--silent')
                common_setflags%silent=.TRUE.
            case('-v','--version')
                !ctime is set by the compiler using c type macro and cpp. Make
                !sure the file is suffixed with capital F or that the compiler
                !is flagged to use the preprocessor.
                ctime=CINFO
                write(*,*)'      Trajman experimental :: ',trim(ctime)
                write(*,*)'      Branch :: ',branch
                write(*,*)'      Revision :: ',revision
                write(*,*)'      Committed ',revdate
                write(*,*)'      Copyright (c) 2010 Jon Kapla :: Contact: jon.kapla@mmk.su.se'
                write(*,*)'      --------------------------------------------------------------'
                write(*,*)'      Trajman comes with NO WARRANTY, to the extent permitted by law.'
                write(*,*)'      You may redistribute copies of Trajman under the terms of the'
                write(*,*)'      GNU General Public License.'
                write(*,*)'      For more information about these matters, see the file named COPYING.'
                write(*,*)

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
        write(*,*)'  -s, --silent'
        write(*,*)'          No output to STDOUT'
        write(*,*)'  -'         
        write(*,*)'          Read from standard input, STDIN'
        write(*,*)
        p=getpid()
        write(pid,*)p
        p=system("cat $(readlink /proc/"//trim(adjustl(pid))//"/exe)_documentation.txt")
        stop
    end subroutine print_help!}}}

    subroutine summary!{{{
        integer(kind=ik) :: i
        character(len=100) :: A,B
        A=''
        A="          Molecules Atoms(mol)  Atoms(tot)"
        if(global_setflags%leaflets_defined)A=trim(A)//'  Upper(mol)  Lower(mol)'
        B=''
        B=stringconv([[(' ',i=1,5)],[('-',i=6,len_trim(A))]])
        write(*,*)
        if(.NOT.global_setflags%whole)then
            write(0,*)"     WARNING: Whole not set!"
            write(0,*)"        If you want to make broken molecules"
            write(0,*)"        whole, add 'set whole' to the input."
        end if
        write(*,*)
        write(*,*)trim(A)
        write(*,*)trim(B)
        do i=1,size(molt)
            write(*,'(5X,A5,A2,I7,2(1A,I11))',advance="no")&
            trim(molt(i)%molname),': ',molt(i)%nmol,' ',molt(i)%natoms,&
            ' ',molt(i)%natoms*molt(i)%nmol!,
            if(global_setflags%leaflets_defined)then
                write(*,'(2(1A,I11))')' ',size(molt(i)%upper),' ',size(molt(i)%lower)
            else
                write(*,*)
            end if
        end do
        write(*,*)trim(B)
        write(*,'(5X,A5,A2,I7,2(1A,I11))',advance="no")&
        "Total",': ',sum(molt(:)%nmol),' ',sum(molt(:)%natoms),' ',size(coor,2)
        if(global_setflags%leaflets_defined)then
            write(*,'(2(1A,I11))')' ',sum([(size(molt(i)%upper),i=1,size(molt))]),' ',sum([(size(molt(i)%lower),i=1,size(molt))])
        else
            write(*,*)
        end if
        write(*,*)
        if(global_setflags%whole)write(*,*)"     Function Whole is active."
        if(global_setflags%apl)write(*,*)"     Area per lipid calculation is enabled."
        !if(.not.global_setflags%distminmax%switch)write(*,*)"     All or some of the &
        !distributions have automatic limits.[DEFAULT]"
        if(global_setflags%distminmax%switch)write(*,*)"     All or some of the &
        distributions have user defined limits.[NOT DEFAULT]"
        if(.not.global_setflags%VSnorm)write(*,*)"     All or some of the VS functions are &
        normalized with binsize and frames.[NOT DEFAULT]"
        !if(global_setflags%VSnorm)write(*,*)"     All or some of the VS functions are &
        !normalized with molecules and frames.[DEFAULT]"
        if(global_setflags%zrdf)write(*,*)'      WARNING! Z-rdf not implemented.'
        if(global_setflags%zrdf)write(*,*)'      The rdftype z setting works only for MBL!'

        if(skipframes/=0)write(*,*)"     Skipping first ",trim(adjustl(intstr(skipframes)))," frames."
        write(*,*)
    end subroutine summary!}}}

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
!        external :: f77_molfile_read_next,f77_molfile_init,f77_molfile_open_read,&
!        f77_molfile_close_read,f77_molfile_finish
        character(kind=1,len=1),allocatable ::charvector(:),arguments(:,:)
        character(kind=1,len=3) :: funcstr
        character(kind=1,len=20) :: arg2,inttest
        character(kind=1,len=2000) :: infile,infilebefore
        integer(kind=ik) ::&
        ios,i,j,aind1,aind2,aind3,aind4,findex,p
        integer(kind=4) :: natm 
        type(instruct) :: trajop
        call getwords(charvector,arguments)
        trajop%findex=0
        trajop%set=global_setflags
        trajop%setapl=.FALSE.
        p=0;funcstr=''
        select case(trim(stringconv(arguments(:,1)))) ! Arg 1
            case('groupres')
                !call reallocate(groupres,merge(size(groupres),0,allocated(groupres))+1)
                allocate(groupres(size(arguments,2)-1))
                do i=1,size(groupres)
                    groupres(i)=trim(stringconv(arguments(:,1+i)))
                    write(*,*)groupres(i)
                enddo
            case('init')
               j=scan(stringconv(arguments(:,2)),'.',BACK=.TRUE.)
               j=merge(scan(stringconv(arguments(:,2)),'/',BACK=.TRUE.),j,j==0)
               select case(trim(stringconv(arguments(j+1:,2))))
                  case('gro')
                     call initgro(arguments(:,2))
                  case(dlpoly3histfile)
                     call inithist(arguments(:,2))
                  case default
                     stop 'INIT: unknown fileformat'
               end select
            case('traj')
                    call reallocate(trajfile,merge(size(trajfile),0,allocated(trajfile))+1)
                    !write(*,*)size(trajfile),'SIZE TRJF'
                    call reallocate(tunit,merge(size(tunit),0,allocated(tunit))+1)
                    !write(*,*)size(tunit),'SIZE TUNIT'
                i=size(trajfile)
                tunit(size(tunit))=tunit(size(tunit-1))+1
                if(size(tunit)==1)tunit(size(tunit))=12
                allocate(trajfile(i)%filename(len_trim(stringconv(arguments(:,2)))))
                trajfile(i)%filename=arguments(1:size(trajfile(i)%filename),2)
                infile=stringconv(trajfile(i)%filename)
                j=scan(infile,'.',BACK=.TRUE.)
                j=merge(scan(infile,'/',BACK=.TRUE.),j,j==0)
                trajtype=trim(infile(j+1:))
                if(i-1>0)then 
                    ! Check that the filetypes are the same, otherwise stop
                    infilebefore=stringconv(trajfile(i-1)%filename)
                    j=scan(infilebefore,'.',BACK=.TRUE.)
                    j=merge(scan(infilebefore,'/',BACK=.TRUE.),j,j==0)
                    if(trajtype/=trim(infilebefore(j+1:j&
                         +min(len_trim(infilebefore(j+1:)),len(trajtype)))))then
                        stop 'TRAJ: MIXED FILETYPES!'
                    end if
                end if
                select case(trim(trajtype))
                    case('gro')
                        open(tunit(size(tunit)),file=stringconv(trajfile(i)%filename),status='old')
                    case('trr','dcd','pdb')!MOLFILEPLUGIN FROM VMD
                        tunit(size(tunit))=-1
                        if(size(tunit)==1)call f77_molfile_init
                        call f77_molfile_open_read(tunit(size(tunit)),natm,stringconv(trajfile(i)%filename),trajtype)
                     !   write(*,*)tunit,'TUNIT'
                    case(dlpoly3histfile(1:min(len(dlpoly3histfile),len(trajtype))))
                        tunit(size(tunit))=-1
                        if(size(tunit)==1)call f77_molfile_init
                        call f77_molfile_open_read(tunit(size(tunit)),natm,stringconv(trajfile(i)%filename),dlpoly3histtag)
                    case default
                        tunit(size(tunit))=-1
                        if(size(tunit)==1)call f77_molfile_init
                        call f77_molfile_open_read(tunit(size(tunit)),natm,stringconv(trajfile(i)%filename),'auto')
                end select

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
                trajop%dpctype=0
                p=3
                funcstr='SV_'
            case('membraneposition','MP')
                if(.NOT.global_setflags%centerofmembrane)&
                stop 'MP: Needs centerofmembrane!'
                trajop%findex=6
                p=2
                funcstr='MP_'
            case('corr','co','CO')
                trajop%findex=7
                p=1
                funcstr='CO_'
            case('dc','DC','dipolarcoupling','d0','d1','d2','d3','d4','d5')
                trajop%findex=8
                p=3
                select case(trim(stringconv(arguments(:,1)))) ! Arg 1
                case('d1')
                    funcstr='D1_'
                    trajop%dpctype=1
                case('d2')
                    funcstr='D2_'
                    trajop%dpctype=2
                case('d3')
                    funcstr='D3_'
                    trajop%dpctype=3
                case('d4')
                    funcstr='D4_'
                    trajop%dpctype=4
                case('d5')
                    funcstr='D5_'
                    trajop%dpctype=5
                case default
                    funcstr='DC_'
                    trajop%dpctype=0
                end select
            case('average','combine')
                if(trim(stringconv(arguments(:,1)))=='combine')then
                    trajop%set%molaverage=.TRUE.
                    funcstr='MA_'
                else
                    funcstr='AV_'
                end if
                trajop%findex=9
                p=2 !faktiskt 2, men andra arg ej molekylnamn
            case('define')
                 trajop%findex=10
                 call define(trajop,arguments)
            case('al','AL','apl')
                if(.NOT.trajop%set%apl)stop 'AL: Needs apl!'
                if(.NOT.trajop%set%leaflets_defined)stop 'AL: Needs leaflets!'
                trajop%findex=11
                inttest=trim(stringconv(arguments(:,size(arguments,2))))
                if (atomindex(trim(stringconv(arguments(:,i+1))),molt(:)%molname,size(molt))&
                    +atomindex(trim(stringconv(arguments(:,i+1))))==0 .and.&
                    inttest(1:1)=='%')then
                    p=size(arguments,2)-1
                else
                    p=size(arguments,2)
                end if
                funcstr='AL_'
            case('MI','mi','moi')
                trajop%findex=12
                funcstr='MI_'
                p=2
            case('RDF','rdf','rf')
                trajop%findex=13
!                if(trajop%set%zrdf)write(*,*)'WARNING! Z-rdf not implemented. Works only &
!                for MBL'
                funcstr='RF_'
                p=3
            case('X','x','Y','y','Z','z')
                trajop%findex=14
                select case(trim(stringconv(arguments(:,1)))) ! Arg 1
                case('X','x')
                    funcstr='XC_'
                case('Y','y')
                    funcstr='YC_'
                case('Z','z')
                    funcstr='ZC_'
                end select
                p=2
            case('GD','gd')
                if(global_setflags%apl)stop 'APL cannot be active'
                trajop%findex=15
                funcstr='GD_'
                p=1
            case('shuffle_GD')
                trajop%findex=16
                funcstr='SH_'
                p=2
            case('boxx','bx','BX','boxy','by','BY','boxz','bz','BZ','BXBY','BXBZ','BYBZ','BXBYBZ')
                trajop%findex=17
                select case(trim(stringconv(arguments(:,1)))) ! Arg 1
                case('boxx','bx','BX')
                    funcstr='BX_'
                case('boxy','by','BY')
                    funcstr='BY_'
                case('boxz','bz','BZ')
                    funcstr='BZ_'
                case('BXBY')
                    funcstr='XY_'
                case('BXBZ')
                    funcstr='XZ_'
                case('BYBZ')
                    funcstr='YZ_'
                case('BXBYBZ')
                    funcstr='BV_'
                end select
                p=1
            case('boxapl')
                trajop%findex=18
                funcstr='BA_'
                p=2
            case('versus','VS','VERSUS')
                trajop%findex=19
                funcstr='VS_'
                p=1
            case('rotcorr','RC','ROTCORR')
                trajop%findex=20
                funcstr='RC_'
                p=1
            case('constant','CS')
                trajop%findex=21
                funcstr='CS_'
                if(.not.trajop%set%const%switch)stop 'CS: To use a constant &
                you need to set it first!'
                p=2
            case('mbl','MBL','MB')
                trajop%findex=22
                funcstr='MB_'
                p=3
            case('tagaverage','tagcombine','tagsum')
                if(trim(stringconv(arguments(:,1)))=='tagcombine')then
                    trajop%set%molaverage=.TRUE.
                    funcstr='MA_'
                elseif(trim(stringconv(arguments(:,1)))=='tagaverage')then
                    funcstr='AV_'
                else
                    trajop%set%instructionsum=.TRUE.
                    funcstr='SU_'
                end if
                trajop%findex=23
                p=2 !function+number of tags
            case('deallocate','dealloc')
                trajop%findex=24
                funcstr=''
                p=1
            case('karplus')
                if(.NOT.global_setflags%karplus)&
                stop 'KP: Needs five karplus_params!&
            <set karplus_params A B C D E>'
                trajop%findex=25
                funcstr='KP_'
                select case(global_setflags%karplus_fnc)
                case(3,4)
                    p=9
                case default
                    p=5
                end select
            case('SP','sp','splay')
                trajop%findex=26
                funcstr='SP_'
                p=3
            case('exit')
                stop

            case default 
                write(*,*)"Not a valid input, ",">",&
                trim(stringconv(arguments(:,1))),"<"
                write(*,*)len(trim(stringconv(arguments(:,1)))),size(arguments(:,1)),arguments(:,1)
                stop
            end select
            trajop%cv%mean=0
            trajop%cv%meandev=0
            trajop%cv%n=0
            !------------------------------------------------------
            if(p>=2)allocate(trajop%atoms(p-1))
            trajop%instructionstring=''
            trajop%set%corrindex=''
            trajop%ref=''
                select case(trajop%findex)
                case(0,10)
                case(7,19,20)
                    select case(size(arguments,2)-1)
                    case(2,3)
                        trajop%set%corrindex(1)=trim(stringconv(arguments(:,2)))
                        trajop%set%corrindex(2)=trim(adjustl(stringconv(arguments(:,3))))
                        if(size(arguments,2)==4)trajop%ref=trim(stringconv(arguments(:,4)))!TAG the function
                        trajop%instructionstring=trim(funcstr)
                    case(6,7)
                        trajop%set%corrindex(1)=trim(stringconv(arguments(:,2)))
                        trajop%set%corrindex(2)=trim(adjustl(stringconv(arguments(:,3))))
                        trajop%set%corrindex(3)=trim(adjustl(stringconv(arguments(:,4))))
                        trajop%set%corrindex(4)=trim(adjustl(stringconv(arguments(:,5))))
                        trajop%set%corrindex(5)=trim(adjustl(stringconv(arguments(:,6))))
                        trajop%set%corrindex(6)=trim(adjustl(stringconv(arguments(:,7))))
                        if(size(arguments,2)==8)trajop%ref=trim(stringconv(arguments(:,8)))!TAG the function
                        trajop%instructionstring=trim(funcstr)
                    case default
                        stop 'CORR/VS: Wrong number of arguments.&
                                Use two CHAR arguments for corr or VS, six for&
                                ROTCORR'
                    end select
                   
                case(9)
                    trajop%instructionstring=funcstr//&
                    trim(concatargs(arguments(:,2:size(arguments,2))))
                    !write(*,*)trajop%instructionstring
                    !stop
                    arg2=trim(stringconv(arguments(:,p)))
                    read(arg2,*,iostat=ios)trajop%average_count
                    if(ios/=0)then
                        write(*,*)'Input is not an integer',&
                        trim(stringconv(arguments(:,p-1)))
                        stop
                    endif
                    if(size(arguments,2)>p)trajop%ref=trim(stringconv(arguments(:,p+1)))
                case(23) !Tagaverage,tagcombine
                    trajop%average_count=readint(trim(stringconv(arguments(:,p))))
                    !read(arg2,*,iostat=ios)trajop%average_count
                    if(.not.allocated(trajop%set%averagetags))&
                        allocate(trajop%set%averagetags(trajop%average_count))
                    do i=1,trajop%average_count
                        trajop%set%averagetags(i)=trim(adjustl(stringconv(arguments(:,p+i))))
                    end do
                    if(size(arguments,2)>(p+trajop%average_count))then
                        trajop%ref=trim(stringconv(arguments(:,p+trajop%average_count+1))) !TAG
                    else
                        write(*,*)"ERROR:TAGAVERAGE:TAGCOMBINE:TAGSUM:"
                        write(*,*)"You need to name the function with a tag!"
                        stop
                    end if
                    trajop%instructionstring=funcstr//'tags_'//trim(trajop%ref)!//':_'//&
                    !    trim(stringconv(arguments(:,p+1)))//'--'//&
                    !    trim(stringconv(arguments(:,p+trajop%average_count)))
                   ! write(*,*)trajop%instructionstring
                   ! stop
                case(24) !DEALLOCATE INSTR
                    if((size(arguments,2)-p)<= 0)then
                        write(*,*)'ERROR INPUT:DEALLOCATE: tags missing?'
                        stop
                    end if
                    allocate(trajop%set%dealloctags(size(arguments,2)-p))
                    do i=1,size(trajop%set%dealloctags)
                        trajop%set%dealloctags(i)=trim(stringconv(arguments(:,p+i)))
                    end do
                case default
                    if(p>=2)then
                        do i=1,p-1
                            trajop%atoms(i)=atomindex(trim(stringconv(arguments(:,i+1))),molt(:)%molname,size(molt))&
                            +atomindex(trim(stringconv(arguments(:,i+1))))
                            if(trajop%atoms(i)==0)then
                                write(*,*)'Input is not an atom or a molecule '&
                                ,trim(stringconv(arguments(:,i+1)))
                                stop
                            endif
                        end do
                        if(size(arguments,2)>p)trajop%ref=trim(stringconv(arguments(:,p+1)))
                        trajop%instructionstring=funcstr//trim(concatargs(arguments(:,2:p)))
                    else
                        trajop%instructionstring=trim(funcstr)
                        if(size(arguments,2)>1)trajop%ref=trim(stringconv(arguments(:,2)))
                    end if
                    if(trajop%findex==12)trajop%atoms(1)=molt(trajop%atoms(1))%firstatom !MOI
                end select
            !------------------------------------------------------
    end subroutine procinp!}}}

    subroutine define(trajop,arguments)!{{{
        type(instruct) :: trajop
        character(kind=1,len=1) :: arguments(:,:)
        character(kind=1,len=size(arguments,1)) :: arg2,arg3,arg4
        integer(kind=ik) :: i,j,p
        if(size(arguments,2)>=2)arg2=trim(stringconv(arguments(:,2)))
        if(size(arguments,2)>=3)arg3=trim(stringconv(arguments(:,3)))
        if(size(arguments,2)>=4)arg4=trim(stringconv(arguments(:,4)))
        select case(arg2)
            case('centerofmembrane')
                trajop%define=1
                if(.NOT. global_setflags%folding)stop 'CENTEROFMEMBRANE: Needs folding'
                if(allocated(common_setflags%membrane_moltypes))deallocate(common_setflags%membrane_moltypes)
                allocate(common_setflags%membrane_moltypes(size(arguments,2)-2))
                do i=1,size(common_setflags%membrane_moltypes)
                    common_setflags%membrane_moltypes(i)=atomindex(trim(stringconv(arguments(:,2+i))),molt(:)%molname,size(molt))
                end do
                call center_of_membrane(common_setflags%membrane_moltypes)
                global_setflags%centerofmembrane=.TRUE.
            case('atom')
                atomsdefined=atomsdefined+1
                allocate(trajop%atoms(1))
                trajop%define=2
                trajop%newatom%atomname=trim(stringconv(arguments(:,3)))
                trajop%newatom%atype=trim(stringconv(arguments(:,4)))
                trajop%newatom%molecule=trim(stringconv(arguments(:,5)))
                if(size(arguments,2)>5)then
                   ! trajop%newatom%angle=readreal(trim(stringconv(arguments(:,6))))
                    do i=6,size(arguments,2)
                        if(i==6)allocate(trajop%newatom%helpers(1))
                        j=size(trajop%newatom%helpers)
                        if(i>6) call reallocate(trajop%newatom%helpers,j+1)
                        trajop%newatom%helpers(size(trajop%newatom%helpers))=atomindex(trim(stringconv(arguments(:,i))))
                        if(trajop%newatom%helpers(size(trajop%newatom%helpers))==0)then
                            write(*,*)'DEFINE: atom: helpers: >',trim(stringconv(arguments(:,i))),'<' 
                            stop
                        end if
                    end do
                end if
                ! Reallocate indexing vectors to add one atom:
                call reallocatechar(atomnames,size(atomnames)+1)
                call reallocate(natoms,size(natoms)+1)
                call reallocate(moltypeofuatom,size(moltypeofuatom)+1)
                call reallocate(masses,size(masses)+1)
                call reallocate(mgratios,size(mgratios)+1)
                ! Handle the indexing for the new atom:
                atomnames(size(atomnames))=trajop%newatom%atomname
                !write(*,*)atomnames(size(atomnames)),' name'

                do i=1,size(atomd)
                    p=len_trim(atomd(i)%aname)
                 !   write(*,*)p,' p',atomd(i)%aname(1:p)==atomnames(size(atomnames))(1:p)
                    if(atomd(i)%aname(1:p)==atomnames(size(atomnames))(1:p))then
                        masses(size(masses))=atomd(i)%mass
                        mgratios(size(mgratios))=atomd(i)%mgratio
                        exit
                    else
                        masses(size(masses))=1._rk
                        mgratios(size(mgratios))=0
                    end if
                end do
                write(*,*)mgratios(size(atomnames)),' mgr'
                i=size(atomnames)
                j=atomindex(trim(trajop%newatom%molecule),molt(:)%molname,size(molt))
                moltypeofuatom(i)=j
                molt(j)%natoms=molt(j)%natoms+1
                natoms(size(atomnames))=molt(j)%natoms
                !atomnames(size(atomnames))=trajop%newatom%atomname
                trajop%atoms(1)=size(atomnames)
                call reallocateatom(atom,size(atom)+1)
                atom(size(atom))%aname=trajop%newatom%atomname
                atom(size(atom))%mname=molt(j)%molname
                atom(size(atom))%moltype=j
                allocate(atom(size(atom))%coor(3,molt(j)%nmol))
            case('leaflet')
                trajop%define=3
                global_setflags%leaflets_defined=.TRUE.
                select case(arg3)
                    case('both')
                        trajop%set%leaflet=0
                        global_setflags%leaflet=0
                    case('lower')
                        trajop%set%leaflet=1
                        global_setflags%leaflet=1
                    case('upper')
                        trajop%set%leaflet=2
                        global_setflags%leaflet=2
                end select
            case('fold')
                select case(arg3)
                case('centerofmembrane')
                    global_setflags%foldcenterofmembrane=.TRUE.
                    trajop%define=4
                    call foldmol('com')
                end select
            case default
                write(*,*)'>',arg3,'< is not a valid for "define"!'
                stop
            end select
    end subroutine define!}}}

    subroutine set(args)!arg2,arg3)!{{{
        character(kind=1,len=1) :: args(:,:)
        character(kind=1,len=size(args,1)) :: arg2,arg3,arg4,arg5,arg6,arg7
        integer(kind=ik) :: i,j,p,ios,imol,kl,ku
        integer(kind=ik),allocatable :: helpers(:)
        real(kind=rk) :: v(3)
        arg2='';arg3='';arg4='';arg5='';arg6='';arg7=''
        if(size(args,2)>=2)arg2=trim(stringconv(args(:,2)))
        if(size(args,2)>=3)arg3=trim(stringconv(args(:,3)))
        if(size(args,2)>=4)arg4=trim(stringconv(args(:,4)))
        if(size(args,2)>=5)arg5=trim(stringconv(args(:,5)))
        if(size(args,2)>=6)arg6=trim(stringconv(args(:,6)))
        if(size(args,2)>=7)arg7=trim(stringconv(args(:,7)))
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
            case('skipframes','sf')
                read(arg3,*,iostat=ios)skipframes
                if (ios/=0)stop 'SET: skipframes should be of type integer'
            case('atomnames')
                open(nunit,file=arg3,action='read',status='old',iostat=ios)
                if(ios/=0)then
                    write(*,*)"Error, cannot open ",&
                    "file ",":",arg3,":"
                    stop
                endif
                do i=1,size(atomnames)
                read(nunit,*,iostat=ios)atomnames(i)
                if(ANY(atomnames(1:i-1)==atomnames(i)))then
                    write(*,*)"ERROR ATOMNAMES: The name ",trim(atomnames(i))," is not unique!"
                    stop
                end if
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
                    !write(*,*)size(atomd),'SIZE ATOMD'
                    i=j
                    !Räkna rader i filen
                    do while(ios/=endf)
                      read(nunit,*,iostat=ios)
                      i=i+1
                    end do
                    call reallocate(atomd,i-1)
                    !write(*,*)i,size(atomd),' SIZE atomdrealloc'
                    call globals
                    i=j+1
                    rewind(nunit)
                    ios=0
                    do while(ios/=endf)
                      read(nunit,*,iostat=ios)atomd(i)%aname,atomd(i)%mass,atomd(i)%mgratio
                      !write(*,*)atomd(i)%aname,atomd(i)%mass,atomd(i)%mgratio,'TESTREAD'
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
                            !write(*,*)'AMASS ',masses(j),mgratios(j),atomd(i)%aname,atomnames(j)(1:p)
                        endif
                     end do
                end do
                do i=1,size(masses)
                write(*,*)'AMASS ',masses(i),mgratios(i),atomnames(i)
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
                    !write(*,*)trim(global_setflags%writeframe(global_setflags%wftot)%outformat)
                else
                    write(*,*)'ERROR:SET:writeframe: Needs an integer'
                    stop
                end if
            case('area_per_lipid','apl','griddiff')
                
                if(size(args,2)>=3)then
                    call apl_atomlist(args(:,3:))
                    global_setflags%apl=.TRUE.
                    global_setflags%gd=.FALSE.
                    if(arg2=='griddiff')global_setflags%apl=.FALSE.
                    if(arg2=='griddiff')global_setflags%gd=.TRUE.
                else
                    write(*,*)'SET: Area per lipid: Requires atomnames'
                end if
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
                if(.NOT. allocated(masses))stop 'FOLDING: Needs atom masses'
                global_setflags%folding=.TRUE.
            case('aplgrid')
                read(arg3,*,iostat=ios)global_setflags%aplgrid(1)
                if(ios/=0)stop 'SET:aplgrid: Arg1'
                read(arg4,*,iostat=ios)global_setflags%aplgrid(2)
                if(ios/=0)stop 'SET:aplgrid: Arg2'
            case('leaflets')
                        if(global_setflags%whole)call whole
                        if(global_setflags%folding)call foldmol
                        if(.NOT.global_setflags%centerofmembrane)then
                            stop 'LEAFLETS: Needs a defined center of membrane'
                        end if
                        do i=1,size(molt)
                        kl=0;ku=0
                            do imol=1,molt(i)%nmol
                            if(.not.global_setflags%foldcenterofmembrane)then
                                if(dot_product(center_of_molecule(i,imol)&
                                -centerofmembrane,director) < 0._rk)then
                                    kl=kl+1
                                    call reallocate(molt(i)%lower,kl)
                                    molt(i)%lower(kl)=imol
                                else if(dot_product(center_of_molecule(i,imol)&
                                -centerofmembrane,director) > 0._rk)then
                                    ku=ku+1
                                    call reallocate(molt(i)%upper,ku)
                                    molt(i)%upper(ku)=imol
                                end if
                            else
                                v=center_of_molecule(i,imol)
                                if(v(3)-centerofmembrane(3) < 0._rk)then
                                    kl=kl+1
                                    call reallocate(molt(i)%lower,kl)
                                    molt(i)%lower(kl)=imol
                                else if(v(3)-centerofmembrane(3) > 0._rk)then
                                    ku=ku+1
                                    call reallocate(molt(i)%upper,ku)
                                    molt(i)%upper(ku)=imol
                                end if
                            end if
                            end do
                            if(.not.allocated(molt(i)%upper))then
                                allocate(molt(i)%upper(0))
                            end if
                            if(.not.allocated(molt(i)%lower))then
                                allocate(molt(i)%lower(0))
                            end if
                        end do
                        if(ku==0 .or. kl==0)then
                            write(*,*)'SET: LEAFLETS: Upper: ',ku,' Lower: ',kl
                            stop
                        end if
            case ('whole')
                global_setflags%whole=.TRUE.
            case('torsion_shift','tshift')
                if(size(args,2)>=2)then
                    global_setflags%tshift=readint(arg3)
                else
                    stop 'SET: Torsion shift needs an integer'
                end if
            case('cscale')
                if(size(args,2)>=2)then
                    common_setflags%traj_cscale=readreal(arg3)
                else
                    stop 'SET: Coordinate scale needs a real number'
                end if
            case('rdf_binsize')
                if(size(args,2)>=2)then
                    global_setflags%rdf_binsize=readreal(arg3)
                else
                    stop 'SET: RDF_Binsize needs a real number'
                end if
            case('rdftype')
                select case(arg3)
                case('xy')
                    global_setflags%xyrdf=.TRUE.
                    global_setflags%zrdf=.FALSE.
                case('z')
                    global_setflags%xyrdf=.FALSE.
                    global_setflags%zrdf=.TRUE.
                case default
                    global_setflags%xyrdf=.FALSE.
                    global_setflags%zrdf=.FALSE.
                end select
            case('shuffle_atoms')
                if(.not.allocated(common_setflags%shuffle_atoms))&
                allocate(common_setflags%shuffle_atoms(size(args,2)-2))
                do i=1,size(common_setflags%shuffle_atoms)
                    common_setflags%shuffle_atoms(i)=atomindex(stringconv(args(:,2+i)))
                end do
            case('constant')
                read(arg3,*,iostat=ios)global_setflags%const%value
                global_setflags%const%switch=.TRUE.
                if (ios/=0)stop 'SET: Constant should be of type real'
            case('VS_norm')
                select case(arg3)
                case('normal')
                    global_setflags%VSnorm=.TRUE.
                case('density')
                    global_setflags%VSnorm=.FALSE.
                case default
                    stop 'SET: VS_norm: normal or density?'
                end select
            case('distriblimits')
                if(size(args,2)<=3)then
                    global_setflags%distminmax%switch=.FALSE.
                else
                    global_setflags%distminmax%switch=.TRUE.
                    read(arg3,*,iostat=ios)global_setflags%distminmax%mi
                    if (ios/=0)stop 'SET: Disttriblimit min should be of type real'
                    read(arg4,*,iostat=ios)global_setflags%distminmax%ma
                    if (ios/=0)stop 'SET: Disttriblimit max should be of type real'
                end if
            case('scaling')
                if(size(args,2)>=3)then
                    global_setflags%scaling%switch=.true.
                    global_setflags%scaling%typ=arg3
                    select case(trim(global_setflags%scaling%typ))
                    case('x','y','z','xy','yx','xz','zx','yz','zy','xyz')
                    case('sub')
                        global_setflags%scaling%val=readreal(arg4)
                    case default
                        stop 'SCALING: Check input!'
                    end select
                else
                    global_setflags%scaling%typ=''
                    global_setflags%scaling%switch=.false.

                end if
            case('slice')
                !global_setflags%slice%switch=.TRUE.
                select case(arg3)
                case('Z_in')
                    global_setflags%slice%typ=arg3
                case('Z_out')
                    global_setflags%slice%typ=arg3
                case('off')
                    global_setflags%slice%switch=.FALSE.
                    global_setflags%slice%switch2=.FALSE.
                case default
                    write(*,*)'SLICE: >',arg3,'< Only slicing in Z (Z_in and Z_out) is implemented!'
                    stop
                end select
                if(size(args,2)==5)then
                    global_setflags%slice%switch=.TRUE.
                    global_setflags%slice%switch2=.FALSE.
                    global_setflags%slice%lower=readreal(arg4)
                    global_setflags%slice%upper=readreal(arg5)
                else if(size(args,2)==7)then
                    global_setflags%slice%switch=.TRUE.
                    global_setflags%slice%switch2=.TRUE.
                    global_setflags%slice%lower=readreal(arg4)
                    global_setflags%slice%upper=readreal(arg5)
                    global_setflags%slice%lower2=readreal(arg6)
                    global_setflags%slice%upper2=readreal(arg7)
                else
                    stop 'SLICE: Wrong number of arguments!'
                end if
            case('ch_bondlength')
                global_setflags%ch_bondlength=readreal(arg3)
            case('karplus_fnc')
                global_setflags%karplus_fnc=readint(arg3)
            case('karplus_params')
                !Acos²(theta+D)+Bcos(theta+D)+Csin(theta+D)+E
                global_setflags%karplus=.TRUE.
                if(global_setflags%karplus_fnc==0) &
                    stop 'karplus_params need karplus_fnc'
                global_setflags%karplus_params(1)=readreal(arg3) !A
                global_setflags%karplus_params(2)=readreal(arg4) !B
                global_setflags%karplus_params(3)=readreal(arg5) !C
                global_setflags%karplus_params(4)=readreal(arg6) !D
                !if(global_setflags%karplus_fnc==1)global_setflags%karplus_params(5)=readreal(arg7) !E
                global_setflags%karplus_params(5)=readreal(arg7) !E
            case('coorsys')
                if(size(args,2)<5)stop 'COORSYS needs 3 helper atoms'
                !if(.not.allocated(global_setflags%coorsys_helpers))&
                !    allocate(global_setflags%coorsys_helpers(size(args,2)-2))
                global_setflags%coorsys=.TRUE.
                global_setflags%coorsys_helpers(1)=atomindex(arg3)
                global_setflags%coorsys_helpers(2)=atomindex(arg4)
                global_setflags%coorsys_helpers(3)=atomindex(arg5)
            case('coorsys_type')
                global_setflags%coorsys_type=readint(arg3)
            case('sph1','sph2','spt1','spt2')
                ! Helper atoms for the splay and tilt function (SP)
                if(.not.allocated(helpers))&
                    allocate(helpers(size(args,2)-2))
                do i=1,size(helpers)
                    helpers(i)=atomindex(trim(stringconv(args(:,i+2))))
                end do
                select case(arg2)
                case('sph1')
                    call reallocate(global_setflags%sph1,size(helpers))
                    global_setflags%sph1=helpers
                case('sph2')
                    call reallocate(global_setflags%sph2,size(helpers))
                    global_setflags%sph2=helpers
                case('spt1')
                    call reallocate(global_setflags%spt1,size(helpers))
                    global_setflags%spt1=helpers
                case('spt2')
                    call reallocate(global_setflags%spt2,size(helpers))
                    global_setflags%spt2=helpers
                end select
                deallocate(helpers)
                
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
