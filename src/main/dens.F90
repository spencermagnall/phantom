!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: densityforce
!
!  DESCRIPTION:
!  This module is the "guts" of the code
!  Calculates density by iteration with smoothing length
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, eos, fastmath, io, io_summary, kdtree,
!    kernel, linklist, mpidens, mpiderivs, mpiutils, nicil, options, part,
!    stack, timestep, timing, viscosity
!+
!--------------------------------------------------------------------------
module densityforce
 use dim,     only:maxdvdx,maxvxyzu,maxp,minpart,maxxpartvecidens,maxrhosum,&
                   maxdusttypes,maxdustlarge
 use part,    only:maxBevol,mhd,dvdx
 use kdtree,  only:iorder,inoderange
 use kernel,  only:cnormk,wab0,gradh0,dphidh0,radkern2
 use mpidens, only:celldens,stackdens
 use timing,  only:getused,printused,print_time



 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 public :: densityiterate,get_neighbour_stats,get_alphaloc

 !--indexing for xpartveci array
 integer, parameter :: &
       ixi  = 1, &
       iyi  = 2, &
       izi  = 3, &
       ivxi = 4, &
       ivyi = 5, &
       ivzi = 6, &
       ieni = 7, &
       iBevolxi = 8, &
       iBevolyi = 9, &
       iBevolzi = 10, &
       ipsi = 11, &
       ifxi = 12, &
       ifyi = 13, &
       ifzi = 14, &
       iradxii = 15

 !--indexing for rhosum array
 integer, parameter :: &
       irhoi            = 1, &
       igradhi          = 2, &
       igradsofti       = 3, &
       idivvi           = 4, &
       idvxdxi          = 5, &
       idvxdyi          = 6, &
       idvxdzi          = 7, &
       idvydxi          = 8, &
       idvydyi          = 9, &
       idvydzi          = 10, &
       idvzdxi          = 11, &
       idvzdyi          = 12, &
       idvzdzi          = 13, &
       idaxdxi          = 14, &
       idaxdyi          = 15, &
       idaxdzi          = 16, &
       idaydxi          = 17, &
       idaydyi          = 18, &
       idaydzi          = 19, &
       idazdxi          = 20, &
       idazdyi          = 21, &
       idazdzi          = 22, &
       irxxi            = 23, &
       irxyi            = 24, &
       irxzi            = 25, &
       iryyi            = 26, &
       iryzi            = 27, &
       irzzi            = 28, &
       idivBi           = 29, &
       idBxdxi          = 30, &
       idBxdyi          = 31, &
       idBxdzi          = 32, &
       idBydxi          = 33, &
       idBydyi          = 34, &
       idBydzi          = 35, &
       idBzdxi          = 36, &
       idBzdyi          = 37, &
       idBzdzi          = 38, &
       irhodusti        = 39, &
       irhodustiend     = 39 + (maxdustlarge - 1), &
       iradfxi          = irhodustiend + 1, &
       iradfyi          = irhodustiend + 2, &
       iradfzi          = irhodustiend + 3


 

 !--kernel related parameters
 !real, parameter    :: cnormk = 1./pi, wab0 = 1., gradh0 = -3.*wab0, radkern2 = 4F.0
 integer, parameter :: isizecellcache = 50000
 integer, parameter :: isizeneighcache = 0

 integer, parameter :: maxdensits = 50

 !--statistics which can be queried later
 integer, private         :: maxneighact,nrelink
 integer(kind=8), private :: nneightry,maxneightry,nneighact,ncalc
 integer(kind=8), private :: nptot = -1

 type denstreestacklocal
  integer :: nodeindex1 = 0 , nodeindex2 = 0
 end type 

 type denstreestack
  integer :: nodeindex1 = 0
  integer :: interactions(2) = 0 
 end type 

 private

contains

!----------------------------------------------------------------
!+
!  this is the main routine for the whole code
!+
!----------------------------------------------------------------
! subroutine densityiterate(icall,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
!                           fxyzu,fext,alphaind,gradh,rad,radprop)
!  use dim,       only:maxp,maxneigh,ndivcurlv,ndivcurlB,maxvxyzu,maxalpha, &
!                      mhd_nonideal,nalpha,use_dust
!  use io,        only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
!  use linklist,  only:ifirstincell,get_neighbour_list,get_hmaxcell,get_cell_list,&
!                      get_cell_location,set_hmaxcell,sync_hmax_mpi,update_hmax_remote,node_is_active
!  use part,      only:mhd,rhoh,dhdrho,rhoanddhdrho,&
!                      ll,get_partinfo,iactive,&
!                      hrho,iphase,igas,idust,iamgas,periodic,&
!                      all_active,dustfrac,Bxyz,set_boundaries_to_active
! #ifdef FINVSQRT
!  use fastmath,  only:finvsqrt
! #endif

!  use mpiutils,  only:reduceall_mpi,barrier_mpi,reduce_mpi,reduceall_mpi
! #ifdef MPI
!  use stack,     only:reserve_stack,swap_stacks
!  use stack,     only:stack_remote => dens_stack_1
!  use stack,     only:stack_waiting => dens_stack_2
!  use stack,     only:stack_redo => dens_stack_3
!  use mpiderivs, only:send_cell,recv_cells,check_send_finished,init_cell_exchange, &
!                      finish_cell_exchange,recv_while_wait,reset_cell_counters
! #endif
!  use timestep,  only:rhomaxnow
!  use part,      only:ngradh
!  use viscosity, only:irealvisc
!  use io_summary,only:summary_variable,iosumhup,iosumhdn
!  integer,      intent(in)    :: icall,npart,nactive
!  real,         intent(inout) :: xyzh(:,:)
!  real,         intent(in)    :: vxyzu(:,:),fxyzu(:,:),fext(:,:)
!  real,         intent(in)    :: Bevol(:,:)
!  real(kind=4), intent(out)   :: divcurlv(:,:)
!  real(kind=4), intent(out)   :: divcurlB(:,:)
!  real(kind=4), intent(out)   :: alphaind(:,:)
!  real(kind=4), intent(out)   :: gradh(:,:)
!  real,         intent(out)   :: stressmax
!  real,         intent(in)    :: rad(:,:)
!  real,         intent(inout) :: radprop(:,:)

!  integer, save :: listneigh(maxneigh)
!  real,   save :: xyzcache(3,isizecellcache)
! !$omp threadprivate(xyzcache,listneigh)

!  integer :: i,icell
!  integer :: nneigh,np
!  integer :: nwarnup,nwarndown,nwarnroundoff

!  logical :: getdv,realviscosity,getdB,converged
!  logical :: iactivei,iamgasi,iamdusti
!  integer :: iamtypei

!  logical :: remote_export(nprocs)

!  real    :: rhomax

!  integer                   :: npcell,istart,iend

!  type(celldens)            :: cell

!  logical                   :: redo_neighbours

! #ifdef MPI
!  integer                   :: j,k,l
!  integer                   :: irequestsend(nprocs),irequestrecv(nprocs)
!  type(celldens)            :: xrecvbuf(nprocs),xsendbuf
!  integer                   :: mpiits,nlocal
!  real                      :: ntotal
!  logical                   :: iterations_finished
!  logical                   :: do_export

!  call init_cell_exchange(xrecvbuf,irequestrecv)
!  stack_waiting%n = 0
!  stack_remote%n = 0
!  stack_redo%n = 0
! #endif

!  if (iverbose >= 3 .and. id==master) &
!     write(iprint,*) ' cell cache =',isizecellcache,' neigh cache = ',isizeneighcache,' icall = ',icall

!  if (icall==0 .or. icall==1) then
!     call reset_neighbour_stats(nneightry,nneighact,maxneightry,maxneighact,ncalc,nrelink)
!     nwarnup       = 0
!     nwarndown     = 0
!     nwarnroundoff = 0
!     np = 0
!  endif
!  !
!  ! flag for whether or not we need to calculate velocity derivatives
!  ! whilst doing the density iterations (needed for viscosity switches
!  ! and for physical viscosity)
!  !
!  realviscosity = (irealvisc > 0)
!  getdv = ((maxalpha==maxp .or. ndivcurlv >= 4) .and. icall <= 1) .or. &
!          (maxdvdx==maxp .and. (use_dust .or. realviscosity))
!  if (getdv .and. ndivcurlv < 1) call fatal('densityiterate','divv not stored but it needs to be')
!  getdB = (mhd .and. (ndivcurlB >= 4 .or. mhd_nonideal))

!  if ( all_active ) stressmax  = 0.   ! condition is required for independent timestepping


! #ifdef MPI
!  ! number of local only cells
!  nlocal = 0
!  call reset_cell_counters

! #endif

!  rhomax = 0.0
!  call get_cell_list(istart,iend)

! !$omp parallel default(none) &
! !$omp shared(icall,istart,iend) &
! !!$omp shared(ncells) &
! !!$omp shared(ll) &
! !$omp shared(ifirstincell) &
! !$omp shared(xyzh) &
! !$omp shared(vxyzu) &
! !$omp shared(fxyzu) &
! !$omp shared(fext) &
! !$omp shared(gradh) &
! !$omp shared(iphase) &
! !$omp shared(Bevol) &
! !$omp shared(divcurlv) &
! !$omp shared(divcurlB) &
! !$omp shared(alphaind) &
! !$omp shared(dustfrac) &
! !$omp shared(Bxyz) &
! !$omp shared(dvdx) &
! !$omp shared(id) &
! !$omp shared(nprocs) &
! !$omp shared(getdB) &
! !$omp shared(getdv) &
! !$omp shared(realviscosity) &
! !$omp shared(iverbose) &
! !$omp shared(iprint) &
! !$omp shared(rad,radprop)&
! #ifdef MPI
! !$omp shared(xrecvbuf) &
! !$omp shared(xsendbuf) &
! !$omp shared(irequestrecv) &
! !$omp shared(irequestsend) &
! !$omp shared(stack_remote) &
! !$omp shared(stack_waiting) &
! !$omp shared(stack_redo) &
! !$omp shared(iterations_finished) &
! !$omp shared(mpiits) &
! !$omp reduction(+:nlocal) &
! !$omp private(do_export) &
! !$omp private(j) &
! !$omp private(k) &
! !$omp private(l) &
! !$omp private(ntotal) &
! #endif
! !$omp private(remote_export) &
! !$omp private(nneigh) &
! !$omp private(npcell) &
! !$omp private(cell) &
! !$omp private(iamgasi) &
! !$omp private(iamtypei) &
! !$omp private(iactivei) &
! !$omp private(iamdusti) &
! !$omp private(converged) &
! !$omp private(redo_neighbours) &
! !$omp reduction(+:ncalc) &
! !$omp reduction(+:np) &
! !$omp reduction(max:maxneighact) &
! !$omp reduction(max:maxneightry) &
! !$omp reduction(+:nneighact) &
! !$omp reduction(+:nneightry) &
! !$omp reduction(+:nrelink) &
! !$omp reduction(+:stressmax) &
! !$omp reduction(max:rhomax) &
! !$omp private(i)

! !$omp do schedule(runtime)
!   over_leaf_nodes: do icell=istart,iend
! ! over_cells: do icell=1,int(ncells)
!     if (ifirstincell(icell) <= 0) cycle over_leaf_nodes
! !    if (.not.node_is_active(icell)) cycle over_leaf_nodes

!     !
!     !--get the neighbour list and fill the cell cache
!     !
!     call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
!                            remote_export=remote_export)
! #ifdef MPI
!     if (any(remote_export)) then
!        do_export = .true.
!     else
!        do_export = .false.
!     endif
! #endif
!     cell%icell                   = icell
! #ifdef MPI
!     cell%owner                   = id
! #endif
!     cell%nits                    = 0
!     cell%nneigh                  = 0
!     cell%remote_export(1:nprocs) = remote_export

!     call start_cell(cell,iphase,xyzh,vxyzu,fxyzu,fext,Bevol,rad)

!     call get_cell_location(icell,cell%xpos,cell%xsizei,cell%rcuti)
!     call get_hmaxcell(icell,cell%hmax)

! #ifdef MPI
! !$omp critical
!     call recv_cells(stack_remote,xrecvbuf,irequestrecv)
! !$omp end critical

!     if (do_export) then
! !$omp critical
!        if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
!        ! make a reservation on the stack
!        call reserve_stack(stack_waiting,cell%waiting_index)
!        ! export the cell: direction 0 for exporting
!        call send_cell(cell,0,irequestsend,xsendbuf)
! !$omp end critical
!     endif
! #endif

!     call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad)

! #ifdef MPI
!     if (do_export) then
!        ! write directly to stack
!        stack_waiting%cells(cell%waiting_index) = cell
!     else
! #endif
!        converged = .false.
!        local_its: do while (.not. converged)
!           call finish_cell(cell,converged)
!           call compute_hmax(cell,redo_neighbours)
!           if (icall == 0) converged = .true.
!           if (.not. converged) then
!              if (redo_neighbours) then
!                 call set_hmaxcell(cell%icell,cell%hmax)
!                 call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
!                                       cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti, &
!                                       remote_export=remote_export)
! #ifdef MPI
!                 cell%remote_export(1:nprocs) = remote_export
!                 if (any(remote_export)) then
!                    do_export = .true.
! !$omp critical
!                    if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
!                    call reserve_stack(stack_waiting,cell%waiting_index)
!                    ! direction export (0)
!                    call send_cell(cell,0,irequestsend,xsendbuf)
! !$omp end critical
!                 endif
! #endif
!                 nrelink = nrelink + 1
!              endif
!              call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad)
! #ifdef MPI
!              if (do_export) then
!                 stack_waiting%cells(cell%waiting_index) = cell
!                 exit local_its
!              endif
! #endif
!           endif
!        enddo local_its
! #ifdef MPI
!        if (.not. do_export) then
! #endif
!           call store_results(icall,cell,getdv,getdB,realviscosity,stressmax,xyzh,gradh,divcurlv,divcurlB,alphaind, &
!                              dvdx,vxyzu,Bxyz,dustfrac,rhomax,nneightry,nneighact,maxneightry,maxneighact,np,ncalc,&
!                              radprop)
! #ifdef MPI
!           nlocal = nlocal + 1
!        endif
!     endif
! #endif
!  enddo over_leaf_nodes
! !$omp enddo

! #ifdef MPI
! !$omp barrier

! !$omp single
!  if (stack_waiting%n > 0) call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
!  call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,irequestsend)
! !$omp end single

! !$omp single
!  if (iverbose>=6) then
!     ntotal = real(nlocal) + real(stack_waiting%n)
!     if (ntotal > 0) then
!        write(iprint,*) id,'local ratio = ',real(nlocal)/ntotal
!     else
!        write(iprint,*) id,'local ratio = 0'
!     endif
!  endif

!  mpiits = 0
!  iterations_finished = .false.
! !$omp end single
!  remote_its: do while(.not. iterations_finished)
! !$omp single
!     mpiits = mpiits + 1
!     call reset_cell_counters
! !$omp end single

!     igot_remote: if (stack_remote%n > 0) then
! !$omp do schedule(runtime)
!        over_remote: do i = 1,stack_remote%n
!           cell = stack_remote%cells(i)

!           ! icell is unused (-1 here)
!           call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
!                                   cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti)

!           call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,&
!                             rad)

!           cell%remote_export(id+1) = .false.

!           ! communication happened while computing contributions to remote cells
! !$omp critical
!           call recv_cells(stack_waiting,xrecvbuf,irequestrecv)
!           call check_send_finished(stack_waiting,irequestsend,irequestrecv,xrecvbuf)
!           ! direction return (1)
!           call send_cell(cell,1,irequestsend,xsendbuf)
! !$omp end critical
!        enddo over_remote
! !$omp enddo
! !$omp barrier
! !$omp single
!        ! reset remote stack
!        stack_remote%n = 0
!        ! ensure send has finished
!        call check_send_finished(stack_waiting,irequestsend,irequestrecv,xrecvbuf)
! !$omp end single
!     endif igot_remote
! !$omp barrier
! !$omp single
!     call recv_while_wait(stack_waiting,xrecvbuf,irequestrecv,irequestsend)
!     call reset_cell_counters
! !$omp end single
!     iam_waiting: if (stack_waiting%n > 0) then
! !$omp do schedule(runtime)
!        over_waiting: do i = 1, stack_waiting%n
!           cell = stack_waiting%cells(i)

!           if (any(cell%remote_export(1:nprocs))) then
!              print*,id,cell%remote_export(1:nprocs)
!              print*,id,'mpiits',mpiits
!              print*,id,'owner',cell%owner
!              print*,id,'icell',cell%icell
!              print*,id,'npcell',cell%npcell
!              print*,id,'xpos',cell%xpos
!              print*,id,'stackpos',i
!              print*,id,'waitindex',cell%waiting_index
!              call fatal('dens', 'not all results returned from remote processor')
!           endif

!           call finish_cell(cell,converged)
!           call compute_hmax(cell,redo_neighbours)

!           ! communication happened while finishing cell
! !$omp critical
!           call recv_cells(stack_remote,xrecvbuf,irequestrecv)
! !$omp end critical
!           if (.not. converged) then
!              call set_hmaxcell(cell%icell,cell%hmax)
!              call get_neighbour_list(-1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
!                                     cell_xpos=cell%xpos,cell_xsizei=cell%xsizei,cell_rcuti=cell%rcuti, &
!                                     remote_export=remote_export)
!              cell%remote_export(1:nprocs) = remote_export
! !$omp critical
!              call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
!              call reserve_stack(stack_redo,cell%waiting_index)
!              ! direction export (0)
!              call send_cell(cell,0,irequestsend,xsendbuf)
! !$omp end critical
!              call compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,&
!              rad)

!              stack_redo%cells(cell%waiting_index) = cell
!           else
!              call store_results(icall,cell,getdv,getdB,realviscosity,stressmax,xyzh,gradh,divcurlv,divcurlB,alphaind, &
!                                 dvdx,vxyzu,Bxyz,dustfrac,rhomax,nneightry,nneighact,maxneightry,maxneighact,np,ncalc,&
!                                 radprop)

!           endif

!        enddo over_waiting
! !$omp enddo
! !$omp barrier
! !$omp single
!        ! reset stacks
!        stack_waiting%n = 0
!        call check_send_finished(stack_remote,irequestsend,irequestrecv,xrecvbuf)
! !$omp end single
!     endif iam_waiting
! !$omp barrier
! !$omp single
!     call recv_while_wait(stack_remote,xrecvbuf,irequestrecv,irequestsend)
! !$omp end single

! !$omp single
!     if (reduceall_mpi('max',stack_redo%n) > 0) then
!        call swap_stacks(stack_waiting, stack_redo)
!     else
!        iterations_finished = .true.
!     endif
!     stack_redo%n = 0
! !$omp end single

!  enddo remote_its

! #endif
! !$omp end parallel

!  !call update_hmax_remote

! #ifdef MPI
!  call finish_cell_exchange(irequestrecv,xsendbuf)
!  call sync_hmax_mpi
! #endif

!  !--reduce values
!  if (realviscosity .and. maxdvdx==maxp) then
!     stressmax = reduceall_mpi('max',stressmax)
!  endif
!  rhomax    = reduceall_mpi('max',rhomax)
!  rhomaxnow = rhomax

!  !--boundary particles are no longer treated as active
!  set_boundaries_to_active = .false.

!  if (realviscosity .and. maxdvdx==maxp .and. stressmax > 0. .and. iverbose > 0 .and. id==master) then
!     call warning('force','applying negative stress correction',var='max',val=-stressmax)
!  endif
! !
! !--warnings
! !
!  if (icall==1) then
!     if (nwarnup   > 0) call summary_variable('hupdn',iosumhup,0,real(nwarnup  ))
!     if (nwarndown > 0) call summary_variable('hupdn',iosumhdn,0,real(nwarndown))
!     if (iverbose  >=1) call reduce_and_print_warnings(nwarnup,nwarndown,nwarnroundoff)
!  endif
! !
! !--diagnostics
! !
!  if (icall==0 .or. icall==1) call reduce_and_print_neighbour_stats(np)

! end subroutine densityiterate

subroutine densityiterate(icall,npart,nactive,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                           fxyzu,fext,alphaind,gradh,rad,radprop)
 !$ use omputils
 use omp_lib
 use dim,       only:maxp,maxneigh,ndivcurlv,ndivcurlB,maxvxyzu,maxalpha, &
                     mhd_nonideal,nalpha,use_dust
 use io,        only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
 use linklist,  only:ifirstincell,ncells,get_neighbour_list,get_hmaxcell,get_cell_list,&
                     get_cell_location,set_hmaxcell,sync_hmax_mpi,update_hmax_remote,node_is_active, node 
 use part,      only:mhd,rhoh,dhdrho,rhoanddhdrho,&
                     ll,get_partinfo,iactive,&
                     hrho,iphase,igas,idust,iamgas,periodic,&
                     all_active,dustfrac,Bxyz,set_boundaries_to_active
#ifdef FINVSQRT
 use fastmath,  only:finvsqrt
#endif

 use mpiutils,  only:reduceall_mpi,barrier_mpi,reduce_mpi,reduceall_mpi
#ifdef MPI
 use stack,     only:reserve_stack,swap_stacks
 use stack,     only:stack_remote => dens_stack_1
 use stack,     only:stack_waiting => dens_stack_2
 use stack,     only:stack_redo => dens_stack_3
 use mpiderivs, only:send_cell,recv_cells,check_send_finished,init_cell_exchange, &
                     finish_cell_exchange,recv_while_wait,reset_cell_counters
#endif
 use timestep,  only:rhomaxnow
 use part,      only:ngradh
 use viscosity, only:irealvisc
 use io_summary,only:summary_variable,iosumhup,iosumhdn
 integer,      intent(in)    :: icall,npart,nactive
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(in)    :: vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real(kind=4), intent(out)   :: divcurlv(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real(kind=4), intent(out)   :: alphaind(:,:)
 real(kind=4), intent(out)   :: gradh(:,:)
 real,         intent(out)   :: stressmax
 real,         intent(in)    :: rad(:,:)
 real,         intent(inout) :: radprop(:,:)

 integer, save :: listneigh(maxneigh)
 real,   save :: xyzcache(3,isizecellcache)
!$omp threadprivate(xyzcache,listneigh)

 integer :: i,icell
 integer :: nneigh,np
 integer :: nwarnup,nwarndown,nwarnroundoff

 logical :: getdv,realviscosity,getdB,converged
 logical :: iactivei,iamgasi,iamdusti
 integer :: iamtypei

 logical :: remote_export(nprocs)

 real    :: rhomax

 integer                   :: npcell,istart,iend

 type(celldens),allocatable            :: celllist(:)
 type(celldens)            :: cell 

 logical                   :: redo_neighbours
 type(denstreestack)       :: stack(1000)
 type(denstreestacklocal)  :: stacklocal(1000)
 integer                   :: top,nodeindex1, nodeindex2, leftindex, rightindex, regnodeindex, splitnodeindex, numthreads
 integer, allocatable      :: istacklocal(:)
 integer                   :: k 
 logical, allocatable      :: threadworking(:)
  integer                   :: numnotconverged

#ifdef MPI
 integer                   :: j,k,l
 integer                   :: irequestsend(nprocs),irequestrecv(nprocs)
 type(celldens)            :: xrecvbuf(nprocs),xsendbuf
 integer                   :: mpiits,nlocal
 real                      :: ntotal
 logical                   :: iterations_finished
 logical                   :: do_export

  
 

 call init_cell_exchange(xrecvbuf,irequestrecv)
 stack_waiting%n = 0
 stack_remote%n = 0
 stack_redo%n = 0
#endif

  if (iverbose >= 3 .and. id==master) &
    write(iprint,*) ' cell cache =',isizecellcache,' neigh cache = ',isizeneighcache,' icall = ',icall

 if (icall==0 .or. icall==1) then
    call reset_neighbour_stats(nneightry,nneighact,maxneightry,maxneighact,ncalc,nrelink)
    nwarnup       = 0
    nwarndown     = 0
    nwarnroundoff = 0
    np = 0
 endif
 !
 ! flag for whether or not we need to calculate velocity derivatives
 ! whilst doing the density iterations (needed for viscosity switches
 ! and for physical viscosity)
 !
 realviscosity = (irealvisc > 0)
 getdv = ((maxalpha==maxp .or. ndivcurlv >= 4) .and. icall <= 1) .or. &
         (maxdvdx==maxp .and. (use_dust .or. realviscosity))
 if (getdv .and. ndivcurlv < 1) call fatal('densityiterate','divv not stored but it needs to be')
 getdB = (mhd .and. (ndivcurlB >= 4 .or. mhd_nonideal))

 if ( all_active ) stressmax  = 0.   ! condition is required for independent timestepping


#ifdef MPI
 ! number of local only cells
 nlocal = 0
 call reset_cell_counters

#endif

 rhomax = 0.0
 call get_cell_list(istart,iend)

 !top = 1
 !stack(top) % nodeindex1 = 1
 !stack(top) % interactions(1) = 1 

!print*,size(node)

!!$omp parallel default(none) shared(numthreads)
!     numthreads = omp_get_num_threads()
!!$omp end parallel

allocate(istacklocal(numthreads))
allocate(threadworking(numthreads))
allocate(celllist(size(node)))

threadworking = .false. 

istacklocal = 0

numnotconverged = 0



!$omp parallel default(none) &
!$omp shared(icall,istart,iend) &
!!$omp shared(ncells) &
!!$omp shared(ll) &
!$omp shared(ifirstincell) &
!$omp shared(xyzh) &
!$omp shared(vxyzu) &
!$omp shared(fxyzu) &
!$omp shared(fext) &
!$omp shared(gradh) &
!$omp shared(iphase) &
!$omp shared(Bevol) &
!$omp shared(divcurlv) &
!$omp shared(divcurlB) &
!$omp shared(alphaind) &
!$omp shared(dustfrac) &
!$omp shared(Bxyz) &
!$omp shared(dvdx) &
!$omp shared(id) &
!$omp shared(nprocs) &
!$omp shared(getdB) &
!$omp shared(getdv) &
!$omp shared(realviscosity) &
!$omp shared(iverbose) &
!$omp shared(iprint) &
!$omp shared(rad,radprop)&
#ifdef MPI
!$omp shared(xrecvbuf) &
!$omp shared(xsendbuf) &
!$omp shared(irequestrecv) &
!$omp shared(irequestsend) &
!$omp shared(stack_remote) &
!$omp shared(stack_waiting) &
!$omp shared(stack_redo) &
!$omp shared(iterations_finished) &
!$omp shared(mpiits) &
!$omp reduction(+:nlocal) &
!$omp private(do_export) &
!$omp private(j) &
!$omp private(k) &
!$omp private(l) &
!$omp private(ntotal) &
#endif
!$omp private(remote_export) &
!$omp private(nneigh) &
!$omp private(npcell) &
!$omp private(cell) &
!$omp private(iamgasi) &
!$omp private(iamtypei) &
!$omp private(iactivei) &
!$omp private(iamdusti) &
!$omp private(converged) &
!$omp private(redo_neighbours) &
!$omp reduction(+:ncalc) &
!$omp reduction(+:np) &
!$omp reduction(max:maxneighact) &
!$omp reduction(max:maxneightry) &
!$omp reduction(+:nneighact) &
!$omp reduction(+:nneightry) &
!$omp reduction(+:nrelink) &
!$omp reduction(+:stressmax) &
!$omp reduction(max:rhomax) &
!$omp private(i,k,stacklocal,istacklocal,nodeindex1,nodeindex2,leftindex,rightindex,splitnodeindex,regnodeindex) &
!$omp shared(stack,top,node,threadworking,ncells,celllist,numnotconverged)

!$omp do 
 do i=1,ncells
 if (ifirstincell(i) > 0) then 
     !print*, i 
      !stop 
      celllist(i)%icell                   = i
      celllist(i)%nits                    = 0
      celllist(i)%nneigh                  = 0
      celllist(i)%remote_export(1:nprocs) = remote_export
      celllist(i)%rhosums = 0
      call start_cell(celllist(i),iphase,xyzh,vxyzu,fxyzu,fext,Bevol,rad)
      call get_cell_location(i,celllist(i)%xpos,celllist(i)%xsizei,celllist(i)%rcuti)
      !print*,"cell xpos ", celllist(i) % xpos 
      !stop 
      call get_hmaxcell(i,celllist(i)%hmax)

 endif 
 enddo 

!$omp end do 

!$omp end parallel 
!print*, "Setup cells correctly"

!!$omp single 
call dual_traversal(celllist, icall,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                           fxyzu,fext,alphaind,gradh,rad,radprop)
!!$omp end single 
!stop 

!!$omp do 


do i=1,ncells
      if (ifirstincell(i) > 0) then 
            converged = .false.

            !print*,"celllist hmax: ", celllist(i) % hmax 
            call finish_cell(celllist(i),converged)
            call compute_hmax(celllist(i),redo_neighbours)
            !print*,"celllist hmax: ", celllist(i) % hmax 
            if (icall == 0) converged = .true.


             if (.not. converged) then 
                  call set_hmaxcell(celllist(i)%icell,celllist(i)%hmax)
                  !print*, celllist(i)%hmax 
            
                  celllist(i)%nneigh                  = 0
                  celllist(i)%remote_export(1:nprocs) = remote_export
                  celllist(i)%rhosums = 0
                 

            endif

      endif 

enddo 
!!$omp end do 

!!$omp single 
 call dual_traversal(celllist, icall,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                           fxyzu,fext,alphaind,gradh,rad,radprop)
!!$omp end single 

 !!$omp do 
      do i=1,ncells
            call finish_cell(celllist(i),converged)
            call compute_hmax(celllist(i),redo_neighbours)

            call store_results(icall,celllist(i),getdv,getdB,realviscosity,stressmax,xyzh,gradh,divcurlv,divcurlB,alphaind, &
                               dvdx,vxyzu,Bxyz,dustfrac,rhomax,nneightry,nneighact,maxneightry,maxneighact,np,ncalc,&
                               radprop)
      enddo 
 !!$omp end do 
! do i=1,ncells
! !print*,"I is: ", i
! if (i==65550) cycle 
! if (ifirstincell(i) > 0) then 


! !print*, "leaf node: ", i
! converged = .false. 
! !if (.not. converged) print*,i
! !do while (.not. converged)
!       !call = 1
!       print*,i

!       print*,"celllist hmax: ", celllist(i) % hmax 
!       call finish_cell(celllist(i),converged)
!       call compute_hmax(celllist(i),redo_neighbours)
!       print*,"celllist hmax: ", celllist(i) % hmax 
!       if (icall == 0) converged = .true.
!       print*, converged
!       print*, redo_neighbours
!       !stop 
      
!       if (.not. converged) then 
!            call set_hmaxcell(celllist(i)%icell,celllist(i)%hmax)
!            print*, celllist(i)%hmax 
!            stop 
!             celllist(i)%nneigh                  = 0
!             celllist(i)%remote_export(1:nprocs) = remote_export
!             celllist(i)%rhosums = 0
!             call dual_traversal(celllist, icall,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
!                            fxyzu,fext,alphaind,gradh,rad,radprop)
!       !      !i = 1
!       endif
! !     z print*," Not converged: ", i
!       !stop 
! !enddo 
! call store_results(icall,celllist(i),getdv,getdB,realviscosity,stressmax,xyzh,gradh,divcurlv,divcurlB,alphaind, &
!                                dvdx,vxyzu,Bxyz,dustfrac,rhomax,nneightry,nneighact,maxneightry,maxneighact,np,ncalc,&
!                                radprop)
! endif 
! enddo 

!  !$omp end do 

!!$omp end parallel 


!deallocate(celllist)


end subroutine densityiterate


subroutine dual_traversal(celllist, icall,xyzh,vxyzu,divcurlv,divcurlB,Bevol,stressmax,&
                           fxyzu,fext,alphaind,gradh,rad,radprop)
use omp_lib
 use dim,       only:maxp,maxneigh,ndivcurlv,ndivcurlB,maxvxyzu,maxalpha, &
                     mhd_nonideal,nalpha,use_dust
 use io,        only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs
 use linklist,  only:ifirstincell,ncells,get_neighbour_list,get_hmaxcell,get_cell_list,&
                     get_cell_location,set_hmaxcell,sync_hmax_mpi,update_hmax_remote,node_is_active, node 
 use part,      only:mhd,rhoh,dhdrho,rhoanddhdrho,&
                     ll,get_partinfo,iactive,&
                     hrho,iphase,igas,idust,iamgas,periodic,&
                     all_active,dustfrac,Bxyz,set_boundaries_to_active
 use viscosity, only:irealvisc
 type(celldens), intent(inout) ::  celllist(:)
 integer,      intent(in)    :: icall
 real,         intent(inout) :: xyzh(:,:)
 real,         intent(in)    :: vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,         intent(in)    :: Bevol(:,:)
 real(kind=4), intent(out)   :: divcurlv(:,:)
 real(kind=4), intent(out)   :: divcurlB(:,:)
 real(kind=4), intent(out)   :: alphaind(:,:)
 real(kind=4), intent(out)   :: gradh(:,:)
 real,         intent(out)   :: stressmax
 real,         intent(in)    :: rad(:,:)
 real,         intent(inout) :: radprop(:,:)
 
 integer, save :: listneigh(maxneigh)
 real,   save :: xyzcache(3,isizecellcache)
!$omp threadprivate(xyzcache,listneigh)

 integer :: i,icell
 integer :: nneigh,np
 integer :: nwarnup,nwarndown,nwarnroundoff

 logical :: getdv,realviscosity,getdB,converged
 logical :: iactivei,iamgasi,iamdusti
 integer :: iamtypei

 logical :: remote_export(nprocs)

 real    :: rhomax

 integer                   :: npcell,istart,iend

 !type(celldens)            :: celllist(size(node))
 type(celldens)            :: cell 

 logical                   :: redo_neighbours
 type(denstreestack)       :: stack(10000)
 type(denstreestacklocal)  :: stacklocal(10000)
 integer                   :: top,nodeindex1, nodeindex2, leftindex, rightindex, regnodeindex, splitnodeindex, numthreads
 integer, allocatable      :: istacklocal(:)
 integer                   :: k 
 logical, allocatable      :: threadworking(:)



!$omp parallel default(none) shared(numthreads)
   numthreads = omp_get_num_threads()
!$omp end parallel
print*, "numthreads: ",numthreads
!numthreads = 1
allocate(istacklocal(numthreads))
allocate(threadworking(numthreads))


!!$omp single 

! flag for whether or not we need to calculate velocity derivatives
 ! whilst doing the density iterations (needed for viscosity switches
 ! and for physical viscosity)
 !
 realviscosity = (irealvisc > 0)
 getdv = ((maxalpha==maxp .or. ndivcurlv >= 4) .and. icall <= 1) .or. &
         (maxdvdx==maxp .and. (use_dust .or. realviscosity))
 if (getdv .and. ndivcurlv < 1) call fatal('densityiterate','divv not stored but it needs to be')
 getdB = (mhd .and. (ndivcurlB >= 4 .or. mhd_nonideal))

 if ( all_active ) stressmax  = 0.   ! condition is required for independent timestepping

call get_cell_list(istart,iend)

 top = 1
 stack(top) % nodeindex1 = 1
 stack(top) % interactions(1) = 1 
 
 !print*, "entered single"
 !print*, "top is: "

!print*,size(node)


threadworking = .false. 

istacklocal = 0
!open(unit=13, file="neighbours.txt", position='append')

!!$omp end single 
!!$ call init_omp()
!print*,"finished single" 
!$omp parallel default(none) &
!$omp shared(icall,istart,iend) &
!!$omp shared(ncells) &
!!$omp shared(ll) &
!$omp shared(ifirstincell) &
!$omp shared(xyzh) &
!$omp shared(vxyzu) &
!$omp shared(fxyzu) &
!$omp shared(fext) &
!$omp shared(gradh) &
!$omp shared(iphase) &
!$omp shared(Bevol) &
!$omp shared(divcurlv) &
!$omp shared(divcurlB) &
!$omp shared(alphaind) &
!$omp shared(dustfrac) &
!$omp shared(Bxyz) &
!$omp shared(dvdx) &
!$omp shared(id) &
!$omp shared(nprocs) &
!$omp shared(getdB) &
!$omp shared(getdv) &
!$omp shared(realviscosity) &
!$omp shared(iverbose) &
!$omp shared(iprint) &
!$omp shared(rad,radprop)&
#ifdef MPI
!$omp shared(xrecvbuf) &
!$omp shared(xsendbuf) &
!$omp shared(irequestrecv) &
!$omp shared(irequestsend) &
!$omp shared(stack_remote) &
!$omp shared(stack_waiting) &
!$omp shared(stack_redo) &
!$omp shared(iterations_finished) &
!$omp shared(mpiits) &
!$omp reduction(+:nlocal) &
!$omp private(do_export) &
!$omp private(j) &
!$omp private(k) &
!$omp private(l) &
!$omp private(ntotal) &
#endif
!$omp private(remote_export) &
!$omp private(nneigh) &
!$omp private(npcell) &
!$omp private(cell) &
!$omp private(iamgasi) &
!$omp private(iamtypei) &
!$omp private(iactivei) &
!$omp private(iamdusti) &
!$omp private(converged) &
!$omp private(redo_neighbours) &
!$omp reduction(+:ncalc) &
!$omp reduction(+:np) &
!$omp reduction(max:maxneighact) &
!$omp reduction(max:maxneightry) &
!$omp reduction(+:nneighact) &
!$omp reduction(+:nneightry) &
!$omp reduction(+:nrelink) &
!$omp reduction(+:stressmax) &
!$omp reduction(max:rhomax) &
!$omp private(i,k,stacklocal,nodeindex1,nodeindex2,leftindex,rightindex,splitnodeindex,regnodeindex) &
!$omp shared(stack,istacklocal,top,node,threadworking,ncells,celllist)


do while (any(istacklocal > 0) .or. top > 0 .or. any(threadworking))

     !$ k=omp_get_thread_num() + 1
     !print*,"k is:", k

    if (istacklocal(k) > 0 ) then ! pop off local stack
        !!$omp critical 
        !print*,"istacklocal", istacklocal(k)
        nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
        nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2
        !print*, "Node indexes : ", nodeindex1, nodeindex2

        ! If the nodes are free then we can do work with them 
        !if (nodes(nodeindex1) % nodefree .and. nodes(nodeindex2) % nodefree) then 
          istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
        !  threadworking(k) = .false.
        !endif 
        !!$omp critical  

    else 
      !$omp critical (stack)
      if (top > 0) then 
        !!$OMP TASK 
       ! nodeindex1 = stack(top) % nodeindex1
        !nodeindex2 = stack(top) % nodeindex2
          ! If the nodes are free then we can do work with them 
        !if (nodes(nodeindex1) % nodefree .and. nodes(nodeindex2) % nodefree) then 
          !top = top - 1
          call global_to_local(stack,stacklocal,k,istacklocal,top)
          ! pop some work of the local stack 
          nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
          nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2
          !print*, "node indexes :", nodeindex1, nodeindex2
          istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
          !threadworking(k) = .false.
          !nodes(nodeindex1) % nodefree  = .false.
          !nodes(nodeindex2) % nodefree  = .false.
      else 
        threadworking(k) = .false.
      endif
      !$omp end critical (stack)
    endif 


  !nneigh = 0 
  !listneigh = 0

 
 !print*, "Node 1: ", nodeindex1, "Node 2: ", nodeindex2
  !print*, "Top: ", top 
  !print*, "istacklocal(k): ", istacklocal(k)

if (threadworking(k)) then 
      
      !print*, "thread: ", k, " working."
        ! if r^2 < rcut 
       if(trial_neighbour(node(nodeindex1), node(nodeindex2))) then
            !if (nodeindex2 == 131080 .and. nodeindex1 == 65550) stop 
            ! IF both leafs 
            if (ifirstincell(nodeindex1) > 0 .AND. ifirstincell(nodeindex2) > 0) then 
            !      if (nodeindex1 /= nodeindex2) then 

                  !!$omp critical (rhosum) 
                  call get_neighbour_list(nodeindex1,listneigh,nneigh,xyzh,xyzcache,isizecellcache,getj=.false., &
                              remote_export=remote_export,currentnodeindex=nodeindex1, interactnodeindex=nodeindex2)
                                                      !print*, "nneigh: ", nneigh
                        
                  
                  ! Put omp lock here 
                  !$omp critical (rhosum)
                  if (nneigh > 0) then 
                  call compute_cell(celllist(nodeindex1),listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext,xyzcache,rad)
                  endif  
                  !$omp end critical (rhosum)

                 

                  !if(nodeindex1 == 65550) then 
                    !do k=1,nneigh
                    !write(13,*) nodeindex2
                  !enddo
                  !endif 
                  !k = 1
                  !if identical split 
            else if (nodeindex1 == nodeindex2) then 
                  ! SPLIT THE IDENTICAL NODE

                  leftindex = node(nodeindex1) % leftchild
                  rightindex = node(nodeindex1) % rightchild 

                 !print*,"leftindex: ", leftindex
                  !print*,"rightindex: ", rightindex

                  !$omp critical (stack)

                  !if (leftindex /= 0 .and. rightindex /= 0) then 
                        top = top + 1
                        stack(top) % nodeindex1 = leftindex
                        stack(top) % interactions(1) = leftindex
                        stack(top) % interactions(2) = rightindex

                        top = top + 1 
                        stack(top) % nodeindex1 = rightindex
                        stack(top) % interactions(1) = rightindex
                        stack(top) % interactions(2) = leftindex
                  !endif 

                  !$omp end critical (stack)

            else 
                   ! FIND LARGER NODE 

                  if (node(nodeindex1) % size > node(nodeindex2) % size) then
                   !splitnode = nodes(nodeindex1)
                        splitnodeindex = nodeindex1
                        !regnode = nodes(nodeindex2)
                        regnodeindex = nodeindex2
                  else 
                        !splitnode = nodes(nodeindex2)
                        splitnodeindex = nodeindex2
                        !regnode = nodes(nodeindex1)
                        regnodeindex = nodeindex1
                  endif 

                  if (ifirstincell(splitnodeindex) <= 0) then 
                        ! split larger node and push interactions to stack 
                        leftindex = node(splitnodeindex) % leftchild
                        rightindex = node(splitnodeindex) % rightchild

                         !print*,"leftindex: ", leftindex
                         !print*,"rightindex: ", rightindex

                        ! If splitting node push to global stack 
                        if (nodeindex1 == splitnodeindex) then 
                              !$omp critical (stack)
                              if (leftindex /= 0) then 
                                    top = top + 1
                                    stack(top) % nodeindex1 = leftindex
                                    stack(top) % interactions(1) = regnodeindex
                              endif 
                              if (rightindex /= 0) then 
                                    top = top + 1
                                    stack(top) % nodeindex1 = rightindex
                                    stack(top) % interactions(1) = regnodeindex
                              endif 
                              !$omp end critical (stack)

                        else 
                  
                              if (leftindex /= 0) then 
                                    istacklocal(k) = istacklocal(k) + 1
                                    stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
                                    stacklocal(istacklocal(k)) % nodeindex2 = leftindex
                              endif

                              if (rightindex /= 0) then 
                                    istacklocal(k) = istacklocal(k) + 1
                                    stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
                                    stacklocal(istacklocal(k)) % nodeindex2 = rightindex
                              endif 

                        endif 

                          !if (nodeindex2 == 131080 .and. nodeindex1 == 32775) then 

                          !print*,"leftindex: ", leftindex
                          !print*,"rightindex: ", rightindex
                          !print*,"splitnodeindex,", splitnodeindex
                          !print*,"regnodeindex: ", regnodeindex
                          !endif 
                  ! IF we have reached a leafnode we need to split the other node 
                  else 
                        leftindex = node(regnodeindex) % leftchild
                        rightindex = node(regnodeindex) % rightchild

                         ! If splitting node push to global stack 
                        if (nodeindex1 == regnodeindex) then 
                              !$omp critical (stack)
                              if (leftindex /= 0) then 
                                    top = top + 1
                                    stack(top) % nodeindex1 = leftindex
                                    stack(top) % interactions(1) = splitnodeindex
                              endif 
                              if (rightindex /= 0) then 
                                    top = top + 1
                                    stack(top) % nodeindex1 = rightindex
                                    stack(top) % interactions(1) = splitnodeindex
                              endif 
                              !$omp end critical (stack)

                        else 
                  
                              if (leftindex /= 0) then 
                                    istacklocal(k) = istacklocal(k) + 1
                                    stacklocal(istacklocal(k)) % nodeindex1 = splitnodeindex
                                    stacklocal(istacklocal(k)) % nodeindex2 = leftindex
                              endif

                              if (rightindex /= 0) then 
                                    istacklocal(k) = istacklocal(k) + 1
                                    stacklocal(istacklocal(k)) % nodeindex1 = splitnodeindex
                                    stacklocal(istacklocal(k)) % nodeindex2 = rightindex
                              endif 

                        endif 

                  endif 


            endif 

     

      endif 

endif 





enddo 
!$omp end parallel

end subroutine dual_traversal


!----------------------------------------------------------------
!+
!  Internal subroutine that computes the contribution to
!  the density sums from a list of neighbours
!
!  MAKE SURE THIS ROUTINE IS INLINED BY THE COMPILER
!+
!----------------------------------------------------------------
 subroutine get_density_sums(i,xpartveci,hi,hi1,hi21,iamtypei,iamgasi,iamdusti,&
                                 listneigh,nneigh,nneighi,dxcache,xyzcache,rhosum,&
                                 ifilledcellcache,ifilledneighcache,getdv,getdB,&
                                 realviscosity,xyzh,vxyzu,Bevol,fxyzu,fext,ignoreself,&
                                 rad)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use kernel,   only:get_kernel,get_kernel_grav1
 use part,     only:iphase,iamgas,iamdust,iamtype,maxphase,ibasetype,igas,idust,rhoh,massoftype,iradxi
 use dim,      only:ndivcurlv,gravity,maxp,nalpha,use_dust,do_radiation
 integer,      intent(in)    :: i
 real,         intent(in)    :: xpartveci(:)
 real(kind=8), intent(in)    :: hi,hi1,hi21
 integer,      intent(in)    :: iamtypei
 logical,      intent(in)    :: iamgasi,iamdusti
 integer,      intent(in)    :: listneigh(:)
 integer,      intent(in)    :: nneigh
 integer,      intent(out)   :: nneighi
 real,         intent(inout) :: dxcache(:,:)
 real,         intent(in)    :: xyzcache(:,:)
 real,         intent(inout)   :: rhosum(:)
 logical,      intent(in)    :: ifilledcellcache,ifilledneighcache
 logical,      intent(in)    :: getdv,realviscosity
 logical,      intent(in)    :: getdB
 real,         intent(in)    :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,         intent(in)    :: Bevol(:,:)
 logical,      intent(in)    :: ignoreself
 real,         intent(in)    :: rad(:,:)
 integer(kind=1)             :: iphasej
 integer                     :: iamtypej
 integer                     :: j,n,iloc
 real                        :: dx,dy,dz,runix,runiy,runiz
 real                        :: rij2,rij,rij1,q2i,qi,q2prev,rij1grkern
 real                        :: wabi,grkerni,dwdhi,dphidhi
 real                        :: projv,dvx,dvy,dvz,dax,day,daz
 real                        :: projdB,dBx,dBy,dBz,fxi,fyi,fzi,fxj,fyj,fzj
 real                        :: rhoi, rhoj
 logical                     :: same_type,gas_gas,iamdustj
 real                        :: dradenij

 !rhosum(:) = 0.
 if (ignoreself) then
    nneighi = 1   ! self
 else
    nneighi = 0
 endif

 ! defaults for type determination
 ! these are determined from iphase if multiple phases are used
 same_type = .true.
 gas_gas   = .true.

 dphidhi   = 0.
 dx = 0. ! to avoid compiler warnings
 dy = 0.
 dz = 0.
 dvx = 0.
 dvy = 0.
 dvz = 0.
 if (nalpha > 1) then
    fxi = xpartveci(ifxi)
    fyi = xpartveci(ifyi)
    fzi = xpartveci(ifzi)
 endif

 loop_over_neigh: do n = 1,nneigh

    j = listneigh(n)
    !--do self contribution separately to avoid problems with 1/sqrt(0.)
    if ((ignoreself) .and. (j==i)) cycle loop_over_neigh

    if (ifilledneighcache .and. n <= isizeneighcache) then
       rij2 = dxcache(1,n)
    else
       if (ifilledcellcache .and. n <= isizecellcache) then
          ! positions from cache are already mod boundary
          dx = xpartveci(ixi) - xyzcache(1,n)
          dy = xpartveci(iyi) - xyzcache(2,n)
          dz = xpartveci(izi) - xyzcache(3,n)
       else
          dx = xpartveci(ixi) - xyzh(1,j)
          dy = xpartveci(iyi) - xyzh(2,j)
          dz = xpartveci(izi) - xyzh(3,j)
       endif
#ifdef PERIODIC
       if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
       if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
       if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
       rij2 = dx*dx + dy*dy + dz*dz
       if (n <= isizeneighcache) dxcache(1,n) = rij2
    endif

    q2i = rij2*hi21

    ! print*, rij2
    ! print*,hi21
    ! print*,q2i
    ! print*, radkern2
    ! print*, "dx,dy,dz"
    ! print*, dx, dy, dz 
    ! stop 
    !print*,i,' neighb',j,' #',n,' of ',nneigh,isizecellcache,' q2i=',dx,dy,dz,hi21
!
!--do interaction if r/h < compact support size
!
    if (q2i < radkern2) then
      !print*, "do"
       if (ifilledneighcache .and. n <= isizeneighcache) then
          q2prev = dxcache(2,n)
          if (q2prev < radkern2) then
             rij = dxcache(3,n)
          else
             rij = sqrt(rij2)
          endif
       else
          rij = sqrt(rij2)
       endif

       qi = rij*hi1
       !--kernel and gradient
       if (gravity) then
          call get_kernel_grav1(q2i,qi,wabi,grkerni,dphidhi)
       else
          call get_kernel(q2i,qi,wabi,grkerni)
       endif

       if (n <= isizeneighcache) then
          !   could possibly ONLY store q2i if q2i>q2prev so that
          !   the maximum number of sqrts are stored
          dxcache(2,n) = q2i ! if h decreasing we don
          dxcache(3,n) = rij
          dxcache(4,n) = grkerni
          !--can ONLY fill this on first pass
          if (.not.ifilledneighcache) then
             dxcache(5,n) = dx
             dxcache(6,n) = dy
             dxcache(7,n) = dz
          endif
       endif

       !
       ! Density, gradh and div v are only computed using
       ! neighbours of the same type
       !
       if (maxphase==maxp) then
          iphasej   = iphase(j)
          iamtypej  = iamtype(iphasej)
          iamdustj  = iamdust(iphasej)
          same_type = ((iamtypei == iamtypej) .or. (ibasetype(iamtypej)==iamtypei))
          gas_gas   = (iamgasi .and. same_type)  ! this ensure that boundary particles are included in gas_gas calculations
       endif

       sametype: if (same_type) then
          dwdhi = (-qi*grkerni - 3.*wabi)
          rhosum(irhoi)      = rhosum(irhoi) + wabi
          rhosum(igradhi)    = rhosum(igradhi) + dwdhi
          rhosum(igradsofti) = rhosum(igradsofti) + dphidhi
          nneighi            = nneighi + 1
          !
          ! calculate things needed for viscosity switches
          ! and real viscosity
          !
          if (getdv .or. getdB .or. do_radiation) then
             rij1 = 1./(rij + epsilon(rij))
             if (ifilledneighcache .and. n <= isizeneighcache) then
                !--dx,dy,dz are either in neighbour cache or have been calculated
                dx = dxcache(5,n)
                dy = dxcache(6,n)
                dz = dxcache(7,n)
             endif
             rij1grkern = rij1*grkerni
             runix = dx*rij1grkern
             runiy = dy*rij1grkern
             runiz = dz*rij1grkern

             if (getdv) then
                !--get dv and den
                dvx = xpartveci(ivxi) - vxyzu(1,j)
                dvy = xpartveci(ivyi) - vxyzu(2,j)
                dvz = xpartveci(ivzi) - vxyzu(3,j)
                projv = dvx*runix + dvy*runiy + dvz*runiz
                rhosum(idivvi) = rhosum(idivvi) + projv

                if (maxdvdx > 0 .or. ndivcurlv > 1 .or. nalpha > 1) then
                   rhosum(idvxdxi) = rhosum(idvxdxi) + dvx*runix
                   rhosum(idvxdyi) = rhosum(idvxdyi) + dvx*runiy
                   rhosum(idvxdzi) = rhosum(idvxdzi) + dvx*runiz
                   rhosum(idvydxi) = rhosum(idvydxi) + dvy*runix
                   rhosum(idvydyi) = rhosum(idvydyi) + dvy*runiy
                   rhosum(idvydzi) = rhosum(idvydzi) + dvy*runiz
                   rhosum(idvzdxi) = rhosum(idvzdxi) + dvz*runix
                   rhosum(idvzdyi) = rhosum(idvzdyi) + dvz*runiy
                   rhosum(idvzdzi) = rhosum(idvzdzi) + dvz*runiz

                   if (nalpha > 1 .and. gas_gas) then
                      !--divergence of acceleration for Cullen & Dehnen switch
                      fxj = fxyzu(1,j) + fext(1,j)
                      fyj = fxyzu(2,j) + fext(2,j)
                      fzj = fxyzu(3,j) + fext(3,j)
                      dax = fxi - fxj
                      day = fyi - fyj
                      daz = fzi - fzj

                      rhosum(idaxdxi) = rhosum(idaxdxi) + dax*runix
                      rhosum(idaxdyi) = rhosum(idaxdyi) + dax*runiy
                      rhosum(idaxdzi) = rhosum(idaxdzi) + dax*runiz
                      rhosum(idaydxi) = rhosum(idaydxi) + day*runix
                      rhosum(idaydyi) = rhosum(idaydyi) + day*runiy
                      rhosum(idaydzi) = rhosum(idaydzi) + day*runiz
                      rhosum(idazdxi) = rhosum(idazdxi) + daz*runix
                      rhosum(idazdyi) = rhosum(idazdyi) + daz*runiy
                      rhosum(idazdzi) = rhosum(idazdzi) + daz*runiz
                   endif
                   rhosum(irxxi) = rhosum(irxxi) - dx*runix
                   rhosum(irxyi) = rhosum(irxyi) - dx*runiy
                   rhosum(irxzi) = rhosum(irxzi) - dx*runiz
                   rhosum(iryyi) = rhosum(iryyi) - dy*runiy
                   rhosum(iryzi) = rhosum(iryzi) - dy*runiz
                   rhosum(irzzi) = rhosum(irzzi) - dz*runiz

                endif
             endif

             if (getdB .and. gas_gas) then
                ! we need B instead of B/rho, so used our estimated h here
                ! either it is close enough to be converged,
                ! or worst case it runs another iteration and re-calculates
                rhoi = rhoh(real(hi), massoftype(igas))
                rhoj = rhoh(xyzh(4,j), massoftype(igas))
                dBx = xpartveci(iBevolxi)*rhoi - Bevol(1,j)*rhoj
                dBy = xpartveci(iBevolyi)*rhoi - Bevol(2,j)*rhoj
                dBz = xpartveci(iBevolzi)*rhoi - Bevol(3,j)*rhoj
                projdB = dBx*runix + dBy*runiy + dBz*runiz

                ! difference operator of divB
                rhosum(idivBi) = rhosum(idivBi) + projdB

                rhosum(idBxdxi) = rhosum(idBxdxi) + dBx*runix
                rhosum(idBxdyi) = rhosum(idBxdyi) + dBx*runiy
                rhosum(idBxdzi) = rhosum(idBxdzi) + dBx*runiz
                rhosum(idBydxi) = rhosum(idBydxi) + dBy*runix
                rhosum(idBydyi) = rhosum(idBydyi) + dBy*runiy
                rhosum(idBydzi) = rhosum(idBydzi) + dBy*runiz
                rhosum(idBzdxi) = rhosum(idBzdxi) + dBz*runix
                rhosum(idBzdyi) = rhosum(idBzdyi) + dBz*runiy
                rhosum(idBzdzi) = rhosum(idBzdzi) + dBz*runiz
             endif

             if (do_radiation .and. gas_gas) then
                rhoi = rhoh(real(hi), massoftype(igas))
                rhoj = rhoh(xyzh(4,j), massoftype(igas))
                dradenij = rad(iradxi,j)*rhoj - xpartveci(iradxii)*rhoi
                rhosum(iradfxi) = rhosum(iradfxi) + dradenij*runix
                rhosum(iradfyi) = rhosum(iradfyi) + dradenij*runiy
                rhosum(iradfzi) = rhosum(iradfzi) + dradenij*runiz
             endif

          endif
       elseif (use_dust .and. (iamgasi  .and. iamdustj)) then
          iloc = irhodusti + iamtypej - idust
          rhosum(iloc) = rhosum(iloc) + wabi
       endif sametype

    elseif (n <= isizeneighcache) then
       ! q2prev > radkern2 from cache indicates rij has NOT been calculated for this pair
       dxcache(2,n) = q2i
       if (.not.ifilledneighcache) then
          dxcache(5,n) = dx
          dxcache(6,n) = dy
          dxcache(7,n) = dz
       endif
    endif
 enddo loop_over_neigh

end subroutine get_density_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract the matrix used in exact linear
!  interpolations from the summations calculated during
!  the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_rmatrix_from_sums(rhosum,denom,rmatrix,idone)
 real,    intent(in)  :: rhosum(:)
 real,    intent(out) :: denom
 real,    intent(out) :: rmatrix(6)
 logical, intent(out) :: idone
 real :: rxxi,rxyi,rxzi,ryyi,ryzi,rzzi

 rxxi = rhosum(irxxi)
 rxyi = rhosum(irxyi)
 rxzi = rhosum(irxzi)
 ryyi = rhosum(iryyi)
 ryzi = rhosum(iryzi)
 rzzi = rhosum(irzzi)

 denom = rxxi*ryyi*rzzi + 2.*rxyi*rxzi*ryzi &
        - rxxi*ryzi*ryzi - ryyi*rxzi*rxzi - rzzi*rxyi*rxyi

 rmatrix(1) = ryyi*rzzi - ryzi*ryzi    ! xx
 rmatrix(2) = rxzi*ryzi - rzzi*rxyi    ! xy
 rmatrix(3) = rxyi*ryzi - rxzi*ryyi    ! xz
 rmatrix(4) = rzzi*rxxi - rxzi*rxzi    ! yy
 rmatrix(5) = rxyi*rxzi - rxxi*ryzi    ! yz
 rmatrix(6) = rxxi*ryyi - rxyi*rxyi    ! zz
 idone = .true.

 return
end subroutine calculate_rmatrix_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract div and curl v from the sums
!  calculated during the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_divcurlv_from_sums(rhosum,termnorm,divcurlvi,xi_limiter,ndivcurlv,denom,rmatrix)
 use part, only:nalpha
 integer, intent(in)  :: ndivcurlv
 real,    intent(in)  :: rhosum(:),denom,rmatrix(6)
 real,    intent(in)  :: termnorm
 real,    intent(out) :: divcurlvi(5),xi_limiter
 real :: div_a
 real :: gradaxdx,gradaxdy,gradaxdz,gradaydx,gradaydy,gradaydz,gradazdx,gradazdy,gradazdz
 real :: ddenom,gradvxdxi,gradvxdyi,gradvxdzi
 real :: gradvydxi,gradvydyi,gradvydzi,gradvzdxi,gradvzdyi,gradvzdzi
 real :: dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi
 real :: fac,traceS,txy,txz,tyz,txx,tyy,tzz!,Ri
 logical, parameter :: use_exact_linear = .true.

 !--divergence of the velocity field
 if (ndivcurlv >= 1) divcurlvi(1) = -rhosum(idivvi)*termnorm

 !--curl of the velocity field
 if (ndivcurlv >= 4) then
    divcurlvi(2) = -(rhosum(idvzdyi) - rhosum(idvydzi))*termnorm
    divcurlvi(3) = -(rhosum(idvxdzi) - rhosum(idvzdxi))*termnorm
    divcurlvi(4) = -(rhosum(idvydxi) - rhosum(idvxdyi))*termnorm
 endif

 !--time derivative of div v, needed for Cullen-Dehnen switch
 if (nalpha >= 2) then
    !--Divvdt For switch
    if (use_exact_linear) then
       ddenom = 1./denom
       call exactlinear(gradaxdx,gradaxdy,gradaxdz,rhosum(idaxdxi),rhosum(idaxdyi),rhosum(idaxdzi),rmatrix,ddenom)
       call exactlinear(gradaydx,gradaydy,gradaydz,rhosum(idaydxi),rhosum(idaydyi),rhosum(idaydzi),rmatrix,ddenom)
       call exactlinear(gradazdx,gradazdy,gradazdz,rhosum(idazdxi),rhosum(idazdyi),rhosum(idazdzi),rmatrix,ddenom)
       div_a = -(gradaxdx + gradaydy + gradazdz)

       call exactlinear(gradvxdxi,gradvxdyi,gradvxdzi, &
                         rhosum(idvxdxi),rhosum(idvxdyi),rhosum(idvxdzi),rmatrix,ddenom)
       call exactlinear(gradvydxi,gradvydyi,gradvydzi, &
                         rhosum(idvydxi),rhosum(idvydyi),rhosum(idvydzi),rmatrix,ddenom)
       call exactlinear(gradvzdxi,gradvzdyi,gradvzdzi, &
                         rhosum(idvzdxi),rhosum(idvzdyi),rhosum(idvzdzi),rmatrix,ddenom)

       dvxdxi = -gradvxdxi
       dvxdyi = -gradvxdyi
       dvxdzi = -gradvxdzi
       dvydxi = -gradvydxi
       dvydyi = -gradvydyi
       dvydzi = -gradvydzi
       dvzdxi = -gradvzdxi
       dvzdyi = -gradvzdyi
       dvzdzi = -gradvzdzi
    else
       div_a = -termnorm*(rhosum(idaxdxi) + rhosum(idaydyi) + rhosum(idazdzi))
       dvxdxi = -termnorm*rhosum(idvxdxi)
       dvxdyi = -termnorm*rhosum(idvxdyi)
       dvxdzi = -termnorm*rhosum(idvxdzi)
       dvydxi = -termnorm*rhosum(idvydxi)
       dvydyi = -termnorm*rhosum(idvydyi)
       dvydzi = -termnorm*rhosum(idvydzi)
       dvzdxi = -termnorm*rhosum(idvzdxi)
       dvzdyi = -termnorm*rhosum(idvzdyi)
       dvzdzi = -termnorm*rhosum(idvzdzi)
    endif
    divcurlvi(5) = div_a - (dvxdxi**2 + dvydyi**2 + dvzdzi**2 + &
                             2.*(dvxdyi*dvydxi + dvxdzi*dvzdxi + dvydzi*dvzdyi))
    !if (divcurlvi(1) < 0.) then
    !   Ri = -1.
    !else
    !   Ri = 1.
    !endif
    txx = dvxdxi - divcurlvi(1)/3.
    tyy = dvydyi - divcurlvi(1)/3.
    tzz = dvzdzi - divcurlvi(1)/3.
    txy = 0.5*(dvxdyi + dvydxi)
    txz = 0.5*(dvxdzi + dvzdxi)
    tyz = 0.5*(dvydzi + dvzdyi)
    fac    = max(-divcurlvi(1),0.)**2 !(2.*(1. - Ri)**4*divcurlvi(1))**2
    !traceS = txx**2 + tyy**2 + tzz**2 + 2.*(txy**2 + txz**2 + tyz**2)
    !traceS = txy**2 + txz**2 + tyz**2
    traceS = (dvzdyi - dvydzi)**2 + (dvxdzi - dvzdxi)**2 + (dvydxi - dvxdyi)**2
    if (fac + traceS > 0.) then
       xi_limiter = fac/(fac + traceS)
    else
       xi_limiter = 1.
    endif
 else
    xi_limiter = 1.
 endif

end subroutine calculate_divcurlv_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract div, curl and grad B from the sums
!  calculated during the density loop
!+
!----------------------------------------------------------------
pure subroutine calculate_divcurlB_from_sums(rhosum,termnorm,divcurlBi,gradBi,ndivcurlB)
 integer, intent(in)  :: ndivcurlB
 real,    intent(in)  :: rhosum(:)
 real,    intent(in)  :: termnorm
 real,    intent(out) :: divcurlBi(ndivcurlB),gradBi

 ! we need these for adaptive resistivity switch
 if (ndivcurlB >= 1) divcurlBi(1) = -rhosum(idivBi)*termnorm
 if (ndivcurlB >= 4) then
    divcurlBi(2) = -(rhosum(idBzdyi) - rhosum(idBydzi))*termnorm
    divcurlBi(3) = -(rhosum(idBxdzi) - rhosum(idBzdxi))*termnorm
    divcurlBi(4) = -(rhosum(idBydxi) - rhosum(idBxdyi))*termnorm
 endif
 gradBi = rhosum(idBxdxi) * rhosum(idBxdxi) &
        + rhosum(idBxdyi) * rhosum(idBxdyi) &
        + rhosum(idBxdzi) * rhosum(idBxdzi) &
        + rhosum(idBydxi) * rhosum(idBydxi) &
        + rhosum(idBydyi) * rhosum(idBydyi) &
        + rhosum(idBydzi) * rhosum(idBydzi) &
        + rhosum(idBzdxi) * rhosum(idBzdxi) &
        + rhosum(idBzdyi) * rhosum(idBzdyi) &
        + rhosum(idBzdzi) * rhosum(idBzdzi)
 gradBi = sqrt(termnorm * termnorm * gradBi)

end subroutine calculate_divcurlB_from_sums

!----------------------------------------------------------------
!+
!  Internal utility to extract velocity gradients from summations
!  calculated during the density loop.
!+
!----------------------------------------------------------------
subroutine calculate_strain_from_sums(rhosum,termnorm,denom,rmatrix,dvdx)
 real, intent(in)  :: rhosum(:)
 real, intent(in)  :: termnorm,denom
 real, intent(in)  :: rmatrix(6)
 real, intent(out) :: dvdx(9)

 real :: ddenom,gradvxdxi,gradvxdyi,gradvxdzi
 real :: gradvydxi,gradvydyi,gradvydzi,gradvzdxi,gradvzdyi,gradvzdzi
 real :: dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi

! if (abs(denom) > tiny(denom)) then ! do exact linear first derivatives
 if (.false.) then ! do exact linear first derivatives
    ddenom = 1./denom
    call exactlinear(gradvxdxi,gradvxdyi,gradvxdzi, &
                     rhosum(idvxdxi),rhosum(idvxdyi),rhosum(idvxdzi),rmatrix,ddenom)
    call exactlinear(gradvydxi,gradvydyi,gradvydzi, &
                     rhosum(idvydxi),rhosum(idvydyi),rhosum(idvydzi),rmatrix,ddenom)
    call exactlinear(gradvzdxi,gradvzdyi,gradvzdzi, &
                     rhosum(idvzdxi),rhosum(idvzdyi),rhosum(idvzdzi),rmatrix,ddenom)

    !print*,'dvxdxi = ',-rhosum(idvxdxi)*termnorm,gradvxdxi
    dvxdxi = -gradvxdxi
    dvxdyi = -gradvxdyi
    dvxdzi = -gradvxdzi
    dvydxi = -gradvydxi
    dvydyi = -gradvydyi
    dvydzi = -gradvydzi
    dvzdxi = -gradvzdxi
    dvzdyi = -gradvzdyi
    dvzdzi = -gradvzdzi
 else

    !--these make rho*dv/dx_i
    dvxdxi = -rhosum(idvxdxi)*termnorm
    dvxdyi = -rhosum(idvxdyi)*termnorm
    dvxdzi = -rhosum(idvxdzi)*termnorm
    dvydxi = -rhosum(idvydxi)*termnorm
    dvydyi = -rhosum(idvydyi)*termnorm
    dvydzi = -rhosum(idvydzi)*termnorm
    dvzdxi = -rhosum(idvzdxi)*termnorm
    dvzdyi = -rhosum(idvzdyi)*termnorm
    dvzdzi = -rhosum(idvzdzi)*termnorm
 endif

 dvdx(:) = (/dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi,dvzdxi,dvzdyi,dvzdzi/)

end subroutine calculate_strain_from_sums

!----------------------------------------------------------------
!+
!  Internal subroutine to get maximum of stress tensor
!  (to avoid tensile instability)
!+
!----------------------------------------------------------------
pure subroutine get_max_stress(dvdx,divvi,rho1i,stressmax,shearvisc,bulkvisc)
 use part, only:strain_from_dvdx
 real, intent(in)    :: dvdx(9), divvi, rho1i, shearvisc, bulkvisc
 real, intent(inout) :: stressmax
 real :: strainmax,stressiso,strain(6)

 strain = strain_from_dvdx(dvdx)

 ! shearvisc = eta/rho, so this is eta/rho**2
 strainmax = -shearvisc*rho1i*maxval(strain) ! 1/rho*L^2/T*1/L*L/T = 1/rho*L^2/T^2
 stressiso = (2./3.*shearvisc - bulkvisc)*divvi*rho1i ! NB: divv -ve at high Mach no.

 ! we initialise stressmax to zero, as for the purpose of preventing the
 ! tensile instability we only care if the total stress is negative
 ! if stress tensor is positive, don't need correction (stressmax=0)
 stressmax = max(stressmax,-(stressiso + strainmax))
 stressmax = 0.

end subroutine get_max_stress

!----------------------------------------------------------------
!+
!  Internal subroutine that inverts the matrix to get an
!  exact linear derivative
!+
!----------------------------------------------------------------
pure subroutine exactlinear(gradAx,gradAy,gradAz,dAx,dAy,dAz,rmatrix,ddenom)
 real, intent(out) :: gradAx,gradAy,gradAz
 real, intent(in)  :: dAx,dAy,dAz
 real, intent(in)  :: rmatrix(6)
 real, intent(in)  :: ddenom
 !
 !--we return the gradient as the following matrix inversion:
 !  gradAx =(dAx*termxx + dAy*termxy + dAz*termxz)*ddenom
 !  gradAy =(dAx*termxy + dAy*termyy + dAz*termyz)*ddenom
 !  gradAz =(dAx*termxz + dAy*termyz + dAz*termzz)*ddenom
 !
 gradAx =(dAx*rmatrix(1) + dAy*rmatrix(2) + dAz*rmatrix(3))*ddenom
 gradAy =(dAx*rmatrix(2) + dAy*rmatrix(4) + dAz*rmatrix(5))*ddenom
 gradAz =(dAx*rmatrix(3) + dAy*rmatrix(5) + dAz*rmatrix(6))*ddenom

end subroutine exactlinear

!-------------------------------------------------------------------------------
!+
!  function to return alphaloc from known values of d(divv)/dt and sound speed
!  for use in Cullen & Dehnen (2010) switch
!+
!-------------------------------------------------------------------------------
pure real function get_alphaloc(divvdti,spsoundi,hi,xi_limiter,alphamin,alphamax)
 !use kernel, only:radkern
 real, intent(in) :: divvdti,spsoundi,hi,xi_limiter,alphamin,alphamax
 real :: source
 real :: temp

 source = 10.*hi**2*xi_limiter*max(-divvdti,0.)
 temp = spsoundi**2 !+ source
 if (temp > epsilon(temp)) then
    get_alphaloc = max(min(source/temp,alphamax),alphamin)
 else
    get_alphaloc = alphamin
 endif

end function get_alphaloc

!----------------------------------------------------------------
!+
!  subroutine to reduce and print warnings across processors
!  related to h-rho iterations
!+fxyzu
!----------------------------------------------------------------
subroutine reduce_and_print_warnings(nwarnup,nwarndown,nwarnroundoff)
 use mpiutils, only:reduce_mpi
 use io,       only:id,master,iprint
 integer, intent(inout) :: nwarnup,nwarndown,nwarnroundoff

 nwarnup       = int(reduce_mpi('+',nwarnup))
 nwarndown     = int(reduce_mpi('+',nwarndown))
 nwarnroundoff = int(reduce_mpi('+',nwarnroundoff))

#ifndef NOWARNRESTRICTEDHJUMP
 if (id==master .and. nwarnup > 0) then
    write(iprint,*) ' WARNING: restricted h jump (up) ',nwarnup,' times'
 endif
 if (id==master .and. nwarndown > 0) then
    write(iprint,*) ' WARNING: restricted h jump (down) ',nwarndown,' times'
 endif
#endif
 if (id==master .and. nwarnroundoff > 0) then
    write(iprint,*) ' WARNING: denom in exact linear gradients zero on ',nwarnroundoff,' particles'
 endif

end subroutine reduce_and_print_warnings

!----------------------------------------------------------------
!+
!  query function to return neighbour statistics
!  (must be called *after* density evaluation
!   and will be correct ONLY on the master thread)
!+
!----------------------------------------------------------------
subroutine get_neighbour_stats(trialmean,actualmean,maxtrial,maxactual,nrhocalc,nactualtot)
 real,            intent(out) :: trialmean,actualmean
 integer,         intent(out) :: maxtrial,maxactual
 integer(kind=8), intent(out) :: nrhocalc,nactualtot

 if (nptot > 0) then
    trialmean    = nneightry/real(nptot)
    actualmean   = nneighact/real(nptot)
    maxtrial     = int(maxneightry)
    maxactual    = maxneighact
    nrhocalc     = ncalc
    nactualtot   = nneighact
 else ! densityforce has not been called
    trialmean = -1; actualmean   = -1
    maxtrial  = -1; maxactual    = -1
    nrhocalc  = -1; nactualtot   = -1
 endif

end subroutine get_neighbour_stats

subroutine reset_neighbour_stats(nneightry,nneighact,maxneightry,maxneighact,ncalc,nrelink)
 integer,         intent(out) :: maxneighact,nrelink
 integer(kind=8), intent(out) :: ncalc,nneightry,nneighact,maxneightry

 nneightry = 0
 nneighact = 0
 maxneightry = 0
 maxneighact = 0
 ncalc = 0_8
 nneighact = 0
 nrelink = 0_8

end subroutine reset_neighbour_stats

!----------------------------------------------------------------
!+
!  function to collate neighbour-finding statistics across
!  processors and print the results
!+
!----------------------------------------------------------------
subroutine reduce_and_print_neighbour_stats(np)
 use mpiutils, only:reduce_mpi,barrier_mpi
 use io,       only:iprint,id,master,iverbose
 integer, intent(in) :: np

 call barrier_mpi()
 nptot = reduce_mpi('+',np)
 call barrier_mpi()
 nneightry   = reduce_mpi('+',nneightry)
 call barrier_mpi()
 nneighact   = reduce_mpi('+',nneighact)
 call barrier_mpi()
 maxneightry = reduce_mpi('max',maxneightry)
 call barrier_mpi()
 maxneighact = int(reduce_mpi('max',maxneighact))
 call barrier_mpi()
 nrelink     = int(reduce_mpi('+',nrelink))
 call barrier_mpi()
 ncalc       = reduce_mpi('+',ncalc)

 if (id==master .and. iverbose >= 2 .and. nptot > 0 .and. nneighact > 0) then
    write(iprint,"(1x,a,f11.2,2(a,f7.2))") 'trial neigh mean  :',nneightry/real(nptot), &
                 ', real neigh mean = ',nneighact/real(nptot), &
                 ' ratio try/act= ',nneightry/real(nneighact)
    write(iprint,"(1x,a,i11,a,i8)")   'trial neigh max   :',maxneightry,', max real neigh = ',maxneighact
    write(iprint,"(1x,a,i11,a,f7.3)") 'n neighbour calls :',nrelink, ', mean per part   = ',nrelink/real(nptot) + 1
    write(iprint,"(1x,a,i11,a,f7.3)") 'n density calcs   :',ncalc,', mean per part   = ',ncalc/real(nptot)
 endif

end subroutine reduce_and_print_neighbour_stats
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine compute_cell(cell,listneigh,nneigh,getdv,getdB,Bevol,xyzh,vxyzu,fxyzu,fext, &
                             xyzcache,rad)
 use dim,         only:maxvxyzu
 use part,        only:get_partinfo,iamgas,mhd,igas,maxphase
 use viscosity,   only:irealvisc
#ifdef MPI
 use io,          only:id
#endif

 type(celldens),  intent(inout)  :: cell

 integer,         intent(in)     :: listneigh(:)
 integer,         intent(in)     :: nneigh
 logical,         intent(in)     :: getdv
 logical,         intent(in)     :: getdB
 real,            intent(in)     :: Bevol(:,:)
 real,            intent(in)     :: xyzh(:,:),vxyzu(:,:),fxyzu(:,:),fext(:,:)
 real,            intent(in)     :: xyzcache(3,isizecellcache)
 real,            intent(in)     :: rad(:,:)

 real                            :: dxcache(7,isizeneighcache)

 real(kind=8)                    :: hi
 real(kind=8)                    :: hi1,hi21,hi31,hi41

 integer                         :: iamtypei
 logical                         :: iactivei,iamgasi,iamdusti

 logical                         :: realviscosity
 logical                         :: ignoreself
 integer                         :: nneighi
 integer                         :: i,lli

 realviscosity = (irealvisc > 0)

 over_parts: do i = 1,cell%npcell
    ! skip particles already converged
    if (cell%converged(i)) cycle over_parts

    lli = iorder(cell%arr_index(i))
    ! note: only active particles have been sent here
    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
    endif

    hi    = cell%h(i)
    hi1   = 1./hi
    hi21  = hi1*hi1
    hi31  = hi1*hi21
    hi41  = hi21*hi21

#ifdef MPI
    if (cell%owner == id) then
       ignoreself = .true.
    else
       ignoreself = .false.
    endif
#else
    ignoreself = .true.
#endif

    call get_density_sums(lli,cell%xpartvec(:,i),hi,hi1,hi21,iamtypei,iamgasi,iamdusti,&
                          listneigh,nneigh,nneighi,dxcache,xyzcache,cell%rhosums(:,i),&
                          .true.,.false.,getdv,getdB,realviscosity,&
                          xyzh,vxyzu,Bevol,fxyzu,fext,ignoreself,rad)

    cell%nneightry = nneigh
    cell%nneigh(i) = nneighi

 enddo over_parts

end subroutine compute_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
pure subroutine compute_hmax(cell,redo_neighbours)
 use kernel, only:radkern
 type(celldens), intent(inout) :: cell
 logical,         intent(out)  :: redo_neighbours
 real                          :: hmax_old,hmax

 redo_neighbours = .false.
 if (cell%npcell > 0) then
    hmax_old = cell%hmax
    hmax     = 1.01*maxval(cell%h(1:cell%npcell))
    if (hmax > hmax_old) redo_neighbours = .true.
    cell%hmax  = hmax
    cell%rcuti = radkern*hmax
 endif
end subroutine compute_hmax
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine start_cell(cell,iphase,xyzh,vxyzu,fxyzu,fext,Bevol,rad)
 use io,          only:fatal
 use dim,         only:maxp,maxvxyzu,do_radiation
 use part,        only:maxphase,get_partinfo,mhd,igas,iamgas,&
                       iamboundary,ibasetype,iradxi

 type(celldens),     intent(inout) :: cell
 integer(kind=1),    intent(in)    :: iphase(:)
 real,               intent(in)    :: xyzh(:,:)
 real,               intent(in)    :: vxyzu(:,:)
 real,               intent(in)    :: fxyzu(:,:)
 real,               intent(in)    :: fext(:,:)
 real,               intent(in)    :: Bevol(:,:)
 real,               intent(in)    :: rad(:,:)

 integer :: i,ip
 integer :: iamtypei
 logical :: iactivei,iamgasi,iamdusti

 cell%npcell = 0
 over_parts: do ip = inoderange(1,cell%icell),inoderange(2,cell%icell)
    if (ip <= 0) exit over_parts
    i = iorder(ip)

    if (i < 0) then
       cycle over_parts
    endif

    if (maxphase==maxp) then
       call get_partinfo(iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamdusti = .false.
       iamgasi  = .true.
    endif
    if (.not.iactivei) then ! skip boundary particles + inactive particles
       cycle over_parts
    endif

    cell%npcell = cell%npcell + 1
    cell%converged(cell%npcell) = .false.

    cell%arr_index(cell%npcell)               = ip
    cell%iphase(cell%npcell)                  = iphase(i)

    cell%xpartvec(ixi,cell%npcell)            = xyzh(1,i)
    cell%xpartvec(iyi,cell%npcell)            = xyzh(2,i)
    cell%xpartvec(izi,cell%npcell)            = xyzh(3,i)

    cell%h(cell%npcell)                       = xyzh(4,i)
    cell%h_old(cell%npcell)                   = xyzh(4,i)

    cell%xpartvec(ivxi,cell%npcell)           = vxyzu(1,i)
    cell%xpartvec(ivyi,cell%npcell)           = vxyzu(2,i)
    cell%xpartvec(ivzi,cell%npcell)           = vxyzu(3,i)

    if (maxvxyzu >= 4) then
       cell%xpartvec(ieni,cell%npcell)        = vxyzu(4,i)
    endif

    cell%xpartvec(ifxi,cell%npcell)           = fxyzu(1,i) + fext(1,i)
    cell%xpartvec(ifyi,cell%npcell)           = fxyzu(2,i) + fext(2,i)
    cell%xpartvec(ifzi,cell%npcell)           = fxyzu(3,i) + fext(3,i)

    if (mhd) then
       if (iamgasi) then
          cell%xpartvec(iBevolxi,cell%npcell) = Bevol(1,i)
          cell%xpartvec(iBevolyi,cell%npcell) = Bevol(2,i)
          cell%xpartvec(iBevolzi,cell%npcell) = Bevol(3,i)
          cell%xpartvec(ipsi,cell%npcell)     = Bevol(4,i)
       else
          cell%xpartvec(iBevolxi:ipsi,cell%npcell)   = 0. ! to avoid compiler warning
       endif
    endif

    if (do_radiation) cell%xpartvec(iradxii,cell%npcell) = rad(iradxi,i)
    !print*,"Cell h: ", cell % h
    !stop 

 enddo over_parts

end subroutine start_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine finish_cell(cell,cell_converged)
 use io,       only:iprint,fatal
 use part,     only:get_partinfo,iamgas,maxphase,massoftype,igas,hrho
 use options,  only:tolh

 type(celldens),  intent(inout) :: cell
 logical,         intent(out)   :: cell_converged
 real                           :: rhosum(maxrhosum)
 real                           :: dhdrhoi,rhohi,omegai
 real                           :: rhoi
 real(kind=8)                   :: gradhi
 real                           :: func,dfdh1,hi,hi_old,hnew
 real                           :: pmassi, xyzh(4)
 integer                        :: i,iamtypei !,nwarnup,nwarndown
 logical                        :: iactivei,iamgasi,iamdusti,converged

 cell_converged = .true.
 over_parts: do i = 1,cell%npcell
    if (cell%converged(i)) cycle over_parts
    cell%nits = cell%nits + 1

    hi = cell%h(i)
    hi_old = cell%h_old(i)
    rhosum = cell%rhosums(:,i)

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
    endif
    !if (.not.iactivei) print*,' ERROR: should be no inactive particles here',iamtypei,iactivei

    pmassi = massoftype(iamtypei)

    call finish_rhosum(rhosum,pmassi,hi,.true.,rhoi=rhoi,rhohi=rhohi,&
                       gradhi=gradhi,dhdrhoi_out=dhdrhoi,omegai_out=omegai)

    func = rhohi - rhoi
    if (omegai > tiny(omegai)) then
       dfdh1 = dhdrhoi/omegai
    else
       dfdh1 = dhdrhoi/abs(omegai + epsilon(omegai))
    endif
    hnew = hi - func*dfdh1
    if (hnew > 1.2*hi) then
       ! nwarnup   = nwarnup + 1
       hnew      = 1.2*hi
    elseif (hnew < 0.8*hi) then
       ! nwarndown = nwarndown + 1
       hnew      = 0.8*hi
    endif

    converged = ((abs(hnew-hi)/hi_old) < tolh .and. omegai > 0. .and. hi > 0.)
    cell%converged(i) = converged
    if (cell_converged) cell_converged = converged

    if ((.not. converged) .and. (cell%nits >= maxdensits)) then
       xyzh(1) = cell%xpartvec(ixi,i)
       xyzh(2) = cell%xpartvec(iyi,i)
       xyzh(3) = cell%xpartvec(izi,i)
       write(iprint,*) 'ERROR: density iteration failed after ',cell%nits,' iterations'
       write(iprint,*) 'hnew = ',hnew,' hi_old = ',hi_old,' nneighi = ',cell%nneigh(i)
       write(iprint,*) 'rhoi = ',rhoi,' gradhi = ',gradhi
       write(iprint,*) 'error = ',abs(hnew-hi)/hi_old,' tolh = ',tolh
       write(iprint,*) 'itype = ',iamtypei
       write(iprint,*) 'x,y,z = ',xyzh(1:3)
       call fatal('densityiterate','could not converge in density',iorder(cell%arr_index(i)),'error',abs(hnew-hi)/hi_old)
    endif

    if (converged) then
       cell%h(i) = hi
    else
       cell%h(i) = hnew
    endif

 enddo over_parts

end subroutine finish_cell
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
pure subroutine finish_rhosum(rhosum,pmassi,hi,iterating,rhoi,rhohi,gradhi,gradsofti,dhdrhoi_out,omegai_out)
 use part,  only:rhoh,dhdrho
 real,          intent(in)              :: rhosum(maxrhosum)
 real,          intent(in)              :: pmassi
 real,          intent(in)              :: hi
 logical,       intent(in)              :: iterating !false for the last bit where we are computing the final result
 real,          intent(out)             :: rhoi
 real(kind=8),  intent(out)             :: gradhi
 real,          intent(out),  optional  :: rhohi
 real(kind=8),  intent(out),  optional  :: gradsofti
 real,          intent(out),  optional  :: dhdrhoi_out
 real,          intent(out),  optional  :: omegai_out

 real           :: omegai,dhdrhoi
 real(kind=8)   :: hi1,hi21,hi31,hi41

 hi1   = 1./hi
 hi21  = hi1*hi1
 hi31  = hi1*hi21
 hi41  = hi21*hi21

 rhoi   = cnormk*pmassi*(rhosum(irhoi) + wab0)*hi31
 gradhi = cnormk*pmassi*(rhosum(igradhi) + gradh0)*hi41

 dhdrhoi = dhdrho(hi,pmassi)
 omegai = 1. - dhdrhoi*gradhi
 gradhi = 1./omegai

 if (iterating) then
    rhohi = rhoh(hi,pmassi)
    dhdrhoi_out = dhdrhoi
    omegai_out = omegai
 else
    gradsofti = pmassi*(rhosum(igradsofti) + dphidh0)*hi21 ! NB: no cnormk in gradsoft
    gradsofti = gradsofti*dhdrhoi
 endif

end subroutine finish_rhosum
!--------------------------------------------------------------------------
!+
!--------------------------------------------------------------------------
subroutine store_results(icall,cell,getdv,getdb,realviscosity,stressmax,xyzh,&
                         gradh,divcurlv,divcurlB,alphaind,dvdx,vxyzu,Bxyz,&
                         dustfrac,rhomax,nneightry,nneighact,maxneightry,&
                         maxneighact,np,ncalc,radprop)
 use part,        only:hrho,get_partinfo,iamgas,&
                       maxphase,massoftype,igas,n_R,n_electronT,&
                       eta_nimhd,iohm,ihall,iambi,ndustlarge,ndustsmall,xyzh_soa,&
                       store_temperature,temperature,maxgradh,idust,&
                       ifluxx,ifluxz,ithick
 use io,          only:fatal,real4
 use eos,         only:get_temperature,get_spsound
 use dim,         only:maxp,ndivcurlv,ndivcurlB,nalpha,mhd_nonideal,use_dust,&
                       do_radiation
 use options,     only:ieos,alpha,alphamax,use_dustfrac
 use viscosity,   only:bulkvisc,shearparam
 use nicil,       only:nicil_get_ion_n,nicil_get_eta,nicil_translate_error
 use linklist,    only:set_hmaxcell,iorder
 use kernel,      only:radkern
 use part,        only:xyzh_soa,store_temperature,temperature
 use kdtree,      only:treecache

 integer,         intent(in)    :: icall
 type(celldens),  intent(in)    :: cell
 logical,         intent(in)    :: getdv
 logical,         intent(in)    :: getdB
 logical,         intent(in)    :: realviscosity
 real,            intent(inout) :: stressmax
 real,            intent(inout) :: xyzh(:,:)
 real(kind=4),    intent(inout) :: gradh(:,:)
 real(kind=4),    intent(inout) :: divcurlv(:,:)
 real(kind=4),    intent(inout) :: divcurlB(:,:)
 real(kind=4),    intent(inout) :: alphaind(:,:)
 real(kind=4),    intent(inout) :: dvdx(:,:)
 real,            intent(in)    :: vxyzu(:,:)
 real,            intent(out)   :: dustfrac(:,:)
 real,            intent(out)   :: Bxyz(:,:)
 real,            intent(inout) :: rhomax
 integer(kind=8), intent(inout) :: nneightry
 integer(kind=8), intent(inout) :: nneighact
 integer(kind=8), intent(inout) :: maxneightry
 integer,         intent(inout) :: maxneighact
 integer,         intent(inout) :: np
 integer(kind=8), intent(inout) :: ncalc
 real,            intent(inout) :: radprop(:,:)

 real         :: rhosum(maxrhosum)

 integer      :: iamtypei,i,lli,ierr,l
 logical      :: iactivei,iamgasi,iamdusti
 logical      :: igotrmatrix,igotspsound
 real         :: hi,hi1,hi21,hi31,hi41
 real         :: pmassi,rhoi
 real(kind=8) :: gradhi,gradsofti
 real         :: psii
 real         :: Bxi,Byi,Bzi,gradBi
 real         :: vxyzui(4)
 real         :: spsoundi,xi_limiter
 real         :: divcurlvi(5),rmatrix(6),dvdxi(9)
 real         :: divcurlBi(ndivcurlB)
 real         :: temperaturei,Bi
 real         :: rho1i,term,denom,rhodusti(maxdustlarge)

 do i = 1,cell%npcell
    lli = iorder(cell%arr_index(i))
    hi = cell%h(i)
    rhosum = cell%rhosums(:,i)

    if (hi < 0.) call fatal('densityiterate','hi < 0 after iterations',lli,var='h',val=hi)

    hi1   = 1./hi
    hi21  = hi1*hi1
    hi31  = hi1*hi21
    hi41  = hi21*hi21

    if (maxphase==maxp) then
       call get_partinfo(cell%iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
    else
       iactivei = .true.
       iamtypei = igas
       iamgasi  = .true.
    endif

    pmassi = massoftype(iamtypei)

    call finish_rhosum(rhosum,pmassi,hi,.false.,rhoi=rhoi,gradhi=gradhi,gradsofti=gradsofti)

    !print*, "rho i", rhoi

    !if (rhosum(1,i) > 0) then 
    ! print*, "rhosum", rhosum
    !endif 
    !
    !--store final results of density iteration
    !
    xyzh(4,lli) = hrho(rhoi,pmassi)
    treecache(4,cell%arr_index(i)) = xyzh(4,lli)
    !xyzh_soa(4,cell%arr_index(i)) = xyzh(4,lli)

    if (xyzh(4,lli) < 0.) call fatal('densityiterate','setting negative h from hrho',i,var='rhoi',val=real(rhoi))

    if (maxgradh==maxp) then
       gradh(1,lli) = real(gradhi,kind=kind(gradh))
#ifdef GRAVITY
       gradh(2,lli) = real(gradsofti,kind=kind(gradh))
#endif
    endif

    rho1i  = 1./rhoi
    rhomax = max(rhomax,real(rhoi))
    if (use_dust .and. .not. use_dustfrac) then
       !
       ! for 2-fluid dust compute dust density on gas particles
       ! and store it in dustfrac as dust-to-gas ratio
       ! so that rho times dustfrac gives dust density
       !
       dustfrac(:,lli) = 0.
       if (iamgasi) then
          do l=1,ndustlarge
             rhodusti(l) = cnormk*massoftype(idust+l-1)*(rhosum(irhodusti+l-1))*hi31
             dustfrac(ndustsmall+l,lli) = rhodusti(l)*rho1i ! dust-to-gas ratio
          enddo
       endif
    endif
    !
    ! store divv and curl v and related quantities
    !
    igotrmatrix = .false.
    igotspsound = .false.

    term = cnormk*pmassi*gradhi*rho1i*hi41
    if (getdv) then
       call calculate_rmatrix_from_sums(rhosum,denom,rmatrix,igotrmatrix)
       call calculate_divcurlv_from_sums(rhosum,term,divcurlvi,xi_limiter,ndivcurlv,denom,rmatrix)
       divcurlv(1:ndivcurlv,lli) = real(divcurlvi(1:ndivcurlv),kind=kind(divcurlv)) ! save to global memory
       !
       ! Cullen & Dehnen (2010) viscosity switch, set alphaloc
       !
       if (nalpha >= 2 .and. iamgasi) then
          igotspsound = .true.
          vxyzui(1) = cell%xpartvec(ivxi,i)
          vxyzui(2) = cell%xpartvec(ivyi,i)
          vxyzui(3) = cell%xpartvec(ivzi,i)
          vxyzui(4) = cell%xpartvec(ieni,i)

          if (store_temperature) then
             spsoundi = get_spsound(ieos,xyzh(:,lli),real(rhoi),vxyzui(:),temperature(lli))
          else
             spsoundi = get_spsound(ieos,xyzh(:,lli),real(rhoi),vxyzui(:))
          endif
          alphaind(2,lli) = real4(get_alphaloc(divcurlvi(5),spsoundi,hi,xi_limiter,alpha,alphamax))
       endif
    else ! we always need div v for h prediction
       if (ndivcurlv >= 1) divcurlv(1,lli) = -real4(rhosum(idivvi)*term)
       if (nalpha >= 2) alphaind(2,lli) = 0.
    endif
    !
    ! store div B, curl B and related quantities
    !
    if (mhd .and. iamgasi) then
       ! construct B from B/rho (conservative to primitive)
       Bxi = cell%xpartvec(iBevolxi,i) * rhoi
       Byi = cell%xpartvec(iBevolyi,i) * rhoi
       Bzi = cell%xpartvec(iBevolzi,i) * rhoi
       psii = cell%xpartvec(ipsi,i)

       ! store primitive variables (if icall < 2)
       if (icall==0 .or. icall==1) then
          Bxyz(1,lli) = Bxi
          Bxyz(2,lli) = Byi
          Bxyz(3,lli) = Bzi
       endif

       if (getdB) then
          term = cnormk*pmassi*gradhi*rho1i*hi41
          call calculate_divcurlB_from_sums(rhosum,term,divcurlBi,gradBi,ndivcurlB)
          divcurlB(:,lli) = real(divcurlBi(:),kind=kind(divcurlB))
       else
          divcurlBi(:) = 0.
       endif
       !
       !--calculate Z_grain, n_electron and non-ideal MHD coefficients
       !
       if (mhd_nonideal) then
          if (.not. igotspsound) then
             vxyzui(1) = cell%rhosums(ivxi,i)
             vxyzui(2) = cell%rhosums(ivyi,i)
             vxyzui(3) = cell%rhosums(ivzi,i)
             vxyzui(4) = cell%rhosums(ieni,i)
          endif
          temperaturei = get_temperature(ieos,cell%xpartvec(ixi:izi,i),real(rhoi),vxyzui(:))
          Bi           = sqrt(Bxi*Bxi + Byi*Byi + Bzi*Bzi)
          call nicil_get_ion_n(real(rhoi),temperaturei,n_R(:,lli),n_electronT(lli),ierr)
          if (ierr/=0) then
             call nicil_translate_error(ierr)
             if (ierr > 0) call fatal('densityiterate','error in Nicil in calculating number densities')
          endif
          call nicil_get_eta(eta_nimhd(iohm,lli),eta_nimhd(ihall,lli),eta_nimhd(iambi,lli),Bi &
                            ,real(rhoi),temperaturei,n_R(:,lli),n_electronT(lli),ierr)
          if (ierr/=0) then ! ierr is reset in the above subroutine
             call nicil_translate_error(ierr)
             if (ierr > 0) call fatal('densityiterate','error in Nicil in calculating eta')
          endif
       endif
    endif
    !
    !--get strain tensor from summations
    !
    if (maxdvdx==maxp .and. getdv) then
       if (.not.igotrmatrix) call calculate_rmatrix_from_sums(cell%rhosums(:,i),denom,rmatrix,igotrmatrix)
       call calculate_strain_from_sums(cell%rhosums(:,i),term,denom,rmatrix,dvdxi)
       ! check for negative stresses to prevent tensile instability
       if (realviscosity) call get_max_stress(dvdxi,divcurlvi(1),rho1i,stressmax,shearparam,bulkvisc)
       ! store strain tensor
       dvdx(:,lli) = real(dvdxi(:),kind=kind(dvdx))
    endif

    if (do_radiation.and.iamgasi) radprop(ifluxx:ifluxz,lli) = cell%rhosums(iradfxi:iradfzi,i)*term

    ! stats
    nneightry = nneightry + cell%nneightry
    nneighact = nneighact + cell%nneigh(i)
    maxneightry = max(int(maxneightry),cell%nneightry)
    maxneighact = max(maxneighact,cell%nneigh(i))
 enddo
 np = np + cell%npcell
 ncalc = ncalc + cell%nits

end subroutine store_results

logical function well_separated(node1,node2) result(bool)
  use dtypekdtree, only:kdnode
  type(kdnode), intent(in) :: node1, node2
  real :: cm1(3), cm2(3), dx(3), zmag
  real :: theta,rmax1,rmax2
  !logical :: bool

  !bool =.true.
  !return 
  theta = 0.1

  cm1 = node1 % xcen
  cm2 = node2 % xcen 

  ! Get the magnitude of the difference between CoM's
  dx = cm1 - cm2
  !print*,dx
  zmag = sqrt(dx(1)*dx(1) + dx(2)*dx(2)+ dx(3)*dx(3))

  ! get the rmax for nodes 1 and 2
  rmax1 = node1 % size
  rmax2 = node2 % size
  !call distance_to_corner(node1,rmax1)
  !call distance_to_corner(node2,rmax2)


  !print*, "Rmax values: "
  !print*, rmax1, rmax2

  !print*,"Zmag: ", zmag
  !print*, "rmax1 + rmax2/theta: ", (rmax1 + rmax2)/theta
  if (zmag > (rmax1+rmax2)/theta) then
    bool = .TRUE.
  else 
    bool = .FALSE.
  endif  
 end function well_separated

 logical function trial_neighbour(node1,node2)
#ifdef PERIODIC
  use boundary, only:dxbound,dybound,dzbound
#endif 
  use dtypekdtree, only:kdnode
  use kernel, only:radkern
  type(kdnode), intent(in) :: node1, node2
  real :: dx, dy, dz
  real :: xsizej,xsizei,hmaxj,hmaxi,rcuti,rcutj,rcut,rcut2,r2 
  real :: xoffset,yoffset,zoffset
  real :: dr,r(3) 
  

  ! PERIODIC STUFF 
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz
  !print*, "periodic"
  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
  !print*, "periodic"
#endif

  dx     = node1%xcen(1) - node2%xcen(1)      ! distance between node centres
  dy     = node1%xcen(2) - node2%xcen(2)
  dz     = node1%xcen(3) - node2%xcen(3)

  !call get_child_nodes(n,il,ir)
  xoffset = 0.
  yoffset = 0.
  zoffset = 0.
#ifdef PERIODIC
    if (abs(dx) > hdlx) then            ! mod distances across boundary if periodic BCs
       xoffset = dxbound*SIGN(1.0,dx)
       dx = dx - xoffset
    endif
    if (abs(dy) > hdly) then
       yoffset = dybound*SIGN(1.0,dy)
       dy = dy - yoffset
    endif
    if (abs(dz) > hdlz) then
       zoffset = dzbound*SIGN(1.0,dz)
       dz = dz - zoffset
    endif
#endif
   r = (/dx,dy,dz/)
   ! find vector mag 
  dr = norm2(r)
  xsizej = node2%size
  hmaxj  = node2%hmax
  xsizei = node1%size
  hmaxi  = node1%hmax
  rcuti  = radkern*hmaxi
  !rcuti = node1%rcuti
  rcutj  = radkern*hmaxj
  rcut   = max(rcuti,rcutj)
  !rcut = rcuti

  rcut2  = (xsizei + xsizej + rcut)**2
  !r2     = dx*dx + dy*dy + dz*dz 
  ! add node size 
  dr = dr ! - xsizei - xsizej
  r2 = dr*dr

  !print*, (xsizei + xsizej + rcut)
  !print*, dx,dy,dz 
  trial_neighbour = r2 < rcut2



 end function trial_neighbour

real function get_r2(node1, node2)
#ifdef PERIODIC
  use boundary, only:dxbound,dybound,dzbound
#endif 
  use dtypekdtree, only:kdnode
  use kernel, only:radkern
  type(kdnode), intent(in) :: node1, node2
  real :: dx, dy, dz
  real :: xsizej,xsizei,hmaxj,hmaxi,rcuti,rcutj,rcut,rcut2,r2 
  real :: xoffset,yoffset,zoffset
  real :: dr, r(3)


  ! PERIODIC STUFF 
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz

  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
  !print*, "periodic"
#endif

  dx     = node1%xcen(1) - node2%xcen(1)      ! distance between node centres
  dy     = node1%xcen(2) - node2%xcen(2)
  dz     = node1%xcen(3) - node2%xcen(3)

  !call get_child_nodes(n,il,ir)
  xoffset = 0.
  yoffset = 0.
  zoffset = 0.
#ifdef PERIODIC
    if (abs(dx) > hdlx) then            ! mod distances across boundary if periodic BCs
       xoffset = dxbound*SIGN(1.0,dx)
       dx = dx - xoffset
    endif
    if (abs(dy) > hdly) then
       yoffset = dybound*SIGN(1.0,dy)
       dy = dy - yoffset
    endif
    if (abs(dz) > hdlz) then
       zoffset = dzbound*SIGN(1.0,dz)
       dz = dz - zoffset
    endif
#endif
  
   r = (/dx,dy,dz/)
   ! find vector mag 
  dr = norm2(r)
  ! add node size 
  
  xsizej = node2%size
  hmaxj  = node2%hmax
  xsizei = node1%size
  hmaxi  = node1%hmax
  rcuti  = radkern*hmaxi
  !rcuti = node1%rcuti
  rcutj  = radkern*hmaxj
  !rcut   = max(rcuti,rcutj)
  rcut = rcuti
  rcut2  = (xsizei + xsizej+ rcut)**2
  dr = dr  - xsizei - xsizej
  !r2     = dx*dx + dy*dy + dz*dz
  r2 = dr*dr 

  get_r2 = r2 

 end function get_r2

real function get_rcut2(node1,node2)
#ifdef PERIODIC
  use boundary, only:dxbound,dybound,dzbound
#endif 
  use dtypekdtree, only:kdnode
  use kernel, only:radkern
  type(kdnode), intent(in) :: node1, node2
  real :: dx, dy, dz
  real :: xsizej,xsizei,hmaxj,hmaxi,rcuti,rcutj,rcut,rcut2,r2 
  real :: xoffset,yoffset,zoffset


  ! PERIODIC STUFF 
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz

  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
  !print*, "periodic"
#endif

  dx     = node1%xcen(1) - node2%xcen(1)      ! distance between node centres
  dy     = node1%xcen(2) - node2%xcen(2)
  dz     = node1%xcen(3) - node2%xcen(3)

  !call get_child_nodes(n,il,ir)
  xoffset = 0.
  yoffset = 0.
  zoffset = 0.
#ifdef PERIODIC
    if (abs(dx) > hdlx) then
        print*, "periodic"            ! mod distances across boundary if periodic BCs
       xoffset = dxbound*SIGN(1.0,dx)
       dx = dx - xoffset
    endif
    if (abs(dy) > hdly) then
      print*, "periodic"
       yoffset = dybound*SIGN(1.0,dy)
       dy = dy - yoffset
    endif
    if (abs(dz) > hdlz) then
      print*, "periodic"
       zoffset = dzbound*SIGN(1.0,dz)
       dz = dz - zoffset
    endif
#endif
  xsizej = node2%size
  hmaxj  = node2%hmax
  xsizei = node1%size
  hmaxi  = node1%hmax
  rcuti  = radkern*hmaxi
  !rcuti = node1%rcuti
  rcutj  = radkern*hmaxj
  !rcut   = max(rcuti,rcutj)
  rcut = rcuti
  rcut2  = (xsizei + xsizej + rcut)**2
  !print*,"rcut: ", (xsizei + xsizej + rcut)
  !print*, "r: ", dx,dy,dz
  r2     = dx*dx + dy*dy + dz*dz

  get_rcut2 = rcut2

 end function get_rcut2



 subroutine global_to_local(stack,stacklocal,threadnumber,istacklocal,top)
  type(denstreestack), intent(inout) :: stack(:)
  type(denstreestacklocal), intent(inout) :: stacklocal(:)
  integer, intent(in) :: threadnumber
  integer, intent(inout) :: istacklocal(:),top 
  integer :: theinteractions(2), thenode,i 

  ! POP FROM GLOBAL AND PUSH TO LOCAL STACK 

  thenode = 0 
  theinteractions = 0

  ! POP FROM GLOBAL
  !!$omp critical (stack)
  thenode = stack(top) % nodeindex1
  theinteractions =  stack(top) % interactions
  ! cleanup interactions
  stack(top) % interactions = 0
  top = top - 1 
  !!$omp end critical (stack)

  do i=1,2
      if (theinteractions(i) /= 0) then 
      !print*,"i is: ", i
            istacklocal(threadnumber) = istacklocal(threadnumber) + 1
            stacklocal(istacklocal(threadnumber)) % nodeindex1 = thenode
            stacklocal(istacklocal(threadnumber)) % nodeindex2 = theinteractions(i)
      endif 



  enddo 



 end subroutine global_to_local


end module densityforce
