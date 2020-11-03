!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module kdtree
!
! This module implements the k-d tree build
!    and associated tree walking routines
!
! :References:
!    Gafton & Rosswog (2011), MNRAS 418, 770-781
!    Benz, Bowers, Cameron & Press (1990), ApJ 348, 647-667
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, balance, boundary, dim, domain, dtypekdtree,
!   fastmath, io, kernel, mpiderivs, mpiutils, part
!
 use dim,         only:maxp,ncellsmax,minpart
 use io,          only:nprocs
 use dtypekdtree, only:kdnode,ndimtree
 use part,        only:ll,iphase,xyzh_soa,iphase_soa,maxphase,dxi

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)
#ifdef MPI
 type(kdnode),     allocatable :: refinementnode(:)
#endif
 integer,          allocatable :: list(:)
 !$omp threadprivate(list)

!
!--tree parameters
!
 integer, public :: irootnode
 character(len=1), parameter, public :: labelax(3) = (/'x','y','z'/)
 integer, parameter :: maxlevelcrazy = 31
!
!--runtime options for this module
!
 real, public :: tree_accuracy = 0.5
 logical, private :: done_init_kdtree = .false.
 logical, private :: already_warned = .false.
 integer, private :: numthreads

 public :: allocate_kdtree, deallocate_kdtree
 public :: maketree, revtree, getneigh, kdnode,get_nodesize_max
#ifdef MPI
 public :: maketreeglobal
#endif
 public :: empty_tree
 public :: compute_fnode, expand_fgrav_in_taylor_series,evaluate,interact,interact_standalone,translate_fgrav_in_taylor_series

 integer, parameter, public :: lenfgrav = 20

 integer, public :: maxlevel_indexed, maxlevel

 type kdbuildstack
    integer :: node
    integer :: parent
    integer :: level
    integer :: npnode
    real    :: xmin(ndimtree)
    real    :: xmax(ndimtree)
 end type


 type denstreestacklocal
  integer :: nodeindex1 =0, nodeindex2 = 0 
 end type 

 type denstreestack
  integer :: nodeindex1 = 0 
  integer :: interactions(2) = 0
 end type 

 type evaluate_stack_data

 integer :: nodeindex = 0
 real :: fnode(lenfgrav) = 0
 real :: z0(3) = 0

 end type 

 private

contains

subroutine allocate_kdtree
 use allocutils, only:allocate_array

 call allocate_array('inoderange', inoderange, 2, ncellsmax+1)
 call allocate_array('inodeparts', inodeparts, maxp)
#ifdef MPI
 call allocate_array('refinementnode', refinementnode, ncellsmax+1)
#endif
 !$omp parallel
 call allocate_array('list', list, maxp)
 !$omp end parallel

end subroutine allocate_kdtree

subroutine deallocate_kdtree
 if (allocated(inoderange)) deallocate(inoderange)
 if (allocated(inodeparts)) deallocate(inodeparts)
#ifdef MPI
 if (allocated(refinementnode)) deallocate(refinementnode)
#endif
 !$omp parallel
 if (allocated(list)) deallocate(list)
 !$omp end parallel

end subroutine deallocate_kdtree


!--------------------------------------------------------------------------------
!+
!  Routine to build the tree from scratch
!
!  Notes/To do:
!  -openMP parallelisation of maketree_stack (done - April 2013)
!  -test centre of mass vs. geometric centre based cell sizes (done - 2013)
!  -need analysis module that times (and checks) set_linklist and tree walk
!   for a given dump
!  -test bottom-up vs. top-down neighbour search
!  -should we try to store tree structure with particle arrays?
!  -need to compute centre of mass and moments for each cell on the fly (c.f. revtree?)
!  -need to implement long-range gravitational interaction (done - May 2013)
!  -implement revtree routine to update tree w/out rebuilding (done - Sep 2015)
!+
!-------------------------------------------------------------------------------
subroutine maketree(node, xyzh, np, ndim, ifirstincell, ncells, refinelevels)
 use io,   only:fatal,warning,iprint,iverbose
!$ use omp_lib
 type(kdnode),    intent(out)   :: node(:) !ncellsmax+1)
 integer,         intent(in)    :: np,ndim
 real,            intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,         intent(out)   :: ifirstincell(:) !ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 integer, optional, intent(out)  :: refinelevels

 integer :: i,npnode,il,ir,istack,nl,nr,mymum
 integer :: nnode,minlevel,level,nqueue
 real :: xmini(ndim),xmaxi(ndim),xminl(ndim),xmaxl(ndim),xminr(ndim),xmaxr(ndim)
 integer, parameter :: istacksize = 512
 type(kdbuildstack), save :: stack(istacksize)
 !$omp threadprivate(stack)
 type(kdbuildstack) :: queue(istacksize)
!$ integer :: threadid
 integer :: npcounter
 logical :: wassplit,finished
 character(len=10) :: string

 irootnode = 1
 ifirstincell = 0

 ir = 0
 il = 0
 nl = 0
 nr = 0
 wassplit = .false.
 finished = .false.

 ! construct root node, i.e. find bounds of all particles
 call construct_root_node(np,npcounter,irootnode,ndim,xmini,xmaxi,ifirstincell,xyzh)
 dxi = xmaxi-xmini

 if (inoderange(1,irootnode)==0 .or. inoderange(2,irootnode)==0 ) then
    call fatal('maketree','no particles or all particles dead/accreted')
 endif

! Put root node on top of stack
 ncells = 1
 maxlevel = 0
 minlevel = 31
 istack = 1

 ! maximum level where 2^k indexing can be used (thus avoiding critical sections)
 ! deeper than this we access cells via a stack as usual
 maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 1
  !maxlevel_indexed = 2
 ! default number of cells is the size of the `indexed' part of the tree
 ! this can be *increased* by building tree beyond indexed levels
 ! and is decreased afterwards according to the maximum depth actually reached
 ncells = 2**(maxlevel_indexed+1) - 1

 ! need to number of particles in node during build to allocate space for local link list
 ! this is counted above to remove dead/accreted particles
 call push_onto_stack(queue(istack),irootnode,0,0,npcounter,xmini,xmaxi,ndim)

 if (.not.done_init_kdtree) then
    ! 1 thread for serial, overwritten when using OpenMP
    numthreads = 1

    ! get number of OpenMPthreads
    !$omp parallel default(none) shared(numthreads)
!$  numthreads = omp_get_num_threads()
    !$omp end parallel
    done_init_kdtree = .true.
 endif

 nqueue = numthreads
 ! build using a queue to build level by level until number of nodes = number of threads
 over_queue: do while (istack  <  nqueue)
    ! if the tree finished while building the queue, then we should just return
    ! only happens for small particle numbers
    if (istack <= 0) then
       finished = .true.
       exit over_queue
    endif
    ! pop off front of queue
    call pop_off_stack(queue(1), istack, nnode, mymum, level, npnode, xmini, xmaxi, ndim)

    ! shuffle queue forward
    do i=1,istack
       queue(i) = queue(i+1)
    enddo

    ! construct node
    call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
            il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, &
            ncells, ifirstincell, minlevel, maxlevel, ndim, xyzh, wassplit, list)

    if (wassplit) then ! add children to back of queue
       if (istack+2 > istacksize) call fatal('maketree',&
                                       'queue size exceeded in tree build, increase istacksize and recompile')

       istack = istack + 1
       call push_onto_stack(queue(istack),il,nnode,level+1,nl,xminl,xmaxl,ndim)
       istack = istack + 1
       call push_onto_stack(queue(istack),ir,nnode,level+1,nr,xminr,xmaxr,ndim)
    endif

 enddo over_queue

 ! fix the indices

 done: if (.not.finished) then

    ! build using a stack which builds depth first
    ! each thread grabs a node from the queue and builds its own subtree

    !$omp parallel default(none) &
    !$omp shared(queue) &
    !$omp shared(ll, ifirstincell) &
    !$omp shared(xyzh) &
    !$omp shared(np, ndim) &
    !$omp shared(node, ncells) &
    !$omp shared(nqueue) &
    !$omp private(istack) &
    !$omp private(nnode, mymum, level, npnode, xmini, xmaxi) &
    !$omp private(ir, il, nl, nr) &
    !$omp private(xminr, xmaxr, xminl, xmaxl) &
    !$omp private(threadid) &
    !$omp private(wassplit) &
    !$omp reduction(min:minlevel) &
    !$omp reduction(max:maxlevel)
    !$omp do schedule(static)
    do i = 1, nqueue

       stack(1) = queue(i)
       istack = 1

       over_stack: do while(istack > 0)

          ! pop node off top of stack
          call pop_off_stack(stack(istack), istack, nnode, mymum, level, npnode, xmini, xmaxi, ndim)

          ! construct node
          call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
              il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, &
              ncells, ifirstincell, minlevel, maxlevel, ndim, xyzh, wassplit, list)

          if (wassplit) then ! add children to top of stack
             if (istack+2 > istacksize) call fatal('maketree',&
                                       'stack size exceeded in tree build, increase istacksize and recompile')

             istack = istack + 1
             call push_onto_stack(stack(istack),il,nnode,level+1,nl,xminl,xmaxl,ndim)
             istack = istack + 1
             call push_onto_stack(stack(istack),ir,nnode,level+1,nr,xminr,xmaxr,ndim)
          endif

       enddo over_stack
    enddo
!$omp enddo
!$omp end parallel

 endif done

 ! decrease number of cells if tree is entirely within 2^k indexing limit
 if (maxlevel < maxlevel_indexed) then
    ncells = 2**(maxlevel+1) - 1
 endif
 if (maxlevel > maxlevel_indexed .and. .not.already_warned) then
    write(string,"(i10)") 2**(maxlevel-maxlevel_indexed)
    if (iverbose > 0) call warning('maketree','maxlevel > max_indexed: will run faster if recompiled with '// &
               'NCELLSMAX='//trim(adjustl(string))//'*maxp,')
 endif

 if (present(refinelevels)) refinelevels = minlevel

 if (iverbose >= 3) then
    write(iprint,"(a,i10,3(a,i2))") ' maketree: nodes = ',ncells,', max level = ',maxlevel,&
       ', min leaf level = ',minlevel,' max level indexed = ',maxlevel_indexed
 endif

end subroutine maketree

!----------------------------
!+
! routine to empty the tree
!+
!----------------------------
subroutine empty_tree(node)
 type(kdnode), intent(out) :: node(:)
 integer :: i

!$omp parallel do private(i)
 do i=1,size(node)
    node(i)%xcen = 0.
    node(i)%size = 0.
    node(i)%hmax = 0.
    node(i)%leftchild = 0
    node(i)%rightchild = 0
    node(i)%parent = 0
#ifdef GRAVITY
    node(i)%mass  = 0.
    node(i)%quads = 0.
#endif
 enddo
!$omp end parallel do

end subroutine empty_tree

!---------------------------------
!+
! routine to construct root node
!+
!---------------------------------
subroutine construct_root_node(np,nproot,irootnode,ndim,xmini,xmaxi,ifirstincell,xyzh)
#ifdef PERIODIC
 use boundary, only:cross_boundary
 use domain,   only:isperiodic
#endif
#ifdef IND_TIMESTEPS
 use part, only:iphase,iactive
#endif
 use part, only:isdead_or_accreted
 integer,         intent(in)  :: np,irootnode,ndim
 integer,         intent(out) :: nproot
 real,            intent(out) :: xmini(ndim), xmaxi(ndim)
 integer,         intent(inout) :: ifirstincell(:)
 real,            intent(inout) :: xyzh(:,:)
 integer :: i,ncross
 real    :: xminpart,yminpart,zminpart,xmaxpart,ymaxpart,zmaxpart
 real    :: xi, yi, zi

 xminpart = xyzh(1,1)
 yminpart = xyzh(2,1)
 zminpart = xyzh(3,1)
 xmaxpart = xminpart
 ymaxpart = yminpart
 zmaxpart = zminpart

 ncross = 0
 nproot = 0
 !$omp parallel default(none) &
 !$omp shared(np,xyzh) &
 !$omp shared(inodeparts,iphase,xyzh_soa,iphase_soa,nproot) &
#ifdef PERIODIC
 !$omp shared(isperiodic) &
 !$omp reduction(+:ncross) &
#endif
 !$omp private(i,xi,yi,zi) &
 !$omp reduction(min:xminpart,yminpart,zminpart) &
 !$omp reduction(max:xmaxpart,ymaxpart,zmaxpart)
 !$omp do schedule(guided,1)
 do i=1,np
    if (.not.isdead_or_accreted(xyzh(4,i))) then
#ifdef PERIODIC
       call cross_boundary(isperiodic,xyzh(:,i),ncross)
#endif
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       xminpart = min(xminpart,xi)
       yminpart = min(yminpart,yi)
       zminpart = min(zminpart,zi)
       xmaxpart = max(xmaxpart,xi)
       ymaxpart = max(ymaxpart,yi)
       zmaxpart = max(zmaxpart,zi)
    endif
 enddo
 !$omp enddo
 !$omp end parallel

 do i=1,np
    isnotdead: if (.not.isdead_or_accreted(xyzh(4,i))) then
       nproot = nproot + 1

#ifdef IND_TIMESTEPS
       if (iactive(iphase(i))) then
          inodeparts(nproot) = i  ! +ve if active
       else
          inodeparts(nproot) = -i ! -ve if inactive
       endif
#else
       inodeparts(nproot) = i
#endif
       xyzh_soa(nproot,:) = xyzh(:,i)
       iphase_soa(nproot) = iphase(i)
    endif isnotdead
 enddo

 if (nproot /= 0) then
    inoderange(1,irootnode) = 1
    inoderange(2,irootnode) = nproot
 else
    inoderange(:,irootnode) = 0
 endif

 if (ndim==2) then
    xmini(:) = (/xminpart,yminpart/)
    xmaxi(:) = (/xmaxpart,ymaxpart/)
 else
    xmini(:) = (/xminpart,yminpart,zminpart/)
    xmaxi(:) = (/xmaxpart,ymaxpart,zmaxpart/)
 endif

end subroutine construct_root_node

! also used for queue push
pure subroutine push_onto_stack(stackentry,node,parent,level,npnode,xmin,xmax,ndim)
 type(kdbuildstack), intent(out) :: stackentry
 integer,            intent(in)  :: node,parent,level,ndim
 integer,            intent(in)  :: npnode
 real,               intent(in)  :: xmin(ndim),xmax(ndim)

 stackentry%node   = node
 stackentry%parent = parent
 stackentry%level  = level
 stackentry%npnode = npnode
 stackentry%xmin   = xmin
 stackentry%xmax   = xmax

end subroutine push_onto_stack

! also used for queue pop
pure subroutine pop_off_stack(stackentry, istack, nnode, mymum, level, npnode, xmini, xmaxi, ndim)
 type(kdbuildstack), intent(in)    :: stackentry
 integer,            intent(inout) :: istack
 integer,            intent(out)   :: nnode, mymum, level, npnode
 integer,            intent(in)    :: ndim
 real,               intent(out)   :: xmini(ndim), xmaxi(ndim)

 nnode  = stackentry%node
 mymum  = stackentry%parent
 level  = stackentry%level
 npnode = stackentry%npnode
 xmini  = stackentry%xmin
 xmaxi  = stackentry%xmax
 istack = istack - 1

end subroutine pop_off_stack

!--------------------------------------------------------------------
!+
!  create all the properties for a given node such as centre of mass,
!  size, max smoothing length, etc
!  will also split the node if necessary, setting wassplit=true
!  returns the left and right child information if split
!+
!--------------------------------------------------------------------
subroutine construct_node(nodeentry, nnode, mymum, level, xmini, xmaxi, npnode, doparallel,&
            il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, &
            ncells, ifirstincell, minlevel, maxlevel, ndim, xyzh, wassplit, list, &
            groupsize)
 use dim,       only:maxtypes
 use part,      only:massoftype,igas,iamtype,maxphase,maxp,npartoftype
 use io,        only:fatal,error
#ifdef MPI
 use mpiderivs, only:get_group_cofm,reduce_group
#endif
 type(kdnode),      intent(out)   :: nodeentry
 integer,           intent(in)    :: nnode, mymum, level
 integer,           intent(in)    :: ndim
 real,              intent(inout) :: xmini(ndim), xmaxi(ndim)
 integer,           intent(in)    :: npnode
 logical,           intent(in)    :: doparallel
 integer, optional, intent(in)    :: groupsize ! used for global node construction
 integer,           intent(out)   :: il, ir, nl, nr
 real,              intent(out)   :: xminl(ndim), xmaxl(ndim), xminr(ndim), xmaxr(ndim)
 integer(kind=8),   intent(inout) :: ncells
 integer,           intent(out)   :: ifirstincell(:)
 integer,           intent(inout) :: maxlevel, minlevel
 real,              intent(in)    :: xyzh(:,:)
 logical,           intent(out)   :: wassplit
 integer,           intent(out)   :: list(:) ! not actually sent out, but to avoid repeated memory allocation/deallocation

 real                           :: xyzcofm(ndim)
 real                           :: totmass_node
#ifdef MPI
 real    :: xyzcofmg(ndim)
 real    :: totmassg
#endif
 integer :: npnodetot

 logical :: nodeisactive
 integer :: i,npcounter,i1
 real    :: xi,yi,zi,hi,dx,dy,dz,dr2
 real    :: r2max, hmax
 real    :: xcofm,ycofm,zcofm,fac,dfac
 real    :: x0(ndimtree)
 integer :: iaxis
 real    :: xpivot
#ifdef GRAVITY
 real    :: quads(6)
#endif
 real    :: pmassi

 nodeisactive = .false.
 if (inoderange(1,nnode) > 0) then
    checkactive: do i = inoderange(1,nnode),inoderange(2,nnode)
       if (inodeparts(i) > 0) then
          nodeisactive = .true.
          exit checkactive
       endif
    enddo checkactive
    npcounter = inoderange(2,nnode) - inoderange(1,nnode) + 1
 else
    npcounter = 0
 endif

 if (npcounter /= npnode) then
    print*,'constructing node ',nnode,': found ',npcounter,' particles, expected:',npnode,' particles for this node'
    call fatal('maketree', 'expected number of particles in node differed from actual number')
 endif

 ! following lines to avoid compiler warnings on intent(out) variables
 ir = 0
 il = 0
 nl = 0
 nr = 0
 if ((.not. present(groupsize)) .and. (npnode  <  1)) return ! node has no particles, just quit

 r2max = 0.
 hmax  = 0.
 xyzcofm(:) = 0.
 xcofm = 0.
 ycofm = 0.
 zcofm = 0.
!
! to avoid round off error from repeated multiplication by pmassi (which is small)
! we compute the centre of mass with a factor relative to gas particles
! but only if gas particles are present
!
 pmassi = massoftype(igas)
 fac    = 1.
 totmass_node = 0.
 if (pmassi > 0.) then
    dfac = 1./pmassi
 else
    pmassi = massoftype(maxloc(npartoftype(2:maxtypes),1)+1)
    if (pmassi > 0.) then
       dfac = 1./pmassi
    else
       dfac = 1.
    endif
 endif

 i1=inoderange(1,nnode)
 ! during initial queue build which is serial, we can parallelise this loop
 if (npnode > 1000 .and. doparallel) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(maxp,maxphase) &
    !$omp shared(npnode,massoftype,dfac) &
    !$omp shared(xyzh_soa,i1,iphase_soa) &
    !$omp private(i,xi,yi,zi,hi) &
    !$omp firstprivate(pmassi,fac) &
    !$omp reduction(+:xcofm,ycofm,zcofm,totmass_node) &
    !$omp reduction(max:hmax)
    do i=i1,i1+npnode-1
       xi = xyzh_soa(i,1)
       yi = xyzh_soa(i,2)
       zi = xyzh_soa(i,3)
       hi = xyzh_soa(i,4)
       hmax  = max(hmax,hi)
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase_soa(i)))
          fac    = pmassi*dfac ! to avoid round-off error
       endif
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = xyzh_soa(i,1)
       yi = xyzh_soa(i,2)
       zi = xyzh_soa(i,3)
       hi = xyzh_soa(i,4)
       hmax  = max(hmax,hi)  
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase_soa(i)))
          fac    = pmassi*dfac ! to avoid round-off error
       endif
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
 endif
 if (ndim==2) then
    xyzcofm(1:2) = (/xcofm,ycofm/)
 else
    xyzcofm = (/xcofm,ycofm,zcofm/)
 endif

 ! if we have no particles, then cofm is zero anyway
 if (totmass_node > 0.) then
    xyzcofm(:)   = xyzcofm(:)/(totmass_node*dfac)
 endif

#ifdef MPI
 ! if this is global node construction
 if (present(groupsize)) then
    call get_group_cofm(xyzcofm,totmass_node,level,xyzcofmg,totmassg)
    xyzcofm = xyzcofmg
    totmass_node = totmassg
 endif
#endif

 ! checks the reduced mass in the case of global maketree
 if (totmass_node<=0.) call fatal('mtree','totmass_node==0',val=totmass_node)

!--for gravity, we need the centre of the node to be the centre of mass
 x0(:) = xyzcofm(:)
 r2max = 0.
#ifdef GRAVITY
 quads(:) = 0.
#endif

 !--compute size of node
 !!$omp parallel do if (npnode > 1000 .and. doparallel) &
 !!$omp default(none) schedule(static) &
 !!$omp shared(npnode,xyzh_soa,x0,i1) &
 !!$omp private(i,xi,yi,zi,dx,dy,dz,dr2,pmassi) &
#ifdef GRAVITY
 !!$omp shared(iphase_soa,massoftype) &
 !!$omp reduction(+:quads) &
#endif
 !!$omp reduction(max:r2max)
 do i=i1,i1+npnode-1
    xi = xyzh_soa(i,1)
    yi = xyzh_soa(i,2)
    zi = xyzh_soa(i,3)
    dx    = xi - x0(1)
    dy    = yi - x0(2)
#ifdef TREEVIZ
    dz    = 0.
#else
    dz    = zi - x0(3)
#endif
    dr2   = dx*dx + dy*dy + dz*dz
    r2max = max(r2max,dr2)
#ifdef GRAVITY
    if (maxphase==maxp) then
       pmassi = massoftype(iamtype(iphase_soa(i)))
    endif
    quads(1) = quads(1) + pmassi*(3.*dx*dx - dr2)  ! Q_xx
    quads(2) = quads(2) + pmassi*(3.*dx*dy)        ! Q_xy = Q_yx
    quads(3) = quads(3) + pmassi*(3.*dx*dz)        ! Q_xz = Q_zx
    quads(4) = quads(4) + pmassi*(3.*dy*dy - dr2)  ! Q_yy
    quads(5) = quads(5) + pmassi*(3.*dy*dz)        ! Q_yz = Q_zy
    quads(6) = quads(6) + pmassi*(3.*dz*dz - dr2)  ! Q_zz
#endif
 enddo
 !!$omp end parallel do

#ifdef MPI
 if (present(groupsize)) then
    r2max    = reduce_group(r2max,'max',level)
    hmax     = reduce_group(hmax,'max',level)
    xmini(1) = reduce_group(xmini(1),'min',level)
    xmini(2) = reduce_group(xmini(2),'min',level)
    xmini(3) = reduce_group(xmini(3),'min',level)
    xmaxi(1) = reduce_group(xmaxi(1),'max',level)
    xmaxi(2) = reduce_group(xmaxi(2),'max',level)
    xmaxi(3) = reduce_group(xmaxi(3),'max',level)
#ifdef GRAVITY
    quads(1) = reduce_group(quads(1),'+',level)
    quads(2) = reduce_group(quads(2),'+',level)
    quads(3) = reduce_group(quads(3),'+',level)
    quads(4) = reduce_group(quads(4),'+',level)
    quads(5) = reduce_group(quads(5),'+',level)
    quads(6) = reduce_group(quads(6),'+',level)
#endif
 endif

 if (present(groupsize)) then
    npnodetot = reduce_group(npnode,'+',level)
 else
#endif
    npnodetot = npnode
#ifdef MPI
 endif
#endif

 ! assign properties to node
 nodeentry%xcen    = x0(:)
 nodeentry%size    = sqrt(r2max) + epsilon(r2max)
 nodeentry%hmax    = hmax
 nodeentry%parent  = mymum
#ifdef GRAVITY
 nodeentry%mass    = totmass_node
 nodeentry%quads   = quads
#endif
#ifdef TREEVIZ
 nodeentry%xmin(:) = xmini(:)
 nodeentry%xmax(:) = xmaxi(:)
#endif

 wassplit = (npnodetot > minpart)

 if (.not. wassplit) then
    nodeentry%leftchild  = 0
    nodeentry%rightchild = 0
    maxlevel = max(level,maxlevel)
    minlevel = min(level,minlevel)
    !--filling link list here is unnecessary (already filled)
    !  except with individual timesteps where we mark leaf node as active/inactive
#ifdef IND_TIMESTEPS
    !
    !--mark leaf node as active (contains some active particles)
    !  or inactive by setting the firstincell to +ve (active) or -ve (inactive)
    !
    if (nodeisactive) then
       ifirstincell(nnode) = 1
    else
       ifirstincell(nnode) = -1
    endif
#else
    ifirstincell(nnode) = 1
#endif
 else ! split this node and add children to stack
    iaxis  = maxloc(xmaxi - xmini,1) ! split along longest axis
    xpivot = xyzcofm(iaxis)          ! split on centre of mass

    ! create two children nodes and point to them from current node
    ! always use G&R indexing for global tree
    if ((level < maxlevel_indexed) .or. (present(groupsize))) then
       il = 2*nnode   ! indexing as per Gafton & Rosswog (2011)
       ir = il + 1
    else
       ! need locks when changing ncells to get node labels correct
       !$omp critical(addncells)
       if (ncells+2 > ncellsmax) call fatal('maketree',&
          'number of nodes exceeds array dimensions, increase ncellsmax and recompile',ival=int(ncellsmax))

       ncells = ncells + 1
       il = int(ncells)
       ncells = ncells + 1
       ir = int(ncells)
       !$omp end critical(addncells)
       !$omp flush(ncells)
    endif
    nodeentry%leftchild  = il
    nodeentry%rightchild = ir

    ifirstincell(nnode) = 0

    if (npnode > 0) then
       !call sort_old(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
       !                         inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts)
       !print*,ir,inodeparts(inoderange(1,il):inoderange(2,il))
       call sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                  inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts)
       !if (il==8) print*,nnode,il,iaxis,nl,nr,xpivot,inodeparts(inoderange(1,il):inoderange(2,il))

       !if (any(xyzh_soa(inoderange(1,il):inoderange(2,il),iaxis) > xpivot)) print*,' ERROR x > xpivot on left'
       !if (any(xyzh_soa(inoderange(1,ir):inoderange(2,ir),iaxis) <= xpivot)) print*,' ERROR x <= xpivot on right'

       if (nr + nl  /=  npnode) then
          call error('maketree','number of left + right != parent while splitting (likely cause: NaNs in position arrays)')
       endif

       ! see if all the particles ended up in one node, if so, arbitrarily build 2 cells
       if ( (.not. present(groupsize)) .and. ((nl==npnode) .or. (nr==npnode)) ) then
          ! no need to move particles because if they all ended up in one node,
          ! then they are still in the original order
          nl = npnode / 2
          inoderange(1,il) = inoderange(1,nnode)
          inoderange(2,il) = inoderange(1,nnode) + nl - 1
          inoderange(1,ir) = inoderange(1,nnode) + nl
          inoderange(2,ir) = inoderange(2,nnode)
          nr = npnode - nl
       endif

       xminl(1) = minval(xyzh_soa(inoderange(1,il):inoderange(2,il),1))
       xminl(2) = minval(xyzh_soa(inoderange(1,il):inoderange(2,il),2))
       xminl(3) = minval(xyzh_soa(inoderange(1,il):inoderange(2,il),3))
       xmaxl(1) = maxval(xyzh_soa(inoderange(1,il):inoderange(2,il),1))
       xmaxl(2) = maxval(xyzh_soa(inoderange(1,il):inoderange(2,il),2))
       xmaxl(3) = maxval(xyzh_soa(inoderange(1,il):inoderange(2,il),3))
       xminr(1) = minval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),1))
       xminr(2) = minval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),2))
       xminr(3) = minval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),3))
       xmaxr(1) = maxval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),1))
       xmaxr(2) = maxval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),2))
       xmaxr(3) = maxval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),3))
    else
       nl = 0
       nr = 0
       xminl = xmini
       xmaxl = xmaxi
       xminr = xmini
       xmaxr = xmaxi
    endif

 endif

end subroutine construct_node

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis
!+
!----------------------------------------------------------------
subroutine sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts)
 integer, intent(in)  :: iaxis,imin,imax
 integer, intent(out) :: min_l,max_l,min_r,max_r,nl,nr
 real, intent(inout)  :: xpivot,xyzh_soa(:,:)
 integer(kind=1), intent(inout) :: iphase_soa(:)
 integer,         intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot
 integer(kind=1) :: iphase_swap
 integer :: inodeparts_swap,i,j
 real :: xyzh_swap(4)

 !print*,'nnode ',imin,imax,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
 j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
 ! k = 0
 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:4)  = xyzh_soa(i,1:4)
          iphase_swap     = iphase_soa(i)

          inodeparts(i)   = inodeparts(j)
          xyzh_soa(i,1:4) = xyzh_soa(j,1:4)
          iphase_soa(i)   = iphase_soa(j)

          inodeparts(j)   = inodeparts_swap
          xyzh_soa(j,1:4) = xyzh_swap(1:4)
          iphase_soa(j)   = iphase_swap

          i = i + 1
          j = j - 1
          i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
          j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
          ! k = k + 1
       endif
    endif
 enddo
 if (.not.i_lt_pivot) i = i - 1
 if (j_lt_pivot)      j = j + 1

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine sort_particles_in_cell

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
subroutine getneigh(node,xpos,xsizei,rcuti,ndim,listneigh,nneigh,xyzh,xyzcache,ixyzcachesize,ifirstincell,&
& get_hj,fnode,remote_export,neighbours_found)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 use io,       only:fatal,id
 use part,     only:gravity
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use kernel,   only:radkern
 type(kdnode), intent(in)           :: node(:) !ncellsmax+1)
 integer, intent(in)                :: ndim,ixyzcachesize
 real,    intent(in)                :: xpos(ndim)
 real,    intent(in)                :: xsizei,rcuti
 integer, intent(out)               :: listneigh(:) !maxneigh)
 integer, intent(out)               :: nneigh
 real,    intent(in)                :: xyzh(:,:)
 real,    intent(out)               :: xyzcache(:,:)
 integer, intent(in)                :: ifirstincell(:)
 logical, intent(in)                :: get_hj
 real,    intent(out),    optional  :: fnode(lenfgrav)
 logical, intent(out),    optional  :: remote_export(:)
 integer, intent(out),    optional ::  neighbours_found
 integer, parameter :: istacksize = 300
 integer :: maxcache
 integer :: nstack(istacksize)
 integer :: ipart,n,istack,il,ir,npnode
 real :: dx,dy,dz,xsizej,rcutj
 real :: rcut,rcut2,r2
 real :: xoffset,yoffset,zoffset,tree_acc2,radius 
 logical :: open_tree_node
 integer :: neighleaf
#ifdef GRAVITY
 real :: quads(6)
 real :: dr,totmass_node
#endif
#ifdef PERIODIC
 real :: hdlx,hdly,hdlz

 hdlx = 0.5*dxbound
 hdly = 0.5*dybound
 hdlz = 0.5*dzbound
#endif
 tree_accuracy = 0.5
 tree_acc2 = tree_accuracy*tree_accuracy
 if (present(fnode)) fnode(:) = 0.
 rcut     = rcuti

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache(1,:))
 else
    maxcache = 0
 endif

 if (present(remote_export)) remote_export = .false.

 neighleaf = 0 
 nneigh = 0
 istack = 1
 nstack(istack) = irootnode
 !open_tree_node = .false.
 !print*,"entered routine"
  ! if (present(neighbours_found)) then 
  ! open(8,file="nodes.txt",position='append')
  ! print*, "neighbours_found: ", neighbours_found
  ! print*, "neighbours found present? ", present(neighbours_found)
!endif 
 over_stack: do while(istack /= 0)
    n = nstack(istack)
    istack = istack - 1
    dx = xpos(1) - node(n)%xcen(1)      ! distance between node centres
    dy = xpos(2) - node(n)%xcen(2)
#ifndef TREEVIZ
    dz = xpos(3) - node(n)%xcen(3)
#endif
    xsizej       = node(n)%size
#ifdef GRAVITY
    totmass_node = node(n)%mass
    quads        = node(n)%quads
#endif
    il      = node(n)%leftchild
    ir      = node(n)%rightchild
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
    r2    = dx*dx + dy*dy + dz*dz !+ xsizei*xsizei
    if (get_hj) then  ! find neighbours within both hi and hj
       rcutj = radkern*node(n)%hmax
       rcut  = max(rcuti,rcutj)
    endif
    !rcut = rcuti
    ! call get_node_radius(node,ifirstincell,n,radius)
    ! rcut = max(radius,rcuti)
    rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
#ifdef GRAVITY
    !open_tree_node = tree_acc2*r2 < xsizej*xsizej    ! tree opening criterion for self-gravity
    open_tree_node = .false.
#endif
    
     !print*,"node is: ", n 
    !print*,"parent: ", node(169) % parent 
    !if (n == 84) then 
     !print*,"r2: ", r2, "rcut2: ",  rcut2
     !call get_node_radius(node,ifirstincell,n,radius)
     !print*, "search radius vs rcut + xsizej"
     !print*, radius
     !print*, xsizej + rcut
    !endif 
    ! if (n == 133204 ) then 
    !             print*, "n: ",n, "r2: ", r2,  "rcut2: ", rcut2 
    !             print*, "radkern1 ",rcuti,rcutj 
    !             !rcut = max(radkern*node(childleaf1(i)) % hmax, radkern*node(childleaf2(j)) % hmax) 
    !             print*, "rcut2 arg ", xsizei,xsizej,rcut
    ! endif 
    if ((r2 < rcut2)) then
       if_leaf: if (ifirstincell(n) /= 0) then ! once we hit a leaf node, retrieve contents into trial neighbour cache
              !neighleaf = neighleaf + 1
              !print*, n
              !if (present(neighbours_found)) write(8,*) n
              ! if (n == 65766 ) then 
              !   print*, "r2: ", r2
              !   print*, "rcut2: ", rcut2 
              ! endif 
          if_global_walk: if (present(remote_export)) then
             ! id is stored in ipart as id + 1
             if (ifirstincell(n) /= (id + 1)) then
                remote_export(ifirstincell(n)) = .true.
             endif
          else
             npnode = inoderange(2,n) - inoderange(1,n) + 1
             if_cache_fits: if (nneigh + npnode <= ixyzcachesize) then
                if (maxcache >= 4) then
                   do ipart=1,npnode
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                      xyzcache(nneigh+ipart,4) = 1./xyzh_soa(inoderange(1,n)+ipart-1,4)
                   enddo
                else
                   do ipart=1,npnode
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                   enddo
                endif
             elseif (nneigh < ixyzcachesize) then
                if (maxcache >= 4) then
                   do ipart=1,ixyzcachesize-nneigh
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                      xyzcache(nneigh+ipart,4) = 1./xyzh_soa(inoderange(1,n)+ipart-1,4)
                   enddo
                else
                   do ipart=1,ixyzcachesize-nneigh
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                   enddo
                endif
                do ipart=ixyzcachesize-nneigh+1,npnode
                   listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                enddo
             else
                do ipart=1,npnode
                   listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                enddo
             endif if_cache_fits
             nneigh = nneigh + npnode
          endif if_global_walk
       else
          if (istack+2 > istacksize) call fatal('getneigh','stack overflow in getneigh')
          if (il /= 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
          if (ir /= 0) then
             istack = istack + 1
             nstack(istack) = ir
          endif
       endif if_leaf
! #ifdef GRAVITY
!     elseif (present(fnode) .and. ((.not. present(remote_export)) .or. n < 2*nprocs-1)) then
! !
! !--long range force on node due to distant node, along node centres
! !  along with derivatives in order to perform series expansion
! !
! #ifdef FINVSQRT
!        dr = finvsqrt(r2)
! #else
!        dr = 1./sqrt(r2)
! #endif
!        call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
! #endif
     endif
  enddo over_stack
  !print*, "noneigh getneigh: ", neighleaf 
  !if (present(neighbours_found)) neighbours_found = neighleaf
 !close(8)
end subroutine getneigh

subroutine interact_standalone(node, ifirstincell,xyzh,fxyzu_dtt,poten_dtt,nonodesfound)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
  use omp_lib
  use kernel, only:radkern
  type(kdnode), intent(inout)        :: node(:)
  integer, intent(in)                :: ifirstincell(:)

  real,               intent(in)   :: xyzh(:,:)
  real,               intent(out) :: fxyzu_dtt(:,:) ! Needed for direct sum between neighbours that are not sph neighbours 
  real,               intent(out) :: poten_dtt(:)
  integer, optional, intent(out) :: nonodesfound(:)
  integer :: maxcache
  !integer :: nstack(istacksize)
  integer, parameter                   :: maxstacksize = 10000
  type(denstreestack)       :: stack(maxstacksize)
  type(denstreestacklocal)  :: stacklocal(maxstacksize)
  integer                   :: top,nodeindex1, nodeindex2, leftindex, rightindex, regnodeindex, splitnodeindex, numthreads
  integer, allocatable      :: istacklocal(:)
  integer                   :: k 
  logical, allocatable      :: threadworking(:)
  integer :: n,istack,il,ir,npnode,i1,i2,j,jmax,nn,i
  real :: dx,dy,dz,xsizei,xsizeimax,xsizejmax,xsizej,rcutj,rcuti
  real :: rcut,rcut2,rcut2max,r2,hmax
  real :: xoffset,yoffset,zoffset,tree_acc2
  logical :: open_tree_node
!#ifdef GRAVITY
  real :: quads(6)
  real :: dr,totmass_node
  real :: fnode(20)
  real :: c0,c1(3),c2(3,3),c3(3,3,3)
  integer :: neighleaf,childleaf1(10000),childleaf2(10000),childleafcounter1,childleafcounter2,specialnode
!#endif
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz

  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
#endif

top = 1
stack(top) % nodeindex1 = 1
stack(top) % interactions(1) = 1 
numthreads = 1 
!specialnode = 85
!$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
!$omp end parallel

allocate(istacklocal(numthreads))
allocate(threadworking(numthreads))
!print*,"allocation started"
!allocate(stack(maxstacksize))
!allocate(stacklocal(maxstacksize))
!print*,"allocation finished"

threadworking = .false. 
istacklocal = 0
neighleaf = 0
nonodesfound = 0 

!$omp parallel default(none) &
!$omp shared(node,ifirstincell,threadworking,top,istacklocal,stack,tree_acc2,fxyzu_dtt,poten_dtt,xyzh,specialnode,nonodesfound) &
!$omp private(stacklocal,nodeindex1,nodeindex2,leftindex,rightindex,splitnodeindex,regnodeindex,k) &
!$omp private(dx,dy,dz,totmass_node,xsizei,xsizeimax,xsizej,xsizejmax,quads,r2,open_tree_node,dr,fnode) &
!$omp private(rcut2,rcut2max,rcut,rcuti,rcutj) & 
!$omp private(childleaf1,childleaf2,childleafcounter1,childleafcounter2,i) &
!$omp private(c0,c1,c2,c3) & 
!$omp reduction(+:neighleaf)

do while (any(istacklocal > 0) .or. top > 0 .or. any(threadworking))

     !$ k=omp_get_thread_num() + 1
     !print*,"k: ", k
     !print*,tree_acc2
     !stop 

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
          call pop_local(nodeindex1,nodeindex2,stacklocal,istacklocal,k)
          ! nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
          ! nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2
          ! !print*, "node indexes :", nodeindex1, nodeindex2
          ! istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
          !threadworking(k) = .false.
          !nodes(nodeindex1) % nodefree  = .false.
          !nodes(nodeindex2) % nodefree  = .false.
      else 
        threadworking(k) = .false.
      endif
      !$omp end critical (stack)
    endif 

  if (threadworking(k)) then 
      fnode = 0.
      dx = 0. 
      dy = 0. 
      dz = 0.
      dr = 0.
      totmass_node = 0.
      xsizej = 0.
      xsizei = 0.
      quads = 0. 
      r2 = 0. 
      open_tree_node = .false. 
      rcut = 0.
      rcuti = 0.
      rcutj = 0. 
      rcut2 = 0.
      !tree_acc2 = 0.
      !print*, tree_acc2
      call get_node_data(node(nodeindex1),node(nodeindex2),dx,dy,dz,xsizei,xsizej,fnode,totmass_node,quads)
      !quads = 0.
      r2    = dx*dx + dy*dy + dz*dz
      rcutj = radkern*node(nodeindex2) % hmax
      rcuti = radkern*node(nodeindex1) % hmax 
      rcut  = max(rcuti,rcutj)
      xsizeimax = node(nodeindex1) % maxsize 
      xsizejmax = node(nodeindex2) % maxsize
      !rcut = rcuti
      rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
      rcut2max = (xsizeimax + xsizejmax + rcut)**2   ! node size + search radius

      ! If nodes are equal 
      if (nodeindex1 == nodeindex2) then 
        if (ifirstincell(nodeindex1) > 0) then 
           !print*, "leaf self interaction"
           ! Don't need to check softening here 
           !call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
           !$omp critical (nodefound)
           nonodesfound(nodeindex1) = nonodesfound(nodeindex1) + 1 
           !$omp end critical (nodefound)
           call direct_sum_softened(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
           !stop 
         
        ! If nodes are equal but not leaf then split 
        elseif (ifirstincell(nodeindex1) <= 0) then 
                  leftindex = node(nodeindex1) % leftchild
                  rightindex = node(nodeindex1) % rightchild 

                  !$omp critical (stack)

                  if (leftindex /= 0 .and. rightindex /= 0) then 
                        ! top = top + 1
                        ! stack(top) % nodeindex1 = leftindex
                        ! stack(top) % interactions(1) = leftindex
                        ! stack(top) % interactions(2) = rightindex

                        call push_global(leftindex,leftindex,rightindex,stack,top)

                        ! top = top + 1 
                        ! stack(top) % nodeindex1 = rightindex
                        ! stack(top) % interactions(1) = rightindex
                        ! stack(top) % interactions(2) = leftindex

                        call push_global(rightindex,rightindex,leftindex,stack,top)
                  endif 

                  !$omp end critical (stack)

        endif

      ! If well separated but not sph neighbour, calculate long range interactions 
      else if (well_separated(node(nodeindex1), node(nodeindex2)) .and. r2 >= rcut2max  &
       .and. ifirstincell(nodeindex1) <= 0 .and. ifirstincell(nodeindex2) <= 0) then  ! &
        !.and. .not. check_child_overlap(node,nodeindex1,nodeindex2).and. .not. check_child_overlap(node,nodeindex2,nodeindex1)) then
        !print*, "Well separated nodes: ", nodeindex1, nodeindex2 
        !stop 
        !write(13,*) nodeindex1,nodeindex2
        !print*, "r2 and rcut2: ", r2, rcut2
        !if (ifirstincell(nodeindex1) > 0) stop 
        !if (ifirstincell(nodeindex2) > 0) stop 
        if ( .not. well_separated(node(nodeindex2),node(nodeindex1))) then 
            print*, "MAC BROKEN"
            stop
        endif  

         !$omp critical (nodefound)
        call get_child_leaf(nodeindex1,node,ifirstincell,childleaf1,childleafcounter1)
        call get_child_leaf(nodeindex2,node,ifirstincell,childleaf2,childleafcounter2)

        do i=1,childleafcounter1
          nonodesfound(childleaf1(i)) = nonodesfound(childleaf1(i)) +  childleafcounter2
        enddo 
        !$omp end critical (nodefound)



#ifdef FINVSQRT
       dr = finvsqrt(r2)
#else
       dr = 1./sqrt(r2)
#endif

#ifdef GRAVITY
       ! This should be a lock but this will do for now 
       !$omp critical (node)
       !print*, "fnode before: ", fnode(1:3)
       !fnode = node(nodeindex1) % fnode 
       fnode = 0. 
       quads = 0. 
       c0 = 0.
       c1 = 0.
       c2 = 0.
       c3 = 0.
       !call compute_coeff(dx,dy,dz,dr,totmass_node,quads,c0,c1,c2,c3)
       !if (ifirstincell(nodeindex1) <= 0 .and. ifirstincell(nodeindex2) <= 0) then 
       ! print*, "well separated not leaf:", nodeindex1, nodeindex2
        !stop 
       !endif 
       !if (nodeindex1 == 11 .and. nodeindex2 == 16 .or. nodeindex1 == 16 .and. nodeindex2 == 11) then 
       call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
     !else if (nodeindex1 == 11 .and. nodeindex2 == 17 .or. nodeindex1 == 17 .and. nodeindex2 == 11) then 
       !call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
       !endif 
       !if (nodeindex1 ==35) then 
       !print*, "nodeindex2: ", nodeindex2
      ! print*, "fnode is: ", fnode(1:3)
      ! print*, "c1 is: ", c1 
      ! print*, "force on node1 is: ", fnode(1)*node(nodeindex1) % mass, &
      !  fnode(2)*node(nodeindex1) % mass, fnode(3)*node(nodeindex1) % mass

       !print*, "force actual", -dr*dr*dr*dx*totmass_node, -dr*dr*dr*dy*totmass_node, -dr*dr*dr*dz*totmass_node  
      !endif 
       !stop 
       ! STORE FNODE 
       node(nodeindex1) % fnode = node(nodeindex1) % fnode + fnode 
       ! call get_node_data(node(nodeindex2),node(nodeindex1),dx,dy,dz,xsizei,xsizej,fnode,totmass_node,quads)
       ! r2    = dx*dx + dy*dy + dz*dz
       ! dr    = 1./sqrt(r2)
       ! fnode = 0. 
       ! quads = 0. 
       ! call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
       !  print*, "force on node2 is: ", fnode(1)*node(nodeindex2) % mass, &
       ! fnode(2)*node(nodeindex2) % mass, fnode(3)*node(nodeindex2) % mass
       !$omp end critical (node)
#endif 
       !endif 


    !  If nodes are not well_separated and not equal then splitting needs to occur 
    else 

      !print*, "sizes: ", node(nodeindex1) % size, node(nodeindex2) % size
      ! Find the splitnode 
      if (node(nodeindex1) % size > node(nodeindex2) % size) then
        !splitnode = nodes(nodeindex1)
        splitnodeindex = nodeindex1
        !regnode = nodes(nodeindex2)
        regnodeindex = nodeindex2
      else if (node(nodeindex2) % size > node(nodeindex1) % size)  then 
        !splitnode = nodes(nodeindex2)
        splitnodeindex = nodeindex2
        !regnode = nodes(nodeindex1)
        regnodeindex = nodeindex1
      else 
        ! If sizes are same take numberically larger node 
        if (nodeindex1 > nodeindex2) then 
          splitnodeindex = nodeindex1
          !regnode = nodes(nodeindex2)
          regnodeindex = nodeindex2
        else 
          splitnodeindex = nodeindex2
          !regnode = nodes(nodeindex2)
          regnodeindex = nodeindex1
        endif 

      endif 

      !split larger node and push interactions to stack 
      leftindex = node(splitnodeindex) % leftchild
      rightindex = node(splitnodeindex) % rightchild


      ! ! If split node is not leaf 
       if (ifirstincell(splitnodeindex) <= 0) then 

      !                    !print*,"leftindex: ", leftindex
      !                    !print*,"rightindex: ", rightindex

                        ! If splitting nodeindex1  push to global stack 
                        if (nodeindex1 == splitnodeindex) then 
                              !$omp critical (stack)
                              if (leftindex /= 0) then 
                                    call push_global(leftindex,regnodeindex,0,stack,top)
                              endif 
                              if (rightindex /= 0) then 
                                    call push_global(rightindex,regnodeindex,0,stack,top)
                              endif 
                              !$omp end critical (stack)

                        ! Otherwise if splitting nodeindex2 keep interactions on local stack 
                        else 
                  
                              if (leftindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = leftindex
                                    call push_local(regnodeindex,leftindex,stacklocal,istacklocal,k)
                              endif

                              if (rightindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = rightindex
                                    call push_local(regnodeindex,rightindex,stacklocal,istacklocal,k)
                              endif 

                        endif
    ! Leaf Node - Node 
    ! If split node is leaf but regnode is not then perform direct sum 
    else if (ifirstincell(splitnodeindex) > 0 .and. ifirstincell(regnodeindex) <= 0) then 
        ! ! Check if sph neighbour 
        ! if (r2 < rcut2) then 
        !   ! Call direct sum with softening 
        ! else
        !   ! Normal unsoftened direct sum  
        !   call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        ! endif 

         !split larger node 

        leftindex = node(regnodeindex) % leftchild
        rightindex = node(regnodeindex) % rightchild

                         ! If splitting node push to global stack 
                        if (nodeindex1 == regnodeindex) then 
                              !$omp critical (stack)
                              if (leftindex /= 0) then 
                                    ! top = top + 1
                                    ! stack(top) % nodeindex1 = leftindex
                                    ! stack(top) % interactions(1) = splitnodeindex
                                    call push_global(leftindex,splitnodeindex,0,stack,top)
                              endif 
                              if (rightindex /= 0) then 
                                    ! top = top + 1
                                    ! stack(top) % nodeindex1 = rightindex
                                    ! stack(top) % interactions(1) = splitnodeindex
                                    call push_global(rightindex,splitnodeindex,0,stack,top)
                              endif 
                              !$omp end critical (stack)

                        else 
                  
                              if (leftindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = splitnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = leftindex
                                    call push_local(splitnodeindex,leftindex,stacklocal,istacklocal,k)
                              endif

                              if (rightindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = splitnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = rightindex
                                    call push_local(splitnodeindex,rightindex,stacklocal,istacklocal,k)
                              endif 

                        endif 


    ! Leaf Node - Leaf Node 
    else  
        !print*, "direct leaf nodes"
        !print*, "ifirstincells:", ifirstincell(nodeindex1),ifirstincell(nodeindex2)
        ! Check if sph neighbour 
        if (r2 < rcut2) then 
          ! Call direct sum with softening 
          !$omp critical (nodefound)
          nonodesfound(nodeindex1) = nonodesfound(nodeindex1) + 1 
          !$omp end critical (nodefound)
          call direct_sum_softened(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        else
          ! Normal unsoftened direct sum 
          !$omp critical(nodefound) 
          nonodesfound(nodeindex1) = nonodesfound(nodeindex1) + 1 
          !$omp end critical(nodefound)
          call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
          ! Fix this sub
          !call direct_sum_softened(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        endif 
        
    endif     
endif
endif  
enddo 

!$omp end parallel 
end subroutine interact_standalone



subroutine interact(node, ifirstincell,xyzh,fxyzu_dtt,poten_dtt,nonodesfound)
use omputils
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
  use omp_lib
  use kernel, only:radkern
  type(kdnode), intent(inout)        :: node(:)
  integer, intent(in)                :: ifirstincell(:)

  real,               intent(in)   :: xyzh(:,:)
  real,               intent(out) :: fxyzu_dtt(:,:) ! Needed for direct sum between neighbours that are not sph neighbours 
  real,               intent(out) :: poten_dtt(:)
  integer, optional, intent(out) :: nonodesfound(:)
  integer :: maxcache
  !integer :: nstack(istacksize)
  integer, parameter                   :: maxstacksize = 100000
  type(denstreestack)       :: stack(maxstacksize)
  type(denstreestacklocal)  :: stacklocal(maxstacksize)
  integer                   :: top,nodeindex1, nodeindex2, leftindex, rightindex, regnodeindex, splitnodeindex, numthreads
  integer, allocatable      :: istacklocal(:),threadinteractions(:),threadcounts(:),idlecounts(:),globalcounts(:)
  integer                   :: k,totalinteractions
  logical, allocatable      :: threadworking(:)
  integer :: n,istack,il,ir,npnode,i1,i2,j,jmax,nn,i
  real :: dx,dy,dz,xsizei,xsizeimax,xsizejmax,xsizej,rcutj,rcuti
  real :: rcut,rcut2,rcut2max,r2,hmax
  real :: xoffset,yoffset,zoffset,tree_acc2
  logical :: open_tree_node
!#ifdef GRAVITY
  real :: quads(6)
  real :: dr,totmass_node
  real :: fnode(20)
  real :: c0,c1(3),c2(3,3),c3(3,3,3)
  ! integer :: neighleaf,childleaf1(10000),childleaf2(10000),childleafcounter1,childleafcounter2,specialnode
!#endif
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz

  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
#endif


!tree_accuracy = 0.0
!tree_acc2 = tree_accuracy*tree_accuracy
!print*, tree_acc2 
!print*, tree_accuracy

!print*, "entered interact"
!print*, "fxyzu :", fxyzu
!stop 
top = 1
stack(top) % nodeindex1 = 1
stack(top) % interactions(1) = 1 
numthreads = 1 
!specialnode = 85
!$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
!$omp end parallel

allocate(istacklocal(numthreads))
allocate(threadworking(numthreads))
allocate(threadinteractions(numthreads))
allocate(threadcounts(numthreads))
allocate(idlecounts(numthreads))
allocate(globalcounts(numthreads))
!print*,"allocation started"
!allocate(stack(maxstacksize))
!allocate(stacklocal(maxstacksize))
!print*,"allocation finished"

threadworking = .false. 
istacklocal = 0
!neighleaf = 0
!nonodesfound = 0 
threadinteractions = 0 
totalinteractions = 0 
threadcounts = 0 
idlecounts = 0 
globalcounts = 0 

! setup locks 
!$ call init_omp()

!$omp parallel default(none) &
!$omp shared(node,ifirstincell,threadworking,top,istacklocal,stack,tree_acc2,fxyzu_dtt,poten_dtt,xyzh) &!specialnode,nonodesfound) &
!$omp private(stacklocal,nodeindex1,nodeindex2,leftindex,rightindex,splitnodeindex,regnodeindex,k) &
!$omp private(dx,dy,dz,totmass_node,xsizei,xsizeimax,xsizej,xsizejmax,quads,r2,open_tree_node,dr,fnode) &
!$omp private(rcut2,rcut2max,rcut,rcuti,rcutj) & 
!!$omp private(childleaf1,childleaf2,childleafcounter1,childleafcounter2,i) &
!$omp private(c0,c1,c2,c3) &
!$omp reduction(+:totalinteractions) &
!$omp shared(threadinteractions,ipart_omp_lock,threadcounts,idlecounts,globalcounts,numthreads)

!!$omp reduction(+:neighleaf)

!!$omp single 
!print*, "entered parallel"
!open(unit=8, file="nodes.txt",position='append')
do while (any(istacklocal > 0) .or. top > 0 .or. any(threadworking))

     !$ k=omp_get_thread_num() + 1
     !print*,"k: ", k
     !print*,tree_acc2
     ! print*, "top: ", top 
     ! print*, "istacklocal: ", istacklocal(k)
     !stop 
     threadcounts(k) = threadcounts(k) + 1 

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
          call pop_local(nodeindex1,nodeindex2,stacklocal,istacklocal,k)
          ! nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
          ! nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2
          !print*, "node indexes :", nodeindex1, nodeindex2
          ! istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
          globalcounts(k) = globalcounts(k) + 1 
          !threadworking(k) = .false.
          !nodes(nodeindex1) % nodefree  = .false.
          !nodes(nodeindex2) % nodefree  = .false.
      else 
        threadworking(k) = .false.
        idlecounts(k) = idlecounts(k) + 1
      endif
      !$omp end critical (stack)
    endif 

  if (threadworking(k)) then 
      !threadinteractions(k) = threadinteractions(k) + 1
      !totalinteractions = totalinteractions + 1 
      !fnode = 0.
      !dx = 0. 
      !dy = 0. 
      !dz = 0.
      !dr = 0.
      !totmass_node = 0.
      !xsizej = 0.
      !xsizei = 0.
      !quads = 0. 
      !r2 = 0. 
      !open_tree_node = .false. 
      !rcut = 0.
      !rcuti = 0.
      !rcutj = 0. 
      !rcut2 = 0.
      !tree_acc2 = 0.
      !print*, tree_acc2
      ! print*, "nodeindex1: ", nodeindex1
      ! print*, "nodeindex2: ", nodeindex2
      ! print*, "istacklocal", istacklocal(k)
      call get_node_data(node(nodeindex1),node(nodeindex2),dx,dy,dz,xsizei,xsizej)!,fnode,totmass_node,quads)
      !quads = 0.
      r2    = dx*dx + dy*dy + dz*dz
      rcutj = radkern*node(nodeindex2) % hmax
      rcuti = radkern*node(nodeindex1) % hmax 
      rcut  = max(rcuti,rcutj)
      xsizeimax = node(nodeindex1) % maxsize 
      xsizejmax = node(nodeindex2) % maxsize
      !rcut = rcuti
      rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
      rcut2max = (xsizeimax + xsizejmax + rcut)**2   ! node size + search radius
      !!$omp critical (check)
      
      ! if (r2 /= 0. .and. abs(r2 - rcut2max) <= 1.e-15) then 
      !   print*, "r2: ",r2
      !   print*, "rcut2max: ", rcut2max
      !   print*, "r2 - rcut2max: ", r2 - rcut2max
      !   print*, "condition: ",abs(r2 - rcut2max) <= 1.e-15
      !   print*, "second val: ", 1.e-15
      !   stop 
      ! endif 
      !!$omp end critical (check)
      ! If nodes are equal 
      if (nodeindex1 == nodeindex2) then 
        ! if (ifirstincell(nodeindex1) > 0) then 
        !    !print*, "leaf self interaction"
        !    !call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        !    !stop 
        ! endif 
        ! If nodes are equal but not leaf then split 
        if (ifirstincell(nodeindex1) <= 0) then 
                  leftindex = node(nodeindex1) % leftchild
                  rightindex = node(nodeindex1) % rightchild 

                  ! print*, "leftindex: ", leftindex
                  ! print*, "rightindex: ", rightindex
                  

                  if (leftindex /= 0 .and. rightindex /= 0) then 

                        ! If local stack sufficently full push to global stack
                        !if (nodeindex1 <  64) then
                          
                          if (top < numthreads) then
                            !$omp critical (stack) 
                            call push_global(leftindex,leftindex,rightindex,stack,top)
                            call push_global(rightindex,rightindex,leftindex,stack,top)
                            !$omp end critical (stack)
                          else 
                            call push_local(leftindex,leftindex,stacklocal,istacklocal,k)
                            call push_local(leftindex,rightindex,stacklocal,istacklocal,k)
                            call push_local(rightindex,rightindex,stacklocal,istacklocal,k)
                            call push_local(rightindex,leftindex,stacklocal,istacklocal,k)
                          endif 
                          

                        ! else 
                        !   call push_local(leftindex,leftindex,stacklocal,istacklocal,k)
                        !   call push_local(leftindex,rightindex,stacklocal,istacklocal,k)
                        !   call push_local(rightindex,rightindex,stacklocal,istacklocal,k)
                        !   call push_local(rightindex,leftindex,stacklocal,istacklocal,k)
                         
                        ! endif   
                        !endif 
                        ! if (istacklocal(k) < 100000) then 
                        !   call push_local(leftindex,leftindex,stacklocal,istacklocal,k)
                        !   call push_local(leftindex,rightindex,stacklocal,istacklocal,k)
                        !   call push_local(rightindex,rightindex,stacklocal,istacklocal,k)
                        !   call push_local(rightindex,leftindex,stacklocal,istacklocal,k)
                        ! else


                  endif 

                  

        endif

      ! If well separated but not sph neighbour, calculate long range interactions 
      else if (well_separated(node(nodeindex1), node(nodeindex2)) .and. r2 >= rcut2max)  then   ! &
       !.and. ifirstincell(nodeindex1) <= 0 .and. ifirstincell(nodeindex2) <= 0) then  ! &
        !.and. .not. check_child_overlap(node,nodeindex1,nodeindex2).and. .not. check_child_overlap(node,nodeindex2,nodeindex1)) then
        !print*, "Well separated nodes: ", nodeindex1, nodeindex2 
        !stop 
        !write(13,*) nodeindex1,nodeindex2
        !print*, "r2 and rcut2: ", r2, rcut2
        !if (ifirstincell(nodeindex1) > 0) stop 
        !if (ifirstincell(nodeindex2) > 0) stop 
        ! if ( .not. well_separated(node(nodeindex2),node(nodeindex1))) then 
        !     print*, "MAC BROKEN"
        !     stop
        ! endif  

        !!$omp critical (nodefound)
        ! call get_child_leaf(nodeindex1,node,ifirstincell,childleaf1,childleafcounter1)
        ! call get_child_leaf(nodeindex2,node,ifirstincell,childleaf2,childleafcounter2)

        ! do i=1,childleafcounter1
        !   nonodesfound(childleaf1(i)) = nonodesfound(childleaf1(i)) +  childleafcounter2
        !   if (childleaf1(i) == 66568) then 
        !     do j=1,childleafcounter2
        !       write(8,*) childleaf2(j)
        !        if (childleaf2(j) == 133204) then 
        !         ! print*, "r2 >= rcut2 ", r2 >= rcut2max
        !         ! print*,"r2: ", r2
        !         ! print*,"rcut2: ",rcut2max 
        !         ! print*, "rcut delta", r2 - rcut2max
        !         ! print*, "nodeindex1: ", nodeindex1
        !         ! print*, "nodeindex2: ", nodeindex2 
        !         ! dx = node(childleaf1(i)) % xcen(1) - node(childleaf2(j)) % xcen(1)
        !         ! dy = node(childleaf1(i)) % xcen(2) -  node(childleaf2(j)) % xcen(2)
        !         ! dz = node(childleaf1(i)) % xcen(3) -  node(childleaf2(j)) % xcen(3)
        !         ! r2 = dx*dx + dy*dy + dz*dz  
        !         ! print*, "r2 child ", r2 
        !         ! print*, "radkern1 ",radkern*node(childleaf1(i)) % hmax,radkern*node(childleaf2(j)) % hmax
        !         ! rcut = max(radkern*node(childleaf1(i)) % hmax, radkern*node(childleaf2(j)) % hmax) 
        !         ! print*, "rcut2 arg ", node(childleaf1(i))%maxsize, node(childleaf2(j))%maxsize, rcut
        !         ! print*, "rcut2 arg ", node(childleaf1(i))%size, node(childleaf2(j))%size, rcut
                
        !         ! print*, "rcut2 ", (node(childleaf1(i))%maxsize + node(childleaf2(j))%maxsize + rcut)**2
        !          !stop 
        !        endif 
        !     enddo 

        !   endif 
        ! enddo 
        ! !$omp end critical (nodefound)

       ! print*,"well_separated"
#ifdef FINVSQRT
       dr = finvsqrt(r2)
#else
       dr = 1./sqrt(r2)
#endif

#ifdef GRAVITY
       ! This should be a lock but this will do for now 
       !!$omp critical (node)
       !print*, "fnode before: ", fnode(1:3)
       !fnode = node(nodeindex1) % fnode 
       fnode = 0. 
       quads = node(nodeindex2) % quads
       totmass_node = node(nodeindex2) % mass
       ! c0 = 0.
       ! c1 = 0.
       ! c2 = 0.
       ! c3 = 0.
       !call compute_coeff(dx,dy,dz,dr,totmass_node,quads,c0,c1,c2,c3)
       !if (ifirstincell(nodeindex1) <= 0 .and. ifirstincell(nodeindex2) <= 0) then 
        !print*, "well separated not leaf:", nodeindex1, nodeindex2
        !stop 
       !endif 
       !if (nodeindex1 == 11 .and. nodeindex2 == 16 .or. nodeindex1 == 16 .and. nodeindex2 == 11) then 
       call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
       !!$omp critical (node)
       !$ call omp_set_lock(ipart_omp_lock(nodeindex1/10))
       node(nodeindex1) % fnode = node(nodeindex1) % fnode + fnode 
       !$ call omp_unset_lock(ipart_omp_lock(nodeindex1/10))
     !else if (nodeindex1 == 11 .and. nodeindex2 == 17 .or. nodeindex1 == 17 .and. nodeindex2 == 11) then 
       !call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
       !endif 
       !if (nodeindex1 ==35) then 
       !print*, "nodeindex2: ", nodeindex2
      ! print*, "fnode is: ", fnode(1:3)
      ! print*, "c1 is: ", c1 
      ! print*, "force on node1 is: ", fnode(1)*node(nodeindex1) % mass, &
      !  fnode(2)*node(nodeindex1) % mass, fnode(3)*node(nodeindex1) % mass

       !print*, "force actual", -dr*dr*dr*dx*totmass_node, -dr*dr*dr*dy*totmass_node, -dr*dr*dr*dz*totmass_node  
      !endif 
       !stop 
       ! STORE FNODE 
       ! !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(1) = node(nodeindex1) % fnode(1) + fnode(1)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(2) = node(nodeindex1) % fnode(2) + fnode(2) 
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(3) = node(nodeindex1) % fnode(3) + fnode(3)
       !   !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(4) = node(nodeindex1) % fnode(4) + fnode(4)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(5) = node(nodeindex1) % fnode(5) + fnode(5) 
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(6) = node(nodeindex1) % fnode(6) + fnode(6)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(7) = node(nodeindex1) % fnode(7) + fnode(7)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(8) = node(nodeindex1) % fnode(8) + fnode(8) 
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(9) = node(nodeindex1) % fnode(9) + fnode(9)
       !   !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(10) = node(nodeindex1) % fnode(10) + fnode(10)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(11) = node(nodeindex1) % fnode(11) + fnode(11) 
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(12) = node(nodeindex1) % fnode(12) + fnode(12)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(13) = node(nodeindex1) % fnode(13) + fnode(1)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(14) = node(nodeindex1) % fnode(14) + fnode(14) 
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(15) = node(nodeindex1) % fnode(15) + fnode(15)
       !   !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(16) = node(nodeindex1) % fnode(16) + fnode(16)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(17) = node(nodeindex1) % fnode(17) + fnode(17) 
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(18) = node(nodeindex1) % fnode(18) + fnode(18)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(19) = node(nodeindex1) % fnode(19) + fnode(19)
       !  !$omp ATOMIC UPDATE 
       !  node(nodeindex1) % fnode(20) = node(nodeindex1) % fnode(20) + fnode(20) 
       ! call get_node_data(node(nodeindex2),node(nodeindex1),dx,dy,dz,xsizei,xsizej,fnode,totmass_node,quads)
       ! r2    = dx*dx + dy*dy + dz*dz
       ! dr    = 1./sqrt(r2)
       ! fnode = 0. 
       ! !quads = 0. 
       ! call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
      !   print*, "force on node2 is: ", fnode(1)*node(nodeindex2) % mass, &
      ! fnode(2)*node(nodeindex2) % mass, fnode(3)*node(nodeindex2) % mass
      !  stop 
       !!$omp end critical (node)
#endif 
       !endif 


    !  If nodes are not well_separated and not equal then splitting needs to occur 
    else 
      ! print*, "split 1 "
      !print*, "sizes: ", node(nodeindex1) % size, node(nodeindex2) % size
      ! Find the splitnode 
      if (node(nodeindex1) % size > node(nodeindex2) % size) then
        !splitnode = nodes(nodeindex1)
        splitnodeindex = nodeindex1
        !regnode = nodes(nodeindex2)
        regnodeindex = nodeindex2
      else if (node(nodeindex2) % size > node(nodeindex1) % size)  then 
        !splitnode = nodes(nodeindex2)
        splitnodeindex = nodeindex2
        !regnode = nodes(nodeindex1)
        regnodeindex = nodeindex1
      else 
        ! If sizes are same take numberically larger node 
        if (nodeindex1 > nodeindex2) then 
          splitnodeindex = nodeindex1
          !regnode = nodes(nodeindex2)
          regnodeindex = nodeindex2
        else 
          splitnodeindex = nodeindex2
          !regnode = nodes(nodeindex2)
          regnodeindex = nodeindex1
        endif 

      endif 

      !split larger node and push interactions to stack 
      leftindex = node(splitnodeindex) % leftchild
      rightindex = node(splitnodeindex) % rightchild


      ! ! If split node is not leaf 
       if (ifirstincell(splitnodeindex) <= 0) then 

      !                    !print*,"leftindex: ", leftindex
      !                    !print*,"rightindex: ", rightindex

                        ! If splitting nodeindex1  push to global stack 
                        if (nodeindex1 == splitnodeindex) then 
                              ! print*, "push global"
                              ! print*, leftindex, rightindex
                              ! print*, "top: ", top 
                              
                              if (leftindex /= 0) then 
                                    if (top < numthreads) then 
                                      !$omp critical (stack)
                                      call push_global(leftindex,regnodeindex,0,stack,top)
                                      !$omp end critical (stack)
                                    else 
                                     call push_local(leftindex,regnodeindex,stacklocal,istacklocal,k)
                                    endif 
                              endif 
                              
                              !!$omp end critical (stack)
                              if (rightindex /= 0) then 
                                    if (top < numthreads) then 
                                      !$omp critical (stack)
                                      call push_global(rightindex,regnodeindex,0,stack,top)
                                      !$omp end critical (stack)
                                    else
                                      call push_local(rightindex,regnodeindex,stacklocal,istacklocal,k)
                                    endif 
                              endif 
                              

                        ! Otherwise if splitting nodeindex2 keep interactions on local stack 
                        else 
                          ! print*, "push local"
                  
                              if (leftindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = leftindex
                                    call push_local(regnodeindex,leftindex,stacklocal,istacklocal,k)
                              endif

                              if (rightindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = rightindex
                                    call push_local(regnodeindex,rightindex,stacklocal,istacklocal,k)
                              endif 

                        endif
    ! Leaf Node - Node 
    ! If split node is leaf but regnode is not then perform direct sum 
    else if (ifirstincell(splitnodeindex) > 0 .and. ifirstincell(regnodeindex) <= 0) then 
        
         ! print*, "Leaf Node - Node"
        ! If we are not sph trial neihgbours
         ! if (r2<rcut2) then 
         !  if (.not. check_parent_overlap(node,nodeindex1,nodeindex2)) then 
         !       call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
         !   endif 
         ! else 
         !     call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
         ! endif 

        ! find all the leaf nodes of the regnode 
        !call directsum_not_leaf(nodeindex1, nodeindex2,splitnodeindex,regnodeindex,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        !split larger node 

        leftindex = node(regnodeindex) % leftchild
        rightindex = node(regnodeindex) % rightchild

                         ! If splitting node push to global stack 
                        if (nodeindex1 == regnodeindex) then 
                              
                              if (leftindex /= 0) then 
                                    ! top = top + 1
                                    ! stack(top) % nodeindex1 = leftindex
                                    ! stack(top) % interactions(1) = splitnodeindex
                                    if (top < numthreads) then 
                                      !$omp critical (stack)
                                      call push_global(leftindex,splitnodeindex,0,stack,top)
                                      !$omp end critical (stack)
                                    else
                                      call push_local(leftindex,splitnodeindex,stacklocal,istacklocal,k)
                                    endif 
                              endif 
                              

                              if (rightindex /= 0) then 
                                    ! top = top + 1
                                    ! stack(top) % nodeindex1 = rightindex
                                    ! stack(top) % interactions(1) = splitnodeindex
                                    if (top < numthreads) then 
                                      !$omp critical (stack)
                                      call push_global(rightindex,splitnodeindex,0,stack,top)
                                      !$omp end critical (stack)
                                    else 
                                      call push_local(rightindex,splitnodeindex,stacklocal,istacklocal,k)
                                    endif 
                              endif 
                              

                        else 
                  
                              if (leftindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = splitnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = leftindex
                                    call push_local(splitnodeindex,leftindex,stacklocal,istacklocal,k)
                              endif

                              if (rightindex /= 0) then 
                                    ! istacklocal(k) = istacklocal(k) + 1
                                    ! stacklocal(istacklocal(k)) % nodeindex1 = splitnodeindex
                                    ! stacklocal(istacklocal(k)) % nodeindex2 = rightindex
                                    call push_local(splitnodeindex,rightindex,stacklocal,istacklocal,k)
                              endif 

                        endif 
    ! Leaf Node - Leaf Node 
    else  
        ! print*, "direct leaf nodes"
        !print*, "ifirstincells:", ifirstincell(nodeindex1),ifirstincell(nodeindex2)
        
         if (r2<rcut2) then 
          if (.not. check_parent_overlap(node,nodeindex1,nodeindex2)) then
               !if (nodeindex1 == specialnode) neighleaf = neighleaf + 1
               ! !$omp critical(nodefound) 
               ! nonodesfound(nodeindex1) = nonodesfound(nodeindex1) + 1
               ! !$omp end critical(nodefound)
               ! if (nodeindex1 == 66568) then 
               ! print*, "direct neighbour" 
               ! if (nodeindex2 == 133204) print*, "problem in direct"

               !write(8,*) nodeindex2
              !endif 
               call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
           endif 
         else 
             !if (nodeindex1 == specialnode ) neighleaf = neighleaf + 1
             ! !$omp critical(nodefound)
             ! nonodesfound(nodeindex1) = nonodesfound(nodeindex1) + 1
             ! !$omp end critical(nodefound)
             ! if (nodeindex1 == 66568) then 
             !  print*, "direct not neigh"
             ! if (nodeindex2 == 133204) print*, "problem in direct"
             
             ! !write(8,*) nodeindex2
            
             ! endif 
             call direct_sum_not_neighbour(nodeindex1, nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
         endif 
          

      endif 

  endif 
endif 
enddo 

! !$omp critical (counts)
! print*, "thread number: ", k 
! print*, "idlecounts: ", idlecounts(k)
! print*, "globalcounts: ", globalcounts(k)
! print*, "total counts: ", threadcounts(k)
! !$omp end critical (counts)

!!$omp end single 


!$omp end parallel 

! print*, "Total interactions: ", totalinteractions
! do i=1, numthreads
!   print*, "Number of interactions performed by thread: ", threadinteractions(i)
!   print*, "Percentage of total interactions performed: ", real(threadinteractions(i))/real(totalinteractions)*100.
! enddo 
!print*, "neighleaf: ",neighleaf
!print*,"interact finished"
!close(8)
!stop 


end subroutine interact

subroutine evaluate(node, ifirstincell)
  use omp_lib
  type(kdnode), intent(inout)        :: node(:)
  integer, intent(in)                :: ifirstincell(:)
   integer :: stacksize, top, toplocal 
   ! Openmp disables heap allocation so this should be allocatable 
   type(evaluate_stack_data) :: stack(1000)
   real :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),z1(3)
   real :: bodaccel(3),xbod(3),fnode(20),bodpot
   integer :: bodyindex,i,j,iter,numthreads,toptemp,k
   integer :: currentnodeindex, childnodeindex
   real    :: dx,dy,dz
   logical, allocatable :: threadworking(:)
   integer, allocatable :: istacklocal(:)

   top = 1
   iter = 1

   ! Push the root node 
   stack(1) % nodeindex = 1
   stack(1) % fnode = 0.
   stack(1) % z0 = 0. 
   currentnodeindex = 1 
   numthreads = 1
   ! get number of OpenMPthreads
   !$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
   !$omp end parallel
   allocate(threadworking(numthreads))
   !allocate(istacklocal(numthreads))
   threadworking = .true.


  !$omp parallel default(none) &
  !$omp shared(stack,top,node,numthreads,threadworking) &!,istacklocal) &
  !$omp private(z0,z1,c0,c1,c2,c3) &
  !$omp private(currentnodeindex,bodyindex,xbod,bodaccel,childnodeindex,dx,dy,dz,k,fnode)
  
   !$ k=omp_get_thread_num() + 1
   ! local stack is currently empty 
   !istacklocal(k) = 0 
   ! Just for compiler errors 
   !currentnode = node(1)
   !print*,"Thread is: ", k
   ! here top is the top of the global stack 
   ! Push the root node 
   
  over_stack: do while (top > 0)  
   !do j=1, top
    !print*, "Thread number is: ", j
    !print*, "The current iteration is: ",iter
    !print*, "Top is: ", top
    !print*, "local stack top: ", istacklocal(k)
    !iter = iter + 1 
    ! POP ITEM FROM STACK 
    !print*,"istacklocal: ",istacklocal(k)
    ! pop of local stack 
    !if (istacklocal(k) > 0 ) then 
    !    currentnode = stack(istacklocal(k)) % node
    !    z0 = stack(istacklocal(k)) % z0
    !    c0 = stack(istacklocal(k)) % c0
     !   c1 = stack(istacklocal(k)) % c1
     !   c2 = stack(istacklocal(k)) % c2 
     !   c3 = stack(istacklocal(k)) % c3 

     !   istacklocal(k) = istacklocal(k) - 1 
      !  threadworking(k) = .true.
    !else 
    !$omp critical(globalstack) 
        if (top > 0) then 
        currentnodeindex = stack(top) % nodeindex
        z0 = stack(top) % z0
        fnode = stack(top) % fnode 

        top = top - 1 
        threadworking(k) = .true. 
        !print*,"working on thread",k
        else
        threadworking(k) = .false.
        
        endif 
    !$omp end critical(globalstack)
    !endif 
    !if top
  ! If thread has work to do, do it 
  if (threadworking(k)) then 
    
     ! print*, "data"
     ! print*,"current node", currentnodeindex
     ! print*,"force ",fnode(1:3)
     ! print*,"poten", fnode(20)
     ! print*,"xcen prev: ",z0


    ! Perform evaluate on popped item 

   z1 = node(currentnodeindex) % xcen 
  

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  !print*, "Translating expansion: "
  !call translate_expansion_center(z0,z1,c0,c1,c2,c3)
  dx = z1(1)-z0(1)
  dy = z1(2)-z0(2)
  dz = z1(3)-z0(3)

  ! WILL ONLY WORK FOR A SECOND ORDER EXPANSION
  ! this is disgusting but needed for compilation 
#ifdef GRAVITY 
  !!$omp critical (node)
  ! c0 = fnode(20)
  ! c1(1) =  fnode(1)
  ! c1(2) =  fnode(2)
  ! c1(3) =  fnode(3)
  ! c2(1,1) = fnode(4)
  ! c2(1,2) = fnode(5)
  ! c2(1,3) = fnode(6)
  ! c2(2,1) = fnode(5)
  ! c2(2,2) = fnode(7)
  ! c2(2,3) = fnode(8)
  ! c2(3,1) = fnode(6)
  ! c2(3,2) = fnode(8)
  ! c2(3,3) = fnode(9)
  !if (fnode(20) /= 0.) then 
  !print*, fnode(20)
  !stop 
  !print*, "fnode before: ", fnode
  !call expand_fgrav_in_taylor_series(fnode,dx,dy,dz,c1(1),c1(2),c1(3),c0,c2(1),c2(2),c2(3),c2(4),c2(5),c2(6))
  call translate_fgrav_in_taylor_series(fnode,dx,dy,dz)
 !call translate_expansion_center(z0,z1,c0,c1,c2,c3)
  ! print*,"value comparision"
  ! print*, fnode(1),fnode(2),fnode(3)
  ! print*, c1(1),c1(2),c1(3)
  !!$omp critical (node)
  ! !node(currentnodeindex) % fnode    = 0. 
  ! ! ! c1
  !  node(currentnodeindex) % fnode(1) = node(currentnodeindex) % fnode(1) +  c1(1) 
  !  node(currentnodeindex) % fnode(2) = node(currentnodeindex) % fnode(2) + c1(2)
  !  node(currentnodeindex) % fnode(3) = node(currentnodeindex) % fnode(3) + c1(3)
  ! ! ! ! c2 
  ! !  node(currentnodeindex) % fnode(4) = node(currentnodeindex) % fnode(4) + c2(1)
  ! !  node(currentnodeindex) % fnode(5) = node(currentnodeindex) % fnode(5) + c2(2)
  ! !  node(currentnodeindex) % fnode(6) = node(currentnodeindex) % fnode(6) + c2(3)
  ! !  node(currentnodeindex) % fnode(7) = node(currentnodeindex) % fnode(7) + c2(4)
  ! !  node(currentnodeindex) % fnode(8) = node(currentnodeindex) % fnode(8) + c2(5)
  ! !  node(currentnodeindex) % fnode(9) = node(currentnodeindex) % fnode(9) + c2(6)

  
  !  node(currentnodeindex) % fnode(4:20) = node(currentnodeindex) % fnode(4:20) + fnode(4:20)
  ! node(currentnodeindex) % fnode(20) = node(currentnodeindex) % fnode(20) + c0 
  !!$omp end critical (node)
  !node(currentnodeindex) % fnode(1) = 2.
  !print*, "fnode trans: ", fnode
  !!$omp critical (node) 
  node(currentnodeindex) % fnode = node(currentnodeindex) % fnode + fnode 
  !!$omp end critical (node)
  !print*, "fnode node: ", node(currentnodeindex) % fnode 
  !endif 
  fnode = node(currentnodeindex) % fnode
  !print*, fnode
  
#endif 
  ! FOR CHILDREN OF A 
  ! change center of mass
  z0 = z1
   if (node(currentnodeindex) % leftchild /= 0) then  
    !print*, "Childnode index: ", currentnode % children(i)
    childnodeindex = node(currentnodeindex) % leftchild
    ! if threads are waiting push to global stack 
      !$omp critical(globalstack)
      top = top + 1 
      ! Push Children onto stack 
      stack(top) % nodeindex = childnodeindex
      stack(top) % z0 = z0
      stack(top) % fnode = fnode
      !$omp end critical(globalstack)


    endif 

   if (node(currentnodeindex) % rightchild /= 0) then  
    !print*, "Childnode index: ", currentnode % children(i)
    childnodeindex = node(currentnodeindex) % rightchild
    ! if threads are waiting push to global stack 
      !$omp critical(globalstack)
      top = top + 1 
      ! Push Children onto stack 
      stack(top) % nodeindex = childnodeindex
      stack(top) % z0 = z0
      stack(top) % fnode = fnode
      !$omp end critical(globalstack)


   endif 

 endif 



    
  !print*, "top is: ", top 
enddo over_stack
!stop 
!$omp end parallel


end subroutine evaluate

subroutine directsum_not_leaf(nodeindex1, nodeindex2,splitnodeindex,regnodeindex,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
  !use part, only :massoftype,maxphase,iamtype,iphase
  !use dim,  only :maxp
  use kernel, only:radkern
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
  integer, intent(in) :: nodeindex1, nodeindex2,splitnodeindex,regnodeindex, ifirstincell(:)
  type(kdnode), intent(in) :: node(:)
  real, intent(inout) :: fxyzu_dtt(:,:), poten_dtt(:)
  real, intent(in)    :: xyzh(:,:)
  integer :: n, il,ir,istack
  integer, parameter :: stacksize = 300
  integer :: nstack(stacksize)
  real :: xpos(3), dx,dy,dz,rcut,rcuti,rcutj,rcut2
  real :: xsizei,xsizej,r2
  real :: xoffset,yoffset,zoffset
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz
  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
#endif

  xpos = node(splitnodeindex) % xcen
  rcuti = radkern*node(splitnodeindex)%hmax
  xsizei = node(splitnodeindex) % size
  nstack = 0
  istack = 1 
  nstack(istack) = regnodeindex 

    ! Find the children of regnode 
    do while(istack /=0)
      n = nstack(istack)
      istack = istack - 1 

      dx = xpos(1) - node(n)%xcen(1)      ! distance between node centres
      dy = xpos(2) - node(n)%xcen(2)
#ifndef TREEVIZ
      dz = xpos(3) - node(n)%xcen(3)
#endif
      xsizej       = node(n)%size
      il      = node(n)%leftchild
      ir      = node(n)%rightchild
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
    r2    = dx*dx + dy*dy + dz*dz !+ xsizei*xsizei
    !if (get_hj) then  ! find neighbours within both hi and hj
    rcutj = radkern*node(n)%hmax
    rcut  = max(rcuti,rcutj)
    !endif
    !rcut = rcuti
    ! call get_node_radius(node,ifirstincell,n,radius)
    ! rcut = max(radius,rcuti)
    rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius

    if (r2<rcut2) then 
      ! If we have a leafnode then 
      if (ifirstincell(n) > 0 ) then
        if (nodeindex1 == splitnodeindex) then 
          call direct_sum_not_neighbour(nodeindex1,n,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        else 
          call direct_sum_not_neighbour(n,nodeindex2,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
        endif 

      ! Push children to stack 
      else

        if (il /= 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
        if (ir /= 0) then
             istack = istack + 1
             nstack(istack) = ir
        endif

      endif 

    endif 


    enddo 
  

end subroutine directsum_not_leaf

!-----------------------------------------------------------
!+
!  Compute the gravitational force between the node centres
!  along with the derivatives and second derivatives
!  required for the Taylor series expansions.
!+
!-----------------------------------------------------------
pure subroutine compute_fnode(dx,dy,dz,dr,totmass,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr,totmass
 real, intent(in)    :: quads(6)
 real, intent(inout) :: fnode(lenfgrav)
 real :: dr3,dr4,dr5,dr6,dr3m,rx,ry,rz,qxx,qxy,qxz,qyy,qyz,qzz
 real :: dr4m3,rijQij,riQix,riQiy,riQiz,fqx,fqy,fqz
 real :: dfxdxq,dfxdyq,dfxdzq,dfydyq,dfydzq,dfzdzq
 real :: d2fxxxq,d2fxxyq,d2fxxzq,d2fxyyq,d2fxyzq
 real :: d2fxzzq,d2fyyyq,d2fyyzq,d2fyzzq,d2fzzzq

 ! note: dr == 1/sqrt(r2)
 dr3  = dr*dr*dr
 dr4  = dr*dr3
 dr5  = dr*dr4
 dr6  = dr*dr5
 dr3m  = totmass*dr3
 dr4m3 = 3.*totmass*dr4
 rx  = dx*dr
 ry  = dy*dr
 rz  = dz*dr
 qxx = quads(1)
 qxy = quads(2)
 qxz = quads(3)
 qyy = quads(4)
 qyz = quads(5)
 qzz = quads(6)
 rijQij = (rx*rx*qxx + ry*ry*qyy + rz*rz*qzz + 2.*(rx*ry*qxy + rx*rz*qxz + ry*rz*qyz))
 riQix = (rx*qxx + ry*qxy + rz*qxz)
 riQiy = (rx*qxy + ry*qyy + rz*qyz)
 riQiz = (rx*qxz + ry*qyz + rz*qzz)
 fqx = dr4*(riQix - 2.5*rx*rijQij)
 fqy = dr4*(riQiy - 2.5*ry*rijQij)
 fqz = dr4*(riQiz - 2.5*rz*rijQij)
 ! dfxdxq = dr5*(qxx - 10.*rx*riQix - 2.5*rijQij   + 17.5*rx*rx*rijQij)
 ! dfxdyq = dr5*(qxy -  5.*ry*riQix - 5.0*rx*riQiy + 17.5*rx*ry*rijQij)
 ! dfxdzq = dr5*(qxz -  5.*rx*riQiz - 5.0*rz*riQix + 17.5*rx*rz*rijQij)
 ! dfydyq = dr5*(qyy - 10.*ry*riQiy - 2.5*rijQij   + 17.5*ry*ry*rijQij)
 ! dfydzq = dr5*(qyz -  5.*ry*riQiz - 5.0*rz*riQiy + 17.5*ry*rz*rijQij)
 ! dfzdzq = dr5*(qzz - 10.*rz*riQiz - 2.5*rijQij   + 17.5*rz*rz*rijQij)
 ! d2fxxxq = dr6*(-15.*qxx*rx + 105.*rx*rx*riQix - 15.*riQix - 157.5*rx*rx*rx*rijQij + 52.5*rx*rijQij)
 ! d2fxxyq = dr6*(35.*rx*rx*riQiy -  5.*qxx*ry - 5.*riQiy + 17.5*ry*rijQij - 157.5*rx*rx*ry*rijQij &
 !              + 70.*rx*ry*riQix - 10.*qxy*rx)
 ! d2fxxzq = dr6*(35.*rx*rx*riQiz -  5.*qxx*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*rx*rx*rz*rijQij &
 !              + 70.*rx*rz*riQix - 10.*qxz*rx)
 ! d2fxyyq = dr6*(70.*rx*ry*riQiy - 10.*qxy*ry - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*ry*ry*rijQij &
 !              + 35.*ry*ry*riQix -  5.*qyy*rx)
 ! d2fxyzq = dr6*(35.*rx*ry*riQiz -  5.*qyz*rx  &
 !              + 35.*ry*rz*riQix -  5.*qxz*ry  &
 !              + 35.*rx*rz*riQiy -  5.*qxy*rz                             - 157.5*rx*ry*rz*rijQij)
 ! d2fxzzq = dr6*(70.*rx*rz*riQiz - 10.*qxz*rz - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*rz*rz*rijQij &
 !              + 35.*rz*rz*riQix -  5.*qzz*rx)
 ! d2fyyyq = dr6*(-15.*qyy*ry + 105.*ry*ry*riQiy - 15.*riQiy - 157.5*ry*ry*ry*rijQij + 52.5*ry*rijQij)
 ! d2fyyzq = dr6*(35.*ry*ry*riQiz -  5.*qyy*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*ry*ry*rz*rijQij &
 !              + 70.*ry*rz*riQiy - 10.*qyz*ry)
 ! d2fyzzq = dr6*(70.*ry*rz*riQiz - 10.*qyz*rz - 5.*riQiy + 17.5*ry*rijQij - 157.5*ry*rz*rz*rijQij &
 !              + 35.*rz*rz*riQiy -  5.*qzz*ry)
 ! d2fzzzq = dr6*(-15.*qzz*rz + 105.*rz*rz*riQiz - 15.*riQiz - 157.5*rz*rz*rz*rijQij + 52.5*rz*rijQij)

 fnode( 1) = fnode( 1) - dx*dr3m + fqx ! fx
 fnode( 2) = fnode( 2) - dy*dr3m + fqy ! fy
 fnode( 3) = fnode( 3) - dz*dr3m + fqz ! fz
 fnode( 4) = fnode( 4) + dr3m*(3.*rx*rx - 1.) !+ dfxdxq ! dfx/dx
 fnode( 5) = fnode( 5) + dr3m*(3.*rx*ry)      !+ dfxdyq ! dfx/dy = dfy/dx
 fnode( 6) = fnode( 6) + dr3m*(3.*rx*rz)      !+ dfxdzq ! dfx/dz = dfz/dx
 fnode( 7) = fnode( 7) + dr3m*(3.*ry*ry - 1.) !+ dfydyq ! dfy/dy
 fnode( 8) = fnode( 8) + dr3m*(3.*ry*rz)      !+ dfydzq ! dfy/dz = dfz/dy
 fnode( 9) = fnode( 9) + dr3m*(3.*rz*rz - 1.) !+ dfzdzq ! dfz/dz
 fnode(10) = fnode(10) - dr4m3*(5.*rx*rx*rx - 3.*rx) !+ d2fxxxq ! d2fxdxdx
 fnode(11) = fnode(11) - dr4m3*(5.*rx*rx*ry - ry)    !+ d2fxxyq ! d2fxdxdy
 fnode(12) = fnode(12) - dr4m3*(5.*rx*rx*rz - rz)    !+ d2fxxzq ! d2fxdxdz
 fnode(13) = fnode(13) - dr4m3*(5.*rx*ry*ry - rx)    !+ d2fxyyq ! d2fxdydy
 fnode(14) = fnode(14) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 fnode(15) = fnode(15) - dr4m3*(5.*rx*rz*rz - rx)    !+ d2fxzzq ! d2fxdzdz
 fnode(16) = fnode(16) - dr4m3*(5.*ry*ry*ry - 3.*ry) !+ d2fyyyq ! d2fydydy
 fnode(17) = fnode(17) - dr4m3*(5.*ry*ry*rz - rz)    !+ d2fyyzq ! d2fydydz
 fnode(18) = fnode(18) - dr4m3*(5.*ry*rz*rz - ry)    !+ d2fyzzq ! d2fydzdz
 fnode(19) = fnode(19) - dr4m3*(5.*rz*rz*rz - 3.*rz) !+ d2fzzzq ! d2fzdzdz
 fnode(20) = fnode(20) - totmass*dr - 0.5*rijQij*dr3   ! potential

end subroutine compute_fnode

subroutine compute_coeff(dx,dy,dz,dr,totmass,quads,c0,c1,c2,c3)
 real, intent(in) :: dx,dy,dz,dr,totmass
 real, intent(in) :: quads(6)
 !real, intent(inout) :: coeff(4)
 real, intent(inout) :: c0,c1(3),c2(3,3),c3(3,3,3)
real :: dr3,dr4,dr5,dr6,dr3m,rx,ry,rz,qxx,qxy,qxz,qyy,qyz,qzz
 real :: dr4m3,rijQij,riQix,riQiy,riQiz,fqx,fqy,fqz
 real :: dfxdxq,dfxdyq,dfxdzq,dfydyq,dfydzq,dfzdzq
 real :: d2fxxxq,d2fxxyq,d2fxxzq,d2fxyyq,d2fxyzq
 real :: d2fxzzq,d2fyyyq,d2fyyzq,d2fyzzq,d2fzzzq

 real :: d0,d1
 real :: d2,d3
 !real :: d1arry(3)
 real :: rarry(3)
 real :: r,dr2
 integer :: i, j, k


 r = sqrt(dx**2 + dy**2 + dz**2)

 rarry = (/dx,dy,dz/)

 rx  = dx*dr
 ry  = dy*dr
 rz  = dz*dr

 ! Note dr = 1/r
 !print*, "dr"
 d0 = dr
 dr2 = dr*dr
 dr3 = dr2*dr
 dr4 = dr3*dr
 dr5 = dr4*dr
 !print*, d0
 d1 = -dr3
 !print*, 'd1'
 !print*, d1*totmass
 ! Why is this 3 not 2?????
 d2 = 3.*dr3*dr2
 !print*, d2
 d3 = -5.*d2*dr2

 dr4m3 = 3.*totmass*dr4

 !quads = 0.
 qxx = quads(1)
 qxy = quads(2)
 qxz = quads(3)
 qyy = quads(4)
 qyz = quads(5)
 qzz = quads(6)
 rijQij = (rx*rx*qxx + ry*ry*qyy + rz*rz*qzz + 2.*(rx*ry*qxy + rx*rz*qxz + ry*rz*qyz))
 riQix = (rx*qxx + ry*qxy + rz*qxz)
 riQiy = (rx*qxy + ry*qyy + rz*qyz)
 riQiz = (rx*qxz + ry*qyz + rz*qzz)
 fqx = dr4*(riQix - 2.5*rx*rijQij)
 fqy = dr4*(riQiy - 2.5*ry*rijQij)
 fqz = dr4*(riQiz - 2.5*rz*rijQij)
 dfxdxq = dr5*(qxx - 10.*rx*riQix - 2.5*rijQij   + 17.5*rx*rx*rijQij)
 dfxdyq = dr5*(qxy -  5.*ry*riQix - 5.0*rx*riQiy + 17.5*rx*ry*rijQij)
 dfxdzq = dr5*(qxz -  5.*rx*riQiz - 5.0*rz*riQix + 17.5*rx*rz*rijQij)
 dfydyq = dr5*(qyy - 10.*ry*riQiy - 2.5*rijQij   + 17.5*ry*ry*rijQij)
 dfydzq = dr5*(qyz -  5.*ry*riQiz - 5.0*rz*riQiy + 17.5*ry*rz*rijQij)
 dfzdzq = dr5*(qzz - 10.*rz*riQiz - 2.5*rijQij   + 17.5*rz*rz*rijQij)
 d2fxxxq = dr6*(-15.*qxx*rx + 105.*rx*rx*riQix - 15.*riQix - 157.5*rx*rx*rx*rijQij + 52.5*rx*rijQij)
 d2fxxyq = dr6*(35.*rx*rx*riQiy -  5.*qxx*ry - 5.*riQiy + 17.5*ry*rijQij - 157.5*rx*rx*ry*rijQij &
              + 70.*rx*ry*riQix - 10.*qxy*rx)
 d2fxxzq = dr6*(35.*rx*rx*riQiz -  5.*qxx*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*rx*rx*rz*rijQij &
              + 70.*rx*rz*riQix - 10.*qxz*rx)
 d2fxyyq = dr6*(70.*rx*ry*riQiy - 10.*qxy*ry - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*ry*ry*rijQij &
              + 35.*ry*ry*riQix -  5.*qyy*rx)
 d2fxyzq = dr6*(35.*rx*ry*riQiz -  5.*qyz*rx  &
              + 35.*ry*rz*riQix -  5.*qxz*ry  &
              + 35.*rx*rz*riQiy -  5.*qxy*rz                             - 157.5*rx*ry*rz*rijQij)
 d2fxzzq = dr6*(70.*rx*rz*riQiz - 10.*qxz*rz - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*rz*rz*rijQij &
              + 35.*rz*rz*riQix -  5.*qzz*rx)
 d2fyyyq = dr6*(-15.*qyy*ry + 105.*ry*ry*riQiy - 15.*riQiy - 157.5*ry*ry*ry*rijQij + 52.5*ry*rijQij)
 d2fyyzq = dr6*(35.*ry*ry*riQiz -  5.*qyy*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*ry*ry*rz*rijQij &
              + 70.*ry*rz*riQiy - 10.*qyz*ry)
 d2fyzzq = dr6*(70.*ry*rz*riQiz - 10.*qyz*rz - 5.*riQiy + 17.5*ry*rijQij - 157.5*ry*rz*rz*rijQij &
              + 35.*rz*rz*riQiy -  5.*qzz*ry)
 d2fzzzq = dr6*(-15.*qzz*rz + 105.*rz*rz*riQiz - 15.*riQiz - 157.5*rz*rz*rz*rijQij + 52.5*rz*rijQij)
 ! C0 = totmass * Greens function
 ! scalar 
 c0 = c0 + totmass * dr
 
 ! C1 = MB*Ri*D1
 ! Should be a vector
 !do i=1, 3
 !   c1(i) = c1(i) + totmass*rarry(i)*d1
 !enddo

 ! UNROLL LOOPS
 c1(1) = c1(1) + totmass*rarry(1)*d1 !+ fqx
 c1(2) = c1(2) + totmass*rarry(2)*d1 !+ fqy
 c1(3) = c1(3) + totmass*rarry(3)*d1 !+ fqz

 ! C2 = MB kronecker ij D1 + MB Ri Rj D2
 ! rank 2 tensor 
 !do j=1,3
 !  do i=1,3
 !     c2(i,j) = c2(i,j) + totmass*delta(i,j)*d1 +  totmass*rarry(i)*rarry(j)*d2
 !  enddo 
 !enddo 

 ! UNROLLLLLL

 c2(1,1) = c2(1,1) + totmass*d1  - totmass*3.*rx*rx*d1 !+ dfxdxq 
 c2(1,2) = c2(1,2) - totmass*3.*rx*ry*d1 !+ dfxdyq ! dfx/dy = dfy/dx
 c2(1,3) = c2(1,3) - totmass*3.*rx*rz*d1 !+ dfxdzq ! dfx/dz = dfz/dx
 c2(2,1) = c2(2,1) - totmass*3.*ry*rx*d1 !+ dfxdyq ! dfx/dy = dfy/dx
 c2(2,2) = c2(2,2) + totmass*d1  - totmass*3.*d1*ry*ry !+ dfydyq ! dfy/dy
 c2(2,3) = c2(2,3) - totmass*3.*ry*rz*d1 !+ dfydzq ! dfy/dz = dfz/dy
 c2(3,1) = c2(3,1) - totmass*3.*rz*rx*d1 !+ dfxdzq ! dfx/dz = dfz/dx
 c2(3,2) = c2(3,2) - totmass*3.*rz*ry*d1 !+ dfydzq ! dfy/dz = dfz/dy
 c2(3,3) = c2(3,3) + totmass*d1 - totmass*3.*rz*rz*d1 !+ dfzdzq 
 
 ! C3
 ! rank3 tensor
 !do k=1,3
 !   do j=1,3
 !      do i=1,3
  !       c3(i,j,k) =  c3(i,j,k) +  totmass*delta(i,j)*rarry(k)*d2 + totmass*delta(j,k)*rarry(i)*d2 &
  !        + totmass * delta(k,i)*rarry(j)*d2 + totmass*rarry(i)*rarry(j)*rarry(k)*d3
  !    enddo 
   ! enddo
 !enddo 

 ! UNROLLLLLLLLLLLLLLLLLLL

 c3(1,1,1) = c3(1,1,1) -  dr4m3*(5.*rx*rx*rx - 3.*rx) !+ d2fxxxq ! d2fxdxdx
 c3(1,1,2) = c3(1,1,2) - dr4m3*(5.*rx*rx*ry - ry)    !+ d2fxxyq ! d2fxdxdy
 c3(1,1,3) = c3(1,1,3) - dr4m3*(5.*rx*rx*rz - rz)    !+ d2fxxzq ! d2fxdxdz
 c3(1,2,1) = c3(1,2,1) - dr4m3*(5.*rx*rx*ry - ry)    !+ d2fxxyq ! d2fxdxdy
 c3(1,2,2) = c3(1,2,2) - dr4m3*(5.*rx*ry*ry - rx)    !+ d2fxyyq ! d2fxdydy
 c3(1,2,3) = c3(1,2,3) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 c3(1,3,1) = c3(1,3,1) - dr4m3*(5.*rx*rx*rz - rz)    !+ d2fxxzq ! d2fxdxdz
 c3(1,3,2) = c3(1,3,2) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 c3(1,3,3) = c3(1,3,3) - dr4m3*(5.*rx*rz*rz - rx)    !+ d2fxzzq ! d2fxdzdz
 c3(2,1,1) = c3(2,1,1) - dr4m3*(5.*rx*rx*ry - ry)    !+ d2fxxyq ! d2fxdxdy
 c3(2,1,2) = c3(2,1,2) - dr4m3*(5.*rx*ry*ry - rx)    !+ d2fxyyq ! d2fxdydy
 c3(2,1,3) = c3(2,1,3) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 c3(2,2,1) = c3(2,2,1) - dr4m3*(5.*rx*ry*ry - rx)    !+ d2fxyyq ! d2fxdydy
 c3(2,2,2) = c3(2,2,2) - dr4m3*(5.*ry*ry*ry - 3.*ry) !+ d2fyyyq ! d2fydydy
 c3(2,2,3) = c3(2,2,3) - dr4m3*(5.*ry*ry*rz - rz)    !+ d2fyyzq ! d2fydydz
 c3(2,3,1) = c3(2,3,1) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 c3(2,3,2) = c3(2,3,2) - dr4m3*(5.*ry*ry*rz - rz)    !+ d2fyyzq ! d2fydydz
 c3(2,3,3) = c3(2,3,3) - dr4m3*(5.*ry*rz*rz - ry)    !+ d2fyzzq ! d2fydzdz
 c3(3,1,1) = c3(3,1,1) - dr4m3*(5.*rx*rx*rz - rz)    !+ d2fxxzq ! d2fxdxdz
 c3(3,1,2) = c3(3,1,2) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 c3(3,1,3) = c3(3,1,3) - dr4m3*(5.*rz*rx*rz - rx)    !+ d2fxzzq ! d2fxdxdz
 c3(3,2,1) = c3(3,2,1) - dr4m3*(5.*rx*ry*rz)         !+ d2fxyzq ! d2fxdydz
 c3(3,2,2) = c3(3,2,2) - dr4m3*(5.*ry*ry*rz - rz)    !+ d2fyyzq ! d2fydydz
 c3(3,2,3) = c3(3,2,3) - dr4m3*(5.*ry*rz*rz - ry)    !+ d2fyzzq ! d2fydzdz
 c3(3,3,1) = c3(3,3,1) - dr4m3*(5.*rx*rz*rz - rx)    !+ d2fxzzq ! d2fxdzdz
 c3(3,3,2) = c3(3,3,2) - dr4m3*(5.*ry*rz*rz - ry)    !+ d2fyzzq ! d2fydzdz
 c3(3,3,3) = c3(3,3,3) - dr4m3*(5.*rz*rz*rz - 3.*rz) ! + d2fzzzq ! d2fzdzdz


 !print*, "Coeff 0:"
 !print*, c0

 !print*, "Coeff 1:"
 !print*, c1

 !print*, "Coeff 2:"
 !print*, c2

 !print*, "Coeff 3:"
 !print*, c3 

 !print*, "dr is: "
 !print*,"Components are: "
 !print*, totmass*delta(2,1)*rarry(1)*d2
 !print*, totmass*delta(1,1)*rarry(2)*d2
 !print*, totmass*delta(1,2)*rarry(1)*d2
 !print*, dr


end subroutine compute_coeff

!----------------------------------------------------------------
!+
!  Internal subroutine to compute the Taylor-series expansion
!  of the gravitational force, given the force acting on the
!  centre of the node and its derivatives
!
! INPUT:
!   fnode: array containing force on node due to distant nodes
!          and first derivatives of f (i.e. Jacobian matrix)
!          and second derivatives of f (i.e. Hessian matrix)
!   dx,dy,dz: offset of the particle from the node centre of mass
!
! OUTPUT:
!   fxi,fyi,fzi : gravitational force at the new position
!+
!----------------------------------------------------------------
pure subroutine expand_fgrav_in_taylor_series(fnode,dx,dy,dz,fxi,fyi,fzi,poti,dfxxi,dfxyi,dfxzi,dfyyi,dfyzi,dfzzi)
 real, intent(in)  :: fnode(lenfgrav)
 real, intent(in)  :: dx,dy,dz
 real, intent(out) :: fxi,fyi,fzi,poti
 real, intent(out), optional              :: dfxxi,dfxyi,dfxzi,dfyyi,dfyzi,dfzzi
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 fxi = fnode(1)
 fyi = fnode(2)
 fzi = fnode(3)
 dfxx = fnode(4)
 dfxy = fnode(5)
 dfxz = fnode(6)
 dfyy = fnode(7)
 dfyz = fnode(8)
 dfzz = fnode(9)

 

 d2fxxx = fnode(10)
 d2fxxy = fnode(11)
 d2fxxz = fnode(12)
 d2fxyy = fnode(13)
 d2fxyz = fnode(14)
 d2fxzz = fnode(15)
 d2fyyy = fnode(16)
 d2fyyz = fnode(17)
 d2fyzz = fnode(18)
 d2fzzz = fnode(19)
 poti = fnode(20)

 fxi = fxi + dx*(dfxx   + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
           + dy*(dfxy  + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dz*(dfxz  + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
 fyi = fyi + dx*(dfxy  + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dy*(dfyy  + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
           + dz*(dfyz  + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
 fzi = fzi + dx*(dfxz  + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
           + dy*(dfyz  + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
           + dz*(dfzz  + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))

 ! if (present(dfxxi)) then 
 ! ! The c2 field tensor to be translated 
 !  !print*, "translated"
 !  dfxxi = fnode(4)
 !  dfxyi = fnode(5)
 !  dfxzi = fnode(6)
 !  dfyyi = fnode(7)
 !  dfyzi = fnode(8)
 !  dfzzi = fnode(9)
  ! dfxxi = dfxxi + dx*(d2fxxx + d2fxxy + d2fxxz)
  ! dfxyi = dfxyi + dy*(d2fxxy + d2fxyy + d2fxyz)
  ! dfxzi = dfxzi + dz*(d2fxxz + d2fxyz + d2fxzz)
  ! dfyyi = dfyyi + dy*(d2fxyy + d2fyyy + d2fyyz)
  ! dfyzi = dfyzi + dz*(d2fxyz + d2fyyz + d2fyzz)
  ! dfzzi = dfzzi + dz*(d2fxzz + d2fyzz + d2fzzz)
 ! endif 
 poti = poti - (dx*fxi + dy*fyi + dz*fzi)

 return
end subroutine expand_fgrav_in_taylor_series

pure subroutine translate_fgrav_in_taylor_series(fnode,dx,dy,dz)
 real, intent(inout)  :: fnode(lenfgrav)
 real, intent(in)  :: dx,dy,dz
 real :: fnodeold(lenfgrav)
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 ! fxi = fnode(1)
 ! fyi = fnode(2)
 ! fzi = fnode(3)
 dfxx = fnode(4)
 dfxy = fnode(5)
 dfxz = fnode(6)
 dfyy = fnode(7)
 dfyz = fnode(8)
 dfzz = fnode(9)

 

 d2fxxx = fnode(10)
 d2fxxy = fnode(11)
 d2fxxz = fnode(12)
 d2fxyy = fnode(13)
 d2fxyz = fnode(14)
 d2fxzz = fnode(15)
 d2fyyy = fnode(16)
 d2fyyz = fnode(17)
 d2fyzz = fnode(18)
 d2fzzz = fnode(19)
 ! poti = fnode(20)

 fnodeold = fnode 

 fnode(1) = fnodeold(1) + dx*(dfxx  + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
           + dy*(dfxy  + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dz*(dfxz    + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
 fnode(2) = fnodeold(2) + dx*(dfxy   + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dy*(dfyy   + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
           + dz*(dfyz   + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
 fnode(3) = fnodeold(3) + dx*(dfxz   + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
           + dy*(dfyz   + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
           + dz*(dfzz   + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))

 ! ! The c2 field tensor to be translated 

  fnode(4) = fnodeold(4) + dx*(d2fxxx + d2fxxy + d2fxxz)
  fnode(5) = fnodeold(5) + dy*(d2fxxy + d2fxyy + d2fxyz)
  fnode(6) = fnodeold(6) + dz*(d2fxxz + d2fxyz + d2fxzz)
  fnode(7) = fnodeold(7) + dy*(d2fxyy + d2fyyy + d2fyyz)
  fnode(8) = fnodeold(8) + dz*(d2fxyz + d2fyyz + d2fyzz)
  fnode(9) = fnodeold(9) + dz*(d2fxzz + d2fyzz + d2fzzz)
 ! endif 
   fnode(10:20) = fnodeold(10:20)
  !poti = poti - (dx*fxi + dy*fyi + dz*fzi)

 return
end subroutine translate_fgrav_in_taylor_series

subroutine translate_expansion_center(z0,z1,c0,c1,c2,c3)
  real, intent(in) :: z0(3), z1(3)
  real, intent(inout) :: c0, c1(3), c2(3,3), c3(3,3,3)
  real :: c0old, c1old(3), c2old(3,3), c3old(3,3,3)
  real :: sep1(3), sep2(3,3), sep3(3,3,3)
  real :: sep1c2comp(3), sep2c3comp(3), sep1c3comp(3,3)
  
  c0old = c0
  c1old = c1
  c2old = c2
  c3old = c3
  !print*, "c0"
  !print*, c0old
  !print*, "c1"
  !print*, c1old
  !print*, "c2"
  !print*, c2old
  !print*, "c3"
  !print*, c3old(2,1,1)
  sep1 = 0.
  sep1 = z1 - z0
  !print*, "sep1: ", sep1
  !print*, c2 

  sep1c2comp(1) = sep1(1)*c2(1,1) + sep1(2)*c2(1,2) + sep1(3)*c2(1,3)
  sep1c2comp(2) = sep1(1)*c2(2,1) + sep1(2)*c2(2,2) + sep1(3)*c2(2,3)
  sep1c2comp(3) = sep1(1)*c2(3,1) + sep1(2)*c2(3,2) + sep1(3)*c2(3,3)


  sep2c3comp(1) = sep1(1)*(sep1(1)*c3(1,1,1) + sep1(2)*c3(1,1,2) + sep1(3)*c3(1,1,3)) &
                + sep1(2)*(sep1(1)*c3(1,2,1)+ sep1(2)*c3(1,2,2) + sep1(3)*c3(1,2,3)) &
                + sep1(3)*(sep1(1)*c3(1,3,1)+ sep1(2)*c3(1,3,2) + sep1(3)*c3(1,3,3))

  sep2c3comp(2) = sep1(1)*(sep1(1)*c3(2,1,1) + sep1(2)*c3(2,1,2) + sep1(3)*c3(2,1,3)) &
                + sep1(2)*(sep1(1)*c3(2,2,1)+ sep1(2)*c3(2,2,2) + sep1(3)*c3(2,2,3)) &
                + sep1(3)*(sep1(1)*c3(2,3,1)+ sep1(2)*c3(2,3,2) + sep1(3)*c3(2,3,3))

  sep2c3comp(3) = sep1(1)*(sep1(1)*c3(3,1,1) + sep1(2)*c3(3,1,2) + sep1(3)*c3(3,1,3)) &
                + sep1(2)*(sep1(1)*c3(3,2,1)+ sep1(2)*c3(3,2,2) + sep1(3)*c3(3,2,3)) &
                + sep1(3)*(sep1(1)*c3(3,3,1)+ sep1(2)*c3(3,3,2) + sep1(3)*c3(3,3,3))


  sep1c3comp(1,1) = sep1(1)*c3(1,1,1) + sep1(2)*c3(1,1,2) + sep1(3)*c3(1,1,3)
  sep1c3comp(1,2) = sep1(1)*c3(1,2,1) + sep1(2)*c3(1,2,2) + sep1(3)*c3(1,2,3)
  sep1c3comp(1,3) = sep1(1)*c3(1,3,1) + sep1(2)*c3(1,3,2) + sep1(3)*c3(1,3,3)
  sep1c3comp(2,1) = sep1(1)*c3(2,1,1) + sep1(2)*c3(2,1,2) + sep1(3)*c3(2,1,3)
  sep1c3comp(2,2) = sep1(1)*c3(2,2,1) + sep1(2)*c3(2,2,2) + sep1(3)*c3(2,2,3)
  sep1c3comp(2,3) = sep1(1)*c3(2,3,1) + sep1(2)*c3(2,3,2) + sep1(3)*c3(2,3,3)
  sep1c3comp(3,1) = sep1(1)*c3(3,1,1) + sep1(2)*c3(3,1,2) + sep1(3)*c3(3,1,3)
  sep1c3comp(3,2) = sep1(1)*c3(3,2,1) + sep1(2)*c3(3,2,2) + sep1(3)*c3(3,2,3)
  sep1c3comp(3,3) = sep1(1)*c3(3,3,1) + sep1(2)*c3(3,3,2) + sep1(3)*c3(3,3,3)

  !print*, "Second term value: "
  !print*, 0.5*inner_product2(sep2,c2old)
  !print*, 0.5*dot_product(sep2,c2old)
  !print*, "C2 translated"
  !print*, sep1c2comp

  ! The components of these sums should all have the save order as the coefficent i.e  c0 = scalar, c1 = vector
  c0 = c0old + dot_product(c1old,sep1) !+ 0.5*inner_product2(sep2,c2old) + 1./6.*inner_product3(sep3,c3)
  c1 = c1old + sep1c2comp !+ 0.5*sep2c3comp
  !print*, "c1 is: ",c1
  c2 = c2old !+ sep1c3comp !+ 0.5*inner_product31_to_2(sep1,c3old)
  c3 = c3old

 end subroutine translate_expansion_center

!-----------------------------------------------
!+
!  Routine to update a constructed tree
!  Note: current version ONLY works if
!  tree is built to < maxlevel_indexed
!  That is, it relies on the 2^n style tree
!  indexing to sweep out each level
!+
!-----------------------------------------------
subroutine revtree(node, xyzh, ifirstincell, ncells)
 use dim,  only:maxp
 use part, only:maxphase,iphase,igas,massoftype,iamtype
 use io,   only:fatal
 type(kdnode), intent(inout) :: node(:) !ncellsmax+1)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: ifirstincell(:) !ncellsmax+1)
 integer(kind=8), intent(in) :: ncells
 real :: hmax, r2max
 real :: xi, yi, zi, hi
 real :: dx, dy, dz, dr2
#ifdef GRAVITY
 real :: quads(6)
#endif
 integer :: icell, i, level, il, ir
 real :: pmassi, totmass
 real :: x0(3)
 type(kdnode) :: nodel,noder

 pmassi = massoftype(igas)
 do i=1,int(ncells)
#ifdef GRAVITY
    ! cannot update centre of node without gravity
    ! as it is not at the centre of mass
    node(i)%xcen(:) = 0.
#endif
    node(i)%size    = 0.
    node(i)%hmax    = 0.
#ifdef GRAVITY
    node(i)%mass    = 0.
    node(i)%quads(:)= 0.
#endif
 enddo

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(xyzh, ifirstincell, ncells) &
!$omp shared(node, ll, iphase, massoftype, maxlevel) &
!$omp private(hmax, r2max, xi, yi, zi, hi, il, ir, nodel, noder) &
!$omp private(dx, dy, dz, dr2, icell, i, x0) &
#ifdef GRAVITY
!$omp private(quads) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)
!$omp do schedule(guided, 2)
 over_cells: do icell=1,int(ncells)

    i = abs(ifirstincell(icell))
    if (i==0) cycle over_cells
    print*, icell

    ! find centre of mass
    ! this becomes the new node center
    x0 = 0.
    totmass = 0.0
    calc_cofm: do while (i /= 0)
       print*, "currentpart: ", i 
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase(i)))
       endif
       x0(1) = x0(1) + pmassi*xi
       x0(2) = x0(2) + pmassi*yi
       x0(3) = x0(3) + pmassi*zi
       totmass = totmass + pmassi

       i = abs(ll(i))
    enddo calc_cofm

    x0 = x0/totmass
#ifdef GRAVITY
    node(icell)%xcen(1) = x0(1)
    node(icell)%xcen(2) = x0(2)
    node(icell)%xcen(3) = x0(3)
    print*, "xcen = ", node(icell) % xcen 
#endif

    i = abs(ifirstincell(icell))

    ! update cell size, hmax
    r2max = 0.
    hmax = 0.
#ifdef GRAVITY
    quads = 0.
#endif
    over_parts: do while (i /= 0)
       print*, i 
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       hi = xyzh(4,i)
       dx = xi - node(icell)%xcen(1)
       dy = yi - node(icell)%xcen(2)
       dz = zi - node(icell)%xcen(3)
       dr2 = dx*dx + dy*dy + dz*dz
       r2max = max(dr2, r2max)
       print*, "r2max: ", r2max
       hmax  = max(hi, hmax)
#ifdef GRAVITY
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase(i)))
       endif
       quads(1) = quads(1) + pmassi*(3.*dx*dx - dr2)
       quads(2) = quads(2) + pmassi*(3.*dx*dy)
       quads(3) = quads(3) + pmassi*(3.*dx*dz)
       quads(4) = quads(4) + pmassi*(3.*dy*dy - dr2)
       quads(5) = quads(5) + pmassi*(3.*dy*dz)
       quads(6) = quads(6) + pmassi*(3.*dz*dz - dr2)
#endif
       ! move to next particle in list
       i = abs(ll(i))
       print*, "next part is: ", i 
    enddo over_parts

    node(icell)%size = sqrt(r2max) + epsilon(r2max)
    print*, "node size: ",node(icell)%size
    node(icell)%hmax = hmax
#ifdef GRAVITY
    node(icell)%mass = totmass
    node(icell)%quads = quads
#endif
 enddo over_cells
!$omp enddo
!
! propagate information to parent nodes
! here we sweep across each level at a time
! and update each node from its two children
!
 do level=maxlevel-1,0,-1
!$omp do
    do i=2**level,2**(level+1)-1
       ! get child nodes
       il = node(i)%leftchild
       ir = node(i)%rightchild
       if (il > 0 .and. ir > 0) then
          nodel = node(il)
          noder = node(ir)

          ! If leaf set maxsize = size 
          ! if (ifirstincell(il) > 0) then 
          !    node(il) % maxsize = node(il) % size 
          ! endif 
          ! if (ifirstincell(ir) > 0) then 
          !    node(ir) % maxsize = node(ir) % size 
          ! endif 
          call add_child_nodes(nodel,noder,node(i))
       else
          if (il > 0 .or. ir > 0) then
             ! should never happen, should have two children or none
             call fatal('revtree','node with only one child during tree revision',var='ir',ival=ir)
          endif
       endif
    enddo
!$omp enddo
 enddo
!$omp end parallel

end subroutine revtree

!-----------------------------------------------------------------
!+
!  Update parent node from the properties of the two child nodes
!  IN:
!    l, r - two child nodes
!  OUT:
!    nodei - updated parent node
!+
!-----------------------------------------------------------------
subroutine add_child_nodes(l,r,nodei)
 type(kdnode), intent(in)  :: l,r
 type(kdnode), intent(out) :: nodei
 real :: xl(3),sl,hl
 real :: xr(3),sr,hr
#ifdef GRAVITY
 real :: ql(6),qr(6),mr,ml,mnode,dm,distl,distr
#endif
 real :: dx,dy,dz,dr2,dr

 xl = l%xcen
 hl = l%hmax
 sl = l%maxsize
#ifdef GRAVITY
 ml = l%mass
 ql = l%quads
#endif

 xr = r%xcen
 hr = r%hmax
 sr = r%maxsize
#ifdef GRAVITY
 mr = r%mass
 qr = r%quads
 mnode = ml + mr
 dm    = 1./mnode
#endif
 dx = xl(1) - xr(1)
 dy = xl(2) - xr(2)
 dz = xl(3) - xr(3)
 dr2 = dx*dx + dy*dy + dz*dz
 dr  = sqrt(dr2)
#ifdef GRAVITY
 ! centre of mass
 ! nodei%xcen = (xl*ml + xr*mr)*dm
 ! !  ! distance between left child and node centre
 ! dx = xl(1) - nodei%xcen(1)
 ! dy = xl(2) - nodei%xcen(2)
 ! dz = xl(3) - nodei%xcen(3)
 ! dr = sqrt(dx*dx + dy*dy + dz*dz)
 ! nodei%maxsize = dr+sl
 ! ! distance between right child and node centre
 ! dx = xr(1) - nodei%xcen(1)
 ! dy = xr(2) - nodei%xcen(2)
 ! dz = xr(3) - nodei%xcen(3)
 ! dr = sqrt(dx*dx + dy*dy + dz*dz)
 ! nodei%maxsize = max(nodei%size,dr+sr)
 ! size, formula as in Benz et al. 1990
 ! and from thinking about it...
 !nodei%size = max(ml*dm*dr+sr,mr*dm*dr+sl)
 nodei%maxsize = max(ml*dm*dr+sr,mr*dm*dr+sl) ! size is replaced with max size 
 ! max node size is maximum distance from center of mass
#else
 ! distance between left child and node centre
 dx = xl(1) - nodei%xcen(1)
 dy = xl(2) - nodei%xcen(2)
 dz = xl(3) - nodei%xcen(3)
 dr = sqrt(dx*dx + dy*dy + dz*dz)
 nodei%size = dr+sl
 ! distance between right child and node centre
 dx = xr(1) - nodei%xcen(1)
 dy = xr(2) - nodei%xcen(2)
 dz = xr(3) - nodei%xcen(3)
 dr = sqrt(dx*dx + dy*dy + dz*dz)
 nodei%size = max(nodei%size,dr+sr)
#endif
 nodei%hmax = max(hl,hr)
#ifdef GRAVITY
 nodei%mass = mnode
 ! quadrupole moments, see Benz et al. (1990), this is also
 ! the parallel axis theorem
 nodei%quads(1) = ql(1) + qr(1) + ml*mr*(3.*dx*dx - dr2)*dm
 nodei%quads(2) = ql(2) + qr(2) + ml*mr*(3.*dx*dy)*dm
 nodei%quads(3) = ql(3) + qr(3) + ml*mr*(3.*dx*dz)*dm
 nodei%quads(4) = ql(4) + qr(4) + ml*mr*(3.*dy*dy - dr2)*dm
 nodei%quads(5) = ql(5) + qr(5) + ml*mr*(3.*dy*dz)*dm
 nodei%quads(6) = ql(6) + qr(6) + ml*mr*(3.*dz*dz - dr2)*dm
#endif

end subroutine add_child_nodes

!--------------------------------------------------------------------------------
!+
!  Routine to build the global level tree
!+
!-------------------------------------------------------------------------------
#ifdef MPI
subroutine maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,np,ndim,cellatid,ifirstincell,ncells)
 use io,           only:fatal,warning,id,nprocs
 use mpiutils,     only:reduceall_mpi
 use balance,      only:balancedomains
 use mpiderivs,    only:tree_sync,tree_bcast
 use part,         only:isdead_or_accreted,iactive,ibelong

 type(kdnode), intent(out)     :: nodeglobal(:) !ncellsmax+1)
 type(kdnode), intent(out)     :: node(:) !ncellsmax+1)
 integer,      intent(out)     :: nodemap(:) !ncellsmax+1)
 integer,      intent(out)     :: globallevel
 integer,      intent(out)     :: refinelevels
 integer,      intent(inout)   :: np
 integer,      intent(in)      :: ndim
 real,         intent(inout)   :: xyzh(:,:)
 integer,      intent(out)     :: cellatid(:) !ncellsmax+1)
 integer,      intent(out)     :: ifirstincell(:) !ncellsmax+1)
 real                          :: xmini(ndim),xmaxi(ndim)
 real                          :: xminl(ndim),xmaxl(ndim)
 real                          :: xminr(ndim),xmaxr(ndim)

 integer                       :: minlevel, maxlevel

 integer                       :: idleft, idright
 integer                       :: groupsize,ifirstingroup,groupsplit

 integer(kind=8), intent(out)  :: ncells

 type(kdnode)                  :: mynode(1)

 integer                       :: nl, nr
 integer                       :: il, ir, iself, parent
 integer                       :: level
 integer                       :: nnodestart, nnodeend,locstart,locend
 integer                       :: npcounter

 integer                       :: i, k, offset, roffset, roffset_prev, coffset
 integer                       :: inode
 integer                       :: npnode

 logical                       :: wassplit

 irootnode = 1
 parent = 0
 iself = irootnode
 ifirstincell = 0

 ! root is level 0
 globallevel = int(ceiling(log(real(nprocs)) / log(2.0)))

 minlevel = 31
 maxlevel = 0

 levels: do level = 0, globallevel
    groupsize = 2**(globallevel - level)
    ifirstingroup = (id / groupsize) * groupsize
    if (level == 0) then
       call construct_root_node(np,npcounter,irootnode,ndim,xmini,xmaxi,ifirstincell,xyzh)
       dxi = xmaxi-xmini
    else
       npcounter = np
    endif

    call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
            il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, &
            ncells, ifirstincell, minlevel, maxlevel, ndim, xyzh, wassplit, list, &
            groupsize)

    if (.not.wassplit) then
       call fatal('maketreeglobal','insufficient particles for splitting at the global level: '// &
            'use more particles or less MPI threads')
    endif

    ! set which tree child this proc will belong to next
    groupsplit = ifirstingroup + (groupsize / 2)

    ! record parent for next round
    parent = iself

    ! which half of the tree this proc is on
    if (id < groupsplit) then
       ! i for the next node we construct
       iself = il
       ! the left and right procIDs
       idleft = id
       idright = id + 2**(globallevel - level - 1)
       xmini = xminl
       xmaxi = xmaxl
    else
       iself = ir
       idleft = id - 2**(globallevel - level - 1)
       idright = id
       xmini = xminr
       xmaxi = xmaxr
    endif

    if (np > 0) then
       do i = inoderange(1,il), inoderange(2,il)
          ibelong(abs(inodeparts(i))) = idleft
       enddo
       do i = inoderange(1,ir), inoderange(2,ir)
          ibelong(abs(inodeparts(i))) = idright
       enddo
    endif

    ! move particles to where they belong
    call balancedomains(np)

    ! move particles from old array
    ! this is a waste of time, but maintains compatibility
    npnode = 0
    do i=1,np
       npnode = npnode + 1
#ifdef IND_TIMESTEPS
       if (iactive(iphase(i))) then
          inodeparts(npnode) = i
       else
          inodeparts(npnode) = -i
       endif
#else
       inodeparts(npnode) = i
#endif
       xyzh_soa(npnode,:) = xyzh(:,i)
       iphase_soa(npnode) = iphase(i)
    enddo

    ! set all particles to belong to this node
    inoderange(1,iself) = 1
    inoderange(2,iself) = np

    ! range of newly written tree
    nnodestart = 2**level
    nnodeend = 2**(level + 1) - 1

    ! synchronize tree with other owners if this proc is the first in group
    call tree_sync(mynode, 1,nodeglobal(nnodestart:nnodeend), ifirstingroup, groupsize, level)

    ! at level 0, tree_sync already 'broadcasts'
    if (level > 0) then
       ! tree broadcast to non-owners
       call tree_bcast(nodeglobal(nnodestart:nnodeend), nnodeend - nnodestart + 1, level)
    endif

 enddo levels

 ! local tree
 call maketree(node,xyzh,np,ndim,ifirstincell,ncells,refinelevels)

 ! tree refinement
 refinelevels = int(reduceall_mpi('min',refinelevels),kind=kind(refinelevels))
 roffset_prev = 1

 do i = 1,refinelevels
    offset = 2**(globallevel + i)
    roffset = 2**i

    nnodestart = offset
    nnodeend   = 2*nnodestart-1

    locstart   = roffset
    locend     = 2*locstart-1

    ! index shift the node to the global level
    do k = roffset,2*roffset-1
       refinementnode(k) = node(k)
       coffset = refinementnode(k)%parent - roffset_prev

       refinementnode(k)%parent = 2**(globallevel + i - 1) + id * roffset_prev + coffset

       if (i /= refinelevels) then
          refinementnode(k)%leftchild  = 2**(globallevel + i + 1) + 2*id*roffset + 2*(k - roffset)
          refinementnode(k)%rightchild = refinementnode(k)%leftchild + 1
       else
          refinementnode(k)%leftchild = 0
          refinementnode(k)%rightchild = 0
       endif
    enddo

    roffset_prev = roffset
    ! sync, replacing level with globallevel, since all procs will get synced
    ! and deeper comms do not exist
    call tree_sync(refinementnode(locstart:locend),roffset,nodeglobal(nnodestart:nnodeend),id,1,globallevel)

    ! get the mapping from the local tree to the global tree, for future hmax updates
    do inode = locstart,locend
       nodemap(inode) = nnodestart + (id * roffset) + (inode - locstart)
    enddo
 enddo

 ! cellatid is zero by default
 cellatid = 0
 do i = 1,nprocs
    offset = 2**(globallevel+refinelevels)
    roffset = 2**refinelevels
    do k = 1,roffset
       cellatid(offset + (i - 1) * roffset + (k - 1)) = i
    enddo
 enddo

end subroutine maketreeglobal

#endif 

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
  stack(top) % nodeindex1 = 0 
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


 subroutine get_node_data(currentnode,interactnode,dx,dy,dz,xsizei,xsizej,fnode,totmass_node,quads)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
  type(kdnode), intent(in) :: currentnode, interactnode
  real, intent(out) :: dx,dy,dz,xsizej,xsizei
  real, optional, intent(out) :: totmass_node, quads(6),fnode(20)
  real :: xoffset, yoffset, zoffset
#ifdef PERIODIC
  real :: hdlx,hdly,hdlz

  hdlx = 0.5*dxbound
  hdly = 0.5*dybound
  hdlz = 0.5*dzbound
#endif
    dx = currentnode%xcen(1) - interactnode%xcen(1)      ! distance between node centres
    dy = currentnode%xcen(2) - interactnode%xcen(2)
#ifndef TREEVIZ
    dz = currentnode%xcen(3) - interactnode%xcen(3)
#endif
    xsizej       = interactnode%size
    xsizei       = currentnode%size 
    !hmax         = node(n)%hmax
    !il           = node(n)%leftchild
    !ir           = node(n)%rightchild
#ifdef GRAVITY
    if (present(fnode)) then 
      totmass_node =  interactnode%mass
      quads        =  interactnode%quads
      fnode        = currentnode % fnode
      !fnode = 0. 
    endif 
#endif
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
end subroutine get_node_data

logical function well_separated(node1,node2) result(bool)
  use dtypekdtree, only:kdnode
  type(kdnode), intent(in) :: node1, node2
  real :: cm1(3), cm2(3), dx(3), zmag
  real :: theta,rmax1,rmax2,r2
  !logical :: bool

  !bool =.true.
  !return 
  theta = 1.0

  cm1 = node1 % xcen
  cm2 = node2 % xcen 

  ! Get the magnitude of the difference between CoM's
  dx = cm1 - cm2
  !print*,dx
  zmag = dx(1)*dx(1) + dx(2)*dx(2)+ dx(3)*dx(3)

  ! get the rmax for nodes 1 and 2
  ! Change to sizemax 
  rmax1 = node1 % maxsize
  rmax2 = node2 % maxsize

  ! square them 
  r2 = (rmax1 + rmax2)**2 
  !call distance_to_corner(node1,rmax1)
  !call distance_to_corner(node2,rmax2)


  !print*, "Rmax values: "
  !print*, rmax1, rmax2

   !print*,"Zmag: ", zmag
   !print*, "rmax1 + rmax2/theta: ", (rmax1 + rmax2)/theta
   if (theta == 0.0) then 
     bool = .FALSE.
     return 
   endif 
  if (zmag > (r2)/(theta*theta)) then
    bool = .TRUE.
  else 
    bool = .FALSE.
  endif  



end function well_separated

 subroutine direct_sum_not_neighbour(currentnodeindex, interactnodeindex,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
  !$ use omputils 
  use part, only :massoftype,maxphase,iamtype,iphase
  use dim,  only :maxp
  integer, intent(in) :: currentnodeindex, interactnodeindex, ifirstincell(:)
  type(kdnode), intent(in) :: node(:)
  real, intent(inout) :: fxyzu_dtt(:,:), poten_dtt(:)
  real, intent(in)    :: xyzh(:,:)
  integer :: currentnodepart(10), interactnodepart(10)
  integer :: npcurrent, npinteract 
  integer :: i, j, indexi,indexj
  integer :: iamtypej
  real    :: dx,dy,dz,rij2,rij1,fgrav,fgravj,pmassj,phii

  currentnodepart = 0.
  interactnodepart = 0.
  npcurrent = 0.
  npinteract = 0.
  dx = 0 
  dy = 0
  dz = 0
  rij2 = 0.
  rij1 = 0.
  fgrav = 0.
  fgravj = 0.
  pmassj = 0.
  iamtypej = 0
  i = 0
  j = 0
  indexi = 0
  indexj = 0 
  phii = 0.

  !print*, "enetered direct"
  ! find particles in each cell 
  ! DEBUG THIS SUBROUTINE 
  call get_part_node(currentnodeindex,node,ifirstincell,currentnodepart,npcurrent)
  call get_part_node(interactnodeindex,node,ifirstincell,interactnodepart,npinteract)

   !print*, "Got particles: ", currentnodepart(1:npcurrent), interactnodepart(1:npinteract)
  ! ! compute direct sum 
  ! print*, "node index: ", currentnodeindex, interactnodeindex
  !print*, npcurrent, npinteract
  ! FIX THIS 
  do i=1, npcurrent
    indexi = currentnodepart(i)
    !print*, "indexi: ", indexi
    do j=1, npinteract
      indexj = interactnodepart(j)
      !print*, "indexj: ", indexj
       if (indexj == 1) then 
     !print*, interactnodeindex
     !stop 
    endif 
      if (indexi /= indexj) then 
        !print*, "calc"
        !if (indexi == 94) print*, "indexj: ", indexj, "currentnodeindex: ", currentnodeindex, &
        !"interactnodeindex: ", interactnodeindex, &
        ! "ifirstincell",ifirstincell(interactnodeindex)
      dx = xyzh(1,indexi) - xyzh(1,indexj)
      dy = xyzh(2,indexi) - xyzh(2,indexj)
      dz = xyzh(3,indexi) - xyzh(3,indexj)
      ! find rij2 by squaring distance of particles 
      rij2 = dx*dx + dy*dy + dz*dz 
      ! print*, indexi, indexj
      ! !print*,"rij2:", rij2
#ifdef FINVSQRT
       rij1 = finvsqrt(rij2)
#else
       rij1 = 1./sqrt(rij2)
#endif
       fgrav  = rij1*rij1*rij1
       if (maxphase==maxp) then
          iamtypej = iamtype(iphase(indexj))
       endif
       pmassj = massoftype(iamtypej)
       phii   = -rij1
       fgravj = fgrav*pmassj
       ! print*, fgravj
       ! print*, indexi
       ! ! Instead of fsum append directly to fxyzu
       ! print*, fxyzu(1,indexi)
       !if (isnan(fxyzu(1,indexi))) stop 
       ! should be omp lock here 
       
       !!$omp critical (direct)
       !!$omp ATOMIC UPDATE
       !$ call omp_set_lock(ipart_omp_lock(indexi/10))
       fxyzu_dtt(1,indexi) = fxyzu_dtt(1, indexi) - dx*fgravj
       !!$omp ATOMIC UPDATE 
       fxyzu_dtt(2, indexi) = fxyzu_dtt(2, indexi) - dy*fgravj
       !!$omp ATOMIC UPDATE
       fxyzu_dtt(3, indexi) = fxyzu_dtt(3, indexi) - dz*fgravj
       !!$omp ATOMIC UPDATE 
       poten_dtt(indexi) = poten_dtt(indexi) + pmassj*phii
       !$ call omp_unset_lock(ipart_omp_lock(indexi/10))
       !print*, "poten: ", poten_dtt(indexi)
       !!$omp end critical (direct)
     endif 

    end do 
  end do 


 end subroutine direct_sum_not_neighbour

 subroutine direct_sum_softened(currentnodeindex, interactnodeindex,fxyzu_dtt,poten_dtt,xyzh,node,ifirstincell)
  use part, only :massoftype,maxphase,iamtype,iphase
  use dim,  only :maxp
  use kernel,    only:grkern,kernel_softening,radkern2,cnormk
  use part,      only:igas,iamtype,maxphase,maxp,iphase, &
                     iactive,isdead_or_accreted,massoftype,maxgradh,gradh
  use io,        only:error
  integer, intent(in) :: currentnodeindex, interactnodeindex, ifirstincell(:)
  type(kdnode), intent(in) :: node(:)
  real, intent(inout) :: fxyzu_dtt(:,:), poten_dtt(:)
  real, intent(in)    :: xyzh(:,:)
  !real(kind=4), intent(in)    :: gradh(:,:)
  integer :: currentnodepart(500), interactnodepart(500)
  integer :: npcurrent, npinteract,start 
  integer :: i,j,indexi,indexj,iamtypei,iamtypej
  real :: dx(3),dr(3),fgravi(3),fgravj(3),xi(3)
  real :: rij,rij1,rij2,rij21,pmassi,pmassj
  real :: gradhi,gradsofti,grkerni,grkernj,dsofti,dsoftj
  real :: phii,phij,phiterm,fm,fmi,fmj,phitemp,potensoft0,qi,qj
  real :: hi,hj,hi1,hj1,hi21,hj21,hi41,hj41,q2i,q2j
  logical :: iactivei,iactivej

  currentnodepart = 0
  interactnodepart = 0
  npcurrent = 0
  npinteract = 0
  dx = 0 
  ! dy = 0
  ! dz = 0
  rij2 = 0.
  rij1 = 0.
  fgravi = 0.
  fgravj = 0.
  pmassj = 0.
  iamtypej = 0
  i = 0
  j = 0
  indexi = 0
  indexj = 0 
  phii = 0.

  call get_part_node(currentnodeindex,node,ifirstincell,currentnodepart,npcurrent)
  call get_part_node(interactnodeindex,node,ifirstincell,interactnodepart,npinteract)



  !
!--reset potential (but not force) initially
!
! phi = 0.
 !phitot = 0.

 iactivei = .true.
 iactivej = .true.
 iamtypei = igas
 iamtypej = igas
 pmassi = massoftype(iamtypei)
 pmassj = massoftype(iamtypej)

 call kernel_softening(0.,0.,potensoft0,fmi)
 if (size(gradh(:,1)) < 2) then
    call error('directsum','cannot do direct sum with ngradh < 2')
    return
 endif
 if (maxgradh /= maxp) then
    call error('directsum','gradh not stored (maxgradh /= maxp in part.F90)')
    return
 endif
!
!--calculate gravitational force by direct summation on all particles
!
 start = 1 
 ! print*, "npcurrent: ", npcurrent
 ! print*,"npinteract: ", npinteract
 overi: do i=start,npcurrent
    ! get particle index i
    indexi = currentnodepart(indexi)
    xi(1:3) = xyzh(1:3,indexi)
    hi      = xyzh(4,indexi)
    !print*, "doesnt pass cycle "
    !print*, "dead or accreted: ",isdead_or_accreted(hi)
    !if (isdead_or_accreted(hi)) cycle overi
    !print*, "passes cycle"
    if (maxphase==maxp) then
       iamtypei = iamtype(iphase(indexi))
       iactivei = iactive(iphase(indexi))
       pmassi = massoftype(iamtypei)
    endif

    hi1  = 1./hi
    hi21 = hi1*hi1
    hi41 = hi21*hi21
    gradhi    = gradh(1,indexi)
    gradsofti = gradh(2,indexi)
    fgravi(:) = 0.
    phitemp   = 0.

    overj: do j=start,npinteract
       ! get particle index j 
       indexj = currentnodepart(j)
       !print*,"indexi, indexj: ", indexi, indexj
       dx(1) = xi(1) - xyzh(1,indexj)
       dx(2) = xi(2) - xyzh(2,indexj)
       dx(3) = xi(3) - xyzh(3,indexj)
       hj    = xyzh(4,indexj)
       if (isdead_or_accreted(hj)) cycle overj
       hj1   = 1./hj
       rij2  = dot_product(dx,dx)
       rij   = sqrt(rij2)
       rij1  = 1./rij
       rij21 = rij1*rij1
       dr(:) = dx(:)*rij1
       hj21  = hj1*hj1
       hj41  = hj21*hj21
       if (maxphase==maxp) then
          iamtypej = iamtype(iphase(indexj))
          iactivej = iactive(iphase(indexj))
          pmassj = massoftype(iamtypej)
       endif
       fgravj(:) = 0.
       q2i = rij2*hi21
       q2j = rij2*hj21
       if (q2i < radkern2) then
          qi = sqrt(q2i)
          grkerni = grkern(q2i,qi)
          call kernel_softening(q2i,qi,phii,fmi)
          phii       = phii*hi1
          fmi        = fmi*hi21
          grkerni    = cnormk*grkerni*hi41*gradhi
          dsofti     = 0.5*grkerni*gradsofti
          fgravi(:)  = fgravi(:) - pmassj*dsofti*dr(:)
          fgravj(:)  = fgravj(:) + pmassi*dsofti*dr(:)
       else
          phii = -rij1
          fmi  = rij21
       endif
       if (q2j < radkern2) then
          qj = sqrt(q2j)
          grkernj = grkern(q2j,qj)
          call kernel_softening(q2j,qj,phij,fmj)
          phij       = phij*hj1
          fmj        = fmj*hj21
          grkernj    = cnormk*grkernj*hj41*gradh(1,indexj)
          dsoftj     = 0.5*grkernj*gradh(2,indexj)
          print*, "dsoftj: ", dsoftj
          print*, "pmassi: ", pmassi 
          fgravi(:)  = fgravi(:) - pmassj*dsoftj*dr(:)
          fgravj(:)  = fgravj(:) + pmassi*dsoftj*dr(:)
          print*, "force is: ", fgravi, fgravj
       else
          phij = -rij1
          fmj  = rij21
       endif

       phiterm = 0.5*(phii + phij)
       phitemp = phitemp + pmassj*phiterm
       !phi(j) = phi(j) + pmassi*phiterm
       !phitot = phitot + pmassj*pmassi*phiterm

       fm = 0.5*(fmi + fmj)
       fgravi(1:3) = fgravi(1:3) - pmassj*dr(1:3)*fm
       if (iactivej) then
          !$omp critical(direct)
          fxyzu_dtt(1:3,indexj) = fxyzu_dtt(1:3,indexj) + pmassi*dr(1:3)*fm + fgravj(1:3)
          !$omp end critical(direct)
       endif
    enddo overj
!
!--add self contribution to potential
!
    if (iactivei) then
       !$omp critical (direct)
       fxyzu_dtt(1:3,indexi) = fxyzu_dtt(1:3,indexi) + fgravi(1:3)
       !$omp end critical (direct)
    endif
    !phi(i) = phitemp + pmassi*potensoft0*hi1
    !phitot = phitot + pmassi*phitemp + pmassi*pmassi*potensoft0*hi1

    start = start + 1 
 enddo overi

! phitot = 0.
! do i=1,ntot
!    phitot = phitot + 0.5*pmassi*phi(i)
! enddo
 !phitot = 0.5*phitot

 return

 end subroutine direct_sum_softened

 subroutine get_part_node(nodeindex,node,ifirstincell,part,nopart)
  type(kdnode), intent(in) :: node(:)
  integer,      intent(in) :: nodeindex, ifirstincell(:)
  integer,     intent(out) :: nopart,part(:)
  integer, parameter :: maxstacksize = 100
  integer :: nstack(maxstacksize), istack,n,npnode,currentindex
  integer :: il,ir,ipart,i1,i2,j


  npnode = 0.
  nopart = 0.
  nstack = 0.
  istack = 1
  part = 0.
  ! Push nodeindex to stack 
  nstack(istack) = nodeindex

  do while (istack /= 0)
    !print*, "istack is: ", istack 
    n = nstack(istack)
    istack = istack - 1 
    il      = node(n)%leftchild
    ir      = node(n)%rightchild
    if (ifirstincell(n) /= 0) then
      ! leafnode, add to particles 
      i1=inoderange(1,n)
      i2=inoderange(2,n)
      npnode = i2 - i1 + 1
      do j=1,npnode
        part(nopart+j) = abs(inodeparts(inoderange(1,n)+j-1))
      enddo
      nopart = nopart + npnode 
    else 
        if (il /= 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
          if (ir /= 0) then
             istack = istack + 1
             nstack(istack) = ir
          endif 

    endif 

  enddo 


 end subroutine get_part_node
logical function check_parent_overlap(node,currentnodeindex,interactnodeindex)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 use kernel, only: radkern
 type(kdnode), intent(in) :: node(:)
 integer, intent(in) :: currentnodeindex, interactnodeindex
 integer :: n 
 real :: dx,dy,dz,r2,rcut,rcut2,rcuti,rcutj,xpos(3),xoffset,yoffset,zoffset
 real :: xsizei,xsizej
#ifdef PERIODIC
 real :: hdlx,hdly,hdlz

 hdlx = 0.5*dxbound
 hdly = 0.5*dybound
 hdlz = 0.5*dzbound
#endif
    ! Setup data for current node index 
    xpos = node(currentnodeindex) % xcen
    xsizei = node(currentnodeindex) % size 
    rcuti = radkern*node(currentnodeindex)%hmax

    ! start with the parent of the node we are interacting with 
    n = node(interactnodeindex) % parent 

    do while(n /= 1) ! Traverse till root 
      dx = xpos(1) - node(n)%xcen(1)      ! distance between node centres
      dy = xpos(2) - node(n)%xcen(2)
#ifndef TREEVIZ
      dz = xpos(3) - node(n)%xcen(3)
#endif
      xsizej       = node(n)%size
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
      r2    = dx*dx + dy*dy + dz*dz !+ xsizei*xsizei
      !if (get_hj) then  ! find neighbours within both hi and hj
      rcutj = radkern*node(n)%hmax
      rcut  = max(rcuti,rcutj)
    !endif


      rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius

      ! If parents dont overlap this interaction must be added 
      ! return false 
      if (r2>=rcut2) then 
        check_parent_overlap =  .false.
        return 
      endif 

      ! otherwise check overlap up to root 
      n = node(n) % parent 

    enddo 

    ! if root has been reached then nodes overlap
    check_parent_overlap = .true.


end function check_parent_overlap




subroutine get_child_leaf(nodeindex,node,ifirstincell,childleaf,childleafcounter)
  integer, intent(in) :: nodeindex, ifirstincell(:)
  type(kdnode), intent(in) :: node(:)
  integer, intent(out) :: childleaf(10000), childleafcounter
  integer :: n, istack, nstack(100)
  integer :: il,ir 

  childleafcounter = 0
  childleaf = 0
  istack = 1
  nstack(istack) = nodeindex 

  do while (istack /= 0)
    n = nstack(istack)
    istack = istack - 1
    !print*, n
    il = node(n) % leftchild
    ir = node(n) % rightchild
    if (ifirstincell(n) > 0) then ! If Leaf, add to childleaf 
      childleafcounter = childleafcounter + 1 
      childleaf(childleafcounter) = n 
      

    else 
       if (il > 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
          if (ir > 0) then
             istack = istack + 1
             nstack(istack) = ir
          endif
    endif 


  enddo 



end subroutine get_child_leaf

subroutine push_local(nodeindex1,nodeindex2,stacklocal,istacklocal,k)
  integer, intent(in) :: nodeindex1, nodeindex2,k
  integer, intent(inout) :: istacklocal(:)
  type(denstreestacklocal), intent(inout) :: stacklocal(:)

  ! Increment stack pointer
  istacklocal(k) = istacklocal(k) + 1 
  ! push nodes 
  stacklocal(istacklocal(k)) % nodeindex1 = nodeindex1
  stacklocal(istacklocal(k)) % nodeindex2 = nodeindex2


end subroutine push_local

subroutine pop_local(nodeindex1,nodeindex2, stacklocal, istacklocal,k)
  integer, intent(out) :: nodeindex1, nodeindex2
  integer, intent(inout) :: istacklocal(:)
  type(denstreestacklocal), intent(inout) :: stacklocal(:)
  integer, intent(in) :: k 

  ! pop nodes 
  nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1 
  nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2 

  ! Decrement stack pointer 
  istacklocal(k) = istacklocal(k) - 1 

end subroutine pop_local

subroutine push_global(thenode, interaction1, interaction2,stack,top)
  integer, intent(in) :: thenode, interaction1, interaction2
  integer, intent(inout) :: top 
  type(denstreestack), intent(inout) :: stack(:)

  top = top + 1 
  stack(top) % nodeindex1 = thenode
  stack(top) % interactions(1) = interaction1
  stack(top) % interactions(2) = interaction2

end subroutine push_global

logical function not_sph_neighbour(node,nodeindex1,nodeindex2,r2,rcut2)
  type(kdnode), intent(in) :: node(:)
  integer, intent(in) :: nodeindex1,nodeindex2
  real, intent(in) :: r2,rcut2

  if (r2 < rcut2 .and. (node(nodeindex1) % leftchild == 0 .and. node(nodeindex1) % rightchild == 0)) then 
    not_sph_neighbour = .not. check_parent_overlap(node,nodeindex1,nodeindex2)
  else if ( r2 >= rcut2) then 
    not_sph_neighbour =  .true.
  endif 

end function not_sph_neighbour

logical function check_leaf_overlap(node,currentnodeindex,interactnodeindex)
 type(kdnode), intent(in) :: node
 integer, intent(in) :: currentnodeindex, interactnodeindex
 integer :: childleaf(100), childleafcounter

 !call get_child_leaf(nodeindex,node,ifirstincell,childleaf,childleafcounter) 



end function check_leaf_overlap

logical function check_child_overlap(node,currentnodeindex,interactnodeindex) result(overlap)
 use kernel, only: radkern
 type(kdnode), intent(in) :: node(:)
 integer, intent(in) :: currentnodeindex, interactnodeindex
 integer :: leftchildindex, rightchildindex
 real    :: r2,rcut2,dx,dy,dz,rcuti,rcutj,rcut,xsizei,xsizej
 logical :: overlapleft, overlapright 

 leftchildindex = node(currentnodeindex) % leftchild
 rightchildindex = node(currentnodeindex) % rightchild

 ! Check overlap of left child first 

 call get_node_data(node(leftchildindex),node(interactnodeindex),dx,dy,dz,xsizei,xsizej)

 r2 = dx*dx + dy*dy + dz*dz 

 rcutj = radkern*node(interactnodeindex) % hmax
 rcuti = radkern*node(leftchildindex) % hmax 
 rcut  = max(rcuti,rcutj)
 !rcut = rcuti
 rcut2 = (xsizei + xsizej + rcut)**2 

 overlapleft = r2 < rcut2 

 call get_node_data(node(rightchildindex),node(interactnodeindex),dx,dy,dz,xsizei,xsizej)

 r2 = dx*dx + dy*dy + dz*dz 

 rcutj = radkern*node(interactnodeindex) % hmax
 rcuti = radkern*node(rightchildindex) % hmax 
 rcut  = max(rcuti,rcutj)
 !rcut = rcuti
 rcut2 = (xsizei + xsizej + rcut)**2 

 overlapright = r2 < rcut2 

 !print*, "overlapright: ", overlapright
 !print*, "overlapleft: ", overlapleft
 overlap = overlapleft .or. overlapright


end function check_child_overlap


subroutine get_nodesize_max(node, xyzh, ifirstincell, ncells)
 use dim,  only:maxp
 use part, only:maxphase,iphase,igas,massoftype,iamtype
 use io,   only:fatal
 type(kdnode), intent(inout) :: node(:) !ncellsmax+1)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: ifirstincell(:) !ncellsmax+1)
 integer(kind=8), intent(in) :: ncells
 real :: hmax, r2max
 real :: xi, yi, zi, hi
 real :: dx, dy, dz, dr2
#ifdef GRAVITY
 real :: quads(6)
#endif
 integer :: icell, i, level, il, ir
 real :: pmassi, totmass
 real :: x0(3)
 type(kdnode) :: nodel,noder

 pmassi = massoftype(igas)

! Test for get level routine 

if (node_level(node,irootnode) /= 0) then
 print*, "get level failed" 
 stop 
endif 

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(xyzh, ifirstincell, ncells) &
!$omp shared(node, ll, iphase, massoftype, maxlevel,maxlevel_indexed) &
!$omp private(hmax, r2max, xi, yi, zi, hi, il, ir, nodel, noder) &
!$omp private(dx, dy, dz, dr2, icell, i, x0) &
#ifdef GRAVITY
!$omp private(quads) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)

! Loop over each leaf node and set maxsize = size 
!$omp do 
do i=1, ncells

  if (ifirstincell(i) <= 0) cycle 

  node(i) % maxsize = node(i) % size

enddo 
!$omp enddo 

!$omp end parallel 
!
! propagate information to parent nodes
! here we sweep across each level at a time
! and update each node from its two children
!
  !print*,"maxlevel: ", maxlevel
  !print*,"maxlevel_indexed: ", maxlevel_indexed
 do level=maxlevel-1,maxlevel_indexed,-1
   !print*, "level: ", level 
   !print*, 2**maxlevel_indexed, ncells
   !print*,"in loop"
  !!$omp do 
  do i=2**maxlevel_indexed, ncells
    !print*,"i is: ", i 
    if (node_level(node,i) == level) then
      !print*, "node is: ", i  
      ! get child nodes
      il = node(i)%leftchild
      ir = node(i)%rightchild
      ! if (ifirstincell(il) > 0) then 
      !   node(il) % maxsize = node(il) % size 
      ! endif 
      ! if (ifirstincell(ir) > 0) then 
      !   node(ir) % maxsize = node(ir) % size 
      ! endif 
      !print*, "leftchild: ", il 
      if (il > 0 .and. ir > 0) then
          nodel = node(il)
          noder = node(ir)
          call add_child_nodes(nodel,noder,node(i))
      else
          if (il > 0 .or. ir > 0) then
             ! should never happen, should have two children or none
             call fatal('revtree','node with only one child during tree revision',var='ir',ival=ir)
          endif
       endif

    endif 

  enddo 
  !!$omp enddo 

 enddo 
 do level=min(maxlevel_indexed-1,maxlevel-1),0,-1
!do level=maxlevel-1,0,-1
   !print*, "level: ", level 
   !print*, 2**maxlevel_indexed, nnodes
   !print*, 2**level,2**(level+1)-1
!!$omp do
    ! 1 to nnodes for level > maxlevel 
    do i=2**level,2**(level+1)-1
       if(i==252107) print*, "level of node is: ", level
       ! get child nodes
       il = node(i)%leftchild
       ir = node(i)%rightchild
      !  if (ifirstincell(il) > 0) then 
      !   node(il) % maxsize = node(il) % size 
      ! endif 
      ! if (ifirstincell(ir) > 0) then 
      !   node(ir) % maxsize = node(ir) % size 
      ! endif 
       if (il > 0 .and. ir > 0) then
          nodel = node(il)
          noder = node(ir)
          call add_child_nodes(nodel,noder,node(i))
       else
          if (il > 0 .or. ir > 0) then
             ! should never happen, should have two children or none
             call fatal('revtree','node with only one child during tree revision',var='ir',ival=ir)
          endif
       endif
    enddo
!!$omp enddo
 enddo
!!$omp end parallel



end subroutine get_nodesize_max

integer function node_level(node,nodeindex)
  type(kdnode), intent(in) :: node(:)
  integer, intent(in) :: nodeindex
  integer :: level,currentnodeindex 

  level  = 0 
  currentnodeindex = nodeindex
  !print*, "cell is: ", nodeindex
  ! Loop until we reach root 
  do while (currentnodeindex /= irootnode)
    !print*, "node index is: ", currentnodeindex
    level = level + 1 
    !print*, "level is: ", level 
    currentnodeindex = node(currentnodeindex) % parent
    !print*, "new index is: ", currentnodeindex
    if (currentnodeindex == 0 ) then 
      node_level = 0 
      return 
    endif 
  enddo 

  
  !print*,'level is: ', level 
  node_level = level 


end function node_level


 ! recursive subroutine get_node_radius(node,ifirstincell,n,radius)
 !  use kernel, only: radkern
 !  type(kdnode), intent(in) :: node(:)
 !  integer, intent(inout) :: n
 !  integer, intent(in) :: ifirstincell(:)

 !  real, intent(inout) :: radius 
 !  real :: sizel, sizer

 !    if (ifirstincell(n) /= 0) then 
 !      radius = node(n) % size + node(n) % hmax*radkern
 !    else 
 !      ! get search radius of left child 
 !      call get_node_radius(node,ifirstincell,node(n)%leftchild,sizel)
 !      ! get search radius of right child 
 !      call get_node_radius(node,ifirstincell,node(n)%rightchild,sizer)

 !      radius = max(sizel,sizer)

 !    endif 

 ! end subroutine get_node_radius
end module kdtree
