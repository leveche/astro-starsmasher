      subroutine uvdots
c     evaluate right-hand sides of equations of motion:
c=======================================================================
c  balvdots creates arrays of the accelerations of each particle, using
c  balsara's form of the artificial viscosity pi_{ij}.
c  this incarnation of balvdots implements the form in which only the pi_{ij}
c     term is multiplied by the gradient of a symmetrized w_{ij}.
c     the p/rho^2 terms are each multiplied by the gradient of the kernel
c     evaluated with the smoothing length of that particle only.
c     see springel and hernquist's mnras article (astro-ph/0111016)
c     and monaghan's compressible turbulence article (astro-ph/0204118).
c=======================================================================
      include 'starsmasher.h'
      include 'cooling/cooling.h'
      include 'mpif.h'
      double precision curlvxi,curlvyi,curlvzi
      real*8 csij,udbij,fi,ci,r2,
     $     abscurli,divvi,vxi,vyi,vzi,ci2,ami,por2i,h5,
     $     h2,hpi,csiju,csijgc
      integer itab,j,in,i,ierr
      real*8 uijmax(nmax)
      common/uijmax/ uijmax
      real*8 dwijin,pijin,sxijin,syijin,szijin,divin, diwkin
      real*8 sxijingc,syijingc,szijingc
      logical switchedav
      double precision divv(nmax)
      common/commdivv/divv
      integer mylength, mygravlength
      integer comm_worker
      common/gravworkers/comm_worker
      real*8 myvxdot(ntot),myvydot(ntot),myvzdot(ntot)
      real*8 myvxdotgc(ntot),myvydotgc(ntot),myvzdotgc(ntot)
      integer irank
      real*8 invh, invh2, invh5
      real*8 dx, dy, dz
      real*8 dvx, dvy, dvz
      real*8 invbonnet
      real*8 alpha1, beta1
      integer offset, nni
c     variables used for radiative cooling portion of the code:
      real*8 uunit,udotunit,rhocgs,ucgs,temperature,kappa
      real*8 opacity,pgas,prad
      real*8 gradponrho,scaleheight,columnden
      real*8 kappar
      real*8 ueqfunction,useeostable,uupperlimit!,ulowerlimit
      external ueqfunction
      real*8 zeroin,mucgs
      common/ueqstuff/rhocgs,teq,mucgs
      real*8 grpottot(nmax),uraddoti
      real*8 hpitilde,h2tilde,invh2tilde
      integer itabtilde
      real*8 dudt
      
      if(nav.eq.0) then
         if(nintvar.eq.1)then
            do i=1,n
               udot(i)=0.d0
            enddo
         else
            error stop 'turn off av by setting alpha=beta=0'
         endif
      endif

      switchedav=.false.
      if(nrelax.eq.1
     $     .and. t.lt.treloff .and. alpha.eq.0.d0
     $     .and. beta.eq.0.d0) then
         switchedav=.true.
         alpha=1.d0
         beta=2.d0
         if(myrank.eq.0)
     $        write(69,*)'vdots: during relaxation alpha,beta=',alpha,beta
      endif

c  initialize dv/dt before accumulating:
      do i=1,n
         myvxdot(i)=0.d0
         myvydot(i)=0.d0
         myvzdot(i)=0.d0
         myVXDOTgc(I)=0.d0
         myVYDOTgc(I)=0.d0
         myVZDOTgc(I)=0.d0
         uijmax(i)=0.d0
         udot(i)=0.d0
      end do

      if (nrelax.eq.1 .and.t.lt.treloff) then
         if(myrank.eq.0)write(69,*)
     $        'baludots is not using pi_{ij} right now'
      endif

c     the factor of 2d0 in the following lines is because our pi_{ij} below
c     is going to involve only p_i/rho_i^2 and not, as is often done,
c     (p_i/rho_i^2 + p_j/rho_j^2).  to keep alpha and beta values at the same
c     strength as in other implementations of sph, we multiply by 2 here.
      alpha1 = -2.d0 * alpha
      beta1  =  2.d0 * beta

c     for each particle:
      do i=n_lower,n_upper
         hpi=hp(i)
         hpitilde=hpi - hfloor
         h2=hpi**2
         h2tilde=hpitilde**2
         h5=hpi*h2**2
         invh  = 1.0/hpi
         invh2 = 1.0/h2
         invh2tilde = 1.0/h2tilde
         invh5 = 1.0/h5
         ami=am(i)
         por2i=por2(i)
         ci2=gam*por2i*rho(i)
         ci=sqrt(ci2)
         vxi=vx(i)
         vyi=vy(i)
         vzi=vz(i)
         call getderivs(i,curlvxi,curlvyi,curlvzi,divvi)
         divv(i)=divvi
         abscurli=sqrt(curlvxi**2+curlvyi**2+curlvzi**2)

c     calculate gradwij's to all neighbors:
         offset = first(i)
         nni    = nn(i)
         do in=1,nni
            j=list(offset+in)
       
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            dvx = vxi - vx(j)
            dvy = vyi - vy(j)
            dvz = vzi - vz(j)
            r2 = dx*dx + dy*dy + dz*dz 
            itab=int(ctab*r2*invh2)+1

            if(r2.lt.4d0*h2tilde) then
               itabtilde=int(ctab*r2*invh2tilde)+1
               diwkin = dgtab(itabtilde) * invh2tilde
            else
               diwkin = 0d0
            endif

            dwijin = dwtab(itab) * invh5

            if(u(j).ne.0.d0 .and. u(i).ne.0) then
              divin = dx*dvx + dy*dvy + dz*dvz
               if (divin.lt.0.d0) then
c                  r2=(x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0
c     $                 +(z(i)-z(j))**2.d0
c     cj=sqrt(5.d0/3.d0*por2(j)*rho(j))
                  
c     calculate dinshaw balsara's (db) \mu_{ij}: see equation (a12) of http://adsabs.harvard.edu/abs/2012apj...745...71p
                  fi=abs(divvi)/(abs(divvi)+abscurli+
     $                 0.00001d0*ci*invh)
c                 fj=abs(divv(j))/(abs(divv(j))+
c     $                    sqrt(curlvx(j)**2+curlvy(j)**2+curlvz(j)**2)+
c     $                    0.00001d0*cj/hp(j))
c     modified from very old version:
                  if(r2.gt.0) then
                     if(nav.eq.3) then
                        udbij=divin/(ci*r2**0.5d0)*fi
                     else
                        udbij=hpi*divin/(ci*r2)*fi
                     endif
                  else
                     udbij=0d0
                  endif

c     old version:
c     udbij=div(in)/((ci+cj)*r2**0.5d0)*(fi+fj)
c     udbij=divin/(ci*r2**0.5d0)*fi
                  if(abs(udbij).gt.1.d8) then
                     write(69,*)'c warning large udbij...',i,j, udbij,divin,ci,
     $                    am(i),am(j),u(i),u(j)
                     if(dt.lt.1d-17)then
                        call mpi_finalize(ierr)
                        error stop 'unreasonably small dt. abort run.'
                     endif
                  endif
c                               a     +    b * c   (multiply-add, cpu likes these)
                  pijin=por2i*(alpha1 + beta1*udbij)*udbij
               else
                  pijin=0.d0
               endif
c     hydro correction part.  here we are implementing corrections as shown
c     in equations (a11) and (a12) of gaburov et al. (2010):
               invbonnet = 1.0/(bonet_0mega(i)*am(j))

c     csij will not contain AV (pijin) contributions so that we can
c     later use the acceleration derived from it to calculate scale height.
c     All AV contributions to the acceleration will be accounted for with
c     the help of csijgc below. (gc=gravity correction... this is for terms
c     associated with the gravity correction *and* with AV)
               csij = por2i * (-bonet_omega(i) * invbonnet * diwkin
     $              + dwijin)
               csijgc= 0.5d0 * pijin * dwijin

c     hydro part
               if (nrelax.eq.1 .and.t.lt.treloff) then
                  csiju= csij ! No heating from AV during relaxation
               else
                  csiju= csij + csijgc
               endif

               if(nintvar.eq.2) then
                  udot(i)=udot(i)+am(j)*csiju*divin
               elseif (nrelax.eq.0 .or. t.ge.treloff) then
                  udot(i)=udot(i)+am(j)*pijin*dwijin*divin
               endif

               uijmax(i)=max(uijmax(i),
     $              (por2i + 0.5d0 * pijin)*rho(i))
c               uijmax(j)=max(uijmax(j),
c     $              (por2i + 0.5d0 * pijin)*rho(i))

               if(ngr.ne.0) then
c     gravity correction part
c                  csij = csij - 0.5d0*diwkin * bonet_psi(i) * invbonnet
                  csijgc = csijgc - 0.5d0*diwkin * bonet_Psi(i) * invBonnet
               endif
            else
               csij=0.d0
               csijgc=0.d0
            endif
            
            sxijin = csij * dx !(x(i) - x(j))
            syijin = csij * dy !(y(i) - y(j))
            szijin = csij * dz !(z(i) - z(j))
            sxijingc = csijgc * dx !(x(i) - x(j))
            syijingc = csijgc * dy !(y(i) - y(j))
            szijingc = csijgc * dz !(z(i) - z(j))

c     finally compute each component of dv/dt:

c     note: sxijin has been calculated using pi/rhoi^2.  so in the "gather"
c     part of the sum, everything is as you would expect.  but the gather part
c     of the sum only gives half of the contributions for any given particle
c     i.  for the sake of concreteness, let's say i=5.  the other half of the
c     contributions to the acceleration of particle 5 comes from the scatter
c     portion of the code with j=5 (note: this is *j*, not i).  in the scatter
c     portion, the index i will be a number other than 5, and so the p/rho^2
c     being used inside of sxijin will come from that particle with index
c     different than 5.  the gather contributions to particle 5 all happen
c     together in one big group.  the scatter contributions to particle 5 will
c     happen haphazardly as the loop over particle index i gradually covers
c     all the particles that have 5 as a neighbor.

c     "gather" part of the sum:
            myvxdot(i)=myvxdot(i)-am(j)*sxijin
            myvydot(i)=myvydot(i)-am(j)*syijin
            myvzdot(i)=myvzdot(i)-am(j)*szijin
c     "scatter" part of the sum:
            myvxdot(j)=myvxdot(j)+ami*sxijin
            myvydot(j)=myvydot(j)+ami*syijin
            myvzdot(j)=myvzdot(j)+ami*szijin
c     "Gather" part of the sum:
            myvxdotgc(i)=myvxdotgc(i)-am(j)*sxijingc
            myvydotgc(i)=myvydotgc(i)-am(j)*syijingc
            myvzdotgc(i)=myvzdotgc(i)-am(j)*szijingc
c     "Scatter" part of the sum:
            myvxdotgc(j)=myvxdotgc(j)+ami*sxijingc
            myvydotgc(j)=myvydotgc(j)+ami*syijingc
            myvzdotgc(j)=myvzdotgc(j)+ami*szijingc
         enddo

c     p=(gam-1)*rho*u=a*rho^gam, so a=(gam-1)*rho^(1-gam)*u
         if(nintvar.eq.1 .and. u(i).ne.0.d0)
     $        udot(i)=udot(i)*0.5d0*(gam-1.d0)*rho(i)**(1.d0-gam)
      enddo
c      write(6,'(a)')'hydrompi'

       mylength=n_upper-n_lower+1
c      do irank=0,nprocs-1
c         if(myrank.ne.irank)then
c            call mpi_gatherv(udot(n_lower), mylength, mpi_double_precision,
c     $           udot, recvcounts,displs, mpi_double_precision, irank,
c     $           mpi_comm_world, ierr)
c         else
c            call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
c     $           udot, recvcounts,displs, mpi_double_precision, irank,
c     $           mpi_comm_world, ierr)
c         endif
c      enddo
      call mpi_allgatherv(mpi_in_place, mylength, mpi_double_precision,
     $     udot, recvcounts,displs, mpi_double_precision, 
     $     mpi_comm_world, ierr)

      uunit=gravconst*munit/runit
      udotunit=(gravconst*munit)**1.5d0/runit**2.5d0
      
      if(ncooling.eq.1) then
c     to get the pressure scale height, we will need the v{x,y,z}dot arrays *without* the
c     gravitational acceleration included.  note: the myv{x,y,z}dot arrays are not
c     affected by the following mpi calls.  instead, they are just being summed up to give
c     the v{x,y,z}dot arrays.
c         if(ncooling.eq.1)then
            call mpi_allreduce(myvxdot,vxdot,n,mpi_double_precision,
     $           mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(myvydot,vydot,n,mpi_double_precision,
     $           mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(myvzdot,vzdot,n,mpi_double_precision,
     $           mpi_sum,mpi_comm_world,ierr)
c         end if
         
         punit=gravconst*(munit/runit**2)**2
c     if(myrank.eq.0) write(69,*) 'uunit,udotunit=',uunit,udotunit

         if(neos.eq.1)then
            do i=n_lower,n_upper
               if(u(i).ne.0.d0) then
                  rhocgs=rho(i)*munit/runit**3.d0
                  if(nintvar.eq.1) then
                     ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                    *uunit
                  else
                     ucgs=u(i)*uunit
                  endif
                  temperature = temperatures(i)
                  pgas=rhocgs*boltz*temperature/meanmolecular(i)/
     $                 punit
                  prad=arad*temperature**4/3.d0/punit
                  
                  gradponrho=(vxdot(i)**2+vydot(i)**2+
     $                 vzdot(i)**2)**0.5d0
                  scaleheight=(pgas+prad)/(gradponrho*rho(i))
                  columnden=rho(i)*scaleheight*1.06d0*munit/runit**2
c     kappa is local opacity; kappar is pseudo-mean opacity
                  call usetable(rhocgs,temperature,kappa,kappar)
                  uraddoti=4*sigma*(teq**4-temperature**4)/
     $                 (columnden**2*kappar + 1/kappa)
     $                 /udotunit
                  ueq(i)=(arad*teq**4/rhocgs
     $                 +1.5d0*boltz*teq/meanmolecular(i))
     $                 /uunit
                  
                  tthermal(i)=(ueq(i)-u(i))/uraddoti

c     approximate the equilibrium temperature... not really correct...
c     teq=(10.d0**4+udot(i)*(columnden**2*kappa+ 1/kappa)
c     $           *udotunit/(4*sigma))**0.25d0
                  
               else
                  tthermal(i)=1.d30
                  ueq(i)=0.d0
               endif
            enddo
         else if(neos.eq.2)then
            do i=n_lower,n_upper
               if(u(i).ne.0.d0) then
                  rhocgs=rho(i)*munit/runit**3.d0
                  mucgs=meanmolecular(i)
                  if(nintvar.eq.1) then
                     ucgs=u(i)*rho(i)**(gam-1.d0)/(gam-1.d0)
     $                    *uunit
                  else
                     ucgs=u(i)*uunit
                  endif
                  temperature = temperatures(i)
                  gradponrho=(vxdot(i)**2+vydot(i)**2+
     $                 vzdot(i)**2)**0.5d0
                  scaleheight=useeostable(ucgs,rhocgs,
     $                 meanmolecular(i),3)/(punit*gradponrho*rho(i))
                  columnden=rho(i)*scaleheight*1.06d0*munit/runit**2
c     kappa is local opacity; kappar is pseudo-mean opacity
                  call usetable(rhocgs,temperature,kappa,kappar)
                  uraddoti=4*sigma*(teq**4-temperature**4)/
     $                 (columnden**2*kappar + 1/kappa)
     $                 /udotunit
                  
c     using eos table, solve for u_eq such that temperature t=teq 
                  uupperlimit=1.1d0*(1.5d0*boltz*teq/meanmolecular(i)
     $                 +arad*teq**4/rhocgs)
                  ucgs=zeroin(0d0,uupperlimit,ueqfunction,1.d-11)
                  if(ueqfunction(ucgs).gt.1d-6) then
                     uupperlimit=2*uupperlimit
                     write(69,*)'Having trouble with u_eq for particle'
     $                    ,i,'at time',t,ueqfunction(ucgs)
                     ucgs=zeroin(0d0,uupperlimit,ueqfunction,1.d-11)
                     write(69,*)'Found',ucgs,'for'
     $                    ,i,'at time',t,ueqfunction(ucgs)
                  endif
                  ueq(i)=ucgs/uunit

                  tthermal(i)=(ueq(i)-u(i))/uraddoti
               else
                  tthermal(i)=1.d30
                  ueq(i)=0
               endif
            enddo
         else
            error stop 'radiative cooling requires neos=1 or 2'
         endif

         mylength=n_upper-n_lower+1
c         do irank=0,nprocs-1
c            if(myrank.ne.irank)then
c               call mpi_gatherv(tthermal(n_lower), mylength, mpi_double_precision,
c     $              tthermal, recvcounts,displs, mpi_double_precision, irank,
c     $              mpi_comm_world, ierr)
c               call mpi_gatherv(ueq(n_lower), mylength, mpi_double_precision,
c     $              ueq, recvcounts,displs, mpi_double_precision, irank,
c     $              mpi_comm_world, ierr)
c            else
c               call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
c     $              tthermal, recvcounts,displs, mpi_double_precision, irank,
c     $              mpi_comm_world, ierr)
c               call mpi_gatherv(mpi_in_place, mylength, mpi_double_precision,
c     $              ueq, recvcounts,displs, mpi_double_precision, irank,
c     $              mpi_comm_world, ierr)
c            endif
c         enddo
         call mpi_allgatherv(mpi_in_place, mylength, mpi_double_precision,                                   
     $        tthermal, recvcounts,displs, mpi_double_precision,
     $        mpi_comm_world, ierr)
         call mpi_allgatherv(mpi_in_place, mylength, mpi_double_precision,                                   
     $        ueq, recvcounts,displs, mpi_double_precision,
     $        mpi_comm_world, ierr)


      else if(ncooling.eq.2) then
c     Update the cooling frequency such that no particle ever gets u <= 0

c     Try an initial increase in the frequency so that we do cooling as
c     few times as possible.
         cooling_update_frequency = max_cooling_update_frequency
         do i=1, n
            if (u(i) .ne. 0.d0) then
               dudt = udot(i)+uraddot(i)
               if ( u(i)+dudt*dt .le. 0.d0 ) then
                  write(*,*) "Particle",i,"will cool to u <= 0"
                  write(*,*) "du/dt = ",dudt
                  write(*,*) "u, udot, uraddot = ",u(i),udot(i),
     $                 uraddot(i)
                  write(*,*) "dt = ",dt
                  error stop "balAV3.f"
               end if
               
c     Decrease frequency until the particle will no longer get u <= 0
               do while ( u(i)+dudt*cooling_update_frequency .lt. 0.d0
     $              .and. cooling_update_frequency .gt. dt )
                  cooling_update_frequency = max(
     $                 0.99d0*cooling_update_frequency,dt)
               end do
               
c     Don't need to check the other particles because we are already
c     at the minimum possible time step
               if ( cooling_update_frequency .eq. dt ) exit
            end if
         end do

c     Do the cooling
         if (abs(t-tcool) .ge. cooling_update_frequency) then
            tcool = t           ! Reset the timer

            if (neos.ne.1.and.neos.ne.2) then
               error stop 'radiative cooling requires neos=1 or 2'
            endif

            call get_ncooling2
         end if
      end if

c      write(6,'(a)')'hydrompi_complete'
      if(ngr.ne.0 .and. myrank.lt.ngravprocs)then
         if(myrank.eq.ngravprocs-1) call cpu_time(time0)
         if(nusegpus.eq.1)then
            call lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
         else
            call get_gravity_using_cpus
         endif   
         if(myrank.eq.ngravprocs-1) call cpu_time(time1)
         if(ngravprocs.gt.1) then
            if(nusegpus.eq.1)then
               mygravlength=ngrav_upper-ngrav_lower+1
c               if(myrank.ne.0) then
c                  call mpi_gatherv(grpot(ngrav_lower), mygravlength, mpi_double_precision,
c     $                 grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
c     $                 comm_worker, ierr)
c               else
c                  call mpi_gatherv(mpi_in_place, mygravlength, mpi_double_precision,
c     $                 grpot, gravrecvcounts, gravdispls, mpi_double_precision, 0,
c     $                 comm_worker, ierr)
c               endif
               call mpi_allgatherv(mpi_in_place, mygravlength, mpi_double_precision,                         
     $              grpot, gravrecvcounts,gravdispls, mpi_double_precision,                                  
     $              comm_worker, ierr)
            else
               call mpi_reduce(grpot, grpottot, n, mpi_double_precision,
     $              mpi_sum, 0, mpi_comm_world, ierr)
               if(myrank.eq.0) then
                  do i=1,n
                     grpot(i)=grpottot(i)
                  enddo
               endif
            endif
         endif
         if(myrank.eq.ngravprocs-1)then
            write(6,'(a,f6.3,a,g10.3,a,i4)')'2ndhalf:',time1-time0
         endif
         if(nusegpus.eq.1)then
            do i=ngrav_lower,ngrav_upper
               myvxdot(i)=myvxdot(i)+gx(i)
               myvydot(i)=myvydot(i)+gy(i)
               myvzdot(i)=myvzdot(i)+gz(i)
            enddo
         else
            do i=1,n
               myvxdot(i)=myvxdot(i)+gx(i)
               myvydot(i)=myvydot(i)+gy(i)
               myvzdot(i)=myvzdot(i)+gz(i)
            enddo
         endif
      endif

ccccccccccccccccccccccccccccccccccccccccc
c     reassemble acceleration arrays here
ccccccccccccccccccccccccccccccccccccccccc

c     use mpi_allreduce to collect subtotals, sum them, and redistribute
c     to all processes

c     Now that scaleheight has been calculated, we can include all
c     relevant terms in the calculation of the acceleration:
      do i=1,n
         myvxdot(i)=myvxdot(i)+myvxdotgc(i)
         myvydot(i)=myvydot(i)+myvydotgc(i)
         myvzdot(i)=myvzdot(i)+myvzdotgc(i)
      enddo

      call mpi_allreduce(myvxdot,vxdot,n,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierr)
      call mpi_allreduce(myvydot,vydot,n,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierr)
      call mpi_allreduce(myvzdot,vzdot,n,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierr)

cc     if this is a relaxation calculation, add centrifugal and drag forces:
c      if (nrelax.ne.0 .and. .not. gonedynamic ) then
          call relax        ! add drag and fictitious forces if necessary
c     endif

      if(switchedav)then
         alpha=0.d0
         beta=0.d0
      endif

c      write(69,*)'balav3: finished' 

      return
      end
