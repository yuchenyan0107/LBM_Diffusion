!
!  MIXLBM.f90 ( MRT + FORCE + LOOP 3 )
!
!  FUNCTIONS:
!
!	MIXLBM						- Main
!	UpdateLatticeData			- Update local distribution function f in P(i,j)
!	HydrodynamicMoments			- Compute hydrodynamic moments of f
!	EquilibriumDistribution		- Compute the equilibrium distribution function feq
!	writeDAT					- write results to files
!
!	D(M)						- Fick diffusivity
!	B(M1,M2)					- Stefan-Maxwell resistivity
!	NU(s,rho,rhosigma,MMsigma)	- Kinematic viscosity
!	XI(s,rho,rhosigma,MMsigma)	- Bulk viscosity
!
!	MRT(->A)					- MRT Collisional matrix A
!	TRANS(->G)					- Matrix for change of variables G=(I+theta A)^-1
!	FtoM(f->m)					- Compute the discrete moments
!	MtoF(m->f)					- Compute the distribution function
!
!****************************************************************************
!
!  PROGRAM: MIXLBM
!
!  PURPOSE: To compare different implementations of species coupling 
!			(Fick model, Maxwell-Stefan) in LBM framework
!
!****************************************************************************

program MIXLBM

	implicit none

! Local variables: tempVors
	real(8),dimension(:,:,:,:),allocatable	:: fd,fd_new
	real(8),dimension(:,:),allocatable		:: f
	real(8),dimension(:,:,:),allocatable	:: rsigma,psigma,uxsigma,uysigma
	real(8),dimension(:,:),allocatable		:: r,p,ux,uy
	real(8),dimension(:),allocatable		:: phi
	real(8),dimension(:,:),allocatable		:: BCsigma
	real(8),dimension(:,:),allocatable		:: Output
	real(8),dimension(:),allocatable		:: MM
	real(8),dimension(0:8)					:: m

	real(8),dimension(0:8) :: feq
	real(8) :: pi = 3.141592653589793

	CHARACTER (len=2) :: nameN, nameB

! Internal Functions
	real(8) :: mean

! Local variables: scalars
	real(8) :: rho_new,ux_new,uy_new,deltap,shift

	integer :: nx,ny,nt
	integer :: species,md
	integer :: i,j,k,s,t,nN,nB

! Values
	species = 3;
	nx = 60;
	ny = 3;
	nt = 120000;! 1 hour

! Physical model
! 0- Ideal Fick model
! 1- Fick model
! 2- Maxwell-Stefan model
	md = 2;

! Allocate memory
	allocate( fd(1:nx,1:ny,0:8,1:species) )
	allocate( fd_new(1:nx,1:ny,0:8,1:species) )
	allocate( f(0:8,1:species) )
	allocate( rsigma(1:nx,1:ny,1:species) )
	allocate( psigma(1:nx,1:ny,1:species) )
	allocate( uxsigma(1:nx,1:ny,1:species) )
	allocate( uysigma(1:nx,1:ny,1:species) )
	allocate( r(1:nx,1:ny) )
	allocate( p(1:nx,1:ny) )
	allocate( ux(1:nx,1:ny) )
	allocate( uy(1:nx,1:ny) )
	allocate( phi(1:species) )
	allocate( BCsigma(1:4,1:species) )
	allocate( MM(1:species) )
	allocate( Output(1:nx,1:ny) )

do nN = 20,20,1
do nB = 15,15,1

! Molecular weight
	MM(1) = 1.0d0;
	MM(2) = 2.0d0;
	MM(3) = 3.0d0;

! R0*T = 2.0 in lattice units, where R0 universal gas constant
	phi(1) = 1.0d0/MM(1);
	phi(2) = 1.0d0/MM(2);
	phi(3) = 1.0d0/MM(3);

	deltap = 0.01d0;
	shift  = 0.01d0;

	do i = 1,nx,1
		do j = 1,ny,1

!============================================================================
! Initial conditions
!============================================================================

			psigma(i,j,1)  = 0.319*( (1.0d0-tanh(real(i)-real(nx)/2.0d0))/2.0d0 )+1e-4;
			psigma(i,j,2)  = 0.528*( (1.0d0-tanh(real(i)-real(nx)/2.0d0))/2.0d0 )+1e-4;
			psigma(i,j,3)  = 0.847*( (1.0d0+tanh(real(i)-real(nx)/2.0d0))/2.0d0 )+0.153;

!			psigma(i,j,1,0)  = shift+deltap*(1.d0+sin(2*pi*real(i)/real(nx)));
!			psigma(i,j,2,0)  = shift+deltap*(1.d0+cos(2*pi*real(i)/real(nx)));
!			p(i,j)			 = 2.d0*(shift+deltap)+deltap*( sin(2*pi*real(i)/real(nx))+cos(2*pi*real(i)/real(nx)) );
!			psigma(i,j,3,0)  = 1.0d0-p(i,j);

!	Simple initialization of distribution functions (it may be inaccurate)
			do s = 1,species,1

				rsigma(i,j,s)  = psigma(i,j,s)*3.0d0/phi(s);
				uxsigma(i,j,s) = 0.0d0;
				uysigma(i,j,s) = 0.0d0;

				call EquilibriumDistribution(rsigma(i,j,s),phi(s),uxsigma(i,j,s),uysigma(i,j,s),feq)
				fd(i,j,:,s)		= feq(:);
				fd_new(i,j,:,s) = feq(:);

			enddo

		enddo
	enddo

!============================================================================
! Boundary conditions
!============================================================================

	BCsigma(:,1) = (/psigma(nx,1,1),	-1.0d0,	psigma(1,1,1),	-1.0d0/);
	BCsigma(:,2) = (/psigma(nx,1,2),	-1.0d0, psigma(1,1,2), 	-1.0d0/);
	BCsigma(:,3) = (/psigma(nx,1,3),	-1.0d0,	psigma(1,1,3),	-1.0d0/);

!	BCsigma(:,1) = (/		   -2.0d0,	-1.0d0,	         -2.0d0,	-1.0d0/);
!	BCsigma(:,2) = (/		   -2.0d0,	-1.0d0,	         -2.0d0,	-1.0d0/);
!	BCsigma(:,3) = (/		   -2.0d0,	-1.0d0,	         -2.0d0,	-1.0d0/);

!============================================================================
! Solver
!============================================================================
	write(*,'("Iter: "I"  Variable:"ES)') 0, psigma(30,1,3)

	do t = 1,nt,1
		do i = 1,nx,1
			do j = 1,ny,1

				call UpdateLatticeData(md,nx,ny,species,fd,i,j,phi,MM,BCsigma,f,nN,nB);

				do s = 1,species,1

					fd_new(i,j,:,s) = f(:,s);
					
					call HydrodynamicMoments(f(:,s),rsigma(i,j,s),uxsigma(i,j,s),uysigma(i,j,s))
					psigma(i,j,s) = phi(s)*rsigma(i,j,s)/3.0d0;

				enddo
			enddo
		enddo

! ATTENTION
		do j = 1,ny,1
			i = 1;
			do s=1,species,1
				f(:,s) = fd_new(i+1,j,:,s);
				call FtoM(f(:,s),m)
				m(3) = m(3) + (3.0d0*BCsigma(3,s)/phi(s)-m(0))*phi(s)/3.0d0;
				m(4) = m(4) + (3.0d0*BCsigma(3,s)/phi(s)-m(0))*phi(s)/3.0d0;
				m(0) = m(0) + (3.0d0*BCsigma(3,s)/phi(s)-m(0));
				call MtoF(m,f(:,s))
				fd_new(i,j,:,s) = f(:,s);
				call HydrodynamicMoments(f(:,s),rsigma(i,j,s),uxsigma(i,j,s),uysigma(i,j,s))
				psigma(i,j,s) = phi(s)*rsigma(i,j,s)/3.0d0;
			enddo
		enddo

		do j = 1,ny,1
			i = nx;
			do s=1,species,1
				f(:,s) = fd_new(i-1,j,:,s);
				call FtoM(f(:,s),m)
				m(3) = m(3) + (3.0d0*BCsigma(1,s)/phi(s)-m(0))*phi(s)/3.0d0;
				m(4) = m(4) + (3.0d0*BCsigma(1,s)/phi(s)-m(0))*phi(s)/3.0d0;
				m(0) = m(0) + (3.0d0*BCsigma(1,s)/phi(s)-m(0));
				call MtoF(m,f(:,s))
				fd_new(i,j,:,s) = f(:,s);
				call HydrodynamicMoments(f(:,s),rsigma(i,j,s),uxsigma(i,j,s),uysigma(i,j,s))
				psigma(i,j,s) = phi(s)*rsigma(i,j,s)/3.0d0;
			enddo
		enddo

		fd(:,:,:,:) = fd_new(:,:,:,:);

		write(*,'("Iter: "I"  Variable:"ES)') t, psigma(30,1,3)

	enddo

!============================================================================
! Output
!============================================================================

	write(*,'("nN: "I" nB: "I)') nN, nB

	write (nameN, '(I2)') nN
	write (nameB, '(I2)') nB

	call writeDAT(nx,ny,psigma(:,:,1),'p1-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,psigma(:,:,1,nt),'p1.csv');
	call writeDAT(nx,ny,psigma(:,:,2),'p2-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,psigma(:,:,2,nt),'p2.csv');
	call writeDAT(nx,ny,psigma(:,:,3),'p3-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,psigma(:,:,3,nt),'p3.csv');

	call writeDAT(nx,ny,rsigma(:,:,1),'r1-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,rsigma(:,:,1,nt),'r1.csv');
	call writeDAT(nx,ny,rsigma(:,:,2),'r2-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,rsigma(:,:,2,nt),'r2.csv');
	call writeDAT(nx,ny,rsigma(:,:,3),'r3-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,rsigma(:,:,3,nt),'r3.csv');

	call writeDAT(nx,ny,uxsigma(:,:,1),'u1-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,uxsigma(:,:,1,nt),'u1.csv');
	call writeDAT(nx,ny,uxsigma(:,:,2),'u2-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,uxsigma(:,:,2,nt),'u2.csv');
	call writeDAT(nx,ny,uxsigma(:,:,3),'u3-' // nameN // '-' // nameB // '.csv');
!call writeDAT(nx,ny,uxsigma(:,:,3,nt),'u3.csv');

enddo
enddo

	deallocate( fd )
	deallocate( fd_new )
	deallocate( f )
	deallocate( rsigma )
	deallocate( psigma )
	deallocate( uxsigma )
	deallocate( uysigma )
	deallocate( r )
	deallocate( p )
	deallocate( ux )
	deallocate( uy )
	deallocate( phi )
	deallocate( BCsigma )
	deallocate( MM )
	deallocate( Output )

end program MIXLBM

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!============================================================================
! Subroutines / Functions
!============================================================================
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine HydrodynamicMoments(f,rho,ux,uy)

	implicit none

! Inlet / Outlet
	real(8),dimension(0:8),intent(in) :: f
	real(8),intent(out) :: rho,ux,uy

! Local variables
	real(8),dimension(0:8) :: zitax,zitay
	real(8) :: qx,qy
	integer :: i

! zitax of D2Q9 lattice
    zitax(:) = (/0.0d0, 1.0d0, 0.0d0,-1.0d0, 0.0d0, 1.0d0,-1.0d0,-1.0d0, 1.0d0/)
! zitay of D2Q9 lattice
    zitay(:) = (/0.0d0, 0.0d0, 1.0d0, 0.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0,-1.0d0/)

	rho = 0.0d0;
	qx  = 0.0d0;
	qy  = 0.0d0;

	do i = 0,8,1
		rho = rho + f(i);
		qx  = qx  + zitax(i)*f(i);
		qy  = qy  + zitay(i)*f(i);
	enddo	

	ux = qx/rho;
	uy = qy/rho;

end subroutine HydrodynamicMoments

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine EquilibriumDistribution(rho,phi,ux,uy,feq)
!
! D2Q9 lattice

	implicit none

! Inlet / Outlet
	real(8),intent(in) :: rho,phi,ux,uy
	real(8),dimension(0:8),intent(out) :: feq

	feq(0) = 4.0d0/9.0d0*rho*(	(9.0d0-5.0d0*phi)/4.0d0							-3.0d0/2.0d0*(ux**2+uy**2) );

	feq(1) = 1.0d0/9.0d0*rho*(	phi	+3.0d0*(ux)		+9.0d0/2.0d0*(ux)**2		-3.0d0/2.0d0*(ux**2+uy**2) );
	feq(2) = 1.0d0/9.0d0*rho*(	phi	+3.0d0*(uy)		+9.0d0/2.0d0*(uy)**2		-3.0d0/2.0d0*(ux**2+uy**2) );
	feq(3) = 1.0d0/9.0d0*rho*(  phi	+3.0d0*(-ux)	+9.0d0/2.0d0*(-ux)**2		-3.0d0/2.0d0*(ux**2+uy**2) );
	feq(4) = 1.0d0/9.0d0*rho*(	phi	+3.0d0*(-uy)	+9.0d0/2.0d0*(-uy)**2		-3.0d0/2.0d0*(ux**2+uy**2) );

	feq(5) = 1.0d0/36.0d0*rho*(	phi+3.0d0*(ux+uy)	+9.0d0/2.0d0*(ux+uy)**2		-3.0d0/2.0d0*(ux**2+uy**2) );
	feq(6) = 1.0d0/36.0d0*rho*(	phi+3.0d0*(-ux+uy)	+9.0d0/2.0d0*(-ux+uy)**2	-3.0d0/2.0d0*(ux**2+uy**2) );
	feq(7) = 1.0d0/36.0d0*rho*(	phi+3.0d0*(-ux-uy)	+9.0d0/2.0d0*(-ux-uy)**2	-3.0d0/2.0d0*(ux**2+uy**2) );
	feq(8) = 1.0d0/36.0d0*rho*(	phi+3.0d0*(ux-uy)	+9.0d0/2.0d0*(ux-uy)**2		-3.0d0/2.0d0*(ux**2+uy**2) );

end subroutine EquilibriumDistribution

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine UpdateLatticeData(model,nx,ny,species,fd,i,j,phi,MMs,BCs,newfdPs,nN,nB)
!
! update distribution function of all the species in P(i,j)
! * for every direction k it identifies I(i+...,j+...)
! * if I is a point out of the domain, which border ?
! * if I is bulk, read f(k,s) in I and compute rsigma, ux sigma, uysigma in I
! * if I is border, then;
!	- read f(k,s) in P and assume rsigma, ux sigma, uysigma in I (wall / source)
!	- compute f(k,s) in I as feq(k,s)
! * compute barycentric r,p,ux,uy (Fick model)
! * compute barycentric rstar,pstar,uxstar,uystar (Maxwell-Stefan model)
! * compute local equilibrium feq(k) for species s
! * Collision Step: compute new value at the new time step f^+(k,s) in I
! * Streaming Step: take component BB(k) of f^+(k,s) in I and stream it in P
!	BB(k) is an operator which produce the direction opposite to direction k

	implicit none

! Internal Functions
	real(8) :: B,D,NU,XI

! Inlet / Outlet
	integer,intent(in) :: model,nx,ny,species,nN,nB
! distribution function fd(x,y,v)
	real(8),dimension(1:nx,1:ny,0:8,1:species),intent(in) :: fd
! generic point P(i,j)
	integer,intent(in) :: i,j
! rekionship between pressure and density
	real(8),dimension(1:species),intent(in) :: phi
! molecular weights
	real(8),dimension(1:species),intent(in) :: MMs
! boundary conditions:
! BC(i) = -2 Periodic
! BC(i) = -1 Neumann,	ux = 0 & uy = 0, wall
! BC(i) >  0 Dirichlet, rho = BC(i), imposed pressure
! 1, 0, -1, 0
! 0, 1, 0, -1
	real(8),dimension(1:4,1:species),intent(in) :: BCs

! new distribution function in P(i,j)
	real(8),dimension(0:8,1:species),intent(out) :: newfdPs

! Local variables
	real(8),dimension(0:8)					:: m,mI,mII,mP,mB,fI,fII,fP
	real(8),dimension(1:species)			:: rsigma,psigma,uxsigma,uysigma
	real(8),dimension(1:species)			:: rcxsig,rcysig,rcxsigP,rcysigP,gxrcxsig,gyrcysig
	real(8),dimension(1:species)			:: uxstar,uystar
	real(8),dimension(0:8,1:species)		:: f,ddsig,ddsigP
	real(8),dimension(1:species)			:: lambda,lamd,lamn,lamb
	real(8),dimension(0:8,0:8,1:species)	:: AA,iIpthAA
	real(8),dimension(1:species,1:species)	:: A,CHI
	real(8),dimension(1:species)			:: gjxsigma,gjysigma,jxsigma,jysigma,CHIsigma
	real(8),dimension(0:8)					:: feq,fcoll
	integer,dimension(0:8,1:2)				:: Incr
	integer,dimension(0:8)					:: BB
	integer,dimension(6)					:: IPARAM
	real(8),dimension(5)					:: RPARAM
	real(8),dimension(0:8)					:: tempV,tempV1,tempV2
	real(8),dimension(0:8,0:8,1:species)	:: tempM
	real(8)									:: temp

	real(8) :: r,p,ux,uy,vx,vy
	real(8) :: MM,PHImix
	integer :: k,ik,iik,iiik,s,vs,iI,jI,iII,jII
	integer :: BCdirection, IPATH = 1
	real(8) :: cx,cy,auxrcxsig,auxrcysig

	real(8) :: theta = 0.5d0

	Incr(:,1) = (/0, 1, 0, -1, 0, 1, -1, -1, 1/) ! x-axis
	Incr(:,2) = (/0, 0, 1, 0, -1, 1, 1, -1, -1/) ! y-axis

	BB(:) = (/0,3,4,1,2,7,8,5,6/) ! bounce-back

	do k=0,8,1
		iI = i + Incr(k,1);
		jI = j + Incr(k,2);

		do s=1,species,1

			rsigma(s)  = 0.0d0;
			psigma(s)  = 0.0d0;
			uxsigma(s) = 0.0d0;
			uysigma(s) = 0.0d0;

			! Boundary Conditions
			BCdirection = 0;		
 			if (jI>ny) then		! BC-2 
				BCdirection = 2;
			elseif (jI<1) then	! BC-4
				BCdirection = 4;
			elseif (iI>nx) then	! BC-1
				BCdirection = 1;
			elseif (iI<1) then	! BC-3
				BCdirection = 3;
			endif

			if (BCdirection == 0) then	! bulk domain
				f(:,s) = fd(iI,jI,:,s);
				call HydrodynamicMoments(f(:,s),rsigma(s),uxsigma(s),uysigma(s))
				psigma(s)  = rsigma(s)*phi(s)/3.0d0;
			else
				if (int(BCs(BCdirection,s))==-2) then		! periodic
					if (iI>nx) then 
						iI = 1;  endif
 					if (jI>ny) then 
						jI = 1;  endif
					if (iI<1)  then 
						iI = nx; endif
					if (jI<1)  then 
						jI = ny; endif
					f(:,s) = fd(iI,jI,:,s);
					call HydrodynamicMoments(f(:,s),rsigma(s),uxsigma(s),uysigma(s))		
				elseif (int(BCs(BCdirection,s))==-1) then	! wall: specify ux, uy
					f(:,s) = fd(i,j,:,s);
					call FtoM(f(:,s),m)				
					m(1) = 0.d0;
					m(2) = 0.d0;
					call MtoF(m,f(:,s))
					call HydrodynamicMoments(f(:,s),rsigma(s),uxsigma(s),uysigma(s))		
				else										! source: specify rho
					f(:,s) = fd(i,j,:,s);
					call FtoM(f(:,s),m)				
					psigma(s) = BCs(BCdirection,s);
					rsigma(s) = 3.0d0*psigma(s)/phi(s);
! ATTENTION
					m(3) = m(3) + (rsigma(s)-m(0))*phi(s)/3.0d0;
					m(4) = m(4) + (rsigma(s)-m(0))*phi(s)/3.0d0;
					m(0) = m(0) + (rsigma(s)-m(0));

					call MtoF(m,f(:,s))
					call HydrodynamicMoments(f(:,s),rsigma(s),uxsigma(s),uysigma(s))		
				endif
			endif

! FORCE
!	point I
			call FORCE(s,iI,jI,cx,cy)
			fI(:) = f(:,s);
			call FtoM(fI,mI)
			rcxsig(s) = mI(0)*cx;
			rcysig(s) = mI(0)*cy;

			gxrcxsig(s) = 0.0d0;
			gyrcysig(s) = 0.0d0;

			m = (/0.d0,rcxsig(s),rcysig(s),gxrcxsig(s),gyrcysig(s),0.d0,0.d0,0.d0,0.d0/);
			call MtoF(m,ddsig(:,s))

			if (k==0) then
				ddsigP(:,s) = ddsig(:,s);
				rcxsigP(s)  = rcxsig(s);
				rcysigP(s)  = rcysig(s);
			endif
! FORCE
		enddo

		r  = 0.0d0;
		p  = 0.0d0;
		do s=1,species,1
			r = r + rsigma(s);
			p = p + psigma(s);
		enddo

! Barycentric quantities
		ux = 0.0d0;
		uy = 0.0d0;
		do s=1,species,1
			ux = ux + ( rsigma(s)/r )*uxsigma(s);
			uy = uy + ( rsigma(s)/r )*uysigma(s);
		enddo

! Molar quantities
		vx = 0.0d0;
		vy = 0.0d0;
		do s=1,species,1
			vx = vx + ( psigma(s)/p )*uxsigma(s);
			vy = vy + ( psigma(s)/p )*uysigma(s);
		enddo

! Molecular weight for mixture
		temp = 0.0d0;
		do s=1,species,1
			temp = temp + ( rsigma(s)/r )/MMs(s);
		enddo
		MM = 1/temp;

! Modified quantities (Maxwell-Stefan)
		do s=1,species,1
			uxstar(s) = uxsigma(s);
			uystar(s) = uysigma(s);
			do vs=1,species,1
				CHI(s,vs) = (MM**2)/(MMs(s)*MMs(vs)) * B(MMs(s),MMs(vs),nB)/B(MMs(s),MMs(s),nB);
				uxstar(s) = uxstar(s) + CHI(s,vs)* ( rsigma(vs)/r )*(uxsigma(vs)-uxsigma(s));
				uystar(s) = uystar(s) + CHI(s,vs)* ( rsigma(vs)/r )*(uysigma(vs)-uysigma(s)); 	
			enddo
		enddo

		do s=1,species,1

! Select your model
			if (model==0) then
			! Ideal Fick model
				call EquilibriumDistribution(rsigma(s),phi(s),0.0d0,0.0d0,feq)
				lambda(s) = phi(s)/( 3.0d0*D(MMs(s)) );
				lamd(s) = lambda(s);
				lamn(s) = 1.d0/( 3.0d0*NU(species,r,MMs,rsigma,nN) );
				lamb(s) = ( 2.d0-phi(s) )/( 3.0d0*XI(species,r,MMs,rsigma) );
				call MRT(lamd(s),lamn(s),lamb(s),AA(:,:,s))

			elseif (model==1) then
			! Fick model
				call EquilibriumDistribution(rsigma(s),phi(s),vx,vy,feq)
				lambda(s) = phi(s)/( 3.0d0*D(MMs(s)) );
				lamd(s) = lambda(s);
				lamn(s) = 1.d0/( 3.0d0*NU(species,r,MMs,rsigma,nN) );
				lamb(s) = ( 2.d0-phi(s) )/( 3.0d0*XI(species,r,MMs,rsigma) );
				call MRT(lamd(s),lamn(s),lamb(s),AA(:,:,s))

			elseif (model==2) then
			! Maxwell-Stefan model
				call EquilibriumDistribution(rsigma(s),phi(s),uxstar(s),uystar(s),feq)
				lambda(s) = p*B(MMs(s),MMs(s),nB)/r;
				lamd(s) = lambda(s);
				lamn(s) = 1.d0/( 3.0d0*NU(species,r,MMs,rsigma,nN) );
				lamb(s) = ( 2.d0-phi(s) )/( 3.0d0*XI(species,r,MMs,rsigma) );
				call MRT(lamd(s),lamn(s),lamb(s),AA(:,:,s))

			endif

! TRANSFORMATION f(:,s) -> g(:,s)

			do ik=0,8,1
				tempV(ik) = 0.d0;
				do iik=0,8,1
					tempV(ik) = tempV(ik)-theta*AA(ik,iik,s)*( feq(iik)-f(iik,s) )
				enddo
			enddo
			f(:,s) = f(:,s)+tempV(:)-theta*ddsig(:,s);

			call TRANS(lamd(s),lamn(s),lamb(s),theta,iIpthAA(:,:,s))
			do ik=0,8,1
				do iik=0,8,1
					tempM(ik,iik,s) = 0.d0;
					do iiik=0,8,1
						tempM(ik,iik,s) = tempM(ik,iik,s)+AA(ik,iiik,s)*iIpthAA(iiik,iik,s);
					enddo
				enddo
			enddo
			AA(:,:,s) = tempM(:,:,s);

			do ik=0,8,1
				tempV(ik) = 0.d0;
				do iik=0,8,1
					tempV(ik) = tempV(ik)+iIpthAA(ik,iik,s)*ddsig(iik,s);
				enddo
			enddo
			ddsig(:,s) = tempV(:);

! Collision Step
			do ik=0,8,1
				tempV(ik) = 0.d0;
				do iik=0,8,1
					tempV(ik) = tempV(ik)+AA(ik,iik,s)*( feq(iik)-f(iik,s) );
				enddo
			enddo
			fcoll(:) = f(:,s)+tempV(:)+ddsig(:,s);

! Streaming Step
! only one component in BB(k) direction is streamed
! BB(k) is an operator which produce the direction
! opposite to direction k
			newfdPs(BB(k),s) = fcoll(BB(k));

		enddo

	enddo

	ddsig(:,:) = ddsigP(:,:);
	rcxsig(:)  = rcxsigP(:);
	rcysig(:)  = rcysigP(:);

! TRANSFORMATION g^+(:,s) -> f^+(:,s)
! moments of g^+(:,s)
	do s=1,species,1
		call HydrodynamicMoments(newfdPs(:,s),rsigma(s),uxsigma(s),uysigma(s))
		psigma(s)   = rsigma(s)*phi(s)/3.0d0;
! FORCE: using values at the previous time step
		gjxsigma(s) = rsigma(s)*uxsigma(s)+theta*rcxsig(s);
		gjysigma(s) = rsigma(s)*uysigma(s)+theta*rcysig(s);
	enddo

	r  = 0.0d0;
	p  = 0.0d0;
	do s=1,species,1
		r = r + rsigma(s);
		p = p + psigma(s);
	enddo	

! Barycentric quantities (conserved moments)
	ux = 0.0d0;
	uy = 0.0d0;
	do s=1,species,1
		ux = ux + ( rsigma(s)/r )*uxsigma(s);
		uy = uy + ( rsigma(s)/r )*uysigma(s);
	enddo

! Molecular weight for mixture
	temp = 0.0d0;
	do s=1,species,1
		temp = temp + ( rsigma(s)/r )/MMs(s);
	enddo
	MM = 1/temp;

! Modified quantities (Maxwell-Stefan)
	do s=1,species,1
		do vs=1,species,1
			CHI(s,vs) = (MM**2)/(MMs(s)*MMs(vs)) * B(MMs(s),MMs(vs),nB)/B(MMs(s),MMs(s),nB);
		enddo
	enddo

! Select your model (relaxation FREQUENCY)
	do s=1,species,1
		if (model==0) then
		! Ideal Fick model
			lambda(s) = phi(s)/( 3.0d0*D(MMs(s)) );
			lamd(s) = lambda(s);
			lamn(s) = 1.d0/( 3.0d0*NU(species,r,MMs,rsigma,nN) );
			lamb(s) = ( 2.d0-phi(s) )/( 3.0d0*XI(species,r,MMs,rsigma) );
			call MRT(lamd(s),lamn(s),lamb(s),AA(:,:,s))

		elseif (model==1) then
		! Fick model
			lambda(s) = phi(s)/( 3.0d0*D(MMs(s)) );
			lamd(s) = lambda(s);
			lamn(s) = 1.d0/( 3.0d0*NU(species,r,MMs,rsigma,nN) );
			lamb(s) = ( 2.d0-phi(s) )/( 3.0d0*XI(species,r,MMs,rsigma) );
			call MRT(lamd(s),lamn(s),lamb(s),AA(:,:,s))

		elseif (model==2) then
		! Maxwell-Stefan model
			lambda(s) = p*B(MMs(s),MMs(s),nB)/r;
			lamd(s) = lambda(s);
			lamn(s) = 1.d0/( 3.0d0*NU(species,r,MMs,rsigma,nN) );
			lamb(s) = ( 2.d0-phi(s) )/( 3.0d0*XI(species,r,MMs,rsigma) );
			call MRT(lamd(s),lamn(s),lamb(s),AA(:,:,s))

		endif

	enddo

! Select your model (EQUILIBRIUM)
	if (model==0) then
	! Ideal Fick model
		do s=1,species,1
			call EquilibriumDistribution(rsigma(s),phi(s),0.0d0,0.0d0,feq)

! MOD Ideal Fick model
			do ik=0,8,1
				tempV(ik) = 0.d0;
				do iik=0,8,1
					tempV(ik) = tempV(ik)+theta*AA(ik,iik,s)*feq(iik);
				enddo
			enddo
			tempV1(:) = newfdPs(:,s)+tempV(:)+theta*ddsig(:,s);

			do ik=0,8,1
				tempV2(ik) = 0.d0;
				do iik=0,8,1
					tempV2(ik) = tempV2(ik)+iIpthAA(ik,iik,s)*tempV1(iik);
				enddo
			enddo
			newfdPs(:,s) = tempV2(:); 
!			newfdPs(:,s) = ( newfdPs(:,s)+theta*lambda(s)*feq(:)+theta*ddsig(:,s) )/( 1.0d0+theta*lambda(s) );

		enddo
	elseif (model==1) then
	! Fick model
	! moments of f^+(:,s) -> uxsigma(s), uysigma(s)

		temp = 0.0d0;
		do s=1,species,1
			temp = temp+phi(s)*( rsigma(s)/r );
		enddo
		PHImix = temp;
		
		do s=1,species,1
			do vs=1,species,1
				if (s==vs) then
					A(s,vs) = 1.0d0+theta*lambda(s);
				else
					A(s,vs) = 0.0d0;
				endif
				A(s,vs) = A(s,vs)-theta*lambda(s)*( rsigma(s)/r )*( phi(s)/PHImix );
			enddo
		enddo

		call my_LSLXG(species,A,gjxsigma,IPATH,IPARAM,RPARAM,jxsigma)
		call my_LSLXG(species,A,gjysigma,IPATH,IPARAM,RPARAM,jysigma)
		do s=1,species,1
			uxsigma(s) = jxsigma(s)/rsigma(s);
			uysigma(s) = jysigma(s)/rsigma(s);
		enddo

	! Molar quantities
		vx = 0.0d0;
		vy = 0.0d0;
		do s=1,species,1
			vx = vx + ( psigma(s)/p )*uxsigma(s);
			vy = vy + ( psigma(s)/p )*uysigma(s);
		enddo

		do s=1,species,1
			call EquilibriumDistribution(rsigma(s),phi(s),vx,vy,feq)
			call TRANS(lamd(s),lamn(s),lamb(s),theta,iIpthAA(:,:,s))

! MOD Fick model
			do ik=0,8,1
				tempV(ik) = 0.d0;
				do iik=0,8,1
					tempV(ik) = tempV(ik)+theta*AA(ik,iik,s)*feq(iik);
				enddo
			enddo
			tempV1(:) = newfdPs(:,s)+tempV(:)+theta*ddsig(:,s);

			do ik=0,8,1
				tempV2(ik) = 0.d0;
				do iik=0,8,1
					tempV2(ik) = tempV2(ik)+iIpthAA(ik,iik,s)*tempV1(iik);
				enddo
			enddo
			newfdPs(:,s) = tempV2(:); 

		enddo

	elseif (model==2) then
	! Maxwell-Stefan model
	! moments of f^+(:,s) -> uxsigma(s), uysigma(s)

		do s=1,species,1
			temp = 0.0d0;
			do vs=1,species,1
				temp = temp+CHI(s,vs)*( rsigma(vs)/r );
			enddo
			CHIsigma(s) = temp;
		enddo
		
		do s=1,species,1
			do vs=1,species,1
				if (s==vs) then
					A(s,vs) = 1.0d0+theta*lambda(s)*CHIsigma(s);
				else
					A(s,vs) = 0.0d0;
				endif
				A(s,vs) = A(s,vs)-theta*lambda(s)*( rsigma(s)/r )*CHI(s,vs);
			enddo
		enddo

		call my_LSLXG(species,A,gjxsigma,IPATH,IPARAM,RPARAM,jxsigma)
		call my_LSLXG(species,A,gjysigma,IPATH,IPARAM,RPARAM,jysigma)
		do s=1,species,1
			uxsigma(s) = jxsigma(s)/rsigma(s);
			uysigma(s) = jysigma(s)/rsigma(s);
		enddo

	! Modified quantities (Maxwell-Stefan)
		do s=1,species,1
			uxstar(s) = uxsigma(s);
			uystar(s) = uysigma(s);

			do vs=1,species,1
				uxstar(s) = uxstar(s) + CHI(s,vs)*( rsigma(vs)/r )*(uxsigma(vs)-uxsigma(s));
				uystar(s) = uystar(s) + CHI(s,vs)*( rsigma(vs)/r )*(uysigma(vs)-uysigma(s)); 	
			enddo
		enddo

		do s=1,species,1
			call EquilibriumDistribution(rsigma(s),phi(s),uxstar(s),uystar(s),feq)
			call TRANS(lamd(s),lamn(s),lamb(s),theta,iIpthAA(:,:,s))

! MOD Maxwell-Stefan model
			do ik=0,8,1
				tempV(ik) = 0.d0;
				do iik=0,8,1
					tempV(ik) = tempV(ik)+theta*AA(ik,iik,s)*feq(iik);
				enddo
			enddo
			tempV1(:) = newfdPs(:,s)+tempV(:)+theta*ddsig(:,s);

			do ik=0,8,1
				tempV2(ik) = 0.d0;
				do iik=0,8,1
					tempV2(ik) = tempV2(ik)+iIpthAA(ik,iik,s)*tempV1(iik);
				enddo
			enddo
			newfdPs(:,s) = tempV2(:); 

		enddo
			
	endif

end subroutine UpdateLatticeData

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

real(8) function NU(species,rho,MMsigma,rhosigma,nN)
!
! Kinematic viscosity
! for ( sigma ) species

	implicit none

! Inlet/Outlet
	integer,intent(in) :: species,nN
	real(8),intent(in) :: rho
	real(8),dimension(1:species),intent(in) :: MMsigma
	real(8),dimension(1:species),intent(in) :: rhosigma

	real(8),dimension(1:20) :: k

	k(:)=(/1.0000,&
	    1.1708,&
		1.3707,&
		1.6048,&
		1.8789,&
		2.1998,&
		2.5754,&
		3.0153,&
		3.5302,&
		4.1331,&
		4.8390,&
		5.6654,&
		6.6329,&
		7.7657,&
		9.0919,&
	   10.6446,&
	   12.4625,&
	   14.5908,&
	   17.0826,&
	   20.0000/);

	NU = k(nN);

	return

end function NU

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

real(8) function B(m1,m2,nB)
!
! Binary resistivity of Maxwell-Stefan
! for couple (sigma , varsigma) species
! ** it depends on both species !! **

	implicit none

! Inlet/Outlet
	real(8),intent(in) :: m1,m2
	integer,intent(in) :: nB
	real(8) :: D
	real(8),dimension(1:20) :: k

	k(:)=(/0.500000000000000,&
		0.601274801898418,&
		0.723062774795962,&
		0.869518853351124,&
		1.045639552591273,&
		1.257433429682934,&
		1.512126072666108,&
		1.818406609575491,&
		2.186724147886553,&
		2.629644257653945,&
		3.162277660168372,&
		3.802795747331057,&
		4.573050519273251,&
		5.499320090094956,&
		6.613205195495662,&
		7.952707287670479,&
		9.563524997900334,&
		11.500613197126169,&
		13.830057843624722,&
		16.631330580338215/);

	B = k(nB)*10.0d0*(1.0d0/m1+1.0d0/m2)**(-0.5d0);

! Equivalent to Fick
!	B = 1.0d0/D(m1);

	return

end function B

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

real(8) function D(m)
!
! Diffusivity of Fick
! for ( sigma ) species
! ** it depends on one species !! **

	implicit none

! Inlet/Outlet
	real(8),intent(in) :: m
	real(8),dimension(1:20) :: k

	k(:)=(/0.010000000000000,&
	   0.013645831365889,&
	   0.018620871366629,&
	   0.025409727055493,&
	   0.034673685045253,&
	   0.047315125896148,&
	   0.064565422903466,&
	   0.088104887300801,&
	   0.120226443461741,&
	   0.164058977319954,&
	   0.223872113856834,&
	   0.305492111321552,&
	   0.416869383470336,&
	   0.568852930843842,&
	   0.776247116628692,&
	   1.059253725177290,&
	   1.445439770745929,&
	   1.972422736114855,&
	   2.691534803926918,&
	   3.672823004980850/);

	D = k(5)*0.2d0/m;

	return

end function D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

real(8) function XI(species,rho,MMsigma,rhosigma)
!
! Bulk viscosity
! for ( sigma ) species

	implicit none

! Inlet/Outlet
	integer,intent(in) :: species
	real(8),intent(in) :: rho
	real(8),dimension(1:species),intent(in) :: MMsigma
	real(8),dimension(1:species),intent(in) :: rhosigma

	XI = 0.4d0;

	return

end function XI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine FORCE(s,i,j,cx,cy)

!  SUBROUTINE: MFORCE
!
!  PURPOSE:	Forcing term

!	Ingressi / uscite
	implicit none  

    integer, intent(in) :: s,i,j
    real(8), intent(out) :: cx,cy
	real(8) :: ax,ay,bx,by

	ax = 0.d0;
	ay = 0.d0;

	bx = 0.d0;
	by = 0.d0;

	cx = ax+bx;
	cy = ay+by;

end subroutine FORCE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine writeDAT(nx,ny,M,name)

	implicit none

! Inlet
	integer,intent(in) :: nx,ny
	real(8),dimension(1:nx,1:ny),intent(in) :: M
	character(12),intent(in) :: name

! Local variables
	integer :: i,j
	integer :: IO_Stat_Number1=-1

	open(UNIT=1,FILE=name,IOSTAT=IO_Stat_Number1)

	if (IO_Stat_Number1==0) then
!		print *,'Saving ',name,' ...'

		do j=ny,1,-1
			do i=1,nx,1
				if (i<nx) then
					write(1,'($F",")') M(i,j)
				else
					write(1,'($F)') M(i,j)
				endif
			enddo
			write(1,*)
		enddo

!		print *, '...done' 
	endif
	
	close(UNIT=1)

end subroutine writeDAT

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine my_LSLXG(n,A,b,IPATH,IPARAM,RPARAM,x)

!  SUBROUTINE: My_LSLXG
!
!  PURPOSE:	Versione riscritta e non ottimizzata della subroutine 
!			LSLXG della libreria IMSL
	
!	Ingressi / uscite
	implicit none

	integer, intent(in) :: n
	real(8),dimension(1:n,1:n),intent(in) :: A
	real(8),dimension(1:n),intent(in) :: b

	integer,intent(in) :: IPATH
	integer,dimension(6),intent(in) :: IPARAM
	real(8),dimension(5),intent(in) :: RPARAM

	real(8),dimension(1:n),intent(out) :: x

!	Variabili locali
    real(8), dimension (n) :: s,b_mod
    integer, dimension(n) :: l
	integer :: i,j

!	Calcolo
!	Fattorizzazone LU
    call gauss(n,A,l,s)

!	print*, "Matrice dopo Gauss"
!   do i=1,n
!		print*, (A(i,j), j=1,n)
!	end do
!	pause

!	Risoluzione sistema lineare
	b_mod=b
    call solve(n,A,l,b_mod,x)
!	print*, "Soluzione x(i) ", (x(i), i=1,n)
!	pause

end subroutine my_LSLXG

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine gauss(n,a,l,s)

!  SUBROUTINE: Gauss
!
!  PURPOSE:	Fattorizzazione della matrice mediante metodo di eliminazione
!			di Gauss 
	    
!	Ingressi / uscite
	implicit none
	integer, intent(in) :: n
    real(8), dimension(1:n,1:n), intent(inout) :: a
    real(8), dimension(n), intent(out) :: s
    integer, dimension(n), intent(out) :: l
    
!	Variabili locali
	integer :: i,j,k
	real(8) :: smax, rmax, xmult, r, lk, sum

!	Corpo
	do  i = 1,n
		l(i) = i  
		smax = 0.0
		do j = 1,n
			smax = dmax1(smax,abs(a(i,j)))
        end do
		s(i) = smax 
    end do 
    
	do  k = 1,n-1
		rmax = -1.0
        do i = k,n
			r = abs(a(l(i),k))/s(l(i))  
			if(r <= rmax)  exit    
			j = i   
			rmax = r
		end do 
		lk = l(j) 
        l(j) = l(k) 
        l(k) = lk  
        do i = k+1,n      
			xmult = a(l(i),k)/a(lk,k)   
			do j = k+1,n    
				a(l(i),j) = a(l(i),j) - xmult*a(lk,j) 
			end do 
			a(l(i),k) = xmult 
        end do 
     end do   

end subroutine gauss 
  
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine solve(n,a,l,b,x)

!  SUBROUTINE: Solve
!
!  PURPOSE:	Risoluzione del sistema lineare

!	Ingressi / uscite
	implicit none  

    integer, intent(in) :: n
    real(8), dimension (1:n,1:n), intent(in) :: a
    integer, dimension(n), intent(in) :: l
    real(8), dimension (1:n), intent(inout) :: b
    real(8), dimension (n), intent(out) :: x

!	Variabili locali
    integer :: k, i, j
	real(8) :: sum

!	Corpo
	do k = 1,n-1
         do i = k+1,n      
			b(l(i)) = b(l(i)) - a(l(i),k)*b(l(k)) 
         end do 
	end do 
	x(n) = b(l(n))/a(l(n),n)

    do i = n-1,1,-1     
		sum = b(l(i))       
        do j = i+1,n      
			sum = sum - a(l(i),j)*x(j)  
        end do 
	x(i) = sum/a(l(i),i)
    end do 

end subroutine solve

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine MRT(lamd,lamn,lamb,A)

!  SUBROUTINE: MRT
!
!  PURPOSE:	Collisional matrix

!	Ingressi / uscite
	implicit none  

    real(8), intent(in) :: lamd,lamn,lamb
    real(8), dimension (0:8,0:8), intent(out) :: A
 
	real(8) :: lam0, lam34
! ATTENTION
	lam0  = 0.d0;
	lam34 = 1.d0;
	
	A(0,0) = lam0;      A(0,1) = lam0-lamb;      A(0,2) = lam0-lamb;
	A(0,3) = lam0-lamb;      A(0,4) = lam0-lamb;      A(0,5) = lam0-2.d0*lamb+lam34;
	A(0,6) = lam0-2.d0*lamb+lam34;      A(0,7) = lam0-2.d0*lamb+lam34;      A(0,8) = lam0-2.d0*lamb+lam34;

	A(1,0) = 0.d0;      A(1,1) = lamd/2.d0+lamb/4.d0+lamn/4.d0;      A(1,2) = lamb/4.d0-lamn/4.d0;
	A(1,3) = -lamd/2.d0+lamb/4.d0+lamn/4.d0;      A(1,4) = lamb/4.d0-lamn/4.d0;      A(1,5) = lamd/2.d0+lamb/2.d0-lam34;
	A(1,6) = -lamd/2.d0+lamb/2.d0;      A(1,7) = -lamd/2.d0+lamb/2.d0;      A(1,8) = lamd/2.d0+lamb/2.d0-lam34;

	A(2,0) = 0.d0;      A(2,1) = lamb/4.d0-lamn/4.d0;      A(2,2) = lamd/2.d0+lamb/4.d0+lamn/4.d0;
	A(2,3) = lamb/4.d0-lamn/4.d0;      A(2,4) = -lamd/2.d0+lamb/4.d0+lamn/4.d0;      A(2,5) = lamd/2.d0+lamb/2.d0-lam34;
	A(2,6) = lamd/2.d0+lamb/2.d0-lam34;      A(2,7) = -lamd/2.d0+lamb/2.d0;      A(2,8) = -lamd/2.d0+lamb/2.d0;

	A(3,0) = 0.d0;      A(3,1) = -lamd/2.d0+lamb/4.d0+lamn/4.d0;      A(3,2) = lamb/4.d0-lamn/4.d0;
	A(3,3) = lamd/2.d0+lamb/4.d0+lamn/4.d0;      A(3,4) = lamb/4.d0-lamn/4.d0;      A(3,5) = -lamd/2.d0+lamb/2.d0;
	A(3,6) = lamd/2.d0+lamb/2.d0-lam34;      A(3,7) = lamd/2.d0+lamb/2.d0-lam34;      A(3,8) = -lamd/2.d0+lamb/2.d0;

	A(4,0) = 0.d0;      A(4,1) = lamb/4.d0-lamn/4.d0;      A(4,2) = -lamd/2.d0+lamb/4.d0+lamn/4.d0;
	A(4,3) = lamb/4.d0-lamn/4.d0;      A(4,4) = lamd/2.d0+lamb/4.d0+lamn/4.d0;      A(4,5) = -lamd/2.d0+lamb/2.d0;
	A(4,6) = -lamd/2.d0+lamb/2.d0;      A(4,7) = lamd/2.d0+lamb/2.d0-lam34;      A(4,8) = lamd/2.d0+lamb/2.d0-lam34;

	A(5,0) = 0.d0;      A(5,1) = 0.d0;      A(5,2) = 0.d0;
	A(5,3) = 0.d0;      A(5,4) = 0.d0;      A(5,5) = lamn/4.d0+3.d0/4.d0*lam34;
	A(5,6) = -lamn/4.d0+lam34/4.d0;      A(5,7) = lamn/4.d0-lam34/4.d0;      A(5,8) = -lamn/4.d0+lam34/4.d0;

	A(6,0) = 0.d0;      A(6,1) = 0.d0;      A(6,2) = 0.d0;
	A(6,3) = 0.d0;      A(6,4) = 0.d0;      A(6,5) = -lamn/4.d0+lam34/4.d0;
	A(6,6) = lamn/4.d0+3.d0/4.d0*lam34;      A(6,7) = -lamn/4.d0+lam34/4.d0;      A(6,8) = lamn/4.d0-lam34/4.d0;

	A(7,0) = 0.d0;      A(7,1) = 0.d0;      A(7,2) = 0.d0;
	A(7,3) = 0.d0;      A(7,4) = 0.d0;      A(7,5) = lamn/4.d0-lam34/4.d0;
	A(7,6) = -lamn/4.d0+lam34/4.d0;      A(7,7) = lamn/4.d0+3.d0/4.d0*lam34;      A(7,8) = -lamn/4.d0+lam34/4.d0;

	A(8,0) = 0.d0;      A(8,1) = 0.d0;      A(8,2) = 0.d0;
	A(8,3) = 0.d0;      A(8,4) = 0.d0;      A(8,5) = -lamn/4.d0+lam34/4.d0;
	A(8,6) = lamn/4.d0-lam34/4.d0;      A(8,7) = -lamn/4.d0+lam34/4.d0;      A(8,8) = lamn/4.d0+3.d0/4.d0*lam34;

end subroutine MRT

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine TRANS(lamd,lamn,lamb,theta,G)

!  SUBROUTINE: TRANS
!
!  PURPOSE:	(I+theta A)^-1

!	Ingressi / uscite
	implicit none  

    real(8), intent(in) :: lamd,lamn,lamb,theta
    real(8), dimension (0:8,0:8), intent(out) :: G
	real(8) :: lam0, lam34
! ATTENTION
	lam0  = 0.d0;
	lam34 = 1.d0;
 
	G(0,0) = 1/(1.d0+theta*lam0);
	G(0,1) = -theta*(lam0-lamb)/(theta*lamb+1.d0)/(1.d0+theta*lam0);
	G(0,2) = -theta*(lam0-lamb)/(theta*lamb+1.d0)/(1.d0+theta*lam0);
	G(0,3) = -theta*(lam0-lamb)/(theta*lamb+1.d0)/(1.d0+theta*lam0);
	G(0,4) = -theta*(lam0-lamb)/(theta*lamb+1.d0)/(1.d0+theta*lam0);
	G(0,5) = -(-theta*lamb*lam34+2.d0*theta*lam0*lam34-theta*lam0*lamb+lam0-2.d0*lamb+lam34)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(1.d0+theta*lam0);
	G(0,6) = -(-theta*lamb*lam34+2.d0*theta*lam0*lam34-theta*lam0*lamb+lam0-2.d0*lamb+lam34)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(1.d0+theta*lam0);
	G(0,7) = -(-theta*lamb*lam34+2.d0*theta*lam0*lam34-theta*lam0*lamb+lam0-2.d0*lamb+lam34)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(1.d0+theta*lam0);
	G(0,8) = -(-theta*lamb*lam34+2.d0*theta*lam0*lam34-theta*lam0*lamb+lam0-2.d0*lamb+lam34)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(1.d0+theta*lam0);

	G(1,0) = 0.d0;
	G(1,1) = (theta*theta*lamd*lamn+theta*theta*lamb*lamd+2.d0*theta*theta*lamn*lamb+2.d0*theta*lamd+3.d0*theta*lamn+3.d0*theta*lamb+4.d0)/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(1,2) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(1,3) = (theta*lamd*lamn+theta*lamb*lamd-2.d0*theta*lamn*lamb-lamn-lamb+2.d0*lamd)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(1,4) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(1,5) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;      G(1,6) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(1,7) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(1,8) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;

	G(2,0) = 0.d0;
	G(2,1) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(2,2) = (theta*theta*lamd*lamn+theta*theta*lamb*lamd+2.d0*theta*theta*lamn*lamb+2.d0*theta*lamd+3.d0*theta*lamn+3.d0*theta*lamb+4.d0)/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(2,3) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(2,4) = (theta*lamd*lamn+theta*lamb*lamd-2.d0*theta*lamn*lamb-lamn-lamb+2.d0*lamd)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(2,5) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(2,6) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(2,7) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(2,8) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;

	G(3,0) = 0.d0;
	G(3,1) = (theta*lamd*lamn+theta*lamb*lamd-2.d0*theta*lamn*lamb-lamn-lamb+2.d0*lamd)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(3,2) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(3,3) = (theta*theta*lamd*lamn+theta*theta*lamb*lamd+2.d0*theta*theta*lamn*lamb+2.d0*theta*lamd+3.d0*theta*lamn+3.d0*theta*lamb+4.d0)/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(3,4) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(3,5) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(3,6) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(3,7) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(3,8) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;

	G(4,0) = 0.d0;
	G(4,1) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(4,2) = (theta*lamd*lamn+theta*lamb*lamd-2.d0*theta*lamn*lamb-lamn-lamb+2.d0*lamd)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(4,3) = -(lamb-lamn)*theta/(theta*lamn+1.d0)/(theta*lamb+1.d0)/4.d0;
	G(4,4) = (theta*theta*lamd*lamn+theta*theta*lamb*lamd+2.d0*theta*theta*lamn*lamb+2.d0*theta*lamd+3.d0*theta*lamn+3.d0*theta*lamb+4.d0)/(theta*lamn+1.d0)/(theta*lamb+1.d0)/(theta*lamd+1.d0)/4.d0;
	G(4,5) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(4,6) = (lamd-lamb)*theta/(theta*lamb+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(4,7) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;
	G(4,8) = (theta*lamd*lam34-2.d0*theta*lamb*lamd+theta*lamb*lam34-lamd+2.d0*lam34-lamb)*theta/(theta*lamb+1.d0)/(theta*lam34+1.d0)/(theta*lamd+1.d0)/2.d0;

	G(5,0) = 0.d0;
	G(5,1) = 0.d0;
	G(5,2) = 0.d0;
	G(5,3) = 0.d0;
	G(5,4) = 0.d0;
	G(5,5) = (theta*lam34+3.d0*theta*lamn+4.d0)/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(5,6) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(5,7) = (lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(5,8) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;

	G(6,0) = 0.d0;
	G(6,1) = 0.d0;
	G(6,2) = 0.d0;
	G(6,3) = 0.d0;
	G(6,4) = 0.d0;
	G(6,5) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(6,6) = (theta*lam34+3.d0*theta*lamn+4.d0)/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(6,7) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(6,8) = (lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;

	G(7,0) = 0.d0;
	G(7,1) = 0.d0;
	G(7,2) = 0.d0;
	G(7,3) = 0.d0;
	G(7,4) = 0.d0;
	G(7,5) = (lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(7,6) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(7,7) = (theta*lam34+3.d0*theta*lamn+4.d0)/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(7,8) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;

	G(8,0) = 0.d0;
	G(8,1) = 0.d0;
	G(8,2) = 0.d0;
	G(8,3) = 0.d0;
	G(8,4) = 0.d0;
	G(8,5) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(8,6) = (lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(8,7) = -(lam34-lamn)*theta/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.d0;
	G(8,8) = (theta*lam34+3.d0*theta*lamn+4.d0)/(theta*lamn+1.d0)/(theta*lam34+1.d0)/4.0;

end subroutine TRANS

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine FtoM(fin,mout)

!  SUBROUTINE: FtoM
!
!  PURPOSE:	calculation of moments

!	Ingressi / uscite
	implicit none  

    real(8), dimension (0:8), intent(in) :: fin
    real(8), dimension (0:8), intent(out) :: mout
	real(8), dimension (0:8,0:8) :: M
	integer :: i,j

	M(0,0) = 1.d0;      M(0,1) = 1.d0;      M(0,2) = 1.d0;      M(0,3) = 1.d0;      M(0,4) = 1.d0;
	M(0,5) = 1.d0;      M(0,6) = 1.d0;      M(0,7) = 1.d0;      M(0,8) = 1.d0;

	M(1,0) = 0.d0;      M(1,1) = 1.d0;      M(1,2) = 0.d0;      M(1,3) = -1.d0;      M(1,4) = 0.d0;
	M(1,5) = 1.d0;      M(1,6) = -1.d0;      M(1,7) = -1.d0;      M(1,8) = 1.d0;

	M(2,0) = 0.d0;      M(2,1) = 0.d0;      M(2,2) = 1.d0;      M(2,3) = 0.d0;      M(2,4) = -1.d0;
	M(2,5) = 1.d0;      M(2,6) = 1.d0;      M(2,7) = -1.d0;      M(2,8) = -1.d0;

	M(3,0) = 0.d0;      M(3,1) = 1.d0;      M(3,2) = 0.d0;      M(3,3) = 1.d0;      M(3,4) = 0.d0;
	M(3,5) = 1.d0;      M(3,6) = 1.d0;      M(3,7) = 1.d0;      M(3,8) = 1.d0;

	M(4,0) = 0.d0;      M(4,1) = 0.d0;      M(4,2) = 1.d0;      M(4,3) = 0.d0;      M(4,4) = 1.d0;
	M(4,5) = 1.d0;      M(4,6) = 1.d0;      M(4,7) = 1.d0;      M(4,8) = 1.d0;

	M(5,0) = 0.d0;      M(5,1) = 0.d0;      M(5,2) = 0.d0;      M(5,3) = 0.d0;      M(5,4) = 0.d0;
	M(5,5) = 1.d0;      M(5,6) = -1.d0;      M(5,7) = 1.d0;      M(5,8) = -1.d0;

	M(6,0) = 0.d0;      M(6,1) = 0.d0;      M(6,2) = 0.d0;      M(6,3) = 0.d0;      M(6,4) = 0.d0;
	M(6,5) = 1.d0;      M(6,6) = -1.d0;      M(6,7) = -1.d0;      M(6,8) = 1.d0;

	M(7,0) = 0.d0;      M(7,1) = 0.d0;      M(7,2) = 0.d0;      M(7,3) = 0.d0;      M(7,4) = 0.d0;
	M(7,5) = 1.d0;      M(7,6) = 1.d0;      M(7,7) = -1.d0;      M(7,8) = -1.d0;

	M(8,0) = 0.d0;      M(8,1) = 0.d0;      M(8,2) = 0.d0;      M(8,3) = 0.d0;      M(8,4) = 0.d0;
	M(8,5) = 1.d0;      M(8,6) = 1.d0;      M(8,7) = 1.d0;      M(8,8) = 1.d0;

	do i = 0,8,1
		mout(i) = 0.d0;
		do j = 0,8,1
			mout(i) = mout(i) + M(i,j)*fin(j);
		enddo
	enddo

end subroutine FtoM

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine MtoF(minl,fout)

!  SUBROUTINE: MtoF
!
!  PURPOSE:	calculation of distribution

!	Ingressi / uscite
	implicit none  

    real(8), dimension (0:8), intent(in) :: minl
    real(8), dimension (0:8), intent(out) :: fout
	real(8), dimension (0:8,0:8) :: iM
	integer :: i,j

	iM(0,0) = 1.d0;      iM(0,1) = 0.d0;      iM(0,2) = 0.d0;      iM(0,3) = -1.d0;      iM(0,4) = -1.d0;
	iM(0,5) = 0.d0;      iM(0,6) = 0.d0;      iM(0,7) = 0.d0;      iM(0,8) = 1.d0;

	iM(1,0) = 0.d0;      iM(1,1) = 1.d0/2.d0;      iM(1,2) = 0.d0;      iM(1,3) = 1.d0/2.d0;      iM(1,4) = 0.d0;
	iM(1,5) = 0.d0;      iM(1,6) = -1.d0/2.d0;      iM(1,7) = 0.d0;      iM(1,8) = -1.d0/2.d0;

	iM(2,0) = 0.d0;      iM(2,1) = 0.d0;      iM(2,2) = 1.d0/2.d0;      iM(2,3) = 0.d0;      iM(2,4) = 1.d0/2.d0;
	iM(2,5) = 0.d0;      iM(2,6) = 0.d0;      iM(2,7) = -1.d0/2.d0;      iM(2,8) = -1.d0/2.d0;

	iM(3,0) = 0.d0;      iM(3,1) = -1.d0/2.d0;      iM(3,2) = 0.d0;      iM(3,3) = 1.d0/2.d0;      iM(3,4) = 0.d0;
	iM(3,5) = 0.d0;      iM(3,6) = 1.d0/2.d0;      iM(3,7) = 0.d0;      iM(3,8) = -1.d0/2.d0;

	iM(4,0) = 0.d0;      iM(4,1) = 0.d0;      iM(4,2) = -1.d0/2.d0;      iM(4,3) = 0.d0;      iM(4,4) = 1.d0/2.d0;
	iM(4,5) = 0.d0;      iM(4,6) = 0.d0;      iM(4,7) = 1.d0/2.d0;      iM(4,8) = -1.d0/2.d0;

	iM(5,0) = 0.d0;      iM(5,1) = 0.d0;      iM(5,2) = 0.d0;      iM(5,3) = 0.d0;      iM(5,4) = 0.d0;
	iM(5,5) = 1.d0/4.d0;      iM(5,6) = 1.d0/4.d0;      iM(5,7) = 1.d0/4.d0;      iM(5,8) = 1.d0/4.d0;

	iM(6,0) = 0.d0;      iM(6,1) = 0.d0;      iM(6,2) = 0.d0;      iM(6,3) = 0.d0;      iM(6,4) = 0.d0;
	iM(6,5) = -1.d0/4.d0;      iM(6,6) = -1.d0/4.d0;      iM(6,7) = 1.d0/4.d0;      iM(6,8) = 1.d0/4.d0;

	iM(7,0) = 0.d0;      iM(7,1) = 0.d0;      iM(7,2) = 0.d0;      iM(7,3) = 0.d0;      iM(7,4) = 0.d0;
	iM(7,5) = 1.d0/4.d0;      iM(7,6) = -1.d0/4.d0;      iM(7,7) = -1.d0/4.d0;      iM(7,8) = 1.d0/4.d0;

	iM(8,0) = 0.d0;      iM(8,1) = 0.d0;      iM(8,2) = 0.d0;      iM(8,3) = 0.d0;      iM(8,4) = 0.d0;
	iM(8,5) = -1.d0/4.d0;      iM(8,6) = 1.d0/4.d0;      iM(8,7) = -1.d0/4.d0;      iM(8,8) = 1.d0/4.d0;

	do i = 0,8,1
		fout(i) = 0.d0;
		do j = 0,8,1
			fout(i) = fout(i) + iM(i,j)*minl(j);
		enddo
	enddo

end subroutine MtoF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

real(8) function mean(vect,n)
!
! Kinematic viscosity
! for ( sigma ) species

	implicit none

! Inlet/Outlet
	integer,intent(in) :: n
	real(8),dimension(1:n),intent(in) :: vect
	real(8) :: temp
	integer :: i

	temp = 0.0d0;
	do i=1,n,1
		temp = temp + vect(i);
	enddo
	mean = temp/real(n);

	return

end function mean