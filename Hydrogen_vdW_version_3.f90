!================================================================
!     Written by Peter Szabo
!     14. 05. 2019.
!     Luxembourg University of Technology
!     peter88szabo@gmail.com      
!================================================================
      module mainvar

      double precision,parameter :: pi=4.0*atan(1.d0)

!     degree to radian conversion: 
      double precision,parameter :: d2r=pi/180.d0


!     ------------------------------------------
!     Constants and unit conversion factors:
!     [bohr]*c1=[Angstrom]
      double precision,parameter :: c1=0.5291772d0

!     [Hartree]*c2=[kcal/mol]
      double precision,parameter :: c2=627.51d0

!     [g/mol]*c3=[electron mass unit]
      double precision,parameter :: c3=1838.6836605d0

!     [Hartree]*c4=[eV]
      double precision,parameter :: c4=27.2114

!     [Hartree]*c5=[cm-1]
      double precision,parameter :: c5=219474.d0

!     [femto-sec]*c6=[time in au]
      double precision,parameter :: c6=41.341105

!     [Hartree]*c7=[kJ/mol]
      double precision,parameter :: c7=2625.5d0


!     Gas constant:
      double precision, parameter :: Rgas=8.3144598/1000.d0/c7 !Hartree/K
!     ------------------------------------------

      double precision, allocatable :: TRP_bnd(:,:),TRP_free(:,:)


      endmodule
!==============================================================



      program asdasd
      use mainvar
      implicit none
      integer :: i,j,N,Lb,Lf,iscal
      integer :: iE,nE,ir,nr,Nmax
      double precision :: x,rmin,rmax,dr,Z
      double precision :: Emin,Emax,E,dE
      double precision :: Eprint,umass
      double precision :: RINIT,RFINAL,H_bound_radial
      double precision :: Psi_bnd_gr,Psi_bnd_exc,Psi_free,TRM
      double precision :: EvdW_bound,dpol_bnd,EvdW_free,dpol_free
      double precision, allocatable :: func(:),r(:)

      read(5,*) Nmax, Z, Lb, Lf
      read(5,*) rmin,rmax,dr 
      read(5,*) Emin,Emax,dE, Eprint
      read(5,*) umass, iscal


      if(Nmax > 32) stop "too big Nmax, factorial can not be calculated"
      if(iscal == 1) umass=umass*c3

      write(*,*) "umass =",umass

      nr=nint((rmax-rmin)/dr)+1
      nE=nint((Emax-Emin)/dE)+1

      write(6,*) "nr = ",nr, "nE = ",nE


      allocate(func(nr))
      allocate(r(nr))
      allocate(TRP_free(Nmax,nE))
      allocate(TRP_bnd(Nmax,Nmax))

      do ir=1,nr
         r(ir)=rmin+(ir-1)*dr
      enddo


      open(99,file="Wavefunctions_only_bounds.dat")
      open(133,file="Wavefunctions_bound_free.dat")
      open(144,file="Transition_probabilities_bound_free.dat")
      open(122,file="Transition_probabilities_bound_bound.dat")

     !---------------------------------------------------------------------------
     ! Bound-to-bound transitions
     !---------------------------------------------------------------------------
     ! Only L --> L+1 transitions are considered
     !---------------------------------------------------------------------------
      do n=2,Nmax
         do ir=1,nr
            Psi_bnd_gr=H_bound_radial(1,Z,0,r(ir))*r(ir)

            Psi_bnd_exc=H_bound_radial(n,Z,1,r(ir))*r(ir)
            func(ir)=Psi_bnd_gr*Z*r(ir)*Psi_bnd_exc
            write(99,*) r(ir),Psi_bnd_gr,Psi_bnd_exc,n
         enddo

         call simpson(nr,1,nr,dr,func,TRM)

         TRP_bnd(1,n)=TRM*TRM
         write(122,*) n,TRP_bnd(1,n),sum(func)
         call flush(122)
      enddo
     !---------------------------------------------------------------------------



     !---------------------------------------------------------------------------
     ! Free-to-bound transitions
     !---------------------------------------------------------------------------
      do iE=1,nE

         E=Emin+(iE-1)*dE

         do ir=1,nr
            x=Z*r(ir)

            Psi_bnd_gr=H_bound_radial(1,Z,0,r(ir))*r(ir) !to be normalized
           !Psi_bnd=RINIT(N,Lb,x,Z)*r(ir) !to be normalized
            Psi_free=RFINAL(E,Lf,x,Z)*r(ir) !to be normalized

           !<Psi_free|D|Psi_bnd>:
            func(ir)=Psi_free*x*Psi_bnd_gr

            if(dabs(Eprint-E) <= 1d-5) then
               write(133,*) r(ir),func(ir),Psi_bnd_gr,Psi_free
            endif

         enddo

         call simpson(nr,1,nr,dr,func,TRM)

        !Dipole transition probability from Bound to Free states:
         TRP_free(1,iE)=TRM*TRM

         write(144,*) E,TRP_free(1,iE)
         write(*,*) E,TRP_free(1,iE)

      enddo
     !---------------------------------------------------------------------------




      call E_pert_bound(Nmax,Z,1.d0,EvdW_bound,dpol_bnd)
      call E_pert_free(nE,Emin,Emax,dE,Z,umass,1.d0,EvdW_free,dpol_free)

      write(*,*)"================================="
      write(*,"(28x,2a15)") "  dEvdW_bound ","   alpha_bound: "
      write(*,"(28x,2f15.5)") EvdW_bound,dpol_bnd
      write(*,"(28x,2a15)") "  dEvdW_free  ","   alpha_free:  "
      write(*,"(28x,2f15.5)") EvdW_free,dpol_free
      write(*,*) "-----------------------------------------"
      write(*,"(a28,15x,f15.5)") 'total polarizability in au:', dpol_free+dpol_bnd
      write(*,*)"================================="

      end program





!===========================================================================================
      subroutine E_pert_bound(Nmax,Z,rr,EvdW_bound,dpol_bnd)
!===========================================================================================
!     Second order perturbative van der Waals
!     correction to the ground state energy of the system
!     
!     Only bound-to-bound contribution are
!     considered here
!
!     See in Landau-Lifsitz pp. 344
!     the excercise at the end of the chapter:
!     "The interaction of atoms at large distances"
!===========================================================================================
      use mainvar
      implicit none
      integer :: Nmax,i,j,L1,L2
      double precision :: Ebound
      double precision :: Z,rr,rr3,rr6
      double precision :: Cang, denom,summ,EvdW_bound
      double precision :: dpol_bnd,fac,omega

     !Angular correction factor
      Cang=2.d0/3.d0
!
      omega=0.d0

      summ=0.d0
      dpol_bnd=0.d0
      do i=2,Nmax
         do j=2,Nmax

            denom=Ebound(i,Z)+Ebound(j,Z)-2.d0*Ebound(1,Z)
            summ=summ + TRP_bnd(1,i)*TRP_bnd(1,j)/denom

         enddo

        !fac=(E_bnd(i,L1)-E_bnd(1,L1))/((E_bnd(i,L1)-E_bnd(1,L1))**2 + omega**2)
         fac=1.0d0/(Ebound(i,Z)-Ebound(1,Z)-omega)

         dpol_bnd=dpol_bnd+fac*TRP_bnd(1,i)
      enddo

     !due to the px and py degenerate states
     !pz has no contribution, see Landaul-Lifsitz Eq 29.7
      dpol_bnd=2.d0*dpol_bnd*Cang
 
      rr3=rr*rr*rr
      rr6=rr3*rr3
 
      EvdW_bound=-2.d0*summ*Cang/rr6

      end subroutine
!===========================================================================================



!===========================================================================================
      subroutine E_pert_free(nE,Emin,Emax,dE,Z,umass,rr,EvdW_free,dpol_free)
!===========================================================================================
!     Second order perturbative van der Waals
!     correction to the ground state energy of the system
!     
!     Only bound-to-free contribution are
!     considered here
!===========================================================================================
      use mainvar
      implicit none
      integer :: i,j,Lb,Lf,nE
      double precision :: rr,rr3,rr6,Z,Ebound
      double precision :: Cang, denom,summ,EvdW_free
      double precision :: E1,E2,Emin,Emax,dE,detjac1,detjac2
      double precision :: T1,T2,T3,T4,umass,fac,omega,dpol_free
      double precision, dimension(nE,nE) :: func
      double precision, dimension(nE) :: funcpol
      integer n1,n2

     !Angular correction factor
      Cang=1.d0/3.d0

      omega=0.d0

!---------------------------------------------------------
!     Create the function to integrate it
!---------------------------------------------------------
      do i=1,nE
         E1=Emin+(i-1)*dE
         detjac1=sqrt(umass/E1/2.0)

         do j=1,nE
            E2=Emin+(j-1)*dE
            detjac2=sqrt(umass/E2/2.0)
            denom=E1+E2-2.d0*Ebound(1,Z)
 
            if(abs(denom) > 1.d-6) then
               func(i,j)=TRP_free(1,i)*TRP_free(1,j)*detjac1*detjac2/denom
            else
               func(i,j)=0.d0
            endif
 
         enddo
        !fac=(E_bnd(1,Lb)-E1)/((E_bnd(1,Lb)-E1)**2 + omega**2)
        fac=1.d0/(E1-Ebound(1,Z)-omega) !(E1-E_bnd(1,Lb))/((E1-E_bnd(1,Lb))**2 + omega**2)
 
        funcpol(i)=fac*TRP_free(1,i)*detjac1
 
      enddo
!---------------------------------------------------------

      call simpson(nE,1,nE,dE,funcpol,dpol_free)
      dpol_free=dpol_free*Cang
 
      call twoD_trapezoid(nE,nE,dE,dE,Emin,Emax,Emin,Emax,func,summ)
 
      rr3=rr*rr*rr
      rr6=rr3*rr3
 
      EvdW_free=-summ*Cang/rr6

      end subroutine
!===========================================================================================


!===========================================================================================
      double precision function Ebound(N,Z)
!===========================================================================================
      implicit none
      integer :: N
      double precision :: P,Z

      P=Z/float(N) 
      Ebound=-P*P

      end function
!===========================================================================================



!===========================================================================================
      double precision function H_bound_radial(N,Z,L,r)
!===========================================================================================
      implicit none
      integer :: i, N, L 
      integer*16 :: factorial
      double precision :: A,B,P,HG,X,Z,r,prefac

      P=dsqrt(dfloat(factorial(N+L))/dfloat(factorial(N-L-1))/dfloat(2*N))

      prefac=P*(2.d0*Z/dfloat(N))*dsqrt(2.d0*Z/dfloat(N))/dfloat(factorial(2*L+1))

      A=-float((N-L-1))
      B=2.d0*float(L)+2.d0
      X=2.d0*Z*r/float(N)

     !Confluent Hypergeometric Function:
      call CHGM(A,B,X,HG)

      H_bound_radial=prefac*dexp(-X/2.d0)*HG*(X**L)


      end function
!===========================================================================================

!===========================================================================================
      integer*16 function factorial(N)
!===========================================================================================
      implicit none
      integer :: i,N
      integer*16 :: nfac

      if(N == 0) then
         factorial=1
      elseif(N > 0) then

        nfac=1
        do i=1,N
           nfac=nfac*i
        enddo
        factorial=nfac

      else
        stop "Warning negative value in factorial"
      endif

      end function
!===========================================================================================



!========================================================================
      subroutine simpson(ndim,nmin,nmax,step,func,res)
!========================================================================
!     Integration a numerical data series with Simpson's 1/3-rule
!
!     ndim --> lenght of func array
!     nmin --> lower limit of integration
!     nmax --> upper limit of integration
!     step --> equidistant step in variable
!     func --> function for integration
!     res  --> result, value of the integral
!========================================================================
      implicit none
      integer, intent (in) :: nmin,nmax,ndim
      integer :: i,ndata
      double precision, intent (in) :: step
      double precision :: s0,s1,s2
      double precision, intent (out) :: res
      double precision, intent (in), dimension(ndim) :: func

      res=0.d0
      s0=0.d0
      s1=0.d0
      s2=0.d0
      ndata=nmax-nmin+1

      if(nmax > nmin) then
        do i=nmin+1,nmax-1,2
           s1 = s1+func(i-1)
           s0 = s0+func(i)
           s2 = s2+func(i+1)
        enddo
        res = step*(s1+4.d0*s0+s2)/3.d0

!       If n is even, add the last slice separately
        if(mod(ndata,2) .eq. 0) then
           res=res+step*(5.d0*func(nmax)+8.d0*func(nmax-1)-func(nmax-2))/12.d0
        endif
      endif
      end subroutine
!========================================================================



!===========================================================================================
      subroutine twoD_trapezoid(nx,ny,dx,dy,ax,bx,ay,by,func,summa)
!===========================================================================================
!     Numerical 2D trapezoidal integration
!===========================================================================================
      implicit none
      integer :: i,j,nx,ny
      double precision :: ax,bx,ay,by,dx,dy
      double precision :: T1,T2,T3,T4,summa
      double precision, dimension(nx,ny) :: func

      T1=func(1,1)+func(1,ny)+func(nx,1)+func(nx,ny)

      T2=0.d0
      do i=2,ny
         do j=2,nx
            T2=T2+func(j,i)
         enddo
      enddo
      T2=4.d0*T2

      T3=0.d0
      do i=2,ny
         T3=T3 + ( func(1,i)+func(nx,i) )
      enddo
      T3=2.d0*T3

      T4=0.d0
      do i=2,nx
         T4=T4 + ( func(i,1)+func(i,ny) )
      enddo
      T4=2.d0*T4

      summa=(T1+T2+T3+T4)*dx*dy/4.d0

      end subroutine
!===========================================================================================


!===========================================================================================
        SUBROUTINE CHGM(A,B,X,HG)
!===========================================================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: GAMMA for computing â(x)


!******************************************************************
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions jin.ece.uiuc.edu/routines/routines.html" *
!*                                                                *
!*                              F90 Release By J-P Moreau, Paris. *
!*                                     (www.jpmoreau.fr)          *
!******************************************************************

!===========================================================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        A0=A
        A1=A
        X0=X
        HG=0.0D0
        IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
           HG=1.0D+300
        ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
           HG=1.0D0
        ELSE IF (A.EQ.-1.0D0) THEN
           HG=1.0D0-X/B
        ELSE IF (A.EQ.B) THEN
           HG=DEXP(X)
        ELSE IF (A-B.EQ.1.0D0) THEN
           HG=(1.0D0+X/B)*DEXP(X)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           HG=(DEXP(X)-1.0D0)/X
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           R=1.0D0
           HG=1.0D0
           DO 10 K=1,M
              R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10            HG=HG+R
        ENDIF
        IF (HG.NE.0.0D0) RETURN
        IF (X.LT.0.0D0) THEN
           A=B-A
           A0=A
           X=DABS(X)
        ENDIF
        IF (A.LT.2.0D0) NL=0
        IF (A.GE.2.0D0) THEN
           NL=1
           LA=INT(A)
           A=A-LA-1.0D0
        ENDIF
        DO 30 N=0,NL
           IF (A0.GE.2.0D0) A=A+1.0D0
           IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                 RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                 HG=HG+RG
                 IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
15            CONTINUE
           ELSE
              CALL GAMMA(A,TA)
              CALL GAMMA(B,TB)
              XG=B-A
              CALL GAMMA(XG,TBA)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8
                 R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                 SUM1=SUM1+R1
20               SUM2=SUM2+R2
              HG1=TB/TBA*X**(-A)*DCOS(PI*A)*SUM1
              HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2
              HG=HG1+HG2
           ENDIF
25         IF (N.EQ.0) Y0=HG
           IF (N.EQ.1) Y1=HG
30      CONTINUE
        IF (A0.GE.2.0D0) THEN
           DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
35            A=A+1.0D0
        ENDIF
        IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
        A=A1
        X=X0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
!       ==================================================
!       Purpose: Compute gamma function â(x)
!       Input :  x  --- Argument of â(x)
!                       ( x is not equal to 0,-1,-2,úúú)
!       Output:  GA --- â(x)
!       ==================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,  &
                -0.6558780715202538D0, -0.420026350340952D-1, &
                0.1665386113822915D0,-.421977345555443D-1,    &
                -.96219715278770D-2, .72189432466630D-2,      &
                -.11651675918591D-2, -.2152416741149D-3,      &
                .1280502823882D-3, -.201348547807D-4,         &
                -.12504934821D-5, .11330272320D-5,            &
                -.2056338417D-6, .61160950D-8,                &
                .50020075D-8, -.11812746D-8,                  &
                .1043427D-9, .77823D-11,                      & 
                -.36968D-11, .51D-12,                         & 
                -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

!end of file mchgm.f90
!===========================================================================================
