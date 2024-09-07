c ======================================================================
c Module for sharing variables
c ======================================================================
c     Declare variable COMSTATEV to store cohesive element status
c     The length of COMSTATEV is the total number of all types of elements
       MODULE shareVars
       implicit none
       integer, dimension(14400,7) :: cz_table
       integer, dimension(35820) :: COMSTATEV
       END MODULE shareVars     
      
c ============================================================================
c Access node_cz.txt to get the connectivity between node and cohesive element
c ============================================================================      
      
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      use shareVars
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
      integer ir
C
      if (KINC .LE. 0) then
        open(unit = 212,file=
     &   "elecz.txt",
     & status='old',form='formatted', action='read')
        do ir = 1, 14400
          read(212, *) cz_table(ir, :)
        end do
      else
        close(212)
      end if

      RETURN
      END
   
c ====================================================================
c 2D PPR CZM
c ====================================================================
       SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     & PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     & DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     & PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     & NJPRO, PERIOD)
       USE shareVars
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*),
     & SVARS(*), ENERGY(8), COORDS(MCRD, NNODE), U(NDOFEL),
     & DU(MLVARX,*), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(*),
     & JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*),
     & PREDEF(2, NPREDF, NNODE), LFLAGS(*), JPROPS(*)
       DIMENSION Sc(ndofel,ndofel), Fc(ndofel,nrhs),
     & T(mcrd,nrhs), T_d(mcrd,mcrd), U_l(ndofel), R(mcrd, mcrd),
     & Bc(mcrd,ndofel), Bct(ndofel,mcrd), ShapeN(nnode),
     & del(mcrd), GP(2), GP_w(2), tmp(ndofel,mcrd)
       
       PARAMETER(PI=3.1415926D0)
 
       DOUBLE PRECISION Gn, Gt, Tn_m, Tt_m, alph, beta, ln, lt, th,
     & dn, dt, m, n, Gam_n, Gam_t, dGnt, dGtn,
     & N1, N2, del1, del2, del3, del4, deln_max, delt_max, el_length
     & bottom, left, right
c -----------------------------------------------------------------------
c Read input data & Initialize
       
       
       Gn   = PROPS(1)
       Gt   = PROPS(2)
       Tn_m = PROPS(3)
       Tt_m = PROPS(4)
       alph = PROPS(5)
       beta = PROPS(6)
       ln   = PROPS(7)
       lt   = PROPS(8)
       th   = PROPS(9)
       n_GP = 2
       data GP   / 0.577350269189626 , -0.577350269189626 /
       data GP_W / 1.0 , 1.0 /
       call k_Matrix_Zero (RHS,ndofel,nrhs)
       call k_Matrix_Zero (AMATRX,ndofel,ndofel)
              
c Determine the PPR parameters
c   From Equation (20):
       m = (alph-1)*alph*ln**2/(1-alph*ln**2)
       n = (beta-1)*beta*lt**2/(1-beta*lt**2)
c   From Equation (21):
       dn = alph*Gn/(m*Tn_m)*(1-ln)**(alph-1)
     &     * (alph/m*ln+1)**(m-1)*(alph+m)*ln
c   From Equation (22):
       dt = beta*Gt/(n*Tt_m)*(1-lt)**(beta-1)
     &     * (beta/n*lt+1)**(n-1)*(beta+n)*lt
c   From Equation (11):
       if (Gt .GT. Gn) then
          dGnt = 0
          dGtn = Gt - Gn
       elseif (Gt .LT. Gn) then
          dGnt = Gn - Gt
          dGtn = 0
       else
          dGnt = 0
          dGtn = 0
       end if
c   From Equations (18) and (19):
       if (Gn .EQ. Gt) then
          Gam_n = -Gn*(alph/m)**m
          Gam_t = (beta/n)**n
       else
          Gam_n = (-Gn)**(dGnt/(Gn-Gt))*(alph/m)**m
          Gam_t = (-Gt)**(dGtn/(Gt-Gn))*(beta/n)**n
       end if
c Change from the global coordinates to the local coordinates
       call k_Coords_Transform (R, el_length, COORDS, U, ndofel,
     & nnode, mcrd)
       do i = 0, nnode-1
          U_l(1+i*mcrd) = R(1,1)*U(1+i*mcrd) + R(1,2)*U(2+i*mcrd)
          U_l(2+i*mcrd) = R(2,1)*U(1+i*mcrd) + R(2,2)*U(2+i*mcrd)
       end do
       del1 = U_l(7) - U_l(1)
       del2 = U_l(8) - U_l(2)
       del3 = U_l(5) - U_l(3)
       del4 = U_l(6) - U_l(4)
c Numerical integration to compute RHS and AMATRX
       do i = 1, n_GP
          N1 = 0.5*(1 - GP(i))
          N2 = 0.5*(1 + GP(i))
          del(1) = N1*del1 + N2*del3
          del(2) = N1*del2 + N2*del4
          delt_max = SVARS(n_GP*(i-1)+1)
          deln_max = SVARS(n_GP*(i-1)+2)
          call k_Cohesive_PPR (T, T_d, Gam_n, Gam_t, alph, beta, m, n,
     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max)
          ShapeN(1) = -N1
          ShapeN(2) = -N2
          ShapeN(3) = N2
          ShapeN(4) = N1
          do j = 1, nnode
             do k = 1, mcrd
               do l = 1, mcrd
                  Bc(k,l+(j-1)*mcrd) = ShapeN(j)*R(k,l)
               end do
             end do
          end do
          call k_Matrix_Transpose (Bc,Bct,mcrd,ndofel)
          call k_Matrix_Multiply (Bct,T_d,tmp,ndofel,mcrd,mcrd)
          call k_Matrix_Multiply (tmp,Bc,Sc,ndofel,mcrd,ndofel)
          call k_Matrix_Multiply (Bct,T,Fc,ndofel,mcrd,nrhs)
          thick = 0.5 * el_length * GP_w(i) * th
          call k_Matrix_PlusScalar (AMATRX,Sc,thick,ndofel,ndofel)
          call k_Matrix_PlusScalar (RHS,-Fc,thick,ndofel,nrhs)
c   Update the state variables: SVARS
          if((delt_max.LT.abs(del(1))).AND.(abs(del(1)).GT.lt*dt)) then
             SVARS(n_GP*(i-1)+1) = abs(del(1))
          end if
          if ((deln_max .LT. del(2)) .AND. (del(2) .GT. ln*dn)) then
             SVARS(n_GP*(i-1)+2) = del(2)
          end if
      end do
       
      if (((abs(del(1)).GT.lt*dt*0.035d0).OR.(del(2).GT.ln*dn*0.035d0)) 
     &   .and. (coords(2,1).gt.10.0d0)) then 
      COMSTATEV(JELEM) = 1
      endif

      
      SVARS(11) = COMSTATEV(JELEM) 

       RETURN
       END
      
c
c ====================================================================
c = Cohesive traction-separation relation of the PPR model ===========
       SUBROUTINE k_Cohesive_PPR (T,T_d,Gam_n,Gam_t,alph,beta,m,n,
     &  dn, dt, dGtn, dGnt, del, deln_max, delt_max)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION T(2,1), T_d(2,2), del(2)
       DOUBLE PRECISION Gam_n, Gam_t, alph, beta, m, n, dn, dt,
     &     dGtn, dGnt, deln_max, delt_max, Tn, Tt, deln, delt, sign_dt
       delt = abs(del(1))
       deln = del(2)
       if (del(1) .GE. 0) then
          sign_dt = 1
       else
          sign_dt = -1
       end if
       Tn = 0
c Equation (12):
c Pre-calculation of the normal cohesive traction, Tn
       if (deln .LT. 0) then
          deln = 0
       elseif ((deln .GE. dn) .OR. (delt .GE. dt))  then
          Tn = 0
       elseif (deln .GE. deln_max) then
          Tn = (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln/dn)**alph*(m/alph+deln/dn)**(m-1)
     &       -alph*(1-deln/dn)**(alph-1)*(m/alph+deln/dn)**m)
       else
          Tn = (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max
       end if
c Pre-calculation of the tangential cohesive traction, Tt
       if ((deln .GE. dn) .OR. (delt .GE. dt))  then
          Tt = 0
       elseif (delt .GE. delt_max) then
          Tt = (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)
     &       -beta*(1-delt/dt)**(beta-1)*(delt/dt+n/beta)**n)
       else
          Tt = (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
       end if
c Algorithm 1 (description)
c Normal cohesive interaction
c   (1) Contact condition
       if (del(2). LT. 0) then
          T_d(2,2) = -Gam_n/dn**2*(m/alph)**(m-1)*(alph+m)*
     & (Gam_t*(n/beta)**n + dGtn)
          T_d(2,1) = 0
          T(2,1) = T_d(2,2)*del(2)
       else if ((deln.LT.dn).AND.(delt.LT.dt).AND.(Tn.GE.-1.0E-5)) then
          T(2,1) = Tn
c   (2) Softening condition
          if (deln .GE. deln_max) then
             T_d(2,2) =
     & (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn**2 *
     & ((1-deln/dn)**(alph-2)*(alph-1)*alph*(deln/dn+m/alph)**m -
     & 2*(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**(m-1)*m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-2)*(m-1)*m)
             T_d(2,1) =
     & Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m)
c   (3) Unloading/reloading condition
          else
             T_d(2,2) =
     & (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*((1-deln_max/dn)**alph*(deln_max/dn+m/alph)**(m-1)*m
     & -(1-deln_max/dn)**(alph-1)*alph*(deln_max/dn+m/alph)**m)
     & / deln_max
             T_d(2,1) =
     & Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max
          end if
c   (4) Complete failure condition
       else
          T(2,1) = 0
          T_d(2,2) = 0
          T_d(2,1) = 0
       endif
c Tangential cohesive interaction
       if ((delt.LT.dt) .AND. (deln.LT.dn) .AND. (Tt.GE.-1.0E-5)) then
          T(1,1) = Tt*sign_dt
c   (1) Softening condition
          if (delt .GE. delt_max) then
             T_d(1,1) =
     & (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt**2 *
     & ((1-delt/dt)**(beta-2)*(beta-1)*beta*(delt/dt+n/beta)**n -
     & 2*(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**(n-1)*n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-2)*(n-1)*n)
             T_d(1,2) =
     & Gam_t/dt*(-(1-delt/dt)**(beta-1)*beta*(delt/dt+n/beta)**n +
     & (1-delt/dt)**beta*(delt/dt+n/beta)**(n-1)*n) * sign_dt *
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m)
c   (2) Unloading/reloading condition
          else
             T_d(1,1) =
     & (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &      -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n)
     & / delt_max
             T_d(1,2) =
     & Gam_n/dn*(-(1-deln/dn)**(alph-1)*alph*(deln/dn+m/alph)**m +
     & (1-deln/dn)**alph*(deln/dn+m/alph)**(m-1)*m) * sign_dt *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/beta)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
          end if
c   (3) Complete failure condition
       else
          T(1,1) = 0
          T_d(1,1) = 0
          T_d(1,2) = 0
          endif
c      if (T_d(1,2) .NE. T_d(2,1)) then
c         T_d(1,2) = 0.5*(T_d(1,2) + T_d(2,1))
c         T_d(2,1) = T_d(1,2)
c      endif
       RETURN
       END
c
c =====================================================================
c = Coordinate Transformation =========================================
c   : Coordinate transformation matrix (R) is obtained on the basis of
c   the deformed configuration
       SUBROUTINE k_Coords_Transform (R, el_length, COORDS, U, ndofel,
     & nnode, mcrd)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION R(mcrd,mcrd), COORDS(mcrd,nnode), U(ndofel)
       DIMENSION Co_de(mcrd,nnode), Co_de_m(2,2)
c Variables used in the k_Coords_Transform subroutine
c   Co_de  : Coord. of a cohesive element in the deformed configuration
c   Co_de_m: Mid-points of a cohesive element to compute the orientation
c   el_length: length of a cohesive element
c
       do i = 1, mcrd
          do j = 1, nnode
             Co_de(i,j) = COORDS(i,j) + U(2*(j-1)+i)
          end do
       end do
       do i = 1, 2
          Co_de_m(i,1) = (Co_de(i,1)+Co_de(i,4))*0.5
          Co_de_m(i,2) = (Co_de(i,2)+Co_de(i,3))*0.5
       end do
c Calculate the directional cosine & the transformation matrix
       d_x = Co_de_m(1,2) - Co_de_m(1,1)
       d_y = Co_de_m(2,2) - Co_de_m(2,1)
       el_length = (d_x**2 + d_y**2)**0.5
       cos_a = d_x / el_length
       sin_a = d_y / el_length
       R(1,1) = cos_a
       R(1,2) = sin_a
       R(2,1) = -sin_a
       R(2,2) = cos_a
       RETURN
       END
c
c =====================================================================
c = Matrix operations =================================================
       SUBROUTINE k_Matrix_Zero (A,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(n,m)
       do i = 1, n
          do j = 1, m
             A(i,j) = 0.0
          end do
       end do
       RETURN
       END
c
       SUBROUTINE k_Matrix_Transpose (A,B,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(n,m), B(m,n)
       call k_Matrix_zero (B,m,n)
       do i = 1, n
          do j = 1, m
             B(j,i) = A(i,j)
          end do
       end do
       RETURN
       END
c
       SUBROUTINE k_Matrix_PlusScalar (A,B,c,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(n,m), B(n,m)
       do i = 1, n
          do j = 1, m
             A(i,j) = A(i,j) + c*B(i,j)
          end do
       end do
       RETURN
       END
c
       SUBROUTINE k_Matrix_Multiply (A,B,C,l,n,m)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(l,n), B(n,m), C(l,m)
       call k_Matrix_zero (C,l,m)
       do i = 1, l
          do j = 1, m
             do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B (k,j)
             end do
          end do
       end do
       RETURN
      END
c =====================================================================
c USDFLD subroutine
c =====================================================================
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      use shareVars
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      integer ir, ic, ids, ele1, ele2, NPT, field_temp, czm_number
      integer, dimension(4616,7) :: cz_table_local

      field_temp = 0      
      do id = 2,4
          if (cz_table(NOEL, id).NE.0) then 
              ele1 = cz_table(NOEL, id)
              field_temp = field_temp + COMSTATEV(ele1)
          end if
      end do
c
      FIELD(1) = field_temp
      RETURN
      END
      
c      IF (EPS.GT.0.01) THEN
c Suppose that we know which cohesive elements are attached to this 
C solid element
      !IF (COMSTATEV(NOEL+1).EQ.2) then
      !    FIELD(1) = COMSTATEV(NOEL+1)
      !ELSE IF (COMSTATEV(NOEL+2).EQ.2) THEN
      !    FIELD(1) = COMSTATEV(NOEL+1)
      !END IF
c      FIELD(1)=2
c      STATEV(10)=2
c     END IF
c      STATEV(10)=FIELD(1)
      
c ======================================================================
c Biotite weathering model UMAT
c ======================================================================      
      

C     WEATHERING  
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      use shareVars
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C FOR CALCULATION
      PARAMETER(TERVAL=31104.0D3,AB=0.15D0,AR=3.0D0,CINCIN=3.3333D-4)
C FOR CONSTANT
      PARAMETER(FTOL=1.0D-10,PI=3.1415926D0,IPHASEN=39)
C FOR CRACKS
      PARAMETER(Xc=1.0D-8,AINICRACK=1.0D-4,RC=3.0D0,RHOC=0.5D0)
      
      COMMON/PARAMETERS/BBM,BSM,VBM,VSM,BETA,ALPHA,BMO,SMO
C	    
      DIMENSION Q(IPHASEN,3,3),SAO(3,3,3,3),ACAP(3,3,3,3),HC(3,3,3,3),
     1 STIFFI(3,3,3,3),STIFFM(3,3,3,3),STIFFC(3,3,3,3),POP(3,3,3,3)
     2 ,PC(IPHASEN,3,3,3,3),PP(IPHASEN,3,3,3,3),PLSTRAIN(IPHASEN,3,3),
     3 SYMIF(3,3,3,3),PAO(IPHASEN,3,3,3,3),PR(IPHASEN),SecIdentity(3,3),
     4 PA(IPHASEN,3,3,3,3),EGSTRAIN(IPHASEN,3,3),AAP(3,3,3,3),
     5 DTENSOR(IPHASEN,IPHASEN,3,3,3,3),ZEROF(3,3,3,3),
     6 PSTRAIN(IPHASEN,3,3),PSTRESS(IPHASEN,3,3),ETA(3,3),WEIGHT(37),
     7 TSTRAIN(3,3),DTSTRAIN(3,3),TSTRESS(3,3),HCMATRIX(6,6),
     8 COMPONENTL(39),PSI(IPHASEN),THETA(IPHASEN),PHI(IPHASEN),
     9 COMPONENTO(39),DIRECTIONVECTOR(37,3),DIRECTIONTENSOR(39,3,3),
     & NCRACK(37),RHO(37),OMEGA(3,3),INDI(37),Jd(3,3,3,3),
     & J1(3,3,3,3),Js(3,3,3,3),NEPS(37),avePIA(3,3),
     & Jm(3,3,3,3),Ja(3,3,3,3),Z(3,3,3,3),TCRACK(37,3,3,3,3),
     & HCRACK(3,3,3,3),GCRACK(37),F(37)
C	 
	  DOUBLE PRECISION DAMAGEN,DALAMBDA,WRCT,CINC,AINIBIO,BCON,
     1 SMM,BMM,SMI,BMI,VSTRAINI,DAMAGE,DAMAGE1,DAMAGE2,DAMAGE3,
     2 NBIO,POISONR,Errf,tref 
C      
      DIMENSION TEMP21(3,3),TEMP22(3,3),TEMP23(3,3),
     1 TEMP25(3,3),TEMP26(3,3),TEMP27(3,3),TEMP28(3,3),TEMP24(3,3),
     2 TEMP41(3,3,3,3),TEMP42(3,3,3,3),TEMP43(3,3,3,3),TEMP44(3,3,3,3),
     3 TEMP45(3,3,3,3),TEMP46(3,3,3,3),TEMP47(3,3,3,3),TEMP48(3,3,3,3),
     4 TEMP49(3,3,3,3),TEMP410(3,3,3,3),TEMP411(3,3,3,3),
     5 TEMP412(3,3,3,3),TEMP413(3,3,3,3),TEMP414(3,3,3,3),
     6 TEMP415(3,3,3,3),TEMP416(3,3,3,3),TEMP417(3,3,3,3),
     7 TEMP418(3,3,3,3),TEMP419(3,3,3,3),TEMP31(39,3,3),VECT11(3),
     8 J0(3,3,3,3)     
       DOUBLE PRECISION TEMP11, TEMP12, TEMP13, TOP, BOTTOM
       
       INTEGER WINDI,ID, CELE
C
C   FOR DAMAGE
      
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/

C----------------------------------------------------------------
C VOLUME RATIO
	  CALL WEIGHT37(WEIGHT)
	  DO I=1,37
        NCRACK(I)=1.0D5
        END DO 

C        ORIENTATION	  
	  CALL ROTATION37(PHI,THETA,PSI)	  
	  CALL DIRECTION37(DIRECTIONVECTOR)
	  DO I=2,IPHASEN-1
	     CALL Ai_Bj(DIRECTIONVECTOR(I-1,:),DIRECTIONVECTOR(I-1,:),TEMP21)
	     DIRECTIONTENSOR(I,:,:)=(TEMP21) 
        END DO
        CINICRACK=Xc*AINICRACK        
C---------------------------------------------------------------------------
C 	                      MATERIAL PARAMETERS
C----------------------------------------------------------------------------
C
C STATE VARIABLE: 1 BIOTITE DIRECTION (THETA), 2 WR, 3 DAMAGE, 4 WRCT, 5 VSTRAINI, 6 BOITITE LENGTH, 7-43 CRACK DENSITY
C 44 BIOTITE WIDTH 45-81 Open and closed cracks indicator 82 BIOTITE DIRECTION (PSI) 83 DAMAGE1 84 DAMAGE2 85 DAMEGE3
	  IF (TIME(2).EQ.0) THEN
          STATEV(1)=0.0D0*PI/2.0D0 
		STATEV(2)=0.0D0
		STATEV(3)=0.0D0
		STATEV(44)=CINCIN
		STATEV(5)=0.0D0
		STATEV(6)=CINCIN*AR		
		DO I=7,43
		STATEV(I)=NCRACK(I-6)*AINICRACK**3
		END DO
		DO I=45,81
		STATEV(I)=0.0D0
		END DO
		DO I=83,85
		STATEV(I)=0.0D0
		END DO
		STATEV(4)=0.0D0
          STATEV(82)=PROPS(7)   
          STATEV(86)=0.0D0
	  END IF
	  WR=STATEV(2)
	  WRCT=STATEV(4)
	  VSTRAINI=STATEV(5)
	  AINIBIO=STATEV(6)	  
        COMPONENTL(1)=STATEV(6)
        DO I=2,38
         COMPONENTL(I)=AINICRACK  
         COMPONENTO(I)=CINICRACK
        END DO
	  COMPONENTL(IPHASEN)=1.00001D0
        COMPONENTO(1)=STATEV(44)
	  COMPONENTO(IPHASEN)=1.0D0
        
        OMEGA(:,:)=0.0D0
        DO I=1,37
            RHO(I)=STATEV(I+6)  
            INDI(I)=STATEV(I+44)
            CALL A_Bij(WEIGHT(I)*RHO(I), 
     &      DIRECTIONTENSOR(I+1,:,:),3,TEMP21)
            CALL Aij_PLUS_Bij(OMEGA,TEMP21,OMEGA)
        END DO
        
	  OMEGA(1,1)=abs(STATEV(83)*3.3333D-8)
	  OMEGA(3,3)=abs(STATEV(84)*3.3333D-8)
	  OMEGA(1,2)=abs(STATEV(85)*3.3333D-8)  
        PHI(1)=0.0D0*PI/2.0D0
        THETA(1)=0.0D0*PI/2.0D0
        PSI(1)=STATEV(82)
        NBIO=AB/(AINIBIO**3.0D0)*(3.0D0/4.0D0/PI*AR) 
        PR(1)=4.0D0/3.0D0*PI*NBIO*COMPONENTL(1)**2.0D0*COMPONENTO(1)
        DO I=2,38
        PR(I)=4.0D0/3.0D0*PI*NCRACK(I-1)*COMPONENTL(I)**2.0D0*
     &        COMPONENTO(I)
        END DO
        PR(IPHASEN)=1-PR(1)
        
        DO I=1,IPHASEN
            CALL TRANSFORMATION(PSI(I),THETA(I),PHI(I),Q(I,:,:))
        END DO
        
C----------------------------------------------------------------
	  BMM=PROPS(1)
	  SMM=PROPS(2)
	  BBM=PROPS(3)
	  BSM=PROPS(4)
	  VBM=PROPS(5)
	  VSM=PROPS(6)
        PSI(1)=PROPS(7)
	  POISONR=(3.0D0*BMM-2.0D0*SMM)/2.0D0/(3.0D0*BMM+2.0D0*SMM)

C-----------------------------------------------------------------------
C                          INPUT INITIATION
C-----------------------------------------------------------------------
C    INPUT FOR EACH TIME STEP, FROM XY TO XZY----Y CHANGED TO Z
		CALL MAT1_MAT2(NTENS,STRAN,TSTRAIN,HALF)
		CALL MAT1_MAT2(NTENS,DSTRAN,DTSTRAIN,HALF)
		CALL MAT1_MAT2(NTENS,STRESS,TSTRESS,ONE)
          CALL Aij_PLUS_Bij(TSTRAIN,DTSTRAIN,TSTRAIN)
		CALL SYMIDENDITYF(SYMIF)
		ZEROF(:,:,:,:)=0.0D0
C----------------------------------------------------------------
C                        WEATHERING RATE
C----------------------------------------------------------------
C    THE VOLUME OF WEATHERED BIOTITE
      WIndi = 0   
      do id = 2,4
          if (cz_table(NOEL, id).NE.0) then 
              CELE = cz_table(NOEL, id)
              WIndi = WIndi + COMSTATEV(CELE)
          end if
      end do   
   
      tref = STATEV(86)
      IF ((WINDI.GT.0).AND.(tref.EQ.0.0D0)) THEN
          STATEV(86) = Time(2)
      ENDIF
    
      IF (tref.GT.0.0D0) THEN   
		WRCT=10**(-12.32D0)*((TIME(2)-tref))**(-0.603)
c    Surface	
          TEMP11 = (AINIBIO**3.2+2*(AINIBIO*CINCIN)**(1.6))
     & **(1.0D0/1.6)*4*PI
C    Volume
          TEMP12 = 1.3963D-9
          WR = (2.10D-4)*TEMP11*(10**(-12.32D0))*2.5189D0*3.2974D4
     & *((TIME(2)-tref)*TERVAL)**(0.397)/TEMP12
      else 
          WR = 0.0D0
      ENDIF
      
c      WRCT = 0.0D0
c      TOP = 50.0D0 + 10.0D0*COS(PI/20.0D0*COORDS(1))
c      BOTTOM = 40.0D0 + 10.0D0*COS(PI/20.0D0*COORDS(1))
c      IF (TIME(1).GE.1.0D0 .AND.
c     &         (COORDS(2) > BOTTOM)) THEN   
c		WRCT=10**(-12.32D0)*((TIME(2)))**(-0.603)*(
c     &      TOP - COORDS(2))/
c     &     (TOP-BOTTOM)
c    Surface	
c          TEMP11 = (AINIBIO**3.2+2*(AINIBIO*CINCIN)**(1.6))
c	1 **(1.0D0/1.6)*4*PI
C    Volume
c          TEMP12 = 1.3963D-9
c          WR = (2.10D-4)*TEMP11*(10**(-12.32D0))*2.5189D0*3.2974D4
c     & *(TIME(2)*TERVAL)**(0.397)/TEMP12*(
c     &      TOP-COORDS(2))/
c     &     (TOP-BOTTOM)
c          END IF  
      
      

C    EIGEN STRAIN
          EGSTRAIN(:,:,:)=0.0D0
	    DO I=1,3
          DO J=1,3
            ETA(I,J)=0.0D0
          END DO
          END DO
		EGSTRAIN(:,:,:)=0.0D0
		ETA(2,2)=0.4D0*WR/1.4D0
		CALL Aik_Bkj(TRANSPOSE(Q(1,:,:)),ETA,TEMP21)            
		CALL Aik_Bkj(TEMP21,Q(1,:,:),EGSTRAIN(1,:,:))

C     WEATHERED MODULUS OF INCLUSION
		SMI=(1-WR)*BSM+WR*VSM
		BMI=(1-WR)*BBM+WR*VBM
C     PHASES STIFFNESS
		CALL STIFFNESS(BMI,SMI,STIFFI)
		CALL STIFFNESS(BMM,SMM,STIFFM)
		STIFFC(:,:,:,:)=0.0D0
		PC(1,:,:,:,:)=STIFFI
		PC(IPHASEN,:,:,:,:)=STIFFM
          
C     Crack's Stiffness (switch between closed and open cracks)
          SecIdentity(:,:)=0.0D0
          DO I=1,3
              SecIdentity(I,I)=1.0D0
          END DO
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      DO L=1,3
                          TEMP41(I,J,K,L)=SecIdentity(I,K)*
     1                          SecIdentity(J,L)
                          TEMP42(I,J,K,L)=SecIdentity(J,K)*
     1                          SecIdentity(I,L)
                          TEMP43(I,J,K,L)=SecIdentity(I,J)*
     1                          SecIdentity(K,L)/3.0D0
                      END DO 
                  END DO
              END DO
          END DO
          
          CALL Aijkl_PLUS_Bijkl(TEMP41,TEMP42,TEMP44)           
          CALL Aijkl_MINUS_Bijkl(TEMP41,TEMP42,TEMP45)    
          
          CALL A_Bijkl(0.50D0,TEMP44,3,TEMP46)
          CALL A_Bijkl(0.50D0,TEMP45,3,TEMP47)
          CALL Aijkl_MINUS_Bijkl(TEMP46,TEMP43,TEMP48)

          TEMP11=3.0D0*BMM
          CALL A_Bijkl(TEMP11,TEMP43,3,STIFFC)
          
          DO I=2,IPHASEN-1
              CALL A_Bijkl(INDI(I-1),STIFFC,3,PC(I,:,:,:,:))
          END DO
              
C     PHASES P TENSOR
        DO I=1,IPHASEN
		   CALL PTENSOR(COMPONENTL(I),COMPONENTO(I),BMM,SMM,POP)
		   CALL TRANSP(TRANSPOSE(Q(I,:,:)),POP,PP(I,:,:,:,:))
        END DO

C CONSENTRATION TENSOR
          DO I=1,IPHASEN
              CALL Aijkl_MINUS_Bijkl(PC(I,:,:,:,:),PC(IPHASEN,:,:,:,:),
     1 TEMP41)
              CALL Aijmn_Bmnkl(PP(I,:,:,:,:),TEMP41,TEMP42)
              CALL Aijkl_PLUS_Bijkl(SYMIF,TEMP42,TEMP43)
              CALL INVERSEFOURTH(TEMP43,PAO(I,:,:,:,:))
          END DO

          SAO(:,:,:,:)=0.0D0
          DO I=1,IPHASEN
              CALL A_Bijkl(PR(I),PAO(I,:,:,:,:),3,TEMP41)
              CALL Aijkl_PLUS_Bijkl(SAO,TEMP41,SAO)
          END DO

               
          CALL INVERSEFOURTH(SAO,TEMP41)
          
          DO I=1,IPHASEN
              CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),TEMP41,PA(I,:,:,:,:))
          END DO
          
C HOMOGENIZED STIFFNESS
          HC(:,:,:,:)=0.0D0
          DO I=1,IPHASEN
              CALL Aijmn_Bmnkl(PC(I,:,:,:,:),PA(I,:,:,:,:),TEMP42)
              CALL A_Bijkl(PR(I),TEMP42,3,TEMP43)
              CALL Aijkl_PLUS_Bijkl(HC,TEMP43,HC)
          END DO     

          
C----------------------------------------------------------------
C TOTAL STRESS
        avePIA(:,:)=0.0D0
          DO I=1,IPHASEN
              CALL Aijkl_Bkl(PC(I,:,:,:,:),EGSTRAIN(I,:,:),TEMP22)
              CALL Aij_Bijkl(TEMP22,PA(I,:,:,:,:),TEMP23)
              CALL A_Bij(PR(I),TEMP23,3,TEMP24)
              CALL Aij_MINUS_BijN(avePIA,TEMP24,3,avePIA)
          END DO
          CALL Aijkl_Bkl(HC,TSTRAIN,TEMP25)
          CALL Aij_PLUS_Bij(TEMP25,avePIA,TSTRESS)

C D TENSOR       
C WITH DIFFERENT P
          AAP(:,:,:,:)=0.0D0
          DO I=1,IPHASEN
              CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),PP(I,:,:,:,:),TEMP41)
              CALL A_Bijkl(PR(I),TEMP41,3,TEMP42)
              CALL Aijkl_PLUS_Bijkl(AAP,TEMP42,AAP)
          END DO
          ACAP(:,:,:,:)=0.0D0
          DO I=1,IPHASEN
              CALL Aijkl_MINUS_Bijkl(HC,PC(I,:,:,:,:),TEMP41)
              CALL Aijmn_Bmnkl(TEMP41,PAO(I,:,:,:,:),TEMP42)
              CALL Aijmn_Bmnkl(TEMP42,PP(I,:,:,:,:),TEMP43)
              CALL A_Bijkl(PR(I),TEMP43,3,TEMP44)
              CALL Aijkl_PLUS_Bijkl(ACAP,TEMP44,ACAP)
          END DO          
          DO I=1,IPHASEN
              DO J=1,IPHASEN
                  IF (I.EQ.J) THEN
                      CALL A_Bijkl(PR(I),PA(I,:,:,:,:),3,TEMP41)
                      CALL Aijkl_MINUS_Bijkl(SYMIF,TEMP41,TEMP42)
                      CALL Aijmn_Bmnkl(TEMP42,PAO(I,:,:,:,:),TEMP43)
                      CALL Aijmn_Bmnkl(TEMP43,PP(I,:,:,:,:),TEMP44)
C
                      CALL Aijmn_Bmnkl(PA(I,:,:,:,:),AAP,TEMP45)
                      CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),PP(I,:,:,:,:)
     1 ,TEMP46)
                      CALL Aijkl_MINUS_Bijkl(TEMP45,TEMP46,TEMP47)
                      CALL INVERSEFOURTH(ACAP,TEMP48)
                      CALL Aijmn_Bmnkl(TEMP47,TEMP48,TEMP49)
                      CALL Aijkl_MINUS_Bijkl(SYMIF,PA(I,:,:,:,:),
     1 TEMP410)
                      CALL TRANSPOSEF(TEMP410,TEMP411)
                      CALL Aijkl_MINUS_Bijkl(HC,PC(I,:,:,:,:),TEMP412)
                      CALL Aijmn_Bmnkl(TEMP412,PAO(I,:,:,:,:),TEMP413)
                      CALL Aijmn_Bmnkl(TEMP413,PP(I,:,:,:,:),TEMP414)
                      CALL Aijkl_PLUS_Bijkl(TEMP411,TEMP414,TEMP415)
                      CALL A_Bijkl(PR(I),TEMP415,3,TEMP416)
                      CALL Aijmn_Bmnkl(TEMP49,TEMP416,TEMP417)
                      CALL Aijkl_PLUS_Bijkl(TEMP44,TEMP417,TEMP418)
                      CALL Aijmn_Bmnkl(TEMP418,PC(I,:,:,:,:),
     1 DTENSOR(I,J,:,:,:,:))
                  END IF
                  IF (I.NE.J) THEN
                      CALL A_Bijkl(PR(J),PA(I,:,:,:,:),3,TEMP41)
                      CALL Aijkl_MINUS_Bijkl(ZEROF,TEMP41,TEMP42)
                      CALL Aijmn_Bmnkl(TEMP42,PAO(J,:,:,:,:),TEMP43)
                      CALL Aijmn_Bmnkl(TEMP43,PP(J,:,:,:,:),TEMP44)
C
                      CALL Aijmn_Bmnkl(PA(I,:,:,:,:),AAP,TEMP45)
                      CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),PP(I,:,:,:,:)
     1 ,TEMP46)
                      CALL Aijkl_MINUS_Bijkl(TEMP45,TEMP46,TEMP47)
                      CALL INVERSEFOURTH(ACAP,TEMP48)
                      CALL Aijmn_Bmnkl(TEMP47,TEMP48,TEMP49)
                      CALL Aijkl_MINUS_Bijkl(SYMIF,PA(J,:,:,:,:),
     1 TEMP410)
C                     OPEN(111,FILE='DEBUG.OUT')
C                     WRITE(111,*) TEMP410
                      CALL TRANSPOSEF(TEMP410,TEMP411)
                      CALL Aijkl_MINUS_Bijkl(HC,PC(J,:,:,:,:),TEMP412)
                      CALL Aijmn_Bmnkl(TEMP412,PAO(J,:,:,:,:),TEMP413)
                      CALL Aijmn_Bmnkl(TEMP413,PP(J,:,:,:,:),TEMP414)
                      CALL Aijkl_PLUS_Bijkl(TEMP411,TEMP414,TEMP415)
                      CALL A_Bijkl(PR(J),TEMP415,3,TEMP416)
                      CALL Aijmn_Bmnkl(TEMP49,TEMP416,TEMP417)
                      CALL Aijkl_PLUS_Bijkl(TEMP44,TEMP417,TEMP418)
                      CALL Aijmn_Bmnkl(TEMP418,PC(J,:,:,:,:),
     1 DTENSOR(I,J,:,:,:,:))
                  END IF
              END DO
          END DO

C PHASE STRESS AND STRAIN
          DO J=1,IPHASEN
              TEMP21(:,:)=0.0D0
              DO I=1,IPHASEN
                  CALL Aijkl_Bkl(DTENSOR(J,I,:,:,:,:),
     1            EGSTRAIN(I,:,:),TEMP22)
                  CALL Aij_PLUS_Bij(TEMP21,TEMP22,TEMP21)
              END DO
              CALL Aijkl_Bkl(PA(J,:,:,:,:),TSTRAIN,TEMP22)
              CALL Aij_PLUS_Bij(TEMP21,TEMP22,PSTRAIN(J,:,:))
              CALL Aij_MINUS_BijN(PSTRAIN(J,:,:),EGSTRAIN(J,:,:),
     1        3,TEMP23)
              CALL Aijkl_Bkl(PC(J,:,:,:,:),TEMP23,PSTRESS(J,:,:))
              CALL Aik_Bkj(Q(J,:,:),PSTRAIN(J,:,:),TEMP21)            
              CALL Aik_Bkj(TEMP21,TRANSPOSE(Q(J,:,:)),
     1        PLSTRAIN(J,:,:))
          END DO 
C----------------------------------------------------------------
C     UPDATE COMPONENTS SIZE
C          COMPONENTO(1)=(1.0D0+PLSTRAIN(1,3,3))*CINCIN  
C          COMPONENTL(1)=SQRT((1.0D0+PLSTRAIN(1,1,1))*
C     1    (1.0D0+PLSTRAIN(1,2,2)))*CINCIN*AR
C          CALL TRACE(PSTRAIN(1,:,:),VSTRAINI)
C          DO I=2,IPHASEN-1
C              COMPONENTO(I)=(1.0D0+PLSTRAIN(I,3,3))*CINICRACK
C              IF (0.0D0.GT.COMPONENTO(I)) THEN
c                  WRITE(6,*) 'lstrain',PLSTRAIN(I,:,:)	
c                WRITE(6,*) 'COMPONENTO',COMPONENTO(I),I	   
C                  COMPONENTO(I)=0.001D0
C              END IF	
C              COMPONENTL(I)=SQRT((1.0D0+PLSTRAIN(I,1,1))*
C     1        (1.0D0+PLSTRAIN(I,2,2)))*AINICRACK	  	  
C          END DO  
C----------------------------------------------------------------
C     CRACK PROPAGATION
C---------------------------------------------------------------- 
          CALL Aijmn_Bmnkl(PC(1,:,:,:,:),PAO(1,:,:,:,:),TEMP41)
          CALL A_Bijkl(PR(1),TEMP41,3,TEMP42)
          CALL A_Bijkl(PR(IPHASEN),PC(IPHASEN,:,:,:,:),3,TEMP43)
          CALL Aijkl_PLUS_Bijkl(TEMP42,TEMP43,Z)
          
          DO N=1,37
              CALL Aijmn_Bmnkl(PP(N+1,:,:,:,:),PC(IPHASEN,:,:,:,:),
     &        TEMP41)
              CALL Aijkl_MINUS_Bijkl(SYMIF,TEMP41,TEMP42)
              CALL INVERSEFOURTH(TEMP42,TEMP43)
              CALL A_Bijkl(Xc,TEMP43,3,Tcrack(N,:,:,:,:))
              CALL Aijkl_Bkl(PA(N+1,:,:,:,:),TSTRAIN,TEMP22)
              CALL Aij_Bij(TEMP22,DIRECTIONTENSOR(N+1,:,:),
     &        NEPS(N))
              IF (NEPS(N).GT.0.0D0) THEN
                  INDI(N)=0.0D0
              ELSE
                  INDI(N)=1.0D0
              END IF
          END DO
          
          CALL INVERSEFOURTH(Z,TEMP41) 
          CALL INVERSEFOURTH(SAO,TEMP42) 
          
          
          DO N=1,37
              CALL Aijmn_Bmnkl(TCRACK(N,:,:,:,:),TEMP41,TEMP43)
              TEMP11=4.0D0*PI/3.0D0
              CALL A_Bijkl(TEMP11,TEMP43,3,HCRACK)
c             WRITE(6,*) 'HCRACK(:,:,1,1)',HCRACK(:,:,1,1)
c             WRITE(6,*) 'HCRACK(:,:,2,2)',HCRACK(:,:,2,2)
c              WRITE(6,*) 'HCRACK(:,:,3,3)',HCRACK(:,:,3,3)
c              WRITE(6,*) 'HCRACK(:,:,1,2)',HCRACK(:,:,1,2)
              
C             A_partial_rho
              CALL  Aijmn_Bmnkl(PA(1,:,:,:,:), TCRACK(N,:,:,:,:),TEMP44)  
              CALL  Aijmn_Bmnkl(TEMP44,TEMP42,TEMP45)
              TEMP12=-1.0D0*4.0D0*PI/3.0D0
              CALL A_Bijkl(TEMP12,TEMP45,3,TEMP46)
C              WRITE(6,*) 'TEMP46(:,:,1,1)',TEMP46(:,:,1,1)
C              WRITE(6,*) 'TEMP46(:,:,2,2)',TEMP46(:,:,2,2)
C              WRITE(6,*) 'TEMP46(:,:,3,3)',TEMP46(:,:,3,3)
C              WRITE(6,*) 'TEMP46(:,:,1,2)',TEMP46(:,:,1,2)
c             eigenstress 
              CALL Aijkl_Bkl(PC(1,:,:,:,:),EGSTRAIN(1,:,:),TEMP21)
              CALL A_Bij(-1.0D0, TEMP21, 3, TEMP21)
c              WRITE(6,*) 'TEMP21',TEMP21
c              WRITE(6,*) 'TSTRAIN',TSTRAIN
c             second term
              CALL Aij_Bijkl(TEMP21,TEMP46,TEMP22)
              CALL Aij_Bij(TEMP22,TSTRAIN,TEMP12)
C              WRITE(6,*) 'TEMP12',TEMP12
c             tstress - avePIA              
              CALL Aij_MINUS_BijN(TSTRESS,avePIA,3,TEMP24)
C              WRITE(6,*) 'TEMP24',TEMP24
              CALL Aij_Bijkl(TEMP24,HCRACK,TEMP25)
              CALL Aij_Bij(TEMP25,TEMP24,TEMP13)
              GCRACK(N)=0.50D0*TEMP13+PR(1)*TEMP12
C              WRITE(6,*) 'TEMP13',TEMP13
c              WRITE(6,*) 'TEMP12',TEMP12
c              WRITE(6,*) 'GCRACK(N)',GCRACK(N)
              F(N)= GCRACK(N)-(RC*RHO(N)*EXP(1-RHO(N)/RHOC))
     &        /RHOC 
c              IF (TIME(2).GT.1.0D0) THEN
c              F(N)= GCRACK(N)-(RC*RHO(N)*EXP(1-RHO(N)/RHOC))
c     &        /RHOC  
c              END IF
c              WRITE(6,*) 'Criterion',F(N)
              TEMP13=0.0D0
              DO WHILE (F(N).GT.FTOL.AND.TEMP13.LT.100.0D0) 
              TEMP11= (RC*EXP(-(RHO(N)-RHOC)/RHOC)*(RHO(N)-RHOC))
     &        /(RHOC**2)
              TEMP12=-1.0D0/TEMP11*F(N)
              IF (TEMP12.LT.0.0D0) THEN
                  TEMP12=0.0D0
              END IF
              RHO(N)=RHO(N)+TEMP12
              F(N)= GCRACK(N)-(RC*RHO(N)*EXP(1-RHO(N)/RHOC))
     &        /RHOC 
              TEMP13=TEMP13+1.0D0
C              WRITE(6,*) 'Iterations',TEMP13
              END DO
          END DO


C----------------------------------------------------------------
C     OUTPUT
C----------------------------------------------------------------
          CALL MAT2_MAT1(NTENS,TSTRESS,STRESS,ONE)	 
          CALL MAT4_MAT2(HC,HCMATRIX,1)        
	
          DO I = 1 , NTENS
              DO J = 1 , NTENS
                  DDSDDE (I , J) = HCMATRIX (I , J)
              END DO
          END DO
          CALL Aik_Bkj(Q(1,:,:),PSTRAIN(1,:,:),TEMP21)            
          CALL Aik_Bkj(TEMP21,TRANSPOSE(Q(1,:,:)),TEMP22)	  
          CALL TRACE(PSTRAIN(1,:,:),VSTRAINI)
          STATEV(2)=WR
          STATEV(3)=tref
C          CALL TRACE(OMEGA,STATEV(3))
          STATEV(4)=WRCT
          STATEV(5)=VSTRAINI	 
          STATEV(6)=COMPONENTL(1)
          DO I=7,43
              STATEV(I)=RHO(I-6)
          END DO
          STATEV(44)=COMPONENTO(1)
          DO I=45,81
              STATEV(I)=INDI(I-44)
          END DO
          STATEV(82)=PSI(1)
          OMEGA(:,:)=0.0D0
          DO I=1,37
              CALL A_Bij(WEIGHT(I)*RHO(I), 
     &        DIRECTIONTENSOR(I+1,:,:),3,TEMP21)
              CALL Aij_PLUS_Bij(OMEGA,TEMP21,OMEGA)
          END DO 
          STATEV(83)=abs(OMEGA(1,1)/3.3333D-8)
          STATEV(84)=abs(OMEGA(3,3)/3.3333D-8)
          STATEV(85)=abs(OMEGA(1,2))/3.3333D-8
          RETURN
	 END
C
C
C ====================================================================
C 
C                          Aij_Bij  
C
C ====================================================================
C-----------------------------------------------------------------------	 
      SUBROUTINE Aij_Bij(A,B,C)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3)
      DOUBLE PRECISION C
      DATA ZERO /0.0D0/
C
      C=ZERO
      DO I=1,3
          DO J=1,3
              C=C+A(I,J)*B(I,J)
          END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C               TRACE  
C
C ====================================================================
C-----------------------------------------------------------------------	 
      SUBROUTINE TRACE(TENSOR,R)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3)
      R=0.0D0
C
      DO I=1,3
          R=R+TENSOR(I,I)
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------	 
C
      SUBROUTINE MAT2_MAT1(NTENS,TENSOR,VECTOR,FACT)
C
C ====================================================================
C
C =================== MAT2_MAT1: TENSOR TO VECTOR=====================
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,NTENS
          VECTOR(I)=0.0D0
      END DO
C
      VECTOR( 1 ) = TENSOR(1 , 1)
      VECTOR( 2 ) = TENSOR(2 , 2)
      VECTOR( 3 ) = TENSOR(3 , 3)
      VECTOR( 4 ) = TENSOR(1 , 2)*FACT
C
      RETURN
      END
C
C	 
C-----------------------------------------------------------------------	 
      SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : 2D PLAIN STRAIN VECTOR TO 3D TENSOR  
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,3
        DO J=1,3
          TENSOR(I,J)=0.0D0
        END DO
      END DO
C
      TENSOR(1 , 1) = VECTOR( 1 )
      TENSOR(2 , 2) = VECTOR( 2 )
      TENSOR(3 , 3) = VECTOR( 3 )
      TENSOR(1 , 2) = VECTOR( 4 )*FACT
      TENSOR(2 , 1) = VECTOR( 4 )*FACT
C
      RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C                STIFFNESS 
C
C ====================================================================
      SUBROUTINE STIFFNESS(BU,AMU,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
        ALAMDA=BU-2.0D0*AMU/3.0D0
        DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
                    R(I,J,K,L)=0.0D0
                END DO
            END DO
        END DO
        END DO
        R(1,1,1,1)=2.0*amu+alamda
        R(2,2,2,2)=2.0*amu+alamda
        R(3,3,3,3)=2.0*amu+alamda
        R(1,1,2,2)=alamda
        R(1,1,3,3)=alamda
        R(2,2,1,1)=alamda
        R(2,2,3,3)=alamda
        R(3,3,1,1)=alamda
        R(3,3,2,2)=alamda
        R(2,3,2,3)=amu
        R(2,3,3,2)=amu
        R(3,2,3,2)=amu
        R(3,2,2,3)=amu
        R(2,1,2,1)=amu
        R(2,1,1,2)=amu
        R(1,2,2,1)=amu
        R(1,2,1,2)=amu
        R(1,3,1,3)=amu
        R(1,3,3,1)=amu
        R(3,1,1,3)=amu
        R(3,1,3,1)=amu
		RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C                INVERSEFOURTH 
C
C ====================================================================  
      SUBROUTINE INVERSEFOURTH(F,FINV)
            INCLUDE 'ABA_PARAM.INC' 
      DIMENSION F(3,3,3,3),a(6,6),W(6,6),FINV(3,3,3,3)
	  DOUBLE PRECISION CONS
	  INTEGER IX,IY
	
		a(1,1)=F(1,1,1,1)
		a(1,2)=F(1,1,2,2)
		a(1,3)=F(1,1,3,3)
		a(1,4)=F(1,1,1,2)
		a(1,5)=F(1,1,2,3)
		a(1,6)=F(1,1,3,1)
		a(2,1)=F(2,2,1,1)
		a(2,2)=F(2,2,2,2)
		a(2,3)=F(2,2,3,3)
		a(2,4)=F(2,2,1,2)
		a(2,5)=F(2,2,2,3)
		a(2,6)=F(2,2,3,1)
		a(3,1)=F(3,3,1,1)
		a(3,2)=F(3,3,2,2)
		a(3,3)=F(3,3,3,3)
		a(3,4)=F(3,3,1,2)
		a(3,5)=F(3,3,2,3)
		a(3,6)=F(3,3,3,1)
		a(4,1)=F(1,2,1,1)
		a(4,2)=F(1,2,2,2)
		a(4,3)=F(1,2,3,3)
		a(4,4)=F(1,2,1,2)
		a(4,5)=F(1,2,2,3)
		a(4,6)=F(1,2,3,1)
		a(5,1)=F(2,3,1,1)
		a(5,2)=F(2,3,2,2)
		a(5,3)=F(2,3,3,3)
		a(5,4)=F(2,3,1,2)
		a(5,5)=F(2,3,2,3)
		a(5,6)=F(2,3,3,1)
		a(6,1)=F(3,1,1,1)
		a(6,2)=F(3,1,2,2)
		a(6,3)=F(3,1,3,3)
		a(6,4)=F(3,1,1,2)
		a(6,5)=F(3,1,2,3)
		a(6,6)=F(3,1,3,1)
      CALL INVERSE(a,6,6,W)
C	  
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      DO L=1,3
                          FINV(I,J,K,L)=0.0D0
                      END DO
                  END DO
              END DO
          END DO
C		  
          DO I=1,3
              DO J=1,3
				DO K=1,3
					DO L=1,3
						CONS=2.0D0
						IF (I.EQ.1.AND.J.EQ.1) THEN
						IX=1
						END IF
						IF (I.EQ.2.AND.J.EQ.2) THEN
						IX=2
						END IF
						IF (I.EQ.3.AND.J.EQ.3) THEN
						IX=3
						END IF
						IF (I.EQ.1.AND.J.EQ.2) THEN
						IX=4
						END IF
						IF (I.EQ.2.AND.J.EQ.1) THEN
						IX=4
						END IF
						IF (I.EQ.3.AND.J.EQ.2) THEN
						IX=5
						END IF
						IF (I.EQ.2.AND.J.EQ.3) THEN
						IX=5
						END IF
						IF (I.EQ.1.AND.J.EQ.3) THEN
						IX=6
						END IF 
						IF (I.EQ.3.AND.J.EQ.1) THEN
						IX=6
						END IF
						IF (K.EQ.1.AND.L.EQ.1) THEN
						IY=1
						END IF
						IF (K.EQ.2.AND.L.EQ.2) THEN
						IY=2
						END IF
						IF (K.EQ.3.AND.L.EQ.3) THEN
						IY=3
						END IF
						IF (K.EQ.1.AND.L.EQ.2) THEN
						IY=4
						END IF
						IF (K.EQ.2.AND.L.EQ.1) THEN
						IY=4
						END IF
						IF (K.EQ.3.AND.L.EQ.2) THEN
						IY=5
						END IF
						IF (K.EQ.2.AND.L.EQ.3) THEN
						IY=5
						END IF
						IF (K.EQ.1.AND.L.EQ.3) THEN
						IY=6
						END IF 
						IF (K.EQ.3.AND.L.EQ.1) THEN
						IY=6
						END IF
						IF (IX.LE.3.AND.IY.LE.3) THEN
						CONS=1.0D0
						END IF
						IF (IX.GT.3.AND.IY.GT.3) THEN
						CONS=4.0D0
						END IF    						
                        FINV(I,J,K,L)=W(IX,IY)/CONS
					END DO
				END DO
              END DO
          END DO
        return
      END
C-----------------------------------------------------------------------	  
C ====================================================================
C 
C                ESHELBY 
C
C ====================================================================  
      SUBROUTINE ESHELBY(ANU,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      DIJ=0.0D0
                      DKL=0.0D0
                      DIK=0.0D0
                      DJL=0.0D0
                      DIL=0.0D0
                      DJK=0.0D0
              IF(I.EQ.J) THEN
              DIJ=1.0D0
              END IF
              IF(K.EQ.L) THEN
              DKL=1.0D0
              END IF              
              IF(I.EQ.K) THEN
              DIK=1.0D0
              END IF              
              IF(J.EQ.L) THEN
              DJL=1.0D0
              END IF
              IF(I.EQ.L) THEN
              DIL=1.0D0
              END IF
              IF(J.EQ.K) THEN
              DJK=1.0D0
              END IF
              R(I,J,K,L)=(5.0D0*ANU-1.0D0)/15.0D0/(1.0D0-ANU)*DIJ*DKL
     1+(-5.0D0*ANU+4.0D0)/15.0D0/(1.0D0-ANU)*(DIK*DJL+DIL*DJK)
                  END DO
              END DO 
          END DO
      END DO
	  RETURN
      END
C-----------------------------------------------------------------------	  
        SUBROUTINE TRANSFORMATION(PSI,THETA,PHI,QQ)
CC====================================================================================================
C
C                          TRANSFORMATION MATRIX Q	  
C
C=====================================================================================================
C
        INCLUDE 'ABA_PARAM.INC'
C
        DIMENSION QQ(3,3)
        DOUBLE PRECISION PSI,THETA,PHI
C
        DO I=1,3
            DO J=1,3
                QQ(I,J)=0.0D0
            END DO
        END DO
        QQ(1,1)=COS(PSI)*COS(THETA)*COS(PHI)
     1  -SIN(PSI)*SIN(PHI)
        QQ(2,1)=-COS(PHI)*SIN(PSI)
     1  -COS(PSI)*COS(THETA)*SIN(PHI)
        QQ(3,1)=SIN(THETA)*COS(PSI)
        
        QQ(1,2)=COS(PHI)*COS(THETA)*SIN(PSI)
     1  +COS(PSI)*SIN(PHI)
        QQ(2,2)=-SIN(PSI)*COS(THETA)*SIN(PHI)
     1  +COS(PSI)*COS(PHI)
        QQ(3,2)=SIN(PSI)*SIN(THETA)
        
        QQ(1,3)=-COS(PHI)*SIN(THETA)
        QQ(2,3)=SIN(PHI)*SIN(THETA)
        QQ(3,3)=COS(THETA)
C
        RETURN
        END
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          SYMIDENDITYF	  
C
C=====================================================================================================
C
      SUBROUTINE  SYMIDENDITYF(R)
            INCLUDE 'ABA_PARAM.INC'
         DIMENSION R(3,3,3,3)
         R(:,:,:,:)=0.0D0
        DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
        R(1,1,1,1)=1.0D0
        R(2,2,2,2)=1.0D0
        R(3,3,3,3)=1.0D0
        R(2,3,2,3)=0.5D0
        R(2,3,3,2)=0.5D0
        R(3,2,3,2)=0.5D0
        R(3,2,2,3)=0.5D0
        R(2,1,2,1)=0.5D0
        R(2,1,1,2)=0.5D0
        R(1,2,2,1)=0.5D0
        R(1,2,1,2)=0.5D0
        R(1,3,1,3)=0.5D0
        R(1,3,3,1)=0.5D0
        R(3,1,1,3)=0.5D0
        R(3,1,3,1)=0.5D0
                  END DO
              END DO
          END DO
        END DO
      RETURN
      end 
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          TRANSPOSEF	  
C
C=====================================================================================================
C	
      SUBROUTINE TRANSPOSEF(A,B)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3,3,3),B(3,3,3,3)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      B(I,J,K,L)=A(K,L,I,J)
                  END DO
              END DO
          END DO
      END DO
      RETURN
      END 	
C-----------------------------------------------------------------------
C ====================================================================
C                        MAT4_MAT2                                                                  I
C I        THIS PROGRAM TRANSFORMS THE FOURTH ORDER COMPLIANCE       I
C I        TENSOR TO A SECOND ORDER MATRIX                           I
C I                                                                  I
C ====================================================================
C
      SUBROUTINE MAT4_MAT2(TENSOR,DMATRIX,ICOE)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3,3,3),DMATRIX(6,6)
C
      DATA ZERO,TWO /0.0D0,2.0D0/
C
C     D2 = THE SECOND ORDER STIFFNESS MATRIX
C
      DO I=1,6
        DO J=1,6
          DMATRIX(I,J)=0.0D0
        END DO
      END DO

      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C
      DMATRIX(1,1)=TENSOR(1,1,1,1)
      DMATRIX(1,2)=TENSOR(1,1,2,2)
      DMATRIX(1,3)=TENSOR(1,1,3,3)
      DMATRIX(1,4)=TENSOR(1,1,1,2)*COE1
      DMATRIX(1,5)=TENSOR(1,1,2,3)*COE1
      DMATRIX(1,6)=TENSOR(1,1,1,3)*COE1
C
      DMATRIX(2,1)=TENSOR(2,2,1,1)
      DMATRIX(2,2)=TENSOR(2,2,2,2)
      DMATRIX(2,3)=TENSOR(2,2,3,3)
      DMATRIX(2,4)=TENSOR(2,2,1,2)*COE1
      DMATRIX(2,5)=TENSOR(2,2,2,3)*COE1
      DMATRIX(2,6)=TENSOR(2,2,1,3)*COE1
C
      DMATRIX(3,1)=TENSOR(3,3,1,1)
      DMATRIX(3,2)=TENSOR(3,3,2,2)
      DMATRIX(3,3)=TENSOR(3,3,3,3)
      DMATRIX(3,4)=TENSOR(3,3,1,2)*COE1
      DMATRIX(3,5)=TENSOR(3,3,2,3)*COE1
      DMATRIX(3,6)=TENSOR(3,3,1,3)*COE1
C
      DMATRIX(4,1)=TENSOR(1,2,1,1)*COE1
      DMATRIX(4,2)=TENSOR(1,2,2,2)*COE1
      DMATRIX(4,3)=TENSOR(1,2,3,3)*COE1
      DMATRIX(4,4)=TENSOR(1,2,1,2)*COE2
      DMATRIX(4,5)=TENSOR(1,2,2,3)*COE2
      DMATRIX(4,6)=TENSOR(1,2,1,3)*COE2
C
      DMATRIX(5,1)=TENSOR(2,3,1,1)*COE1
      DMATRIX(5,2)=TENSOR(2,3,2,2)*COE1
      DMATRIX(5,3)=TENSOR(2,3,3,3)*COE1
      DMATRIX(5,4)=TENSOR(2,3,1,2)*COE2
      DMATRIX(5,5)=TENSOR(2,3,2,3)*COE2
      DMATRIX(5,6)=TENSOR(2,3,1,3)*COE2
C
      DMATRIX(6,1)=TENSOR(1,3,1,1)*COE1
      DMATRIX(6,2)=TENSOR(1,3,2,2)*COE1
      DMATRIX(6,3)=TENSOR(1,3,3,3)*COE1
      DMATRIX(6,4)=TENSOR(1,3,1,2)*COE2
      DMATRIX(6,5)=TENSOR(1,3,2,3)*COE2
      DMATRIX(6,6)=TENSOR(1,3,1,3)*COE2
C
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------	  
      SUBROUTINE A_Bij(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION B(3,3),C(3,3)
      DO I=1,3
          DO J=1,3
          C(I,J)=A*B(I,J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE A_Bijkl(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION B(N,N,N,N),C(N,N,N,N)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
          C(I,J,K,L)=A*B(I,J,K,L)
                  END DO
              END DO
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bi(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N)
      C=0.0D0
      DO I=1,3
          C=C+A(I)*B(I)
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bj(A,B,C)
      INCLUDE 'ABA_PARAM.INC' 
      DIMENSION A(3),B(3),C(3,3)
      C(:,:)=0.0D0
      DO I=1,3
          DO J=1,3
          C(I,J)=A(I)*B(J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_PLUS_Bi(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N),C(N)
      DO I=1,3
          C(I)=A(I)+B(I)
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_Bj(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N,N)
      DIMENSION B(N),C(N)
      DO I=1,3
          C(I)=0.0D0
          DO J=1,3
          C(I)=C(I)+A(I,J)*B(J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_MINUS_BijN(A,B,N,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_BijN                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
C      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)-B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_PLUS_Bij(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_PLUS_BijN(A,B,N,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_BijN                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
C      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_Bkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO I=1,3
        DO J=1,3
          C(I,J)=ZERO
          DO K=1,3
            DO L=1,3
              C(I,J)=C(I,J)+A(I,J,K,L)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_MINUS_Bijkl(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
          C(I,J,K,L)=A(I,J,K,L)-B(I,J,K,L)
                END DO
            END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_PLUS_Bijkl(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
          C(I,J,K,L)=A(I,J,K,L)+B(I,J,K,L)
                END DO
            END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijmn_Bmnkl(A,B,C)
CC====================================================================================================
C
C                          Aijmn_Bmnkl=Cijkl	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=0.0D0
              DO M=1,3
                DO N=1,3
            C(I,J,K,L)=C(I,J,K,L)+A(I,J,M,N)*B(M,N,K,L)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aik_Bkj(A,B,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bij(A,B,C)
CC====================================================================================================
C
C                          Ai_Bij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3),B(3,3),C(3)
C
C
      DO J=1,3
          C(J)=0.0D0
        DO I=1,3
            C(J)=C(J)+A(I)*B(I,J)
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aik_BkjN(A,B,N,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE Aij_Bijkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aij_Bijkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DOUBLE PRECISION A(3,3),B(3,3,3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO K=1,3
        DO L=1,3
          C(K,L)=0.0D0
          DO I=1,3
            DO J=1,3
              C(K,L)=C(K,L)+A(I,J)*B(I,J,K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INVERSE(A,N,NP,AINV)
C========================================================================
C
C    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
C    A^{-1} = AINV    
C    this subroutine inverses a (n x n) A matrix
C	 following a Gauss-Jordan elimination process
C
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C  
      DIMENSION A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     1 A0(NP,NP),AINV(NP,NP)
C
      DO J=1,N
        IPIV(J)=0
      END DO
C
C
C     storage of the original A matrix
C
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
C
C	find a pivot among the rows of A that have not already been reduced
C
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                   BIG=ABS(A(J,K))
                   IROW=J
                   ICOL=K
                   PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
              END IF
            END DO
          END IF
        END DO
C
      IPIV(ICOL)=IPIV(ICOL)+1
      INDXR(I)=IROW
      INDXC(I)=ICOL
C	  
C     interchange the rows to put the pivot on the diagonal
C
      IF(irow.ne.icol)THEN
        DO L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
        END DO
      END IF
C     reduction of the row of the pivot
C       
      IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
C       
      PIVIN=1./PIV          ! numerical stabilization
C
      A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
C
C     reduction of the column of the pivot
C
        DO LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.    ! numerical stabilization
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            END DO
          END IF
        END DO
      END DO
C
C
C     unscramble the columns to get A-1
C		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
C
C	restitution process of A and Ainv
C
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
C
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
C     
      RETURN
      END
C----------------------------------------------------------------
C----------------------------------------------------------------
C      UNIFORM DISTRIBUTED
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE r8_uniform_01 ( Iseed,RANUNI )
       INCLUDE 'ABA_PARAM.INC'

C      integer i4_huge 
      parameter ( i4_huge = 2147483647 )
C      integer k
      double precision RANUNI
C      integer seed

      k = Iseed / 127773

      Iseed = 16807 * ( Iseed - k * 127773 ) - k * 2836

      if ( Iseed .lt. 0 ) then
        Iseed = Iseed + i4_huge
      end if

      RANUNI = dble ( Iseed ) * 4.656612875D-10

      return
      end
C----------------------------------------------------------------
C----------------------------------------------------------------
C      LOGNORMAL DISTRIBUTED
C----------------------------------------------------------------
C----------------------------------------------------------------	
      SUBROUTINE rlognorm ( alm,ald,iseed,RANLOG )
             INCLUDE 'ABA_PARAM.INC'
      integer iseed
      a=log((alm**2)/(((alm**2)+ald)**(0.5D0)))
      b=(log(ald/(alm**2)+1.0D0))**(0.5D0)
	  CALL r8_normal_ab ( a, b, iseed,R1 )
      RANLOG=exp(R1)
	  return
      end 
C----------------------------------------------------------------
C----------------------------------------------------------------
C      NORMAL DISTRIBUTED
C----------------------------------------------------------------
C----------------------------------------------------------------	
      SUBROUTINE r8_normal_ab ( a, b, iseed ,RANNOR)
       INCLUDE 'ABA_PARAM.INC'
      parameter ( r8_pi = 3.141592653589793D00 )
      integer iseed
      CALL r8_uniform_01 ( iseed ,R1)
      CALL r8_uniform_01 ( iseed,R2 )
      x = sqrt ( -2.0D0 * log ( R1 ) ) 
     1  * cos ( 2.0D0 * r8_pi * R2 )

      RANNOR = a + b * x

      return
      end	 
C----------------------------------------------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
C
C                          P TENSOR
C
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE PTENSOR(A,C,BM,SM,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
	  DOUBLE PRECISION A,C,BM,SM,W,H,PK,PL,PN,PM,PP
	  W=C/A
	  H=W*(ACOS(W)-W*SQRT(1.0D0-W**2))/(1.0D0-W**2)**(1.5D0)
	  PK=((7.0D0*H-2.0D0*W**2-4.0D0*H*W**2)*SM+3.0D0*(H-2.D0
     1 *W**2+2.0D0*H*W**2)*BM)/8.0D0/(1.0D0-W**2)/SM/(
     2 4.0D0*SM+3.0D0*BM)
	  PL=(SM+3.0D0*BM)*(2.0D0*W**2-H-2.0D0*H*W**2)/4.0D0
     1 /(1.0D0-W**2)/SM/(4.0D0*SM+3.0D0*BM)
	  PN=((6.0D0-5.0D0*H-8.0D0*W**2+8*H*W**2)*SM+3.0D0*
     1(H-2.0D0*W**2+2.0D0*H*W**2)*BM)/2.0D0/(1.0D0-W**2)/
     2 SM/(4.0D0*SM+3.0D0*BM)
	  PM=((15.0D0*H-2.0D0*W**2-12.0D0*H*W**2)*SM+3.0D0*
     1 (3.0D0*H-2.0D0*W**2)*BM)/16.0D0/(1.0D0-W**2)/SM/
     2 (4.0D0*SM+3.0D0*BM)
	  PP=(2.0D0*(4.0D0-3.0D0*H-2.0D0*W**2)*SM+3.0D0*(2.0D0-
     1 3.0D0*H+2.0D0*W**2-3.0D0*H*W**2)*BM)/8.0D0/(1.0D0-W**2)
     2 /SM/(4.0D0*SM+3.0D0*BM)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
					R(I,J,K,L)=0.0D0
                  END DO
              END DO 
          END DO
      END DO
	  R(1,1,1,1)=PK+PM
	  R(2,2,2,2)=PK+PM
	  R(3,3,3,3)=PN
	  R(1,1,2,2)=PK-PM
	  R(1,1,3,3)=PL
	  R(2,2,1,1)=PK-PM
	  R(2,2,3,3)=PL	
	  R(3,3,1,1)=PL	
	  R(3,3,2,2)=PL	
	  R(2,3,2,3)=PP
	  R(2,3,3,2)=PP
	  R(3,2,2,3)=PP
	  R(3,2,3,2)=PP
	  R(3,1,1,3)=PP
	  R(3,1,3,1)=PP
	  R(1,3,1,3)=PP
	  R(1,3,3,1)=PP
	  R(1,2,1,2)=PM
	  R(1,2,2,1)=PM
	  R(2,1,1,2)=PM
	  R(2,1,2,1)=PM	  	  
	  RETURN
      END
C----------------------------------------------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
C
C                          TRANSP
C
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE TRANSP(TQ,P,PO)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TQ(3,3),P(3,3,3,3),PO(3,3,3,3)
C
	  PO(:,:,:,:)=0.0D0
	  DO I=1,3
	  DO J=1,3
	  DO K=1,3
	  DO L=1,3
	  DO II=1,3
	  DO JJ=1,3 
	  DO KK=1,3
	  DO LL=1,3
		PO(I,J,K,L)=PO(I,J,K,L)+TQ(I,II)*TQ(J,JJ)*TQ(K,KK)
     1 *TQ(L,LL)*P(II,JJ,KK,LL)
	  END DO
	  END DO
	  END DO
	  END DO
	  END DO
	  END DO
	  END DO 
	  END DO
	  
C
      RETURN
      END	
C----------------------------------------------------------------
C----------------------------------------------------------------
C 
C              ROTATION37
C 
C----------------------------------------------------------------
C----------------------------------------------------------------	  
      SUBROUTINE ROTATION37(PHI,THETA,PSI)
	  INCLUDE 'ABA_PARAM.INC'
	  DIMENSION PSI(39),THETA(39),PHI(39)
        PSI(:)=0.0D0
        PHI(:)=0.0D0
        THETA(:)=0.0D0
        PSI(2)=0.0D0
        THETA(2)=1.570796D0	  
        PSI(3)=0D0
        THETA(3)=1.570796D0	
        PSI(4)=0.0D0
        THETA(4)=0.0D0	  
        PSI(5)=0.785398D0
        THETA(5)=1.570796D0	  
        PSI(6)=-0.785398D0
        THETA(6)=1.570796D0	  
        PSI(7)=0.0D0
        THETA(7)=0.785398D0	  
        PSI(8)=0.0D0
        THETA(8)=2.356194D0	  
        PSI(9)=0.0D0
        THETA(9)=0.785398D0		  
        PSI(10)=0.0D0
        THETA(10)=2.356194D0		  
        PSI(11)=0.314090D0
        THETA(11)=1.570796D0	  
        PSI(12)=-0.314090D0
        THETA(12)=1.570796D0	 	  
        PSI(13)=1.256706D0
        THETA(13)=1.570796D0	 	  
        PSI(14)=-1.256706D0
        THETA(14)=1.570796D0	 	  
        PSI(15)=0.0D0
        THETA(15)=1.256706D0	 	  
        PSI(16)=0.0D0
        THETA(16)=1.884886D0	 	  
        PSI(17)=0.0D0
        THETA(17)=0.314090D0	  
        PSI(18)=0.0D0
        THETA(18)=2.827502D0		  
        PSI(19)=0.0D0
        THETA(19)=1.256706D0	  
        PSI(20)=0.0D0
        THETA(20)=1.884886D0	  	  
        PSI(21)=0.0D0
        THETA(21)=0.314090D0		  
        PSI(22)=0.0D0
        THETA(22)=2.827502D0	  
        PSI(23)=0.785398D0
        THETA(23)=0.493805D0	  
        PSI(24)=0.785398D0
        THETA(24)=2.647787D0	  
        PSI(25)=-0.785398D0
        THETA(25)=0.493805D0	  
        PSI(26)=-0.785398D0
        THETA(26)=2.647787D0	  
        PSI(27)=1.207103D0
        THETA(27)=1.229027D0	  
        PSI(28)=1.207103D0
        THETA(28)=1.912566D0
        PSI(29)=-1.207103D0
        THETA(29)=1.229027D0
        PSI(30)=-1.207103D0
        THETA(30)=1.912566D0
        PSI(31)=0.363694D0
        THETA(31)=1.229027D0
        PSI(32)=0.363694D0
        THETA(32)=1.912566D0
        PSI(33)=-0.363694D0
        THETA(33)=1.229027D0
        PSI(34)=-0.363694D0
        THETA(34)=1.912566D0  
        PSI(35)=0.785398D0
        THETA(35)=0.955317D0
        PSI(36)=0.785398D0
        THETA(36)=2.186276D0
        PSI(37)=-0.785398D0
        THETA(37)=0.955317D0
        PSI(38)=-0.785398D0
        THETA(38)=2.186276D0       
	  RETURN
	  END
C----------------------------------------------------------------
C----------------------------------------------------------------
C 
C              WEIGHTS37
C 
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE WEIGHT37(WEIGHT)
	  INCLUDE 'ABA_PARAM.INC'
	  DIMENSION WEIGHT(37)
      WEIGHT(1)=0.0214477714606D0
      WEIGHT(2)=0.0214477714606D0
      WEIGHT(3)=0.0214477714606D0
      WEIGHT(4)=0.0422832190396D0
      WEIGHT(5)=0.0422832190396D0
	  WEIGHT(6)=0.0422832190396D0
      WEIGHT(7)=0.0422832190396D0	  
      WEIGHT(8)=0.0422832190396D0	  
      WEIGHT(9)=0.0422832190396D0	  
      WEIGHT(10)=0.01071011181674D0	  
      WEIGHT(11)=0.01071011181674D0	  
      WEIGHT(12)=0.01071011181674D0	  
      WEIGHT(13)=0.01071011181674D0	  
      WEIGHT(14)=0.01071011181674D0
      WEIGHT(15)=0.01071011181674D0	  
      WEIGHT(16)=0.01071011181674D0	  
      WEIGHT(17)=0.01071011181674D0	  
      WEIGHT(18)=0.01071011181674D0	  
      WEIGHT(19)=0.01071011181674D0	  
      WEIGHT(20)=0.01071011181674D0	  
      WEIGHT(21)=0.01071011181674D0	  
      WEIGHT(22)=0.0335541818312D0	  
      WEIGHT(23)=0.0335541818312D0	  
      WEIGHT(24)=0.0335541818312D0	
      WEIGHT(25)=0.0335541818312D0	
      WEIGHT(26)=0.0335541818312D0	
      WEIGHT(27)=0.0335541818312D0	
      WEIGHT(28)=0.0335541818312D0	
      WEIGHT(29)=0.0335541818312D0	
      WEIGHT(30)=0.0335541818312D0	
      WEIGHT(31)=0.0335541818312D0
      WEIGHT(32)=0.0335541818312D0
      WEIGHT(33)=0.0335541818312D0
      WEIGHT(34)=0.0376964619016D0
      WEIGHT(35)=0.0376964619016D0
      WEIGHT(36)=0.0376964619016D0
      WEIGHT(37)=0.0376964619016D0	
	  RETURN
	  END
C----------------------------------------------------------------
C----------------------------------------------------------------
C 
C              DIRECTION37
C 
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE DIRECTION37(DIRECTIONVECTOR)
	  INCLUDE 'ABA_PARAM.INC'
	  DIMENSION DIRECTIONVECTOR(37,3)	  
	  DIRECTIONVECTOR(1,:)=(/ 1.0D0, 0.0D0, 0.0D0 /)
	  DIRECTIONVECTOR(2,:)=(/ 0.0D0, 1.0D0, 0.0D0 /)	  
	  DIRECTIONVECTOR(3,:)=(/ 0.0D0, 0.0D0, 1.0D0 /)
	  DIRECTIONVECTOR(4,:)=(/ 0.707106781D0, 0.707106781D0, 0.0D0 /)	  
	  DIRECTIONVECTOR(5,:)=(/ 0.707106781D0, -0.707106781D0, 0.0D0 /)
	  DIRECTIONVECTOR(6,:)=(/ 0.707106781D0, 0.0D0, 0.707106781D0 /)
	  DIRECTIONVECTOR(7,:)=(/ 0.707106781D0, 0.0D0, -0.707106781D0 /)
	  DIRECTIONVECTOR(8,:)=(/ 0.0D0, 0.707106781D0, 0.707106781D0 /)
	  DIRECTIONVECTOR(9,:)=(/ 0.0D0, 0.707106781D0, -0.707106781D0 /)
	  DIRECTIONVECTOR(10,:)=(/ 0.951077869651D0, 0.308951267775D0,
     1 0.0D0 /)
	  DIRECTIONVECTOR(11,:)=(/ 0.951077869651D0, -0.308951267775D0,
     1 0.0D0 /)	 
	  DIRECTIONVECTOR(12,:)=(/ 0.308951267775D0, 0.951077869651D0,
     1 0.0D0 /)
      DIRECTIONVECTOR(13,:)=(/ 0.308951267775D0, -0.951077869651D0,
     1 0.0D0 /)
	  DIRECTIONVECTOR(14,:)=(/ 0.951077869651D0, 0.0D0,
     1 0.308951267775D0 /)
      DIRECTIONVECTOR(15,:)=(/ 0.951077869651D0, 0.0D0,
     1 -0.308951267775D0 /)
	  DIRECTIONVECTOR(16,:)=(/ 0.308951267775D0, 0.0D0,
     1 0.951077869651D0 /)
      DIRECTIONVECTOR(17,:)=(/ 0.308951267775D0, 0.0D0,
     1 -0.951077869651D0 /)
      DIRECTIONVECTOR(18,:)=(/ 0.0D0, 0.951077869651D0,
	1 0.308951267775D0 /)	 
      DIRECTIONVECTOR(19,:)=(/ 0.0D0, 0.951077869651D0,
     1 -0.308951267775D0 /)	
      DIRECTIONVECTOR(20,:)=(/ 0.0D0, 0.308951267775D0,
     1 0.951077869651D0 /)	
      DIRECTIONVECTOR(21,:)=(/ 0.0D0, 0.308951267775D0,
     1 -0.951077869651D0 /)
      DIRECTIONVECTOR(22,:)=(/ 0.335154591939D0, 0.335154591939D0,
     1 0.880535518310D0 /) 
      DIRECTIONVECTOR(23,:)=(/ 0.335154591939D0, 0.335154591939D0,
     1 -0.880535518310D0 /)
	  DIRECTIONVECTOR(24,:)=(/ 0.335154591939D0, -0.335154591939D0,
     1 0.880535518310D0 /)
      DIRECTIONVECTOR(25,:)=(/ 0.335154591939D0, -0.335154591939D0,
     1 -0.880535518310D0 /)
      DIRECTIONVECTOR(26,:)=(/ 0.335154591939D0, 0.880535518310D0,
     1 0.335154591939D0 /)
      DIRECTIONVECTOR(27,:)=(/ 0.335154591939D0, 0.880535518310D0,
     1 -0.335154591939D0 /)
      DIRECTIONVECTOR(28,:)=(/ 0.335154591939D0, -0.880535518310D0,
     1 0.335154591939D0 /)
      DIRECTIONVECTOR(29,:)=(/ 0.335154591939D0, -0.880535518310D0,
     1 -0.335154591939D0 /)	 
      DIRECTIONVECTOR(30,:)=(/ 0.880535518310D0, 0.335154591939D0,
     1 0.335154591939D0 /)	
      DIRECTIONVECTOR(31,:)=(/ 0.880535518310D0, 0.335154591939D0,
     1 -0.335154591939D0 /)		 
      DIRECTIONVECTOR(32,:)=(/ 0.880535518310D0, -0.335154591939D0,
     1 0.335154591939D0 /)	
      DIRECTIONVECTOR(33,:)=(/ 0.880535518310D0, -0.335154591939D0,
     1 -0.335154591939D0 /)	
      DIRECTIONVECTOR(34,:)=(/ 0.577350269190D0, 0.577350269190D0,
     1 0.577350269190D0 /)	
      DIRECTIONVECTOR(35,:)=(/ 0.577350269190D0, 0.577350269190D0,
     1 -0.577350269190D0 /)
      DIRECTIONVECTOR(36,:)=(/ 0.577350269190D0, -0.577350269190D0,
     1 0.577350269190D0 /)
      DIRECTIONVECTOR(37,:)=(/ 0.577350269190D0, -0.577350269190D0,
     1 -0.577350269190D0 /)	 
	 RETURN
      END
      
      
      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2), COORDS (3)
      CHARACTER*80 SNAME
     
c      F = 5.0D0*SIN(0.314159D0*COORDS(1))
       F = 20-0.002*TIME(2)
c       F = -0.0196*COORDS(2)
c10     - TIME(2)*0.1
      WRITE(6,*) 'F', F

      RETURN
      END
           