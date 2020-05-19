!     A SOLVER FOR THE CHEMICAL MASTER EQUATION
!     BY USING THE KRYLOV-FSP-SSA ALGORITHM
MODULE KRYLOVSOLVER
  USE STATESPACE
CONTAINS
  !=====================INTERFACE FOR KRYLOV-FSP-SSA ALGORITHM
  SUBROUTINE CME_SOLVE(MODEL, T, FSP_IN, FSP_OUT, FSPTOL, EXP_TOL, VERBOSITY)
    IMPLICIT NONE
    !------------------------------------------------INPUT/OUTPUT PARAMETERS
    !     THE UNDERLYING MODEL
    TYPE (CME_MODEL), INTENT(IN) :: MODEL
    !     TIME ELAPSED FOR THE CME SOLUTION
    DOUBLE PRECISION, INTENT(IN) :: T
    !     HASH TABLE OF THE FSP
    TYPE(FINITE_STATE_PROJECTION) :: FSP_IN
    !     TOLERANCE OF THE KRYLOV-FSP-SSA ALGORITHM
    DOUBLE PRECISION, INTENT(IN) :: FSPTOL
    !     TOLERANCE OF THE LOCAL EXPONENTIAL APPROXIMATION
    DOUBLE PRECISION, INTENT(IN) :: EXP_TOL
    !     HASH TABLE OF THE FSP
    TYPE(FINITE_STATE_PROJECTION) :: FSP_OUT
    !-----------------------------------------VERBOSITY AND ERROR FLAG
    INTEGER, INTENT(IN), OPTIONAL :: VERBOSITY
    INTEGER :: ITRACE, IFLAG
    !-------------------------------------APPLY THE KRYLOV-FSP-SSA ALGORITHM
    IF (PRESENT(VERBOSITY)) THEN
       ITRACE = VERBOSITY
    ELSE
       ITRACE = 0
    ENDIF

    PRINT*, 'CALLING DGEXPV_FSP'
    CALL DGEXPV_FSP(MODEL, T, FSP_IN%VECTOR, FSP_OUT, FSP_OUT%VECTOR, &
         FSPTOL, EXP_TOL, ITRACE, IFLAG)

  END SUBROUTINE CME_SOLVE


  !===============================================KRYLOV-FSP-SSA ALGORITHM
  SUBROUTINE DGEXPV_FSP(MODEL, T, V, FSP, W, FSPTOL, KRYTOL, ITRACE, IFLAG)
    IMPLICIT NONE

    INTEGER :: NFULL
    PARAMETER(NFULL = NMAX)

    INTEGER :: M_MAX, M_MIN
    PARAMETER(M_MAX = 100, M_MIN = 10)
    !     WORKSPACE
    INTEGER :: LWSP, LIWSP
    PARAMETER(LWSP = NFULL * (M_MAX + 2) + 5 * (M_MAX + 2)**2 + 7, LIWSP = NFULL)
    DOUBLE PRECISION :: WSP(LWSP)
    INTEGER :: IWSP(LIWSP)

    !     THE UNDERLYING MODEL
    TYPE(CME_MODEL) :: MODEL
    !     TIME ELAPSED FOR THE CME SOLUTION
    DOUBLE PRECISION, INTENT(IN) :: T
    !     INPUT PROBABILITY VECTOR
    DOUBLE PRECISION, INTENT(IN) :: V(:)
    !     OUTPUT PROBABILITY VECTOR
    DOUBLE PRECISION, INTENT(OUT) :: W(:)
    !     TOLERANCE OF THE KRYLOV-FSP-SSA ALGORITHM
    DOUBLE PRECISION, INTENT(IN) :: FSPTOL
    !     TOLERANCE OF THE LOCAL KRYLOV ERROR
    DOUBLE PRECISION :: KRYTOL

    TYPE(FINITE_STATE_PROJECTION) :: FSP
    !
    INTEGER :: ITRACE
    !
    INTEGER :: IFLAG


    !-------------------------PARAMETERS FOR THE STEPSIZE CONTROL OF EXPOKIT
    !     MAXIMUM ALLOWABLE NUMBER OF INTEGRATION STEPS;
    !     0 MEANS AN INFINITE NUMBER OF STEPS
    INTEGER, PARAMETER :: MXSTEP = 0

    INTEGER, PARAMETER :: MXREJECT = 0         !     MAXIMUM ALLOWABLE NUMBER OF REJECTIONS AT EACH STEP;
    !     0 MEANS AN INFINITE NUMBER OF REJECTIONS

    INTEGER, PARAMETER :: IDEG = 6         !     THE PADE APPROXIMATION OF TYPE (IDEG,IDEG) IS USED AS AN
    !     APPROXIMATION TO EXP(H); THE VALUE 0 SWITCHES TO THE UNIFORM
    !     RATIONAL CHEBYSHEV APPROXIMATION OF TYPE (14,14)
    DOUBLE PRECISION, PARAMETER :: DELTA = 1.2D0 !     LOCAL TRUNCATION ERROR `SAFETY FACTOR'
    !     STEPSIZE `SHRINKING FACTOR'
    DOUBLE PRECISION, PARAMETER :: GAMMA = 0.9D0


    !--------------------------------------------PARAMETERS OF THE ALGORITHM
    INTEGER :: N, M
    DOUBLE PRECISION :: ANORM, WTMP(NFULL)
    INTEGER :: I, J, K1, MH, MX, IV, IH, J1V, NS, IFREE, LFREE, IEXPH
    INTEGER :: IREJECT, IBRKFLAG, MBRKDWN, NMULT, NREJECT, NEXPH
    INTEGER :: NSCALE, NSTEP
    DOUBLE PRECISION :: SGN, T_OUT, TBRKDWN, STEP_MIN, STEP_MAX, ERR_LOC
    DOUBLE PRECISION :: S_ERROR, X_ERROR, T_NOW, T_NEW, T_STEP, T_OLD
    DOUBLE PRECISION :: XM, BETA, BREAK_TOL, P1, P2, P3, EPS, RNDOFF
    DOUBLE PRECISION :: VNORM, AVNORM, HJ1J, HIJ, HUMP, SQR1
    INTRINSIC :: AINT, ABS, DBLE, LOG10, MAX, MIN, NINT, SIGN
    DOUBLE PRECISION :: DDOT, DNRM2, DASUM

    DOUBLE PRECISION :: HTMP((M_MAX + 2) * (M_MAX + 2))

    ! FOR THE KRYLOV DIMENSION ADAPTIVITY,
    DOUBLE PRECISION :: K_FACTOR, ORDER
    LOGICAL :: M_CHANGED, ORDEROLD, KESTOLD
    INTEGER :: M_NEW, M_OLD, IMREJECT, JOLD, IH_OLD
    REAL :: COST1, COST2
    DOUBLE PRECISION :: OMEGA, OMEGA_OLD, T_OPT
    DOUBLE PRECISION :: HNORM, NOM
    INTEGER :: NNZ, M_START, M_OPT
    INTEGER :: QIOP

    ! FOR THE FSP ADAPTIVITY
    INTEGER :: N_NOW, IEXPAND
    DOUBLE PRECISION :: WSUM, T_RATIO
    LOGICAL :: DROP(NFULL)
    DOUBLE PRECISION :: ERROR, ERROROLD, TAU_OLD, TFSP
    DOUBLE PRECISION :: FSPORDER
    INTEGER :: IREJECTFSP, ISTART
    DOUBLE PRECISION :: T_SSA
    DOUBLE PRECISION :: WSUM_OLD
    DOUBLE PRECISION :: DROPTOL, DSUM


    !------------------------------------CHECK AND INITIALIZE THE PARAMETERS
    !     ASSUME ||A||=1, USUALLY IT IS UNKNOWN
    ANORM = 1.0D0
    CALL MATRIX_STARTER(FSP, MODEL)

    DO I = 1, 5
       CALL ONESTEP_EXTENDER(FSP, MODEL)
    ENDDO
    !     INITIALIZE THE KRYLOV DIMENSION AND ORTHOGONALIZATION LENGTH
    M = M_MIN
    QIOP = 2
    ISTART = 1
    N = NFULL

    !     CHECK THE RESTRICTIONS ON THE INPUT PARAMETERS
    IFLAG = 0
    IF (LWSP<(N * (M + 2) + 5 * (M + 2)**2 + IDEG + 1))  IFLAG = -1
    IF (LIWSP<(M + 2))                       IFLAG = -2
    IF ((M>=N).OR.(M<=0))                 IFLAG = -3
    IF (IFLAG.NE.0) THEN
       PRINT*, 'BAD SIZES (IN INPUT OF DGEXPV), IFLAG = ', IFLAG
       STOP
    END IF

    !     INITIALIZE THE PARAMETERS
    IBRKFLAG = 0
    NMULT = 0
    NREJECT = 0
    NEXPH = 0
    NSCALE = 0
    T_OUT = ABS(T)
    TBRKDWN = 0.0D0
    STEP_MIN = T_OUT
    STEP_MAX = 0.0D0
    NSTEP = 0
    S_ERROR = 0.0D0
    X_ERROR = 0.0D0
    T_NOW = 0.0D0
    T_NEW = 0.0D0
    P1 = 4.0D0 / 3.0D0
1   P2 = P1 - 1.0D0
    P3 = P2 + P2 + P2
    EPS = ABS(P3 - 1.0D0)
    IF (EPS==0.0D0) GO TO 1
    IF (KRYTOL<=EPS)   KRYTOL = SQRT(EPS)
    RNDOFF = EPS * ANORM
    BREAK_TOL = 1.0D-07
    SGN = SIGN(1.0D0, T)

    CALL DCOPY(FSP%SIZE, V, 1, W, 1)
    BETA = DNRM2(FSP%SIZE, W, 1)
    VNORM = BETA
    HUMP = BETA

    !     OBTAIN THE VERY FIRST STEPSIZE
    SQR1 = SQRT(0.1D0)
    XM = 1.0D0 / DBLE(M)
    P1 = KRYTOL * (((M + 1) / 2.72D0)**(M + 1)) * SQRT(2.0D0 * 3.14D0 * (M + 1))
    T_NEW = (1.0D0 / ANORM) * (P1 / (4.0D0 * BETA * ANORM))**XM
    P1 = 10.0D0**(NINT(LOG10(T_NEW) - SQR1) - 1)
    T_NEW = AINT(T_NEW / P1 + 0.55D0) * P1

    !     INITIALIZE THE FSP PARAMETERS
    N_NOW = FSP%SIZE
    IEXPAND = 0
    WSUM_OLD = 1.0D0
    IREJECTFSP = 0
    DROPTOL = 1.0D-8

    NNZ = (MODEL%NREACTIONS + 1) * FSP%SIZE

    !     PARAMETERS FOR THE KRYLOV SUBSPACE ADAPTIVITY,
    !     BASED ON NIESEN AND WRIGHT
    IMREJECT = 0
    JOLD = 1
    M_NEW = M
    ORDEROLD = .TRUE.
    KESTOLD = .TRUE.

100 IF (T_NOW >= T_OUT) GOTO 500
    !     FIND THE REAL STEPSIZE
    T_STEP = MIN(T_OUT - T_NOW, T_NEW)
    !     COMPUTE THE ARNOLDI MATRICES H_M AND V_M
    N = N_NOW
    M = MIN(N - 1, M_NEW)
    MBRKDWN = M
    !     FIND POINTERS INTO THE WORKSPACE
    K1 = 2
    MH = M + 2
    IV = 1
    IH = IV + N * (M + 1) + N
    IFREE = IH + MH * MH
    LFREE = LWSP - IFREE + 1

    NSTEP = NSTEP + 1

    P1 = 1.0D0 / BETA
    DO I = 1, N
       WSP(IV + I - 1) = P1 * W(I)
    ENDDO
    DO I = 1, MH * MH
       WSP(IH + I - 1) = 0.0D0
    ENDDO

    IREJECT = 0
    !     START THE STOPWATCH
    IF (ITRACE.NE.0) THEN
       PRINT*, 'BEGINNING IOP...'
    ENDIF
101 CONTINUE
    J1V = IV + JOLD * N
    DO 200 J = JOLD, M
       NMULT = NMULT + 1
       CALL FMATVEC(WSP(J1V - N), WSP(J1V), FSP%MATRIX)
       IF (QIOP>0)    ISTART = MAX(1, J - QIOP + 1)
       DO I = ISTART, J
          HIJ = DDOT(N, WSP(IV + (I - 1) * N), 1, WSP(J1V), 1)
          CALL DAXPY(N, -HIJ, WSP(IV + (I - 1) * N), 1, WSP(J1V), 1)
          WSP(IH + (J - 1) * MH + I - 1) = HIJ
       ENDDO
       HJ1J = DNRM2(N, WSP(J1V), 1)
       !           IF A `HAPPY BREAKDOWN' OCCURS, GO STRAIGHTFORWARD AT THE END
       IF (HJ1J<=BREAK_TOL) THEN
          K1 = 0
          IBRKFLAG = 1
          MBRKDWN = J
          TBRKDWN = T_NOW
          T_STEP = T_OUT - T_NOW
          GO TO 300
       END IF
       WSP(IH + (J - 1) * MH + J) = HJ1J
       CALL DSCAL(N, 1.0D0 / HJ1J, WSP(J1V), 1)
       J1V = J1V + N
200    CONTINUE
       NMULT = NMULT + 1
       CALL FMATVEC(WSP(J1V - N), WSP(J1V), FSP%MATRIX)
       AVNORM = DNRM2(N, WSP(J1V), 1)
       !     SET 1 FOR THE 2-CORRECTED SCHEME
300    CONTINUE
       WSP(IH + M * MH + M + 1) = 1.0D0
401    CONTINUE

       !     COMPUTE W=BETA*V*EXP(T_STEP*H)*E1...
       NEXPH = NEXPH + 1
       MX = MBRKDWN + K1
       IF (IDEG.NE.0) THEN
          !           BY THE IRREDUCIBLE RATIONAL PADE APPROXIMATION...
          CALL DGPADMNORM(IDEG, MX, SGN * T_STEP, WSP(IH), MH, &
               WSP(IFREE), LFREE, IWSP, IEXPH, NS, IFLAG, HNORM)
          IEXPH = IFREE + IEXPH - 1
          NSCALE = NSCALE + NS
       ELSE
          !           OR THE UNIFORM RATIONAL CHEBYSHEV APPROXIMATION
          IEXPH = IFREE
          DO I = 1, MX
             WSP(IEXPH + I - 1) = 0.0D0
          ENDDO
          WSP(IEXPH) = 1.0D0
          CALL DGCHBV(MX, SGN * T_STEP, WSP(IH), MH, &
               WSP(IEXPH), WSP(IFREE + MX))
       ENDIF
402    CONTINUE
       !     ESTIMATE THE ERROR
       IF (K1==0) THEN
          ERR_LOC = KRYTOL
       ELSE
          P1 = ABS(WSP(IEXPH + M)) * BETA
          P2 = ABS(WSP(IEXPH + M + 1)) * BETA * AVNORM
          IF (P1>10.0D0 * P2) THEN
             ERR_LOC = P2
             XM = 1.0D0 / DBLE(M)
          ELSE IF (P1>P2) THEN
             ERR_LOC = (P1 * P2) / (P1 - P2)
             XM = 1.0D0 / DBLE(M)
          ELSE
             ERR_LOC = P1
             XM = 1.0D0 / DBLE(M - 1)
          END IF
       ENDIF
       !     REDUCE THE STEPSIZE IF THE ERROR IS OUT OF BOUND
       IF (ISNAN(ERR_LOC)) THEN
          T_STEP = T_STEP / 5.0D0
          GO TO 401
       ENDIF

       !     FIND THE RATIO OF ERROR PER UNIT STEP
       OMEGA_OLD = OMEGA
       OMEGA = ERR_LOC / (KRYTOL * T_STEP)
       !     ESTIMATE THE ORDER
       IF ((M==M_OLD).AND.(T_STEP.NE.T_OLD).AND.(IREJECT>=1)) THEN
          ORDER = MAX(1.0D0, LOG(OMEGA / OMEGA_OLD) / LOG(T_STEP / T_OLD))
          ORDEROLD = .FALSE.
       ELSE IF (ORDEROLD.OR.IREJECT==0) THEN
          ORDER = DBLE(M) / 4.0D0
          ORDEROLD = .TRUE.
       ELSE
          ORDEROLD = .TRUE.
       ENDIF
       !     ESTIMATE KAPPA
       IF ((M.NE.M_OLD).AND.(T_STEP==T_OLD).AND.(IREJECT>=1)) THEN
          K_FACTOR = MAX(1.1D0, (OMEGA / OMEGA_OLD)**(1.0D0 / (M_OLD - M)))
          KESTOLD = .FALSE.
       ELSE IF (KESTOLD .OR. IREJECT==0) THEN
          KESTOLD = .TRUE.
          K_FACTOR = 2.0D0
       ELSE
          KESTOLD = .TRUE.
       ENDIF
       !     RECORD THE OLD STEP SIZE AND KRYLOV DIMENSION
       T_OLD = T_STEP
       M_OLD = M
       !     SUGGEST NEW STEP SIZE AND KRYLOV DIMENSION
       IF ((M==M_MAX).AND.(OMEGA>DELTA).OR.(IMREJECT>4)) THEN
          T_NEW = MIN(T_OUT - T_NOW, &
               MAX(T_STEP / 5.0D0, &
               MIN(5.0D0 * T_STEP, &
               GAMMA * T_STEP * OMEGA**(-1.0D0 / ORDER))))
          P1 = 10.0D0**(NINT(LOG10(T_NEW) - SQR1) - 1)
          T_NEW = AINT(T_NEW / P1) * P1
          M_CHANGED = .FALSE.
       ELSE
          !           COMPUTE OPTIONS FOR THE NEW STEP SIZE AND KRYLOV DIMENSION
          T_OPT = MIN(T_OUT - T_NOW, &
               MAX(T_STEP / 5.0D0, &
               MIN(5.0D0 * T_STEP, &
               GAMMA * T_STEP * OMEGA**(-1.0D0 / ORDER))))
          M_OPT = MIN(MAX(M_MIN, 3 * M / 4, &
               M + CEILING(LOG(OMEGA) / LOG(K_FACTOR))), &
               M_MAX, &
               CEILING(4.0D0 * M / 3.0D0) + 1)

          !           ESTIMATE COSTS
          COST1 = KRYLOV_COST(T_NOW, T_OUT, T_OPT, M, N, HNORM)
          COST2 = KRYLOV_COST(T_NOW, T_OUT, T_STEP, M_OPT, N, HNORM)

          IF (COST1<=COST2) THEN
             T_NEW = T_OPT
             P1 = 10.0D0**(NINT(LOG10(T_NEW) - SQR1) - 1)
             T_NEW = AINT(T_NEW / P1) * P1
             M_NEW = M
             M_CHANGED = .FALSE.
          ELSE
             M_NEW = M_OPT
             T_NEW = T_STEP
             M_CHANGED = .TRUE.
          ENDIF
       ENDIF
       !     IF THE KRYLOV STEP ERROR IS NOT ACCEPTABLE, REJECT THE STEPSIZE
       IF ((K1.NE.0).AND.(OMEGA>DELTA).AND.&
            (MXREJECT==0.OR.IREJECT<MXREJECT)) THEN
          IF (.NOT.M_CHANGED) THEN
             !                 CHOOSE TO CHANGE STEP SIZE
             T_STEP = MIN(T_OUT - T_NOW, &
                  MAX(T_STEP / 5.0D0, &
                  MIN(5.0D0 * T_STEP, T_NEW)))
             P1 = 10.0D0**(NINT(LOG10(T_STEP) - SQR1) - 1)
             T_STEP = AINT(T_STEP / P1 + 0.55D0) * P1
             IF (ITRACE.NE.0) THEN
                PRINT*, 'T_STEP =', T_OLD
                PRINT*, 'ERR_LOC =', ERR_LOC
                PRINT*, 'ERR_REQUIRED =', DELTA * T_OLD * KRYTOL
                PRINT*, 'STEPSIZE REJECTED, DOWN TO:', T_STEP
             ENDIF
             IREJECT = IREJECT + 1
             NREJECT = NREJECT + 1
             IF ((MXREJECT.NE.0).AND.(IREJECT>MXREJECT)) THEN
                PRINT*, 'FAILURE IN DGEXPV: ---'
                PRINT*, 'THE REQUESTED TOLERANCE IS TOO HIGH.'
                PRINT*, 'RERUN WITH A SMALLER VALUE.'
                IFLAG = 2
                RETURN
             ENDIF
             GO TO 401
          ELSE
             !                 CHOOSE TO CHANGE DIMENSION
             NREJECT = NREJECT + 1
             IMREJECT = IMREJECT + 1
             M = M_NEW
             HTMP(1:IFREE - IH + 1) = WSP(IH:IFREE)
             IH_OLD = IH
             MBRKDWN = M
             K1 = 2
             MH = M + 2
             IV = 1
             IH = IV + N * (M + 1) + N
             IFREE = IH + MH * MH
             LFREE = LWSP - IFREE + 1
             T_STEP = MIN(T_OUT - T_NOW, T_NEW)
             DO I = 1, MH * MH
                WSP(IH + I - 1) = 0.0D0
             ENDDO
             !                 COPY THE HESSENBERG MATRIX
             DO J = 1, M_OLD
                DO I = 1, J + 1
                   WSP(IH + (J - 1) * (M + 2) + I - 1) = &
                        HTMP(1 + (J - 1) * (M_OLD + 2) + I - 1)
                ENDDO
             ENDDO
             !                 ONLY DO ARNOLDI FROM THE M_OLD-TH COLUMN
             JOLD = M_OLD
             IF (ITRACE.NE.0) THEN
                PRINT*, 'ERR_LOC =', ERR_LOC
                PRINT*, 'ERR_REQUIRED =', DELTA * T_OLD * KRYTOL
                PRINT*, 'DIMENSION CHANGED INTO M =', M
             ENDIF
             GO TO 101
          ENDIF
       ENDIF
       IMREJECT = 0
       JOLD = 1
       IF (ERR_LOC<1.0D-16) T_NEW = MAX(T_NEW, 2.0D0 * T_STEP)
       MX = MBRKDWN + MAX(0, K1 - 1)
       IREJECTFSP = 0

       ! ENSURE THE FSP-LIKE CRITERIA
       DO
          !     FIND W = BETA*V*EXP(T*H)*E1
          CALL DGEMV('N', N, MX, BETA, WSP(IV), N, WSP(IEXPH), 1, 0.0D0, W, 1)

          !     ENSURE FSP NON-NEGATIVITY AND COMPUTE THE PROBABILITY SUM
          DO I = 1, FSP%SIZE
             IF (W(I)<0.0D0)      W(I) = 0.0D0
          ENDDO
          WSUM = DASUM(FSP%SIZE, W, 1)

          PRINT*, "WSUM= ", WSUM

          ERROR = WSUM_OLD - WSUM
          T_RATIO = (T_NOW + T_STEP) / T_OUT

          !     ENSURE FSP-LIKE CRITERIA...
          IF (WSUM >= (1.0D0 - FERRORBOUND(T_NOW + T_STEP))) THEN
             EXIT
          ELSE
             IEXPAND = 1
             IREJECTFSP = IREJECTFSP + 1

             !           IF THE FSP HAS BEEN REJECTED TOO MANY TIMES, JUMPT TO SSA
             !           LOOP TO EXPAND THE STATE SPACE
             IF (IREJECTFSP>=5) THEN
                W(1:N) = BETA * WSP(IV:IV + N - 1)
                NSTEP = NSTEP - 1
                T_SSA = T_NEW
                GO TO 404
             ELSEIF (IREJECTFSP==1) THEN
                FSPORDER = 2
             ELSE
                FSPORDER = LOG(ERROR / ERROROLD) / &
                     LOG(T_STEP / TAU_OLD) - 1.0D0
             ENDIF

             !           REDUCE THE TIME STEP
             TFSP = GAMMA * T_STEP * (FSPTOL * T_STEP / (ERROR * T_OUT))**&
                  (1.0D0 / (FSPORDER))
             ERROROLD = ERROR
             TAU_OLD = T_STEP

             T_STEP = MIN(T_OUT - T_NOW, &
                  MAX(T_STEP / 5.0D0, MIN(0.9D0 * T_STEP, TFSP)))
             P1 = 10.0D0**(NINT(LOG10(T_STEP) - SQR1) - 1)
             T_STEP = AINT(T_STEP / P1 + 0.55D0) * P1

             NEXPH = NEXPH + 1
             CALL DGPADM(IDEG, MX, SGN * T_STEP, WSP(IH), MH, &
                  WSP(IFREE), LFREE, IWSP, IEXPH, NS, IFLAG)
             IEXPH = IFREE + IEXPH - 1
             NSCALE = NSCALE + NS
          ENDIF
       ENDDO

       !     UPDATE THE TIME COVERED
       T_NOW = T_NOW + T_STEP
       WSUM_OLD = WSUM

       !     DISPLAY INFORMATION ABOUT THE LAST SUCCESSFUL TIME STEP
       IF (ITRACE.NE.0) CALL PRINT_STATS

       !     IF THE INTEGRATION IS AT THE END TIME POINT, THEN THE ALGORITHM
       !     IS FINISHED
       IF (T_NOW>=T_OUT) GO TO 500

       !     REDUCE THE STATE SPACE TO ONLY SUBSTANTIAL STATES
       IF (NSTEP > 1 .AND. IEXPAND /= 1) THEN
          DSUM = WSUM - (1.0D0 - FERRORBOUND(T_NOW))
          IF (DSUM>0.0D0) CALL DROP_STATES(W, FSP, MODEL, DSUM, FMATVEC)
       ENDIF

       !-------------------------------------------STATE SPACE EXPANSION BY SSA
       !-----------------------------IN THE CASE WHERE THE STEPSIZE WAS REDUCED
404    CONTINUE

       IF ((IEXPAND==1).AND.(T_NOW<T_OUT)) THEN

          IF (NSTEP==1) T_NEW = T_STEP
          T_SSA = MIN(T_NEW, T_OUT - T_NOW)

          IF (ITRACE.NE.0) THEN
             PRINT*, 'CALLING SSA'
          ENDIF

          ! EXTEND THE STATE SPACE BY SSA AND 1-STEP REACHABILITY
          CALL SSA_EXTENDER(T_SSA, FSP, MODEL)
          CALL ONESTEP_EXTENDER(FSP, MODEL)

          IF (N_NOW>NFULL) N_NOW = NFULL
          IEXPAND = 0

       ENDIF

       !     ESTIMATE THE NUMBER OF NONZERO ENTRIES IN CURRENT FSP MATRIX
       NNZ = (MODEL%NREACTIONS + 1) * FSP%SIZE
       !     ESTIMATE OTHER PARAMETERS
       N_NOW = FSP%SIZE
       BETA = DNRM2(N_NOW, W, 1)
       HUMP = MAX(HUMP, BETA)
       ERR_LOC = MAX(ERR_LOC, RNDOFF)
       STEP_MIN = MIN(STEP_MIN, T_STEP)
       STEP_MAX = MAX(STEP_MAX, T_STEP)
       S_ERROR = S_ERROR + ERR_LOC
       X_ERROR = MAX(X_ERROR, ERR_LOC)
       P1 = 10.0D0**(NINT(LOG10(T_NEW) - SQR1) - 1)
       T_NEW = AINT(T_NEW / P1 + 0.55D0) * P1
       !     RETURN TO THE BEGINNING FOR THE NEXT TIME STEP
       IF ((MXSTEP==0).OR.(NSTEP<MXSTEP)) GO TO 100
       IFLAG = 1

       !-------------------------------------------OUTPUT IMPORTANT INFORMATION
500    CONTINUE
       !     IN IWSP
       IWSP(1) = NMULT
       IWSP(2) = NEXPH
       IWSP(3) = NSCALE
       IWSP(4) = NSTEP
       IWSP(5) = NREJECT
       IWSP(6) = IBRKFLAG
       IWSP(7) = MBRKDWN
       !     IN WSP
       WSP(1) = STEP_MIN
       WSP(2) = STEP_MAX
       WSP(3) = 0.0D0
       WSP(4) = 0.0D0
       WSP(5) = X_ERROR
       WSP(6) = S_ERROR
       WSP(7) = TBRKDWN
       WSP(8) = SGN * T_NOW
       WSP(9) = HUMP / VNORM
       WSP(10) = BETA / VNORM
     CONTAINS

       !===========================COMPUTE THE MATRIX-VECTOR PRODUCT IN THE CME
       SUBROUTINE FMATVEC(X, Y, MATRIX)
         IMPLICIT NONE
         !------------------------------------------------INPUT/OUTPUT PARAMETERS
         !     VECTOR X
         DOUBLE PRECISION :: X(*)
         !     VECTOR Y=A*X, WHERE A IS THE FSP MATRIX
         DOUBLE PRECISION :: Y(*)
         !     CURRENT FSP MATRIX
         TYPE(FSP_MATRIX) :: MATRIX
         !--------------------------------------------PARAMETERS OF THE ALGORITHM
         INTEGER :: I, J, K, N, BW
         !----------------------------------------------------------COMPUTE Y=A*X
         !     FIND THE SIZE OF THE MATRIX
         N = MATRIX%SIZE
         BW = SIZE(MATRIX%ADJ, DIM = 1)
         !     INITIALIZE VECTOR Y
         DO I = 1, N
            Y(I) = 0.0D0
         ENDDO

         !     COMPUTE Y=A*X
         DO I = 1, N
            DO J = 1, BW
               K = MATRIX%ADJ(J, I)
               IF (K>=1) THEN
                  Y(K) = Y(K) + MATRIX%OFFDIAG(J, I) * X(I)
               ENDIF
            ENDDO
            Y(I) = Y(I) - MATRIX%DIAG(I) * X(I)
         ENDDO
       END SUBROUTINE FMATVEC

       DOUBLE PRECISION FUNCTION FERRORBOUND(TX)
         IMPLICIT NONE
         !-------------------------------------------------------INPUT PARAMETERS
         !     A TIME INTERVAL
         DOUBLE PRECISION, INTENT(IN) :: TX
         !-----------------COMPUTE THE FSP ERROR TOLERANCE FOR THIS TIME INTERVAL
         FERRORBOUND = TX * FSPTOL / T_OUT
       END FUNCTION FERRORBOUND

       DOUBLE PRECISION FUNCTION KRYLOV_COST(T_NOW, T_OUT, TAU, M, N, HNORM)
         ! ESTIMATE THE COST OF THE KRYLOV APPROXIMATION ASSOCIATED WITH A PARTICULAR CHOICE OF STEPSIZE AND BASIS SIZE.
         ! T_NOW: CURRENT TIME
         ! T_OUT: FINAL TIME
         ! TAU: STEPSIZE
         ! M: KRYLOV BASIS SIZE
         ! N: SOLUTION VECTOR LENGTH
         ! HNORM: NORM OF THE HESSENBERG MATRIX COMPUTED FROM THE LAST TIME STEP
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: M, N
         DOUBLE PRECISION, INTENT(IN) :: TAU, HNORM, T_NOW, T_OUT
         DOUBLE PRECISION :: NOM

         NOM = 25.0D0 / 3.0D0 + &
              MAX(0, 2 + INT(LOG(TAU * HNORM) / LOG(2.0D0)))
         KRYLOV_COST = NINT((T_OUT - T_NOW) / TAU) * &
              (2 * (M + 1) * NNZ + &
              (5 * M + 4 * QIOP * M + 2 * QIOP - 2 * QIOP * QIOP + 7) * N + &
              2 * NOM * (M + 2) * (M + 2) * (M + 2))

       END FUNCTION KRYLOV_COST

       SUBROUTINE PRINT_STATS
         IMPLICIT NONE

         PRINT*, 'TIMESTEP', NSTEP, '-------------------------------'
         PRINT*, 'FSP SIZE         =', FSP%SIZE
         PRINT*, 'STEP_SIZE        =', T_STEP
         PRINT*, 'NEXT_STEP        =', T_NEW
         PRINT*, 'T_NOW            =', T_NOW
         PRINT*, 'KRYLOV DIMENSION =', M

       END SUBROUTINE PRINT_STATS

     END SUBROUTINE DGEXPV_FSP

   END MODULE KRYLOVSOLVER
