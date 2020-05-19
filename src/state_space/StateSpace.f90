MODULE STATESPACE
    !     CONTAINS GLOBAL VARIABLES, DATA STRUCTURES, TYPES AND DATA-RELATED SUBROUTINES REQUIRED IN THE KRYLOV-FSP-SSA ALGORITHM
    USE HASHTABLE
    USE BIG_INTEGER_MODULE
    USE MODELMODULE

    IMPLICIT NONE
    !     MAXIMUM NUMBER OF STATES ALLOWED IN THE FSP; MUST BE A PRIME
    !     NUMBER; SOME CHOICES: 809, 1319, 9431, 76667, 100009, 500009, 1162951, 6291469,15485867,32452843,23884543
    INTEGER, PARAMETER :: NMAX = 6291469
    INTEGER, PARAMETER :: MAXNUMBERMOLECULES = 10000

    TYPE FSP_MATRIX
        INTEGER :: SIZE
        INTEGER, ALLOCATABLE :: ADJ(:, :)      ! THE ADJACENCY MATRIX FOR THE MARKOV CHAIN
        DOUBLE PRECISION, ALLOCATABLE :: OFFDIAG(:, :), DIAG(:)
    END TYPE FSP_MATRIX

    TYPE :: FINITE_STATE_PROJECTION
        INTEGER :: MAX_SIZE = NMAX
        INTEGER :: SIZE
        INTEGER, ALLOCATABLE :: STATE(:, :)
        TYPE (BIG_INTEGER), ALLOCATABLE :: KEY(:)
        TYPE (FSP_MATRIX) :: MATRIX
        DOUBLE PRECISION, ALLOCATABLE :: VECTOR(:)

        !     THE LENGTH OF THE HASH TABLE
        INTEGER, PRIVATE :: KTLEN = NMAX
        !     VECTOR STORING THE KEYS OF THE STATES IN THE FSP
        TYPE(BIG_INTEGER), ALLOCATABLE :: KEYTAB(:)
        !     VECTOR STORING THE INDICES OF THE STATES IN THE CURRECT FSP
        INTEGER, ALLOCATABLE :: KVTAB(:)

        !     VECTOR STORING THE REACTION KEYS, SUCH THAT IF KEY(X) IS THE KEY
        !     OF STATE X THEN KEY(X+NU_K)=KEY(X)+RKEYSIGN(K)*REACTIONKEY(K)
        TYPE(BIG_INTEGER), ALLOCATABLE, PRIVATE :: REACTIONKEY(:)
        INTEGER, ALLOCATABLE, PRIVATE :: RKEYSIGN(:)

    CONTAINS
        PROCEDURE :: CREATE => CREATE_FSP
        PROCEDURE :: CLEAR => CLEAR_FSP
        PROCEDURE :: PROBABILITY => POINTWISE_FSP
        PROCEDURE :: ADD => ADD_STATE
        PROCEDURE :: INDEX => INDEX_STATE
    END TYPE FINITE_STATE_PROJECTION


    !----------------------------------------------------------------------|
CONTAINS
    !=============================== INITIALIZING AND CLEARING THE FSP
    SUBROUTINE CREATE_FSP(FSP, MODEL, MAX_SIZE_CUSTOM)
        IMPLICIT NONE

        INTEGER, OPTIONAL, INTENT(IN) :: MAX_SIZE_CUSTOM
        TYPE (CME_MODEL), INTENT(IN) :: MODEL
        CLASS (FINITE_STATE_PROJECTION), INTENT(INOUT) :: FSP

        INTEGER :: N, M

        N = MODEL%NSPECIES
        M = MODEL%NREACTIONS

        IF (PRESENT(MAX_SIZE_CUSTOM)) THEN
            FSP%MAX_SIZE = MAX_SIZE_CUSTOM
            FSP%KTLEN = MAX_SIZE_CUSTOM
        END IF

        ALLOCATE(FSP%STATE(N, FSP%MAX_SIZE), &
                FSP%MATRIX%DIAG(FSP%MAX_SIZE), &
                FSP%MATRIX%OFFDIAG(M, FSP%MAX_SIZE), &
                FSP%MATRIX%ADJ(M, FSP%MAX_SIZE), &
                FSP%KEY(FSP%MAX_SIZE), &
                FSP%REACTIONKEY(M), &
                FSP%RKEYSIGN(M), &
                FSP%KEYTAB(FSP%MAX_SIZE), &
                FSP%KVTAB(FSP%MAX_SIZE), &
                FSP%VECTOR(FSP%MAX_SIZE))

        FSP%SIZE = 0
        FSP%KEYTAB = BIG(0)
        CALL COMPUTE_RKEY(FSP%REACTIONKEY, FSP%RKEYSIGN, N, M, MODEL)

    END SUBROUTINE CREATE_FSP

    SUBROUTINE CLEAR_FSP(FSP)
        IMPLICIT NONE

        CLASS (FINITE_STATE_PROJECTION), INTENT(INOUT) :: FSP

        DEALLOCATE(FSP%STATE, FSP%MATRIX%DIAG, FSP%MATRIX%OFFDIAG, &
                FSP%MATRIX%ADJ, FSP%KEY, FSP%REACTIONKEY, FSP%RKEYSIGN, FSP%KEYTAB, FSP%KVTAB, FSP%VECTOR)

    END SUBROUTINE CLEAR_FSP

    !=============================== RETRIEVE THE POINT-WISE PROBABILITY FROM THE FSP
    DOUBLE PRECISION FUNCTION POINTWISE_FSP(FSP, X)
        IMPLICIT NONE
        CLASS (FINITE_STATE_PROJECTION), INTENT(IN) :: FSP
        INTEGER, INTENT(IN) :: X(:)

        TYPE(BIG_INTEGER) :: KEY
        LOGICAL :: FOUND
        INTEGER :: KA

        !           LOOK UP THE STATE IN THE FSP HASH TABLE
        CALL STATE2KEY(KEY, X, SIZE(FSP%STATE, 1), MAXNUMBERMOLECULES)
        CALL HASH(KEY, 1, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)

        IF (FOUND) THEN
            POINTWISE_FSP = FSP%VECTOR(FSP%KVTAB(KA))
        ELSE
            POINTWISE_FSP = 0.0D0
        ENDIF
    END FUNCTION POINTWISE_FSP
    !=============================== RETRIEVE THE STATE'S INDEX IN THE FSP
    INTEGER FUNCTION INDEX_STATE(FSP, X)
        IMPLICIT NONE
        CLASS (FINITE_STATE_PROJECTION), INTENT(IN) :: FSP
        INTEGER, INTENT(IN) :: X(:)

        TYPE(BIG_INTEGER) :: KEY
        LOGICAL :: FOUND
        INTEGER :: KA

        !           LOOK UP THE STATE IN THE FSP HASH TABLE
        CALL STATE2KEY(KEY, X, SIZE(FSP%STATE, 1), MAXNUMBERMOLECULES)
        CALL HASH(KEY, 1, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)

        IF (FOUND) THEN
            INDEX_STATE = FSP%KVTAB(KA)
        ELSE
            INDEX_STATE = 0
        ENDIF
    END FUNCTION INDEX_STATE
    !===========================ADD A STATE IN THE FSP HASH TABLE AND MATRIX
    SUBROUTINE ADD_STATE(FSP, MODEL, STATE, KEYIN)
        IMPLICIT NONE
        !------------------------------------------------INPUT/OUTPUT PARAMETERS
        !     THE CURRENT FSP HASH TABLE TO BE EXTENDED
        CLASS(FINITE_STATE_PROJECTION) :: FSP
        !     THE UNDERLYING MODEL
        TYPE (CME_MODEL), INTENT(IN) :: MODEL
        !     THE STATE TO BE ADDED
        INTEGER :: STATE(:)
        TYPE(BIG_INTEGER), INTENT(IN), OPTIONAL :: KEYIN

        !  PURPOSE:
        !  =======
        !
        !  ADD_STATE ADDS THE NEW STATE TO THE FSP AND UPDATES THE PRINCPAL SUBMATRIX
        !  ACCORDINGLY.
        !
        ! ARGUMENTS:
        ! =========
        !
        ! STATE          : (IN) THE STATE TO BE ADDED TO THE PROJECTION.
        !
        ! FSP            : (IN/OUT) THE FINITE STATE PROJECTION TO BE UPDATED BY
        !                  INCLUDING THE NEW STATE. THIS IS OF DERVIED TYPE FINITE_STATE_PROJECTION.
        !
        ! MODEL          : THE CME_MODEL OBJECT THAT DEFINES THE UNDERLYING STOCHASTIC MODEL
        !
        ! KEY            : (IN/OUT) THE KEY OF THE STATE TO BE ADDED. NOTE: IT HAS TO BE
        !                   ALREADY COMPUTED PRIOR TO THE CALL OF ADD_STATE.
        !
        ! GLOBAL VARIABLES AFFECTED:
        ! =========================
        !
        ! KEYTAB, KVTAB : VARIABLES FOR STORING THE HASH TABLE.
        !
        ! METHOD:
        ! ======
        !
        ! ESSENTIALLY THE SAME AS MATRIX_STARTER (SEE ABOVE), ONLY THAT ONLY ONE STATE IS ADDED.

        !     KEY OF THE STATE IN THE HASH TABLE
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        INTEGER :: I, K, RS(MODEL%NSPECIES), S, LSIZE
        LOGICAL :: FOUND
        INTEGER :: MODE, KA
        DOUBLE PRECISION :: AIJ
        TYPE(BIG_INTEGER) :: KEY

        ! IF THE KEY IS ALREADY COMPUTED, WE WILL USE IT...
        IF (PRESENT(KEYIN)) THEN
            KEY = KEYIN
        ELSE      ! ... ELSE WE WILL COMPUTE IT
            CALL STATE2KEY(KEY, STATE, MODEL%NSPECIES, MAXNUMBERMOLECULES)
        ENDIF
        !----------------------------------------ADD STATE TO THE FSP HASH TABLE
        !     UPDATE THE HASH TABLE
        MODE = 2
        CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
        IF ((.NOT.FOUND) .AND. KA/=0) THEN
            FSP%SIZE = FSP%SIZE + 1
            FSP%STATE(1:MODEL%NSPECIES, FSP%SIZE) = STATE
            FSP%KEY(FSP%SIZE) = KEY
            FSP%KVTAB(KA) = FSP%SIZE
            FSP%VECTOR(FSP%SIZE) = 0.0D0
            !--------------------------------------------------UPDATE THE FSP MATRIX
            LSIZE = FSP%SIZE
            FSP%MATRIX%SIZE = LSIZE
            MODE = 1
            !     UPDATE THE COLUMN IN THE FSP MATRIX
            FSP%MATRIX%DIAG(LSIZE) = 0.0D0
            !     FOR EACH REACTION...
            DO K = 1, MODEL%NREACTIONS
                !           FIND THE PROPENSITY OF THE REACTION
                AIJ = MODEL%PROPENSITY(STATE, K)
                !           UPDATE THE FSP MATRIX
                FSP%MATRIX%DIAG(LSIZE) = FSP%MATRIX%DIAG(LSIZE) + AIJ
                FSP%MATRIX%OFFDIAG(K, LSIZE) = AIJ
                !           FIND THE NEXT STATE ACCORDING TO THE REACTION
                RS = STATE + MODEL%STOICHIOMETRY(:, K)
                !           CHECK IF THE NEXT STATE IS ILLEGAL
                DO S = 1, MODEL%NSPECIES
                    IF (RS(S).LT.0) RS(1) = -1
                END DO
                IF (RS(1).GE.0) THEN
                    !                 FIND THE KEY OF THE NEXT STATE; NOTE THAT THE
                    !                 SUBROUTINE KEY2KEY IS SPECIFICALLY BUILT BY HUY
                    !                 FOR THIS ALGORITHM, BECAUSE IT'S FASTER TO ADD
                    !                 TWO KEYS THAN TO GENERATE A NEW KEY FROM THE
                    !                 STATE
                    CALL KEY2KEY(FSP%KEY(LSIZE), KEY, K, &
                            FSP%REACTIONKEY, FSP%RKEYSIGN, MODEL%NREACTIONS)
                    CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                    !                 UPDATE THE ADJACENCY MATRIX FOR THE MARKOV CHAIN
                    IF (FOUND) THEN
                        FSP%MATRIX%ADJ(K, LSIZE) = FSP%KVTAB(KA)
                    ELSE
                        FSP%MATRIX%ADJ(K, LSIZE) = 0
                    END IF
                ELSE
                    FSP%MATRIX%ADJ(K, LSIZE) = -1 ! ILLEGAL STATE, WE WILL NEVER EXPLORE THIS DIRECTION AGAIN
                END IF
            ENDDO
            !-----UPDATE THE ROWS IN THE ADJACENCY MATRIX FOR THE MARKOV CHAIN
            !     FOR EACH REACTION...
            DO K = 1, MODEL%NREACTIONS
                CALL KEY2KEYBW(FSP%KEY(LSIZE), KEY, K, FSP%REACTIONKEY, FSP%RKEYSIGN, MODEL%NREACTIONS)
                CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                IF (FOUND) FSP%MATRIX%ADJ(K, FSP%KVTAB(KA)) = LSIZE
            ENDDO
        ENDIF
    END SUBROUTINE ADD_STATE
    !===============================INITIALIZE THE FSP HASH TABLE AND MATRIX
    SUBROUTINE MATRIX_STARTER(FSP, MODEL)
        IMPLICIT NONE
        !------------------------------------------------INPUT/OUTPUT PARAMETERS
        !     THE MODEL
        TYPE(CME_MODEL), INTENT(IN) :: MODEL
        !     THE CURRENT FSP HASH TABLE
        TYPE(FINITE_STATE_PROJECTION) :: FSP
        ! PURPOSE:
        ! ========
        !
        ! MATRIX_STARTER GENERATES THE PRINCIPAL SUBMATRIX OF THE CME MATRIX BY
        ! KEEPING ONLY ENTRIES INDEXED BY THE STATES IN THE FSP.
        !
        ! ARGUMENTS:
        ! =========
        !
        ! FSP        : (INPUT/OUTPUT) THE FINITE STATE PROJECTION OF DERIVED TYPE
        !              FINITE_STATE_PROJECTION. THE KEY ATTRIBUTE IS ALSO UPDATED
        !              AFTER THE CALL OF THE SUBROUTINE. (SEE MODULE FSP_MANAGER)
        !
        ! MODEL      : (INPUT) THE MODEL UNDERLYING THE CME
        !
        ! METHOD:
        ! =======
        !
        ! THE OUTER LOOP SCANS THE FSP FROM THE FIRST TO THE LAST STATE.
        ! FOR EACH STATE THE SUBROUTINE COMPUTES THE WHOLE COLUMN OF THE CME
        ! MATRIX AND STORES IT IN FSP%DIAG AND FSP%MATRIX%OFFDIAG. IT ALSO UPDATES THE
        ! VARIABLE FSP%MATRIX%ADJ (ADJ: ADJACENCY) TO 'CONNECT' THE PREVIOUSLY ADDED STATES
        ! WITH THE NEW STATE.
        !
        ! IN OTHER WORDS, THE SUBROUTINE BUILDS THE (SPARSE) PRINCIPAL SUBMATRIX THAT EXTENDS
        ! BY ONE ROW AND ONE COLUMN AFTER EACH ITERATION OF THE OUTER LOOP.
        !
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        INTEGER :: SD, PD, I, K, S, RS(MODEL%NSPECIES), STATE(MODEL%NSPECIES), MODE, KA
        LOGICAL :: FOUND
        TYPE(BIG_INTEGER) :: KEY

        !-----------------------------------UPDATE THE FSP HASH TABLE AND MATRIX
        SD = MODEL%NSPECIES
        PD = MODEL%NREACTIONS
        !     FIND THE CURRENT FSP STATE SPACE SIZE
        FSP%MATRIX%SIZE = FSP%SIZE
        DO I = 1, FSP%SIZE
            !           FIND THE STATE TO ADD IN
            STATE(1:SD) = FSP%STATE(1:SD, I)
            !           ADD THE STATE IN THE FSP HASH TABLE
            MODE = 2
            CALL STATE2KEY(KEY, STATE, SD, MAXNUMBERMOLECULES)
            CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
            FSP%KVTAB(KA) = I
            FSP%KEY(I) = KEY
            !-----------UPDATE THE COLUMN IN THE FSP MATRIX
            MODE = 1
            FSP%MATRIX%DIAG(I) = 0.0D0
            !           FOR EACH REACTION...
            DO K = 1, PD
                !                 FIND THE NEXT STATE ACCORDING TO THE REACTION
                RS = STATE + MODEL%STOICHIOMETRY(:, K)
                !                 CHECK IF THE NEXT STATE IS ILLEGAL
                DO S = 1, SD
                    IF (RS(S).LT.0) RS(1) = -1
                ENDDO
                !                 UPDATE THE COLUMN IN THE FSP MATRIX
                FSP%MATRIX%DIAG(I) = FSP%MATRIX%DIAG(I) + MODEL%PROPENSITY(STATE, K)
                FSP%MATRIX%OFFDIAG(K, I) = MODEL%PROPENSITY(STATE, K)
                !                 UPDATE THE ADJACENCY MATRIX FOR THE MARKOV CHAIN
                IF (RS(1).GE.0) THEN
                    CALL STATE2KEY(KEY, RS, SD, MAXNUMBERMOLECULES)
                    CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                    IF (FOUND) THEN
                        FSP%MATRIX%ADJ(K, I) = FSP%KVTAB(KA)
                    ELSE
                        FSP%MATRIX%ADJ(K, I) = 0
                    ENDIF
                ELSE
                    FSP%MATRIX%ADJ(K, I) = -1 ! THE STATE IS ILLEGAL
                ENDIF
            ENDDO
            !-----------UPDATE THE ROWS IN THE ADJACENCY MATRIX FOR THE MARKOV CHAIN
            !           FOR EACH REACTION...
            DO K = 1, PD
                !                 FIND THE PREVIOUS STATE ACCORDING TO THE REACTION
                RS = STATE - MODEL%STOICHIOMETRY(:, K)
                !                 CHECK IF THE PREVIOUS STATE IS ILLEGAL
                DO S = 1, SD
                    IF (RS(S).LT.0) RS(1) = -1
                ENDDO
                !                 UPDATE THE ADJACENCY MATRIX FOR THE MARKOV CHAIN
                IF (RS(1).GE.0) THEN
                    CALL STATE2KEY(KEY, RS, SD, MAXNUMBERMOLECULES)
                    CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                    IF (FOUND) FSP%MATRIX%ADJ(K, FSP%KVTAB(KA)) = I
                ENDIF
            ENDDO
        ENDDO
    END SUBROUTINE MATRIX_STARTER
    !===================================EXTEND THE FSP STATE SPACE BY 1-STEP
    SUBROUTINE ONESTEP_EXTENDER(FSP, MODEL)
        IMPLICIT NONE
        !------------------------------------------------INPUT/OUTPUT PARAMETERS
        !     THE CURRENT FSP HASH TABLE TO BE EXTENDED
        TYPE(FINITE_STATE_PROJECTION) :: FSP
        !     THE UNDERLYING MODEL
        TYPE(CME_MODEL), INTENT(IN) :: MODEL
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        INTEGER :: SD, PD
        INTEGER :: I, J, K, RS(MODEL%NSPECIES), STATE(MODEL%NSPECIES), LSIZE, LSIZE_COPY
        LOGICAL :: FOUND
        INTEGER :: MODE, KA
        PARAMETER(MODE = 1)
        TYPE(BIG_INTEGER) :: KEY
        !-----------------------------------UPDATE THE FSP HASH TABLE AND MATRIX
        SD = MODEL%NSPECIES
        PD = MODEL%NREACTIONS
        !     FIND THE CURRENT FSP STATE SPACE SIZE
        LSIZE_COPY = FSP%SIZE
        DO J = 1, LSIZE_COPY
            !           FIND THE CURRENT STATE
            STATE = FSP%STATE(:, J)
            !           FOR EACH REACTION...
            DO K = 1, PD
                IF (FSP%MATRIX%ADJ(K, J).EQ.0) THEN
                    !                       FIND THE NEXT STATE ACCORDING TO THE REACTION
                    RS(:) = STATE + MODEL%STOICHIOMETRY(:, K)

                    CALL KEY2KEY(FSP%KEY(J), KEY, &
                            K, FSP%REACTIONKEY, FSP%RKEYSIGN, PD)
                    !                       CALL STATE2KEY(KEY,RS(1:SD,K),SD,MAXNUMBERMOLECULES)
                    CALL HASH(KEY, MODE, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                    !                       IF THE NEXT STATE IS ALREADY IN THE FSP...
                    IF (FOUND) THEN
                        !                             UPDATE THE ADJACENCY MATRIX FOR THE MARKOV
                        !                             CHAIN, OTHERWISE...
                        FSP%MATRIX%ADJ(K, J) = FSP%KVTAB(KA)
                    ELSE
                        !                             ADD THE STATE IN THE FSP HASH TABLE
                        CALL FSP%ADD(MODEL, RS, KEYIN = KEY)
                        !                             EXIT IF THE HASH TABLE LENGTH IS TOO SMALL
                        IF (FSP%SIZE.GE.FSP%KTLEN) THEN
                            STOP 'OVERFLOW ERROR: FSP SIZE EXCEEDS MEMORY LIMIT.'
                            RETURN
                        END IF
                    END IF
                END IF
            END DO
        END DO
    END SUBROUTINE ONESTEP_EXTENDER

    SUBROUTINE FIND_DROPTOL(SD, LSIZE, W, DROPTOL, DSUM)
        IMPLICIT NONE
        !------------------------------------------------INPUT/OUTPUT PARAMETERS
        !     NUMBER OF SPECIES IN THE MODEL
        INTEGER :: SD
        !     SIZE OF THE CURRENT FSP STATE SPACE
        INTEGER :: LSIZE
        !     THE TRUNCATED VECTOR APPROXIMATING THE PROBABILITY VECTOR
        DOUBLE PRECISION :: W(:)
        !     THE TRUNCATION THRESHOLD THAT SATISFIES THE FSP CONDITION
        DOUBLE PRECISION :: DROPTOL
        !     MAXIMUM ALLOWABLE VALUE FOR THE SUM OF STATES MARKED AS
        !     'INSIGNIFICANT'
        DOUBLE PRECISION :: DSUM
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        DOUBLE PRECISION :: SUM1
        INTEGER :: I, J

        DROPTOL = 1.0D-08
        DO
            SUM1 = 0.0D0
            DO I = 1, LSIZE
                IF (W(I).LT.DROPTOL .AND. W(I).GT.0) THEN
                    SUM1 = SUM1 + W(I)
                ENDIF
            ENDDO
            IF (SUM1.LT.DSUM) EXIT
            DROPTOL = DROPTOL / 10.0D0
        ENDDO
    END SUBROUTINE FIND_DROPTOL


    !============DROP STATES IN THE FSP STATE SPACE WITH SMALL PROBABILITIES
    SUBROUTINE DROP_STATES(W, FSP, MODEL, DSUM, FMATVEC)
        !     USE FSP_MANAGER

        IMPLICIT NONE
        !------------------------------------------------INPUT/OUTPUT PARAMETERS
        !     THE TRUNCATED VECTOR APPROXIMATING THE PROBABILITY VECTOR
        DOUBLE PRECISION :: W(:)
        !     THE CURRENT FSP HASH TABLE TO BE EXTENDED
        TYPE(FINITE_STATE_PROJECTION) :: FSP

        DOUBLE PRECISION :: DSUM
        !     THE UNDERLYING MODEL
        TYPE (CME_MODEL), INTENT(IN) :: MODEL
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        INTEGER :: SD, PD
        INTEGER, ALLOCATABLE :: ARS(:, :), NEW_INDEX(:), LIST(:, :)
        DOUBLE PRECISION, ALLOCATABLE :: APROPS(:, :), ADIAG(:), W_COPY(:)
        INTEGER :: I, J, K, Q, MODE, KA, LSIZE
        TYPE(BIG_INTEGER), ALLOCATABLE :: KEYLISTCOPY(:)

        LOGICAL :: FOUND
        TYPE(BIG_INTEGER) :: KEY

        !     VECTOR CONTAINING THE INDICES TO BE DROPPED
        LOGICAL, ALLOCATABLE :: DROP(:)
        DOUBLE PRECISION, ALLOCATABLE :: WTMP(:)
        DOUBLE PRECISION :: DROPTOL

        INTEGER :: DROP_COUNT
        !----------------------------------------------INITIALIZE THE PARAMETERS
        SD = MODEL%NSPECIES
        PD = MODEL%NREACTIONS

        LSIZE = FSP%SIZE

        ALLOCATE(DROP(LSIZE), WTMP(LSIZE))

        !                       FIND THE TRUNCATION THRESHOLD THAT SATISFIES THE
        !                       PROBABILITY SUM CONDITION
        CALL FIND_DROPTOL(MODEL%NSPECIES, FSP%SIZE, W, DROPTOL, DSUM)


        !                       MARK THE STATES WITH PROBABILITIES BELOW THE
        !                       THRESHOLD TO BE DROPPED
        DROP_COUNT = 0
        DO I = 1, FSP%SIZE
            DROP(I) = (W(I).LT.DROPTOL)
            IF (W(I) < DROPTOL) THEN
                DROP(I) = .TRUE.
                DROP_COUNT = DROP_COUNT + 1
            ELSE
                DROP(I) = .FALSE.
            END IF
        ENDDO

        CALL FMATVEC(W, WTMP, FSP%MATRIX)
        !                       MAKE SURE THAT THE STATES WITH BIG POSITIVE
        !                       DERIVATIVES, EVEN WITH SMALL PROBABILITIES, ARE
        !                       NOT DROPPED
        DO I = 1, FSP%SIZE
            IF (WTMP(I).GT.1.0D-8) THEN
                DROP(I) = .FALSE.
                DROP_COUNT = DROP_COUNT - 1
            END IF
        ENDDO

        IF (DROP_COUNT*1.0D0/(LSIZE*1.0D0) > 0.1D0) THEN

            ALLOCATE(ARS(PD, LSIZE), NEW_INDEX(LSIZE), LIST(SD, LSIZE), &
                    APROPS(PD, LSIZE), ADIAG(LSIZE), KEYLISTCOPY(LSIZE), &
                    W_COPY(LSIZE))
            W_COPY = 0D0
            Q = 0
            !---------------------------------DROP THE STATES IN THE FSP STATE SPACE
            !     FOR EACH STATE IN THE FSP STATE SPACE...
            DO J = 1, LSIZE
                !           CHECK IF THE STATE IS TO BE DROPPED; IF NOT...
                IF (.NOT.DROP(J)) THEN
                    !                 STORE THE STATE IN A NEW TEMPORARY HASH TABLE,
                    !                 OTHERWISE...
                    Q = Q + 1
                    W_COPY(Q) = W(J)
                    LIST(1:SD, Q) = FSP%STATE(:, J)
                    ADIAG(Q) = FSP%MATRIX%DIAG(J)
                    APROPS(1:PD, Q) = FSP%MATRIX%OFFDIAG(1:PD, J)
                    ARS(1:PD, Q) = FSP%MATRIX%ADJ(1:PD, J)
                    KEYLISTCOPY(Q) = FSP%KEY(J)
                    NEW_INDEX(J) = Q
                    CALL HASH(FSP%KEY(J), 1, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                    FSP%KVTAB(KA) = Q
                ELSE
                    !                 DELETE THE STATE
                    NEW_INDEX(J) = 0
                    CALL HASH(FSP%KEY(J), 3, KA, FOUND, FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                ENDIF
            ENDDO
            !     CLEAN THE PROBABILITY VECTOR
            W(1:LSIZE) = 0.0D0
            LSIZE = Q
            !     COPY THE TEMPORARY HASH TABLE TO THE FSP HASH TABLE
            FSP%SIZE = LSIZE
            FSP%STATE(1:SD, 1:LSIZE) = LIST(1:SD, 1:LSIZE)
            FSP%KEY(1:LSIZE) = KEYLISTCOPY(1:LSIZE)
            W(1:LSIZE) = W_COPY(1:LSIZE)
            FSP%MATRIX%SIZE = LSIZE
            FSP%MATRIX%ADJ(1:PD, 1:LSIZE) = ARS(1:PD, 1:LSIZE)
            FSP%MATRIX%DIAG(1:LSIZE) = ADIAG(1:LSIZE)
            FSP%MATRIX%OFFDIAG(1:PD, 1:LSIZE) = APROPS(1:PD, 1:LSIZE)
            !     RE-INDEX THE ADJACENCY MATRIX
            DO J = 1, LSIZE
                DO K = 1, PD
                    I = FSP%MATRIX%ADJ(K, J)
                    IF (I>0) FSP%MATRIX%ADJ(K, J) = NEW_INDEX(I)
                ENDDO
            ENDDO
        ENDIF
        ! ARRAYS ARE AUTOMATICALLY DEALLOCATED
    END SUBROUTINE DROP_STATES
    !======================================EXTEND THE FSP STATE SPACE BY SSA
    SUBROUTINE SSA_EXTENDER(TIMESTEP, FSP, MODEL)

        IMPLICIT NONE
        !------------------------------------------------INPUT/OUTPUT PARAMETERS
        DOUBLE PRECISION :: TIMESTEP
        TYPE(FINITE_STATE_PROJECTION) :: FSP
        TYPE (CME_MODEL), INTENT(IN) :: MODEL
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        INTEGER :: SD, PD, I, J, J0, K, S, RS(MODEL%NSPECIES), STATE(MODEL%NSPECIES)
        DOUBLE PRECISION :: T, R1, R2, R2A
        DOUBLE PRECISION :: TMP
        INTEGER :: LSIZE_OLD
        LOGICAL :: FOUND
        INTEGER :: MODE, KA
        TYPE(BIG_INTEGER) :: KEY
        !--------------------------------------EXTEND THE FSP STATE SPACE BY SSA
        SD = MODEL%NSPECIES
        PD = MODEL%NREACTIONS
        MODE = 1
        LSIZE_OLD = FSP%SIZE
        !     FROM EACH STATE IN THE CURRENT FSP STATE SPACE...
        DO J0 = 1, LSIZE_OLD
            J = J0
            STATE(1:SD) = FSP%STATE(1:SD, J)
            T = 0.0D0
            !           GENERATE THE NEXT STATE IN THE SSA PATH
            300    CONTINUE
            CALL RANDOM_NUMBER(R1)
            CALL RANDOM_NUMBER(R2)
            T = MIN(TIMESTEP, T + (-LOG(R1) / FSP%MATRIX%DIAG(J)))
            IF (T<=TIMESTEP) THEN
                TMP = FSP%MATRIX%OFFDIAG(1, J)
                K = 1
                R2A = MIN(R2 * FSP%MATRIX%DIAG(J), FSP%MATRIX%DIAG(J))
                301       IF (TMP<R2A .AND. K<PD) THEN
                    K = K + 1
                    TMP = TMP + FSP%MATRIX%OFFDIAG(K, J)
                    GO TO 301
                ENDIF
                RS = STATE + MODEL%STOICHIOMETRY(:, K)
                DO S = 1, SD
                    IF (RS(S).LT.0) RS(1) = -1
                ENDDO
                !                 IF THE NEXT STATE IS ILLEGAL...
                IF (RS(1).LT.0) THEN
                    !                       IGNORE IT, OTHERWISE...
                    FSP%MATRIX%ADJ(K, J) = -1
                ELSE
                    !                       STORE THE NEXT STATE IN THE FSP STATE SPACE
                    IF (FSP%MATRIX%ADJ(K, J).EQ.0) THEN

                        CALL KEY2KEY(FSP%KEY(J), KEY, K, &
                                FSP%REACTIONKEY, FSP%RKEYSIGN, PD)
                        CALL HASH(KEY, MODE, KA, FOUND, &
                                FSP%KTLEN, FSP%KEYTAB, FSP%KVTAB)
                        !                             CHECK IF KEY IS ALREADY IN HASHTABLE, IF
                        !                             NOT...
                        IF (FOUND) THEN
                            J = FSP%KVTAB(KA)
                        ELSE
                            !                                   ADD THE STATE TO THE LIST AND UPDATE
                            !                                   THE FSP HASH TABLE
                            IF (FSP%SIZE.GE.FSP%KTLEN) THEN
!                                PRINT*, 'OVERFLOW ERROR:' &
!                                 ' LSIZE>=N'
                                RETURN
                            ENDIF
                            CALL FSP%ADD(MODEL, RS, KEYIN = KEY)
                            J = FSP%SIZE
                        ENDIF
                    ELSE
                        I = FSP%MATRIX%ADJ(K, J)
                        J = I
                    ENDIF
                    !                       MOVE TO THE NEXT STEP ON THE SSA PATH
                    STATE(1:SD) = FSP%STATE(1:SD, J)
                    IF ((T.LT.TIMESTEP).AND.(J.GE.J0)) GO TO 300
                ENDIF
            ENDIF
        ENDDO
    END SUBROUTINE SSA_EXTENDER


    !=============================COMPUTING THE ABSOLUTE VALUE OF KEY CHANGE
    !==========================================DUE TO THE CHEMICAL REACTIONS
    SUBROUTINE COMPUTE_RKEY(REACTIONKEY, RKEYSIGN, SD, PD, MODEL)
        USE BIG_INTEGER_MODULE
        IMPLICIT NONE
        !------------------------------------------------------OUTPUT PARAMETERS
        !     VECTOR STORING THE REACTION KEYS, SUCH THAT IF KEY(X) IS THE KEY
        !     OF STATE X THEN KEY(X+NU_K)=KEY(X)+RKEYSIGN(K)*REACTIONKEY(K)
        INTEGER :: SD, PD
        TYPE(BIG_INTEGER) :: REACTIONKEY(PD)
        INTEGER :: RKEYSIGN(PD)
        TYPE(CME_MODEL), INTENT(IN) :: MODEL
        !--------------------------------------------PARAMETERS OF THE ALGORITHM
        INTEGER :: I, J, STATE(1:SD), RS(1:SD), SGN
        TYPE(BIG_INTEGER) :: RKEY, C
        !-----------------------------COMPUTING THE ABSOLUTE VALUE OF KEY CHANGE
        C = BIG(MAXNUMBERMOLECULES)
        STATE(1:SD) = 0
        !     FROM THE ZERO STATE, FOR EACH REACTION...
        DO J = 1, PD
            !           FIND THE NEXT STATE ACCORDING TO THE REACTION
            RS = STATE + MODEL%STOICHIOMETRY(:, J)
            SGN = 1
            RKEY = BIG(0)
            !           FIND THE REACTION KEY
            DO I = 1, SD
                IF (SGN * RS(I).LT.0) THEN
                    SGN = -1 * SGN
                    RKEY = ABS(RS(I)) * (C + 1)**(I - 1) - RKEY
                ELSE
                    RKEY = ABS(RS(I)) * (C + 1)**(I - 1) + RKEY
                ENDIF
            ENDDO
            REACTIONKEY(J) = RKEY
            RKEYSIGN(J) = SGN
        ENDDO
    END SUBROUTINE COMPUTE_RKEY

END MODULE STATESPACE
