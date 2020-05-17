MODULE StateSpace
    !     Contains global variables, data structures, types and data-related subroutines required in the Krylov-FSP-SSA algorithm
    USE HashTable
    USE big_integer_module
    USE ModelModule

    IMPLICIT NONE
    !     maximum number of states allowed in the FSP; must be a prime
    !     number; some choices: 809, 1319, 9431, 76667, 100009, 500009, 1162951, 6291469,15485867,32452843,23884543
    INTEGER, PARAMETER :: nmax = 6291469
    INTEGER, PARAMETER :: MaxNumberMolecules = 10000

    TYPE fsp_matrix
        INTEGER :: size
        INTEGER, ALLOCATABLE :: adj(:, :)      ! the adjacency matrix for the Markov chain
        DOUBLE PRECISION, ALLOCATABLE :: offDiag(:, :), Diag(:)
    END TYPE fsp_matrix

    TYPE :: finite_state_projection
        INTEGER :: max_size = nmax
        INTEGER :: size
        INTEGER, ALLOCATABLE :: state(:, :)
        TYPE (big_integer), ALLOCATABLE :: key(:)
        TYPE (fsp_matrix) :: matrix
        DOUBLE PRECISION, ALLOCATABLE :: vector(:)

        !     The length of the hash table
        INTEGER, PRIVATE :: ktlen = nmax
        !     Vector storing the keys of the states in the FSP
        TYPE(big_integer), ALLOCATABLE :: keytab(:)
        !     Vector storing the indices of the states in the currect FSP
        INTEGER, ALLOCATABLE :: kvtab(:)

        !     Vector storing the reaction keys, such that if key(x) is the key
        !     of state x then key(x+nu_k)=key(x)+rkeysign(k)*reactionkey(k)
        TYPE(big_integer), ALLOCATABLE, PRIVATE :: reactionkey(:)
        INTEGER, ALLOCATABLE, PRIVATE :: rkeysign(:)

    CONTAINS
        PROCEDURE :: init => init_fsp
        PROCEDURE :: clear => clear_fsp
        PROCEDURE :: probability => pointwise_fsp
        PROCEDURE :: add => add_state
        PROCEDURE :: index => index_state
    END TYPE finite_state_projection


    !----------------------------------------------------------------------|
CONTAINS
    !=============================== Initializing and clearing the fsp
    SUBROUTINE init_fsp(fsp, model, max_size_custom)
        IMPLICIT NONE

        INTEGER, OPTIONAL, INTENT(in) :: max_size_custom
        TYPE (cme_model), INTENT(in) :: model
        class (finite_state_projection), INTENT(inout) :: fsp

        INTEGER :: N, M

        N = model%nspecies
        M = model%nreactions

        IF (present(max_size_custom)) THEN
            fsp%max_size = max_size_custom
            fsp%ktlen = max_size_custom
        END IF

        ALLOCATE(fsp%state(N, fsp%max_size), &
                fsp%matrix%diag(fsp%max_size), &
                fsp%matrix%offdiag(M, fsp%max_size), &
                fsp%matrix%adj(M, fsp%max_size), &
                fsp%key(fsp%max_size), &
                fsp%reactionkey(M), &
                fsp%rkeysign(M), &
                fsp%keytab(fsp%max_size), &
                fsp%kvtab(fsp%max_size), &
                fsp%vector(fsp%max_size))

        fsp%size = 0
        fsp%keytab = big(0)
        CALL compute_rkey(fsp%reactionkey, fsp%rkeysign, N, M, model)

    END SUBROUTINE init_fsp

    SUBROUTINE clear_fsp(fsp)
        IMPLICIT NONE

        class (finite_state_projection), INTENT(inout) :: fsp

        DEALLOCATE(fsp%state, fsp%matrix%diag, fsp%matrix%offdiag, &
                fsp%matrix%adj, fsp%key, fsp%reactionkey, fsp%rkeysign, fsp%keytab, fsp%kvtab, fsp%vector)

    END SUBROUTINE clear_fsp

    !=============================== Retrieve the point-wise probability from the fsp
    DOUBLE PRECISION FUNCTION pointwise_fsp(fsp, X)
        IMPLICIT NONE
        class (finite_state_projection), INTENT(in) :: fsp
        INTEGER, INTENT(in) :: X(:)

        TYPE(big_integer) :: key
        LOGICAL :: found
        INTEGER :: ka

        !           Look up the state in the FSP hash table
        CALL state2key(key, X, SIZE(fsp%state, 1), MaxNumberMolecules)
        CALL hash(key, 1, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)

        IF (found) THEN
            pointwise_fsp = fsp%vector(fsp%kvtab(ka))
        ELSE
            pointwise_fsp = 0.0d0
        ENDIF
    END FUNCTION pointwise_fsp
    !=============================== Retrieve the state's index in the FSP
    INTEGER FUNCTION index_state(fsp, X)
        IMPLICIT NONE
        class (finite_state_projection), INTENT(in) :: fsp
        INTEGER, INTENT(in) :: X(:)

        TYPE(big_integer) :: key
        LOGICAL :: found
        INTEGER :: ka

        !           Look up the state in the FSP hash table
        CALL state2key(key, X, SIZE(fsp%state, 1), MaxNumberMolecules)
        CALL hash(key, 1, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)

        IF (found) THEN
            index_state = fsp%kvtab(ka)
        ELSE
            index_state = 0
        ENDIF
    END FUNCTION index_state
    !===========================Add a state in the FSP hash table and matrix
    SUBROUTINE add_state(fsp, model, state, keyin)
        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        !     The current FSP hash table to be extended
        class(finite_state_projection) :: fsp
        !     The underlying model
        TYPE (cme_model), INTENT(in) :: model
        !     The state to be added
        INTEGER :: state(:)
        TYPE(big_integer), INTENT(in), OPTIONAL :: keyin

        !  Purpose:
        !  =======
        !
        !  add_state adds the new state to the FSP and updates the princpal submatrix
        !  accordingly.
        !
        ! Arguments:
        ! =========
        !
        ! state          : (in) the state to be added to the projection.
        !
        ! fsp            : (in/out) the finite state projection to be updated by
        !                  including the new state. This is of dervied type finite_state_projection.
        !
        ! model          : the cme_model object that defines the underlying stochastic model
        !
        ! KEY            : (in/out) the key of the state to be added. NOTE: It has to be
        !                   already computed prior to the call of add_state.
        !
        ! Global variables affected:
        ! =========================
        !
        ! KEYTAB, KVTAB : variables for storing the hash table.
        !
        ! Method:
        ! ======
        !
        ! Essentially the same as matrix_starter (see above), only that only one state is added.

        !     Key of the state in the hash table
        !--------------------------------------------Parameters of the algorithm
        INTEGER :: i, k, rs(model%nspecies), s, lsize
        LOGICAL :: found
        INTEGER :: mode, ka
        DOUBLE PRECISION :: aij
        TYPE(big_integer) :: key

        ! If the key is already computed, we will use it...
        IF (PRESENT(keyin)) THEN
            key = keyin
        ELSE      ! ... Else we will compute it
            CALL state2key(key, state, model%nspecies, MaxNumberMolecules)
        ENDIF
        !----------------------------------------Add state to the FSP hash table
        !     Update the hash table
        mode = 2
        CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
        IF ((.NOT.found) .AND. ka/=0) THEN
            fsp%size = fsp%size + 1
            fsp%state(1:model%nspecies, fsp%size) = state
            fsp%key(fsp%size) = key
            fsp%kvtab(ka) = fsp%size
            fsp%vector(fsp%size) = 0.0d0
            !--------------------------------------------------Update the FSP matrix
            lsize = fsp%size
            fsp%matrix%size = lsize
            mode = 1
            !     Update the column in the FSP matrix
            fsp%matrix%diag(lsize) = 0.0d0
            !     For each reaction...
            DO k = 1, model%nreactions
                !           Find the propensity of the reaction
                aij = model%propensity(state, k)
                !           Update the FSP matrix
                fsp%matrix%diag(lsize) = fsp%matrix%diag(lsize) + aij
                fsp%matrix%offdiag(k, lsize) = aij
                !           Find the next state according to the reaction
                rs = state + model%stoichiometry(:, k)
                !           Check if the next state is illegal
                DO s = 1, model%nspecies
                    IF (rs(s).LT.0) rs(1) = -1
                END DO
                IF (rs(1).GE.0) THEN
                    !                 Find the key of the next state; note that the
                    !                 SUBROUTINE key2key is specifically built by Huy
                    !                 for this algorithm, because it's faster to add
                    !                 two keys than to generate a new key from the
                    !                 state
                    CALL key2key(fsp%key(lsize), key, k, &
                            fsp%reactionkey, fsp%rkeysign, model%nreactions)
                    CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                    !                 Update the adjacency matrix for the Markov chain
                    IF (found) THEN
                        fsp%matrix%adj(k, lsize) = fsp%kvtab(ka)
                    ELSE
                        fsp%matrix%adj(k, lsize) = 0
                    END IF
                ELSE
                    fsp%matrix%adj(k, lsize) = -1 ! illegal state, we will never explore this direction again
                END IF
            ENDDO
            !-----Update the rows in the adjacency matrix for the Markov chain
            !     For each reaction...
            DO k = 1, model%nreactions
                CALL key2keybw(fsp%key(lsize), key, k, fsp%reactionkey, fsp%rkeysign, model%nreactions)
                CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                IF (found) fsp%matrix%adj(k, fsp%kvtab(ka)) = lsize
            ENDDO
        ENDIF
    END SUBROUTINE add_state
    !===============================Initialize the FSP hash table and matrix
    SUBROUTINE matrix_starter(fsp, model)
        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        !     The model
        TYPE(cme_model), INTENT(in) :: model
        !     The current FSP hash table
        TYPE(finite_state_projection) :: fsp
        ! Purpose:
        ! ========
        !
        ! matrix_starter generates the principal submatrix of the CME matrix by
        ! keeping only entries indexed by the states in the FSP.
        !
        ! Arguments:
        ! =========
        !
        ! fsp        : (input/output) the finite state projection of derived type
        !              finite_state_projection. The key attribute is also updated
        !              after the call of the subroutine. (see module fsp_manager)
        !
        ! model      : (input) the model underlying the CME
        !
        ! Method:
        ! =======
        !
        ! The outer loop scans the fsp from the first to the last state.
        ! For each state the subroutine computes the whole column of the CME
        ! matrix and stores it in fsp%diag and fsp%matrix%offdiag. It also updates the
        ! variable fsp%matrix%adj (adj: adjacency) to 'connect' the previously added states
        ! with the new state.
        !
        ! In other words, the subroutine builds the (sparse) principal submatrix that extends
        ! by one row and one column after each iteration of the outer loop.
        !
        !--------------------------------------------Parameters of the algorithm
        INTEGER :: sd, pd, i, k, s, rs(model%nspecies), state(model%nspecies), mode, ka
        LOGICAL :: found
        TYPE(big_integer) :: key

        !-----------------------------------Update the FSP hash table and matrix
        sd = model%nspecies
        pd = model%nreactions
        !     Find the current FSP state space size
        fsp%matrix%size = fsp%size
        DO i = 1, fsp%size
            !           Find the state to add in
            state(1:sd) = fsp%state(1:sd, i)
            !           Add the state in the FSP hash table
            mode = 2
            CALL state2key(key, state, sd, MaxNumberMolecules)
            CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
            fsp%kvtab(ka) = i
            fsp%key(i) = key
            !-----------Update the column in the FSP matrix
            mode = 1
            fsp%matrix%diag(i) = 0.0d0
            !           For each reaction...
            DO k = 1, pd
                !                 Find the next state according to the reaction
                rs = state + model%stoichiometry(:, k)
                !                 Check if the next state is illegal
                DO s = 1, sd
                    IF (rs(s).LT.0) rs(1) = -1
                ENDDO
                !                 Update the column in the FSP matrix
                fsp%matrix%diag(i) = fsp%matrix%diag(i) + model%propensity(state, k)
                fsp%matrix%offdiag(k, i) = model%propensity(state, k)
                !                 Update the adjacency matrix for the Markov chain
                IF (rs(1).GE.0) THEN
                    CALL state2key(key, rs, sd, MaxNumberMolecules)
                    CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                    IF (found) THEN
                        fsp%matrix%adj(k, i) = fsp%kvtab(ka)
                    ELSE
                        fsp%matrix%adj(k, i) = 0
                    ENDIF
                ELSE
                    fsp%matrix%adj(k, i) = -1 ! The state is illegal
                ENDIF
            ENDDO
            !-----------Update the rows in the adjacency matrix for the Markov chain
            !           For each reaction...
            DO k = 1, pd
                !                 Find the previous state according to the reaction
                rs = state - model%stoichiometry(:, k)
                !                 Check if the previous state is illegal
                DO s = 1, sd
                    IF (rs(s).LT.0) rs(1) = -1
                ENDDO
                !                 Update the adjacency matrix for the Markov chain
                IF (rs(1).GE.0) THEN
                    CALL state2key(key, rs, sd, MaxNumberMolecules)
                    CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                    IF (found) fsp%matrix%adj(k, fsp%kvtab(ka)) = i
                ENDIF
            ENDDO
        ENDDO
    END SUBROUTINE matrix_starter
    !===================================Extend the FSP state space by 1-step
    SUBROUTINE onestep_extender(fsp, model)
        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        !     The current FSP hash table to be extended
        TYPE(finite_state_projection) :: fsp
        !     The underlying model
        TYPE(cme_model), INTENT(in) :: model
        !--------------------------------------------Parameters of the algorithm
        INTEGER :: sd, pd
        INTEGER :: i, j, k, rs(model%nspecies), state(model%nspecies), lsize, lsize_copy
        LOGICAL :: found
        INTEGER :: mode, ka
        PARAMETER(mode = 1)
        TYPE(big_integer) :: key
        !-----------------------------------Update the FSP hash table and matrix
        sd = model%nspecies
        pd = model%nreactions
        !     Find the current FSP state space size
        lsize_copy = fsp%size
        DO j = 1, lsize_copy
            !           Find the current state
            state = fsp%state(:, j)
            !           For each reaction...
            DO k = 1, pd
                IF (fsp%matrix%adj(k, j).EQ.0) THEN
                    !                       Find the next state according to the reaction
                    rs(:) = state + model%stoichiometry(:, k)

                    CALL key2key(fsp%key(j), key, &
                            k, fsp%reactionkey, fsp%rkeysign, pd)
                    !                       call state2key(key,rs(1:sd,k),sd,MaxNumberMolecules)
                    CALL hash(key, mode, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                    !                       If the next state is already in the FSP...
                    IF (found) THEN
                        !                             Update the adjacency matrix for the Markov
                        !                             chain, otherwise...
                        fsp%matrix%adj(k, j) = fsp%kvtab(ka)
                    ELSE
                        !                             Add the state in the FSP hash table
                        CALL fsp%add(model, rs, keyin = key)
                        !                             Exit if the hash table length is too small
                        IF (fsp%size.GE.fsp%ktlen) THEN
                            STOP 'overflow error: fsp size exceeds memory limit.'
                            RETURN
                        END IF
                    END IF
                END IF
            END DO
        END DO
    END SUBROUTINE onestep_extender

    SUBROUTINE find_droptol(sd, lsize, w, droptol, dsum)
        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        !     Number of species in the model
        INTEGER :: sd
        !     Size of the current FSP state space
        INTEGER :: lsize
        !     The truncated vector approximating the probability vector
        DOUBLE PRECISION :: w(:)
        !     The truncation threshold that satisfies the FSP condition
        DOUBLE PRECISION :: droptol
        !     Maximum allowable value for the sum of states marked as
        !     'insignificant'
        DOUBLE PRECISION :: dsum
        !--------------------------------------------Parameters of the algorithm
        DOUBLE PRECISION :: sum1
        INTEGER :: i, j

        droptol = 1.0d-08
        DO
            sum1 = 0.0d0
            DO i = 1, lsize
                IF (w(i).LT.droptol .AND. w(i).GT.0) THEN
                    sum1 = sum1 + w(i)
                ENDIF
            ENDDO
            IF (sum1.LT.dsum) EXIT
            droptol = droptol / 10.0d0
        ENDDO
    END SUBROUTINE find_droptol


    !============Drop states in the FSP state space with small probabilities
    SUBROUTINE drop_states(w, fsp, model, dsum, fmatvec)
        !     use fsp_manager

        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        !     The truncated vector approximating the probability vector
        DOUBLE PRECISION :: w(:)
        !     The current FSP hash table to be extended
        TYPE(finite_state_projection) :: fsp

        DOUBLE PRECISION :: dsum
        !     The underlying model
        TYPE (cme_model), INTENT(in) :: model
        !--------------------------------------------Parameters of the algorithm
        INTEGER :: sd, pd
        INTEGER, ALLOCATABLE :: Ars(:, :), new_index(:), list(:, :)
        DOUBLE PRECISION, ALLOCATABLE :: Aprops(:, :), Adiag(:), w_copy(:)
        INTEGER :: i, j, k, q, mode, ka, lsize
        TYPE(big_integer), ALLOCATABLE :: keylistcopy(:)

        LOGICAL :: found
        TYPE(big_integer) :: key

        !     Vector containing the indices to be dropped
        LOGICAL, ALLOCATABLE :: drop(:)
        DOUBLE PRECISION, ALLOCATABLE :: wtmp(:)
        DOUBLE PRECISION :: droptol

        INTEGER :: drop_count
        !----------------------------------------------Initialize the parameters
        sd = model%nspecies
        pd = model%nreactions

        lsize = fsp%size

        ALLOCATE(drop(lsize), wtmp(lsize))

        !                       Find the truncation threshold that satisfies the
        !                       probability sum condition
        CALL find_droptol(model%nspecies, fsp%size, w, droptol, dsum)


        !                       Mark the states with probabilities below the
        !                       threshold to be dropped
        drop_count = 0
        DO i = 1, fsp%size
            drop(i) = (w(i).LT.droptol)
            if (w(i) < droptol) then
                drop(i) = .true.
                drop_count = drop_count + 1
            else
                drop(i) = .false.
            end if
        ENDDO

        CALL fmatvec(w, wtmp, fsp%matrix)
        !                       Make sure that the states with big positive
        !                       derivatives, even with small probabilities, are
        !                       not dropped
        DO i = 1, fsp%size
            IF (wtmp(i).GT.1.0d-8) then
                drop(i) = .FALSE.
                drop_count = drop_count - 1
            END IF
        ENDDO

        if (drop_count*1.0d0/(lsize*1.0d0) > 0.1d0) then

            ALLOCATE(Ars(pd, lsize), new_index(lsize), list(sd, lsize), &
                    Aprops(pd, lsize), Adiag(lsize), keylistcopy(lsize), &
                    w_copy(lsize))
            w_copy = 0d0
            q = 0
            !---------------------------------Drop the states in the FSP state space
            !     For each state in the FSP state space...
            DO j = 1, lsize
                !           Check if the state is to be dropped; if not...
                IF (.NOT.drop(j)) THEN
                    !                 Store the state in a new temporary hash table,
                    !                 otherwise...
                    q = q + 1
                    w_copy(q) = w(j)
                    list(1:sd, q) = fsp%state(:, j)
                    Adiag(q) = fsp%matrix%diag(j)
                    Aprops(1:pd, q) = fsp%matrix%offdiag(1:pd, j)
                    Ars(1:pd, q) = fsp%matrix%adj(1:pd, j)
                    keylistcopy(q) = fsp%key(j)
                    new_index(j) = q
                    CALL hash(fsp%key(j), 1, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                    fsp%kvtab(ka) = q
                ELSE
                    !                 Delete the state
                    new_index(j) = 0
                    CALL hash(fsp%key(j), 3, ka, found, fsp%ktlen, fsp%keytab, fsp%kvtab)
                ENDIF
            ENDDO
            !     Clean the probability vector
            w(1:lsize) = 0.0d0
            lsize = q
            !     Copy the temporary hash table to the FSP hash table
            fsp%size = lsize
            fsp%state(1:sd, 1:lsize) = list(1:sd, 1:lsize)
            fsp%key(1:lsize) = keylistcopy(1:lsize)
            w(1:lsize) = w_copy(1:lsize)
            fsp%matrix%size = lsize
            fsp%matrix%adj(1:pd, 1:lsize) = Ars(1:pd, 1:lsize)
            fsp%matrix%diag(1:lsize) = Adiag(1:lsize)
            fsp%matrix%offdiag(1:pd, 1:lsize) = Aprops(1:pd, 1:lsize)
            !     Re-index the adjacency matrix
            DO j = 1, lsize
                DO k = 1, pd
                    i = fsp%matrix%adj(k, j)
                    IF (i>0) fsp%matrix%adj(k, j) = new_index(i)
                ENDDO
            ENDDO
        endif
        ! arrays are automatically deallocated
    END SUBROUTINE drop_states
    !======================================Extend the FSP state space by SSA
    SUBROUTINE SSA_extender(timestep, fsp, model)

        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        DOUBLE PRECISION :: timestep
        TYPE(finite_state_projection) :: fsp
        TYPE (cme_model), INTENT(in) :: model
        !--------------------------------------------Parameters of the algorithm
        INTEGER :: sd, pd, i, j, j0, k, s, rs(model%nspecies), state(model%nspecies)
        DOUBLE PRECISION :: t, r1, r2, r2a
        DOUBLE PRECISION :: tmp
        INTEGER :: lsize_old
        LOGICAL :: found
        INTEGER :: mode, ka
        TYPE(big_integer) :: key
        !--------------------------------------Extend the FSP state space by SSA
        sd = model%nspecies
        pd = model%nreactions
        mode = 1
        lsize_old = fsp%size
        !     From each state in the current FSP state space...
        DO j0 = 1, lsize_old
            j = j0
            state(1:sd) = fsp%state(1:sd, j)
            t = 0.0d0
            !           Generate the next state in the SSA path
            300    CONTINUE
            CALL RANDOM_NUMBER(r1)
            CALL RANDOM_NUMBER(r2)
            t = MIN(timestep, t + (-LOG(r1) / fsp%matrix%diag(j)))
            IF (t<=timestep) THEN
                tmp = fsp%matrix%offdiag(1, j)
                k = 1
                r2a = MIN(r2 * fsp%matrix%diag(j), fsp%matrix%diag(j))
                301       IF (tmp<r2a .AND. k<pd) THEN
                    k = k + 1
                    tmp = tmp + fsp%matrix%offdiag(k, j)
                    go to 301
                ENDIF
                rs = state + model%stoichiometry(:, k)
                DO s = 1, sd
                    IF (rs(s).LT.0) rs(1) = -1
                ENDDO
                !                 If the next state is illegal...
                IF (rs(1).LT.0) THEN
                    !                       Ignore it, otherwise...
                    fsp%matrix%adj(k, j) = -1
                ELSE
                    !                       Store the next state in the FSP state space
                    IF (fsp%matrix%adj(k, j).EQ.0) THEN

                        CALL key2key(fsp%key(j), key, k, &
                                fsp%reactionkey, fsp%rkeysign, pd)
                        CALL hash(key, mode, ka, found, &
                                fsp%ktlen, fsp%keytab, fsp%kvtab)
                        !                             Check if key is already in hashtable, if
                        !                             not...
                        IF (found) THEN
                            j = fsp%kvtab(ka)
                        ELSE
                            !                                   Add the state to the list and update
                            !                                   the FSP hash table
                            IF (fsp%size.GE.fsp%ktlen) THEN
!                                PRINT*, 'overflow error:' &
!                                 ' lsize>=n'
                                RETURN
                            ENDIF
                            CALL fsp%add(model, rs, keyin = key)
                            j = fsp%size
                        ENDIF
                    ELSE
                        i = fsp%matrix%adj(k, j)
                        j = i
                    ENDIF
                    !                       Move to the next step on the SSA path
                    state(1:sd) = fsp%state(1:sd, j)
                    IF ((t.LT.timestep).AND.(j.GE.j0)) go to 300
                ENDIF
            ENDIF
        ENDDO
    END SUBROUTINE SSA_extender


    !=============================Computing the absolute value of key change
    !==========================================due to the chemical reactions
    SUBROUTINE compute_rkey(reactionkey, rkeysign, sd, pd, model)
        USE big_integer_module
        IMPLICIT NONE
        !------------------------------------------------------Output parameters
        !     Vector storing the reaction keys, such that if key(x) is the key
        !     of state x then key(x+nu_k)=key(x)+rkeysign(k)*reactionkey(k)
        INTEGER :: sd, pd
        TYPE(big_integer) :: reactionkey(pd)
        INTEGER :: rkeysign(pd)
        TYPE(cme_model), INTENT(in) :: model
        !--------------------------------------------Parameters of the algorithm
        INTEGER :: i, j, state(1:sd), rs(1:sd), sgn
        TYPE(big_integer) :: rkey, C
        !-----------------------------Computing the absolute value of key change
        C = big(MaxNumberMolecules)
        state(1:sd) = 0
        !     From the zero state, for each reaction...
        DO j = 1, pd
            !           Find the next state according to the reaction
            rs = state + model%stoichiometry(:, j)
            sgn = 1
            rkey = big(0)
            !           Find the reaction key
            DO i = 1, sd
                IF (sgn * rs(i).LT.0) THEN
                    sgn = -1 * sgn
                    rkey = ABS(rs(i)) * (C + 1)**(i - 1) - rkey
                ELSE
                    rkey = ABS(rs(i)) * (C + 1)**(i - 1) + rkey
                ENDIF
            ENDDO
            reactionkey(j) = rkey
            rkeysign(j) = sgn
        ENDDO
    END SUBROUTINE compute_rkey

END MODULE StateSpace
