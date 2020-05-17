MODULE ModelModule
  ! This module will define the 'class' of CME models in Fortran
  USE FortranParser
  IMPLICIT NONE

  INTERFACE
     DOUBLE PRECISION FUNCTION propfunc(state, reaction, parameters)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: state(:), reaction
       DOUBLE PRECISION, INTENT(in), OPTIONAL :: parameters(:)
     END FUNCTION propfunc
  END INTERFACE

  TYPE :: cme_model
     ! Model needs to be loaded before we do anything else
     LOGICAL :: loaded = .FALSE.
     ! Numbers of species, reactions, and uncertain parameters
     INTEGER :: nSpecies
     INTEGER :: nReactions
     INTEGER :: nParameters

     DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: parameter_val

     ! Stoichiometry matrix, stored in dense format, each column corresponds to one reaction
     INTEGER, DIMENSION (:, :), ALLOCATABLE :: stoichiometry

     CHARACTER (len = 20), ALLOCATABLE, DIMENSION (:) :: species_names, parameter_names
     ! Private attribute, the EquationParser type stor\es the byte code to
     ! compute propensity functions
     TYPE (EquationParser), DIMENSION(:), ALLOCATABLE, PRIVATE :: propParser
     PROCEDURE (propfunc), POINTER, nopass :: customprop

   CONTAINS
     PROCEDURE :: init_mem

     PROCEDURE :: load

     PROCEDURE :: reset_parameters

     PROCEDURE :: propensity => propensity_builtin

  END TYPE cme_model

CONTAINS

  SUBROUTINE init_mem(this, n_species, n_reactions, n_parameters)
    IMPLICIT NONE
    class(cme_model) :: this
    INTEGER :: n_species
    INTEGER :: n_reactions
    INTEGER :: n_parameters

    this%nSpecies = n_species
    this%nReactions = n_reactions
    this%nParameters = n_parameters
    ALLOCATE(this%stoichiometry(n_species, n_reactions), this%parameter_val(n_parameters))
  END SUBROUTINE init_mem

  SUBROUTINE load(this, filename)
    ! Purpose:
    ! ========
    ! Load the model from file
    !
    ! Usage:
    ! ======
    !        call model%load(<filename>)
    !        the variable filename is optional
    IMPLICIT NONE

    CLASS (cme_model) :: this
    CHARACTER (len = *), OPTIONAL :: filename

    CHARACTER (len = 400) :: ichar
    CHARACTER (len = 200) :: fstring
    CHARACTER (len = 20), ALLOCATABLE, DIMENSION (:) :: fvar

    LOGICAL :: nspecies_read = .FALSE., nreactions_read = .FALSE., nparam_read = .FALSE., &
         parameters_read = .FALSE., species_names_read = .FALSE., &
         param_names_read = .FALSE.

    INTEGER :: i, ios

    IF (PRESENT(filename)) THEN
       OPEN(unit = 10, file = filename, iostat = ios, status = "old", action = "read")
       IF (ios /= 0) STOP "Error opening file "
    ELSE
       OPEN(unit = 10, file = 'model.input', iostat = ios, status = "old", action = "read")
       IF (ios /= 0) STOP "Error opening file "
    ENDIF

    file_scan : DO
       READ(10, fmt = *, iostat = ios) ichar
       IF (ios /= 0) EXIT

       IF (ichar == 'nspecies') THEN
          READ(10, fmt = *, iostat = ios) this%nspecies
          nspecies_read = .TRUE.
       ENDIF

       IF (ichar == 'nreactions') THEN
          READ(10, fmt = *, iostat = ios) this%nreactions
          nreactions_read = .TRUE.
       ENDIF

       IF (ichar == 'nparameters') THEN
          READ(10, fmt = *, iostat = ios) this%nparameters
          nparam_read = .TRUE.
       ENDIF

       IF (ichar == 'species') THEN
          ALLOCATE(this%species_names(this%nspecies))
          DO i = 1, this%nspecies
             READ(10, fmt = *, iostat = ios) this%species_names(i)
          ENDDO
          species_names_read = .TRUE.
       ENDIF

       !      Read the reactions
       IF (ichar == 'reactions') THEN
          IF (.NOT.species_names_read) STOP 'Model input error: reactions stated before species names are declared.'
          IF (.NOT.nspecies_read) STOP 'Model input error: number of species not declared.'
          IF (.NOT.nreactions_read) STOP 'Model input error: number of reactions not declared.'
          IF (.NOT.ALLOCATED(this%stoichiometry)) ALLOCATE(this%stoichiometry(this%nspecies, this%nreactions))
          DO i = 1, this%nreactions
             READ(10, fmt = '(A)', iostat = ios) ichar
             CALL stoich_input(this%nspecies, this%stoichiometry(1:this%nspecies, i), ichar, this%species_names)
          ENDDO
       ENDIF

       IF (ichar == 'parameters') THEN
          IF (.NOT.nparam_read) STOP 'Model input error: number of parameters not declared before specifying parameter names.'
          ALLOCATE(this%parameter_names(this%nparameters), this%parameter_val(this%nparameters))
          DO i = 1, this%nparameters
             READ(10, fmt = *, iostat = ios) this%parameter_names(i)
          ENDDO
          param_names_read = .TRUE.
       ENDIF

       !     Read the propensity functions, the strings of mathematical expressions have to be listed in the same order as the reactions
       IF (ichar == 'propensities') THEN
          IF ((species_names_read.EQV..FALSE.).OR.(param_names_read.EQV..FALSE.)) THEN
             STOP 'Model input error: Propensities specified before all species and parameters are named.'
          ELSE
             ALLOCATE(fvar(this%nspecies + this%nparameters))
             ALLOCATE(this%propparser(this%nreactions))
             DO i = 1, this%nspecies
                fvar (i) = TRIM(this%species_names(i))
             ENDDO
             DO i = 1, this%nparameters
                fvar (this%nspecies + i) = TRIM(this%parameter_names(i))
             ENDDO
             DO i = 1, this%nreactions
                READ(10, fmt = '(A)', iostat = ios) fstring
                this%propparser(i) = EquationParser(TRIM(fstring), fvar)
             ENDDO
          ENDIF
       ENDIF
    ENDDO file_scan
    CLOSE (10)
    this%loaded = .TRUE.
  END SUBROUTINE load
  !-----------------------------------------------------------------------------|
  DOUBLE PRECISION FUNCTION propensity_builtin(this, state, reaction)
    ! Purpose:
    ! ========
    ! Interface to compute the propensity based on the string expression
    ! in the model input file.
    !
    ! How to use:
    ! ==========
    ! y = model%propensity(state, reaction)
    !
    ! Arguments:
    ! =========
    ! state         : array storing the state
    ! reaction      : integer storing the reaction index
    !

    IMPLICIT NONE

    class (cme_model), INTENT(in) :: this
    INTEGER, INTENT(in) :: state(:)
    INTEGER, INTENT(in) :: reaction

    DOUBLE PRECISION, DIMENSION(this%nspecies + this%nparameters) :: fval
    INTEGER i, N

    IF (ASSOCIATED(this%customprop)) THEN
       propensity_builtin = this%customprop(state(1:this%nspecies), reaction, this%parameter_val(1:this%nParameters))
    ELSE
       N = this%nspecies
       DO i = 1, this%nspecies
          fval(i) = DBLE(state(i))
       ENDDO
       fval(this%nspecies + 1:this%nspecies + this%nparameters) = this%parameter_val(1:this%nparameters)
       propensity_builtin = this%propparser(reaction)%evaluate(fval)
    ENDIF

  END FUNCTION propensity_builtin
  !-----------------------------------------------------------------------------|
  SUBROUTINE reset_parameters(this, pval)
    ! Purpose:
    ! ========
    ! Reset the values of the uncertain parameters in the stochastic chemical kinetics model.
    !
    ! How to use:
    ! ===========
    ! the command
    !           call model%reset_parameters(X)
    ! will reset the parameter values in the model to the values in the array X
    IMPLICIT NONE
    class (cme_model) :: this
    DOUBLE PRECISION, DIMENSION (:) :: pval

    this%parameter_val(1:this%nparameters) = pval(1:this%nparameters)

  END SUBROUTINE reset_parameters
  !-----------------------------------------------------------------------------|
  SUBROUTINE stoich_input(nspecies, stoich_vec, str, species_names)

    !   Purpose:
    !   =======
    !   Output the stoichiometry vector of the reaction given its string expression
    !
    IMPLICIT NONE
    INTEGER, INTENT (inout) :: stoich_vec(*)
    INTEGER, INTENT(in) :: nspecies
    CHARACTER(len = *), INTENT(in) :: str
    CHARACTER(len = *), DIMENSION(nspecies), INTENT(in) :: species_names


    !   Internal variables
    CHARACTER(len = 200), DIMENSION(100) :: terms       ! supporting up to 100 reactants+products
    CHARACTER(len = 50) :: word = ''
    INTEGER :: dir = 0       ! direction of the arrow , 1 if "->", 2 if "<-", if 0 then unspecified.
    INTEGER :: ell, i, nterms, nleft, coeff, k, j
    LOGICAL :: defined = .FALSE.

    nterms = 0
    nleft = 0

    word = ''

    stoich_vec(1:nspecies) = 0
    ell = LEN(TRIM(str))
    string_scan : DO i = 1, ell
       word = TRIM(word) // str(i:i)
       IF (str(i:i) == ' ' .OR. i == ell) THEN
          word = ADJUSTL(word);
          IF (TRIM(word) == '->') THEN
             dir = 1
             nleft = nterms ! number of terms on the left side of the equation
          ELSEIF (TRIM(word) == '<-') THEN
             dir = 2
             nleft = nterms ! number of terms on the left side of the equation
          ELSEIF (TRIM(word) /= '+') THEN
             nterms = nterms + 1
             terms(nterms) = TRIM(word)
          ENDIF
          word = ''
       ENDIF
    END DO string_scan

    IF (dir==0) THEN
       STOP 'Syntax error in chemical reaction, only one side was written.'
    ENDIF

    ! Processing the terms of the reaction equation, read their coefficients to the stoichiometry vector
    DO i = 1, nterms
       WRITE(*, *) terms(i)
       IF (terms(i) /= '0') THEN
          DO j = 1, nspecies
             word = terms(i)
             k = INDEX(terms(i), TRIM(species_names(j)))
             defined = defined .OR. (k/=0)
             IF (k==0) THEN
                coeff = 0
             ELSEIF (word(k:LEN(TRIM(word)))==TRIM(species_names(j))) THEN
                IF (k > 1) THEN
                   READ(word(1:k - 1), *) coeff
                ELSEIF (k == 1) THEN
                   coeff = 1
                ENDIF
             ENDIF
             IF (i <= nleft) THEN
                stoich_vec(j) = stoich_vec(j) - coeff
             ELSE
                stoich_vec(j) = stoich_vec(j) + coeff
             ENDIF
          ENDDO
          IF (.NOT.defined) WRITE(*, *) 'Warning: species ', terms(i), ' not defined in the model.'
       ENDIF
    END DO

    IF (dir == 2) stoich_vec(1:nspecies) = -1 * stoich_vec(1:nspecies)

  END SUBROUTINE stoich_input

END MODULE ModelModule
