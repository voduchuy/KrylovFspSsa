! Created by ${USER_NAME} on 11/27/18.

PROGRAM solve_repress
  !     Solve CME of the repressilator
  USE StateSpace
  USE KrylovSolver
  IMPLICIT NONE

  INTERFACE
     DOUBLE PRECISION FUNCTION clock()
     END FUNCTION clock
  END INTERFACE

  DOUBLE PRECISION :: t = 10.0d0, exp_tol = 1.0d-14, fsp_tol = 1.0d-4
  TYPE (cme_model) :: model
  TYPE (finite_state_projection) :: fsp_in, fsp
  DOUBLE PRECISION, DIMENSION(nmax) :: p0, p
  DOUBLE PRECISION :: tic
  INTEGER :: i

  CALL RANDOM_SEED()

  CALL model%init_mem(3, 6, 3)
  model%stoichiometry = RESHAPE((/1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1/), (/3,6/))
  model%customprop => propensity
  CALL model%reset_parameters([100.0D0, 25.0D0, 1.0D0])
  model%loaded = .TRUE.

  ! initialize the attributes of the fsp based on the underlying model
  tic = clock()
  write(*,*) 'Allocating memory for the FSP...'
  CALL fsp_in%init(model)
  CALL fsp%init(model)
  WRITE(*, fmt = 200) clock() - tic

  fsp_in%size = 1
  fsp_in%state(:, 1) = [22, 0, 0]
  p0(1) = 1.0
  fsp_in%vector = p0
  fsp = fsp_in
  tic = clock()
  write(*,*) 'Calling the solver...'
  CALL cme_solve(model, t, fsp_in, fsp, fsp_tol, exp_tol, verbosity = 1)
  WRITE(*, fmt = 200) clock() - tic
200 FORMAT ('Elapsed time :', F10.2)

  CALL fsp%clear()
  CALL fsp_in%clear()
CONTAINS
  DOUBLE PRECISION FUNCTION propensity(state, reaction, parameters)
    INTEGER, INTENT(in) :: state(:), reaction
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: parameters(:)

    SELECT CASE (reaction)
    CASE(1)
       propensity = parameters(1)/(1d0 + parameters(2)*state(2)**6.0d0)
    CASE(2)
       propensity = parameters(3)*state(1)
    CASE(3)
      propensity = parameters(1)/(1d0 + parameters(2)*state(3)**6.0d0)
    CASE(4)
      propensity = parameters(3)*state(2)
    CASE(5)
      propensity = parameters(1)/(1d0 + parameters(2)*state(1)**6.0d0)
    CASE(6)
      propensity = parameters(3)*state(3)

    END SELECT
  END FUNCTION propensity
END PROGRAM
