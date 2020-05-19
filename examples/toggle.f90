! Created by ${USER_NAME} on 11/27/18.

PROGRAM solve_toggle
  !     Solve CME of the toggle switch
  USE StateSpace
  USE KrylovSolver
  IMPLICIT NONE

  INTERFACE
     DOUBLE PRECISION FUNCTION clock()
     END FUNCTION clock
  END INTERFACE

  DOUBLE PRECISION :: t = 100.0d0, exp_tol = 1.0d-8, fsp_tol = 1.0d-4
  TYPE (cme_model) :: model
  TYPE (finite_state_projection) :: fsp_in, fsp
  DOUBLE PRECISION, DIMENSION(nmax) :: p0, p
  DOUBLE PRECISION :: tic
  INTEGER :: i

  CALL RANDOM_SEED()

  CALL model%CREATE(2, 4, 6)
  model%stoichiometry = RESHAPE((/1,0,-1,0,0,1,0,-1/), (/2,4/))
  model%customprop => toggle_propensity
  CALL model%reset_parameters([1.0d0, 100.0d0, 1.0d0, 1.0d0, 100.0d0, 1.0d0])
  model%loaded = .TRUE.

  WRITE(*,*) model%stoichiometry
  write(*,*) model%parameter_val
  write(*,*) toggle_propensity([0,0], 1, [1.0d0, 100.0d0, 1.0d0, 1.0d0, 100.0d0, 1.0d0])
  WRITE(*,*) model%propensity([0,0], 1)

  ! initialize the attributes of the fsp based on the underlying model
  tic = clock()
  write(*,*) 'Allocating memory for the FSP...'
  CALL fsp_in%CREATE(model)
  CALL fsp%CREATE(model)
  WRITE(*, fmt = 200) clock() - tic

  fsp_in%size = 1
  fsp_in%state(:, 1) = [0, 0]
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
  DOUBLE PRECISION FUNCTION toggle_propensity(state, reaction, parameters)
    INTEGER, INTENT(in) :: state(:), reaction
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: parameters(:)

    SELECT CASE (reaction)
    CASE(1)
       toggle_propensity = parameters(1) + parameters(2)/(1d0 + state(2)**1.5d0)
    CASE(2)
       toggle_propensity = parameters(3)*state(1)
    CASE(3)
       toggle_propensity = parameters(4) + parameters(5)/(1d0 + state(1)**3.5d0)
    CASE(4)
       toggle_propensity = parameters(6)*state(2)
    END SELECT
  END FUNCTION toggle_propensity
END PROGRAM solve_toggle
