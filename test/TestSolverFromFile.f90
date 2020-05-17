PROGRAM TestSolverFromFile
    !     Solve CME of the toggle switch
    USE StateSpace
    USE KrylovSolver
    IMPLICIT NONE

    INTERFACE
        DOUBLE PRECISION FUNCTION clock()
        END FUNCTION clock
    END INTERFACE

    DOUBLE PRECISION :: t = 1000.0d0, exp_tol = 1.0d-8, fsp_tol = 1.0d-6
    TYPE (cme_model) :: model
    TYPE (finite_state_projection) :: fsp_in, fsp
    DOUBLE PRECISION, DIMENSION(nmax) :: p0, p
    DOUBLE PRECISION :: tic
    INTEGER :: i

    CALL RANDOM_SEED()

    CALL model%load('toggle_model.input')

    ! initialize the attributes of the fsp based on the underlying model
    CALL fsp_in%init(model)
    CALL fsp%init(model)

    fsp_in%size = 1
    fsp_in%state(:, 1) = [0, 0]
    p0(1) = 1.0
    fsp_in%vector = p0
    CALL model%reset_parameters([1.0d0, 100.0d0, 1.0d0, 1.0d0, 100.0d0, 1.0d0])
    fsp = fsp_in
    tic = clock()

    CALL cme_solve(model, t, fsp_in, fsp, 1.0d-4, 1.0d-10, verbosity = 1)

    WRITE(*, fmt = 200) clock() - tic
200 FORMAT ('Elapsed time=', F10.2)

    CALL fsp%clear()
    CALL fsp_in%clear()
END PROGRAM TestSolverFromFile
