! Created by ${USER_NAME} on 11/27/18.

PROGRAM solve_goutsias
    !     Solve CME of the Goutsias example
    USE StateSpace
    USE KrylovSolver
    IMPLICIT NONE

    INTERFACE
        DOUBLE PRECISION FUNCTION clock()
        END FUNCTION clock
    END INTERFACE

    INTEGER::  M = 1, D = 2, RNA = 3, DNA = 4, DNAD = 5, DNA2D = 6

    DOUBLE PRECISION :: t = 300.0d0, exp_tol = 1.0d-8, fsp_tol = 1.0d-6
    TYPE (cme_model) :: model
    TYPE (finite_state_projection) :: fsp_in, fsp
    DOUBLE PRECISION, DIMENSION(nmax) :: p0, p
    DOUBLE PRECISION :: tic
    INTEGER :: i

    double precision, dimension(10) :: true_parameters = (/0.043d0, &
            0.0007d0, &
            0.0715d0, &
            0.0039d0, &
            0.0199264663575241d0, &
            0.4791d0, &
            0.000199264663575241d0, &
            0.8765d0 * 1.0D-11, &
            0.0830269431563506104d0, &
            0.5d0 /)

    CALL RANDOM_SEED()

    CALL model%init_mem(6, 10, 10)
    CALL goutsias_stoich(model%stoichiometry)
    model%customprop => goutsias_propensity
    CALL model%reset_parameters(true_parameters)
    model%loaded = .TRUE.

    ! initialize the attributes of the fsp based on the underlying model
    tic = clock()
    write(*, *) 'Allocating memory for the FSP...'
    CALL fsp_in%init(model)
    CALL fsp%init(model)
    WRITE(*, fmt = 200) clock() - tic

    fsp_in%size = 1
    fsp_in%state(:, 1) = [2, 6, 0, 2, 0, 0]
    p0(1) = 1.0
    fsp_in%vector = p0
    fsp = fsp_in
    tic = clock()
    write(*, *) 'Calling the solver...'
    CALL cme_solve(model, t, fsp_in, fsp, fsp_tol, exp_tol, verbosity = 1)
    WRITE(*, fmt = 200) clock() - tic
    200 FORMAT ('Elapsed time :', F10.2)

    CALL fsp%clear()
    CALL fsp_in%clear()
CONTAINS
    DOUBLE PRECISION FUNCTION goutsias_propensity(state, reaction, parameters)
        INTEGER, INTENT(in) :: state(:), reaction
        DOUBLE PRECISION, INTENT(in), OPTIONAL :: parameters(:)

        select case (reaction)
        case (1)
            goutsias_propensity = parameters(1) * state(RNA)
        case (2)
            goutsias_propensity = parameters(2) * state(M)
        case (3)
            goutsias_propensity = parameters(3) * state(DNAD)
        case (4)
            goutsias_propensity = parameters(4) * state(RNA)
        case (5)
            goutsias_propensity = parameters(5) * state(DNA) * state(D)
        case (6)
            goutsias_propensity = parameters(6) * state(DNAD)
        case (7)
            goutsias_propensity = parameters(7) * state(DNAD) * state(D)
        case (8)
            goutsias_propensity = parameters(8) * state(DNA2D)
        case (9)
            goutsias_propensity = parameters(9) * (state(M) * (state(M) - 1) / 2)
        case (10)
            goutsias_propensity = parameters(10) * state(D)
        end select
    END FUNCTION goutsias_propensity

    subroutine goutsias_stoich(stoich)
        implicit none

        integer, dimension(:, :), intent(inout) :: stoich
        integer :: i, j
        do i = 1,6
            do j = 1,10
                stoich(i,j) = 0
            end do
        end do

        ! The Void -----> M
        stoich(M, 1) = 1
        !     M-------> The Void
        stoich(M, 2) = -1
        !    DNA.D --> RNA+DNA.D
        stoich(RNA, 3) = 1
        !    RNA----> The Void
        stoich(RNA, 4) = -1
        !     DNA+D----->DNA.D
        stoich(DNA, 5) = -1
        stoich(D, 5) = -1
        stoich(DNAD, 5) = 1
        !    DNA.D----> DNA+D
        stoich(DNA, 6) = 1
        stoich(D, 6) = 1
        stoich(DNAD, 6) = -1
        !    DNA.D+D--->DNA.2D
        stoich(DNAD, 7) = -1
        stoich(D, 7) = -1
        stoich(DNA2D, 7) = 1
        !     DNA.2D---->DNA.D+D
        stoich(DNAD, 8) = 1
        stoich(D, 8) = 1
        stoich(DNA2D, 8) = -1
        !     M+M----->D
        stoich(M, 9) = -2
        stoich(D, 9) = 1
        !     D-----> M+M
        stoich(M, 10) = 2
        stoich(D, 10) = -1
    end subroutine goutsias_stoich

END PROGRAM solve_goutsias
