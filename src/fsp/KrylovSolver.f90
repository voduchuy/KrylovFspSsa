!     A solver for the Chemical Master equation
!     by using the Krylov-FSP-SSA algorithm
MODULE KrylovSolver
    USE StateSpace
CONTAINS
    !=====================Interface for Krylov-FSP-SSA algorithm
    SUBROUTINE cme_solve(model, t, fsp_in, fsp_out, FSPTol, exp_tol, verbosity)
        IMPLICIT NONE
        !------------------------------------------------Input/output parameters
        !     The underlying model
        TYPE (cme_model), INTENT(in) :: model
        !     Time elapsed for the CME solution
        DOUBLE PRECISION, INTENT(in) :: t
        !     hash table of the FSP
        TYPE(finite_state_projection) :: fsp_in
        !     Tolerance of the Krylov-FSP-SSA algorithm
        DOUBLE PRECISION, INTENT(in) :: FSPTol
        !     Tolerance of the local exponential approximation
        DOUBLE PRECISION, INTENT(in) :: exp_tol
        !     hash table of the FSP
        TYPE(finite_state_projection) :: fsp_out
        !-----------------------------------------verbosity and error flag
        INTEGER, INTENT(in), OPTIONAL :: verbosity
        INTEGER :: itrace, iflag
        !-------------------------------------Apply the Krylov-FSP-SSA algorithm
        IF (PRESENT(verbosity)) THEN
            itrace = verbosity
        ELSE
            itrace = 0
        ENDIF

        PRINT*, 'calling dgexpv_fsp'
        CALL dgexpv_fsp(model, t, fsp_in%vector, fsp_out, fsp_out%vector, &
                FSPTol, exp_tol, itrace, iflag)

    END SUBROUTINE cme_solve


    !===============================================Krylov-FSP-SSA algorithm
    SUBROUTINE dgexpv_fsp(model, t, v, fsp, w, FSPTol, KryTol, itrace, iflag)
        IMPLICIT NONE

        INTEGER :: nfull
        PARAMETER(nfull = nmax)

        INTEGER :: m_max, m_min
        PARAMETER(m_max = 100, m_min = 10)
        !     Workspace
        INTEGER :: lwsp, liwsp
        PARAMETER(lwsp = nfull * (m_max + 2) + 5 * (m_max + 2)**2 + 7, liwsp = nfull)
        DOUBLE PRECISION :: wsp(lwsp)
        INTEGER :: iwsp(liwsp)

        !     The underlying model
        TYPE(cme_model) :: model
        !     Time elapsed for the CME solution
        DOUBLE PRECISION, INTENT(in) :: t
        !     Input probability vector
        DOUBLE PRECISION, INTENT(in) :: v(:)
        !     Output probability vector
        DOUBLE PRECISION, INTENT(out) :: w(:)
        !     Tolerance of the Krylov-FSP-SSA algorithm
        DOUBLE PRECISION, INTENT(in) :: FSPTol
        !     Tolerance of the local Krylov error
        DOUBLE PRECISION :: KryTol

        TYPE(finite_state_projection) :: fsp
        !
        INTEGER :: itrace
        !
        INTEGER :: iflag


        !-------------------------Parameters for the stepsize control of EXPOKIT
        !     maximum allowable number of integration steps;
        !     0 means an infinite number of steps
        INTEGER, PARAMETER :: mxstep = 0

        INTEGER, PARAMETER :: mxreject = 0         !     maximum allowable number of rejections at each step;
        !     0 means an infinite number of rejections

        INTEGER, PARAMETER :: ideg = 6         !     The Pade approximation of type (ideg,ideg) is used as an
        !     approximation to exp(H); the value 0 switches to the uniform
        !     rational Chebyshev approximation of type (14,14)
        DOUBLE PRECISION, PARAMETER :: delta = 1.2d0 !     Local truncation error `safety factor'
        !     Stepsize `shrinking factor'
        DOUBLE PRECISION, PARAMETER :: gamma = 0.9d0


        !--------------------------------------------Parameters of the algorithm
        INTEGER :: n, m
        DOUBLE PRECISION :: anorm, wtmp(nfull)
        DOUBLE PRECISION :: ent_old = 0.0d0, ent = 0.0d0
        INTEGER :: i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph
        INTEGER :: ireject, ibrkflag, mbrkdwn, nmult, nreject, nexph
        INTEGER :: nscale, nstep
        DOUBLE PRECISION :: sgn, t_out, tbrkdwn, step_min, step_max, err_loc
        DOUBLE PRECISION :: s_error, x_error, t_now, t_new, t_step, t_old
        DOUBLE PRECISION :: xm, beta, break_tol, p1, p2, p3, eps, rndoff
        DOUBLE PRECISION :: vnorm, avnorm, hj1j, hij, hump, sqr1
        INTRINSIC :: aint, abs, dble, log10, max, min, nint, sign
        DOUBLE PRECISION :: ddot, dnrm2, dasum

        DOUBLE PRECISION :: htmp((m_max + 2) * (m_max + 2))

        ! For the Krylov dimension adaptivity,
        DOUBLE PRECISION :: k_factor, order
        LOGICAL :: m_changed, orderold, kestold
        INTEGER :: m_new, m_old, imreject, jold, ih_old
        REAL :: cost1, cost2
        DOUBLE PRECISION :: omega, omega_old, t_opt, m_opt
        DOUBLE PRECISION :: hnorm, nom
        INTEGER :: nnz, m_start
        INTEGER :: qiop

        ! For the FSP adaptivity
        INTEGER :: n_now, iexpand
        DOUBLE PRECISION :: wsum, t_ratio
        LOGICAL :: drop(nfull)
        DOUBLE PRECISION :: error, errorold, tau_old, tfsp
        DOUBLE PRECISION :: fsporder
        INTEGER :: irejectfsp, istart
        DOUBLE PRECISION :: t_ssa
        DOUBLE PRECISION :: wsum_old
        DOUBLE PRECISION :: droptol, dsum


        !------------------------------------Check and initialize the parameters
        !     Assume ||A||=1, usually it is unknown
        anorm = 1.0d0
        CALL matrix_starter(fsp, model)

        DO i = 1, 5
            CALL onestep_extender(fsp, model)
        ENDDO
        !     Initialize the Krylov dimension and orthogonalization length
        m = m_min
        qiop = 2
        istart = 1
        n = nfull

        !     Check the restrictions on the input parameters
        iflag = 0
        IF (lwsp.LT.(n * (m + 2) + 5 * (m + 2)**2 + ideg + 1))  iflag = -1
        IF (liwsp.LT.(m + 2))                       iflag = -2
        IF ((m.GE.n).OR.(m.LE.0))                 iflag = -3
        IF (iflag.NE.0) THEN
            PRINT*, 'Bad sizes (in input of DgeXPV), iflag = ', iflag
            STOP
        END IF

        !     Initialize the parameters
        ibrkflag = 0
        nmult = 0
        nreject = 0
        nexph = 0
        nscale = 0
        t_out = ABS(t)
        tbrkdwn = 0.0d0
        step_min = t_out
        step_max = 0.0d0
        nstep = 0
        s_error = 0.0d0
        x_error = 0.0d0
        t_now = 0.0d0
        t_new = 0.0d0
        p1 = 4.0d0 / 3.0d0
        1   p2 = p1 - 1.0d0
        p3 = p2 + p2 + p2
        eps = ABS(p3 - 1.0d0)
        IF (eps.EQ.0.0d0) go to 1
        IF (KryTol.LE.eps)   KryTol = SQRT(eps)
        rndoff = eps * anorm
        break_tol = 1.0d-07
        sgn = SIGN(1.0d0, t)

        CALL dcopy(fsp%size, v, 1, w, 1)
        beta = dnrm2(fsp%size, w, 1)
        vnorm = beta
        hump = beta

        !     Obtain the very first stepsize
        sqr1 = SQRT(0.1d0)
        xm = 1.0d0 / DBLE(m)
        p1 = KryTol * (((m + 1) / 2.72D0)**(m + 1)) * SQRT(2.0D0 * 3.14D0 * (m + 1))
        t_new = (1.0d0 / anorm) * (p1 / (4.0d0 * beta * anorm))**xm
        p1 = 10.0d0**(NINT(LOG10(t_new) - sqr1) - 1)
        t_new = AINT(t_new / p1 + 0.55d0) * p1

        !     Initialize the FSP parameters
        n_now = fsp%size
        iexpand = 0
        wsum_old = 1.0d0
        irejectfsp = 0
        droptol = 1.0d-8

        nnz = (model%nReactions + 1) * fsp%size

        !     Parameters for the Krylov subspace adaptivity,
        !     based on Niesen and Wright
        imreject = 0
        jold = 1
        m_new = m
        orderold = .TRUE.
        kestold = .TRUE.

        100 IF (t_now >= t_out) GOTO 500
        !     Find the real stepsize
        t_step = MIN(t_out - t_now, t_new)
        !     Compute the Arnoldi matrices H_m and V_m
        n = n_now
        m = MIN(n - 1, m_new)
        mbrkdwn = m
        !     Find pointers into the workspace
        k1 = 2
        mh = m + 2
        iv = 1
        ih = iv + n * (m + 1) + n
        ifree = ih + mh * mh
        lfree = lwsp - ifree + 1

        nstep = nstep + 1

        p1 = 1.0d0 / beta
        DO i = 1, n
            wsp(iv + i - 1) = p1 * w(i)
        ENDDO
        DO i = 1, mh * mh
            wsp(ih + i - 1) = 0.0d0
        ENDDO

        ireject = 0
        !     Start the stopwatch
        IF (itrace.NE.0) THEN
            PRINT*, 'beginning IOP...'
        ENDIF
        101 CONTINUE
        j1v = iv + jold * n
        DO 200 j = jold, m
            nmult = nmult + 1
            CALL fmatvec(wsp(j1v - n), wsp(j1v), fsp%matrix)
            IF (qiop.GT.0)    istart = MAX(1, j - qiop + 1)
            DO i = istart, j
                hij = ddot(n, wsp(iv + (i - 1) * n), 1, wsp(j1v), 1)
                CALL daxpy(n, -hij, wsp(iv + (i - 1) * n), 1, wsp(j1v), 1)
                wsp(ih + (j - 1) * mh + i - 1) = hij
            ENDDO
            hj1j = dnrm2(n, wsp(j1v), 1)
            !           If a `happy breakdown' occurs, go straightforward at the end
            IF (hj1j.LE.break_tol) THEN
                k1 = 0
                ibrkflag = 1
                mbrkdwn = j
                tbrkdwn = t_now
                t_step = t_out - t_now
                go to 300
            END IF
            wsp(ih + (j - 1) * mh + j) = hj1j
            CALL DSCAL(n, 1.0d0 / hj1j, wsp(j1v), 1)
            j1v = j1v + n
        200    CONTINUE
        nmult = nmult + 1
        CALL fmatvec(wsp(j1v - n), wsp(j1v), fsp%matrix)
        avnorm = DNRM2(n, wsp(j1v), 1)
        !     Set 1 for the 2-corrected scheme
        300    CONTINUE
        wsp(ih + m * mh + m + 1) = 1.0d0
        401    CONTINUE

        !     Compute w=beta*V*exp(t_step*H)*e1...
        nexph = nexph + 1
        mx = mbrkdwn + k1
        IF (ideg.NE.0) THEN
            !           by the irreducible rational Pade approximation...
            CALL dgpadmnorm(ideg, mx, sgn * t_step, wsp(ih), mh, &
                    wsp(ifree), lfree, iwsp, iexph, ns, iflag, hnorm)
            iexph = ifree + iexph - 1
            nscale = nscale + ns
        ELSE
            !           or the uniform rational Chebyshev approximation
            iexph = ifree
            DO i = 1, mx
                wsp(iexph + i - 1) = 0.0d0
            ENDDO
            wsp(iexph) = 1.0d0
            CALL dgchbv(mx, sgn * t_step, wsp(ih), mh, &
                    wsp(iexph), wsp(ifree + mx))
        ENDIF
        402    CONTINUE
        !     Estimate the error
        IF (k1.EQ.0) THEN
            err_loc = KryTol
        ELSE
            p1 = ABS(wsp(iexph + m)) * beta
            p2 = ABS(wsp(iexph + m + 1)) * beta * avnorm
            IF (p1.GT.10.0d0 * p2) THEN
                err_loc = p2
                xm = 1.0d0 / DBLE(m)
            ELSE IF (p1.GT.p2) THEN
                err_loc = (p1 * p2) / (p1 - p2)
                xm = 1.0d0 / DBLE(m)
            ELSE
                err_loc = p1
                xm = 1.0d0 / DBLE(m - 1)
            END IF
        ENDIF
        !     Reduce the stepsize if the error is out of bound
        IF (isnan(err_loc)) THEN
            t_step = t_step / 5.0d0
            go to 401
        ENDIF

        !     Find the ratio of error per unit step
        omega_old = omega
        omega = err_loc / (KryTol * t_step)
        !     Estimate the order
        IF ((m.EQ.m_old).AND.(t_step.NE.t_old).AND.(ireject.GE.1)) THEN
            order = MAX(1.0d0, LOG(omega / omega_old) / LOG(t_step / t_old))
            orderold = .FALSE.
        ELSE IF (orderold.OR.ireject.EQ.0) THEN
            order = DBLE(m) / 4.0d0
            orderold = .TRUE.
        ELSE
            orderold = .TRUE.
        ENDIF
        !     Estimate kappa
        IF ((m.NE.m_old).AND.(t_step.EQ.t_old).AND.(ireject.GE.1)) THEN
            k_factor = MAX(1.1d0, (omega / omega_old)**(1.0d0 / (m_old - m)))
            kestold = .FALSE.
        ELSE IF (kestold .OR. ireject.EQ.0) THEN
            kestold = .TRUE.
            k_factor = 2.0d0
        ELSE
            kestold = .TRUE.
        ENDIF
        !     Record the old step size and Krylov dimension
        t_old = t_step
        m_old = m
        !     Suggest new step size and Krylov dimension
        IF ((m.EQ.m_max).AND.(omega.GT.delta).OR.(imreject.GT.4)) THEN
            t_new = MIN(t_out - t_now, &
                    MAX(t_step / 5.0d0, &
                            MIN(5.0d0 * t_step, &
                                    gamma * t_step * omega**(-1.0d0 / order))))
            p1 = 10.0d0**(NINT(LOG10(t_new) - sqr1) - 1)
            t_new = AINT(t_new / p1) * p1
            m_changed = .FALSE.
        ELSE
            !           Compute options for the new step size and Krylov dimension
            t_opt = MIN(t_out - t_now, &
                    MAX(t_step / 5.0d0, &
                            MIN(5.0d0 * t_step, &
                                    gamma * t_step * omega**(-1.0d0 / order))))
            m_opt = MIN(MAX(m_min, 3 * m / 4, &
                    m + CEILING(LOG(omega) / LOG(k_factor))), &
                    m_max, &
                    CEILING(4.0d0 * m / 3.0d0) + 1)
            !           Estimate costs
            nom = 25.0d0 / 3.0d0 + &
                    MAX(0, 2 + INT(LOG(t_opt * hnorm) / LOG(2.0d0)))
            cost1 = NINT((t_out - t_now) / t_opt) * &
                    (2 * (m + 1) * nnz + (5 * m + 4 * qiop * m + &
                            2 * qiop - 2 * qiop * qiop + 7) * n + 2 * nom * (m + 2) * (m + 2) * (m + 2))
            nom = 25.0d0 / 3.0d0 + &
                    MAX(0, 2 + INT(LOG(t_step * hnorm) / LOG(2.0d0)))
            cost2 = NINT((t_out - t_now) / t_step) * &
                    (2 * (m_opt + 1) * nnz + &
                            (5 * m_opt + 4 * qiop * m_opt + 2 * qiop - 2 * qiop * qiop + 7) * n + &
                            2 * nom * (m_opt + 2) * (m_opt + 2) * (m_opt + 2))
            IF (cost1.LE.cost2) THEN
                t_new = t_opt
                p1 = 10.0d0**(NINT(LOG10(t_new) - sqr1) - 1)
                t_new = AINT(t_new / p1) * p1
                m_new = m
                m_changed = .FALSE.
            ELSE
                m_new = m_opt
                t_new = t_step
                m_changed = .TRUE.
            ENDIF
        ENDIF
        !     If the Krylov step error is not acceptable, reject the stepsize
        IF ((k1.NE.0).AND.(omega.GT.delta).AND.&
                (mxreject.EQ.0.OR.ireject.LT.mxreject)) THEN
            IF (.NOT.m_changed) THEN
                !                 Choose to change step size
                t_step = MIN(t_out - t_now, &
                        MAX(t_step / 5.0d0, &
                                MIN(5.0d0 * t_step, t_new)))
                p1 = 10.0d0**(NINT(LOG10(t_step) - sqr1) - 1)
                t_step = AINT(t_step / p1 + 0.55d0) * p1
                IF (itrace.NE.0) THEN
                    PRINT*, 't_step =', t_old
                    PRINT*, 'err_loc =', err_loc
                    PRINT*, 'err_required =', delta * t_old * KryTol
                    PRINT*, 'stepsize rejected, down to:', t_step
                ENDIF
                ireject = ireject + 1
                nreject = nreject + 1
                IF ((mxreject.NE.0).AND.(ireject.GT.mxreject)) THEN
                    PRINT*, 'Failure in DGEXPV: ---'
                    PRINT*, 'The requested tolerance is too high.'
                    PRINT*, 'Rerun with a smaller value.'
                    iflag = 2
                    RETURN
                ENDIF
                go to 401
            ELSE
                !                 Choose to change dimension
                nreject = nreject + 1
                imreject = imreject + 1
                m = m_new
                htmp(1:ifree - ih + 1) = wsp(ih:ifree)
                ih_old = ih
                mbrkdwn = m
                k1 = 2
                mh = m + 2
                iv = 1
                ih = iv + n * (m + 1) + n
                ifree = ih + mh * mh
                lfree = lwsp - ifree + 1
                t_step = MIN(t_out - t_now, t_new)
                DO i = 1, mh * mh
                    wsp(ih + i - 1) = 0.0d0
                ENDDO
                !                 Copy the Hessenberg matrix
                DO j = 1, m_old
                    DO i = 1, j + 1
                        wsp(ih + (j - 1) * (m + 2) + i - 1) = &
                                htmp(1 + (j - 1) * (m_old + 2) + i - 1)
                    ENDDO
                ENDDO
                !                 Only do Arnoldi from the m_old-th column
                jold = m_old
                IF (itrace.NE.0) THEN
                    PRINT*, 'err_loc =', err_loc
                    PRINT*, 'err_required =', delta * t_old * KryTol
                    PRINT*, 'Dimension changed into m =', m
                ENDIF
                go to 101
            ENDIF
        ENDIF
        imreject = 0
        jold = 1
        IF (err_loc.LT.1.0d-16) t_new = MAX(t_new, 2.0d0 * t_step)
        mx = mbrkdwn + MAX(0, k1 - 1)
        irejectfsp = 0

        ! Ensure the FSP-like criteria
        DO
            !     Find w = beta*V*exp(t*H)*e1
            CALL dgemv('n', n, mx, beta, wsp(iv), n, wsp(iexph), 1, 0.0d0, w, 1)

            !     Ensure FSP non-negativity and compute the probability sum
            DO i = 1, fsp%size
                IF (w(i).LT.0.0d0)      w(i) = 0.0d0
            ENDDO
            wsum = dasum(fsp%size, w, 1)

            PRINT*, "wsum= ", wsum

            error = wsum_old - wsum
            t_ratio = (t_now + t_step) / t_out

            !     Ensure FSP-like criteria...
            IF (wsum >= (1.0d0 - ferrorbound(t_now + t_step))) THEN
                EXIT
            ELSE
                iexpand = 1
                irejectfsp = irejectfsp + 1

                !           If the FSP has been rejected too many times, jumpt to SSA
                !           loop to expand the state space
                IF (irejectfsp.GE.5) THEN
                    w(1:n) = beta * wsp(iv:iv + n - 1)
                    nstep = nstep - 1
                    t_ssa = t_new
                    go to 404
                ELSEIF (irejectfsp.EQ.1) THEN
                    fsporder = 2
                ELSE
                    fsporder = LOG(error / errorold) / &
                            LOG(t_step / tau_old) - 1.0d0
                ENDIF

                !           Reduce the time step
                tfsp = gamma * t_step * (FSPTol * t_step / (error * t_out))**&
                        (1.0d0 / (fsporder))
                errorold = error
                tau_old = t_step

                t_step = MIN(t_out - t_now, &
                        MAX(t_step / 5.0d0, MIN(0.9d0 * t_step, tfsp)))
                p1 = 10.0d0**(NINT(LOG10(t_step) - sqr1) - 1)
                t_step = AINT(t_step / p1 + 0.55d0) * p1

                nexph = nexph + 1
                CALL dgpadm(ideg, mx, sgn * t_step, wsp(ih), mh, &
                        wsp(ifree), lfree, iwsp, iexph, ns, iflag)
                iexph = ifree + iexph - 1
                nscale = nscale + ns
            ENDIF
        ENDDO

        !     Update the time covered
        t_now = t_now + t_step
        wsum_old = wsum

        !     Display information about the last successful time step
        IF (itrace.NE.0) CALL print_stats

        !     If the integration is at the end time point, then the algorithm
        !     is finished
        IF (t_now.GE.t_out) go to 500

        !     Reduce the state space to only substantial states
        IF (nstep > 1 .AND. iexpand /= 1) THEN
            dsum = wsum - (1.0d0 - ferrorbound(t_now))
            IF (dsum.GT.0.0d0) CALL drop_states(w, fsp, model, dsum, fmatvec)
        ENDIF

        !-------------------------------------------State space expansion by SSA
        !-----------------------------in the case where the stepsize was reduced
        404    CONTINUE

        IF ((iexpand.EQ.1).AND.(t_now.LT.t_out)) THEN

            IF (nstep.EQ.1) t_new = t_step
            t_ssa = MIN(t_new, t_out - t_now)

            IF (itrace.NE.0) THEN
                PRINT*, 'calling SSA'
            ENDIF

            ! Extend the state space by SSA and 1-step reachability
            CALL SSA_extender(t_ssa, fsp, model)
            CALL onestep_extender(fsp, model)

            IF (n_now.GT.nfull) n_now = nfull
            iexpand = 0

        ENDIF

        !     Estimate the number of nonzero entries in current FSP matrix
        nnz = (model%nreactions + 1) * fsp%size
        !     Estimate other parameters
        n_now = fsp%size
        beta = dnrm2(n_now, w, 1)
        hump = MAX(hump, beta)
        err_loc = MAX(err_loc, rndoff)
        step_min = MIN(step_min, t_step)
        step_max = MAX(step_max, t_step)
        s_error = s_error + err_loc
        x_error = MAX(x_error, err_loc)
        p1 = 10.0d0**(NINT(LOG10(t_new) - sqr1) - 1)
        t_new = AINT(t_new / p1 + 0.55d0) * p1
        !     Return to the beginning for the next time step
        IF ((mxstep.EQ.0).OR.(nstep.LT.mxstep)) go to 100
        iflag = 1

        !-------------------------------------------Output important information
        500    CONTINUE
        !     In iwsp
        iwsp(1) = nmult
        iwsp(2) = nexph
        iwsp(3) = nscale
        iwsp(4) = nstep
        iwsp(5) = nreject
        iwsp(6) = ibrkflag
        iwsp(7) = mbrkdwn
        !     In wsp
        wsp(1) = step_min
        wsp(2) = step_max
        wsp(3) = 0.0d0
        wsp(4) = 0.0d0
        wsp(5) = x_error
        wsp(6) = s_error
        wsp(7) = tbrkdwn
        wsp(8) = sgn * t_now
        wsp(9) = hump / vnorm
        wsp(10) = beta / vnorm
    CONTAINS

        !===========================Compute the matrix-vector product in the CME
        SUBROUTINE fmatvec(x, y, matrix)
            IMPLICIT NONE
            !------------------------------------------------Input/output parameters
            !     Vector x
            DOUBLE PRECISION :: x(*)
            !     Vector y=A*x, where A is the FSP matrix
            DOUBLE PRECISION :: y(*)
            !     Current FSP matrix
            TYPE(fsp_matrix) :: matrix
            !--------------------------------------------Parameters of the algorithm
            INTEGER :: i, j, k, n, bw
            !----------------------------------------------------------Compute y=A*x
            !     Find the size of the matrix
            n = matrix%size
            bw = SIZE(matrix%adj, dim = 1)
            !     Initialize vector y
            DO i = 1, n
                y(i) = 0.0d0
            ENDDO

            !     Compute y=A*x
            DO i = 1, n
                DO j = 1, bw
                    k = matrix%adj(j, i)
                    IF (k.GE.1) THEN
                        y(k) = y(k) + matrix%offdiag(j, i) * x(i)
                    ENDIF
                ENDDO
                y(i) = y(i) - matrix%diag(i) * x(i)
            ENDDO
        END SUBROUTINE fmatvec

        DOUBLE PRECISION FUNCTION ferrorbound(tx)
            IMPLICIT NONE
            !-------------------------------------------------------Input parameters
            !     A time interval
            DOUBLE PRECISION, INTENT(in) :: tx
            !-----------------Compute the FSP error tolerance for this time interval
            ferrorbound = tx * FSPTol / t_out
            !     ferrorbound = FSPTol*(tx/t_out)**(2.0d0)
            !     ferrorbound = FSPTol*(tx/t_out)**(10.0d0)
        END FUNCTION ferrorbound

        SUBROUTINE print_stats
            IMPLICIT NONE

            PRINT*, 'integration', nstep, '-------------------------------'
            PRINT*, 'FSP size         =', fsp%size
            PRINT*, 'scale-square     =', ns
            PRINT*, 'step_size        =', t_step
            PRINT*, 'err_loc          =', err_loc
            PRINT*, 'next_step        =', t_new
            PRINT*, 't_now            =', t_now
            PRINT*, 'Krylov dimension =', m
            PRINT*, 'qiop             =', qiop
            PRINT*, 'entropy          =', ent
            PRINT*, 'norm(Hm)         =', hnorm

        END SUBROUTINE print_stats

    END SUBROUTINE dgexpv_fsp

END MODULE KrylovSolver
