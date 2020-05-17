PROGRAM TestModelParser

  ! Test inputting the model from file using model_manager module

  USE ModelModule
  IMPLICIT NONE

  TYPE (cme_model) :: model, model2

  INTEGER :: i, j, r
  DOUBLE PRECISION :: p1, p2, error

  CALL model%load('models/toggle_test_model.input')

  CALL model%reset_parameters((/5000.0d0, 1600.d0, 1.0d0, 1.0d0/))

  WRITE (*,*) 'Model loaded.'
  WRITE (*,*) 'Number of species: ',model%nspecies
  WRITE (*,*) 'Species names: ', model%species_names
  WRITE (*,*) 'Number of parameters:', model%nparameters
  WRITE (*,*) 'Parameter values:', model%parameter_val
  WRITE (*,*) 'Stoichiometry:'
  DO i = 1,model%nspecies
     DO j = 1, model%nreactions
        WRITE (*,fmt = '(I2 4X)',advance='no') model%stoichiometry(i,j)
     ENDDO
     WRITE (*,*) ''
  END DO


  error = 0.0d0
  ! Now test the accuracy of the read-from-file propensities
  DO i = 1, 50
     DO j = 1,50
        DO r = 1, model%nreactions
           p1 = prop((/ i, j /), r, (/0d0 , 0d0/))
           p2 = model%propensity((/i, j/), r)
           error = error + ABS(p1-p2)
           WRITE (*, 201) i, j, r, p1, p2, ABS(p1-p2)
201        FORMAT (1X, 'ie=', I4,2X, 'j=', I4, 2X, 'r=',I4,2X, 'p1=', E9.2, 2X, 'p2=', E9.2, 2X, '|p2-p1|=',E9.2)
        ENDDO
     ENDDO
  ENDDO

  PRINT*,'Propensity testing completed for model 1. Accumulated error = ',error


  CALL model2%load('models/toggle_test_model.input')
  CALL model2%reset_parameters((/5000.0d0, 1600.d0, 1.0d0, 1.0d0/))

  WRITE (*,*) 'Model loaded.'
  WRITE (*,*) 'Number of species: ',model2%nspecies
  WRITE (*,*) 'Species names: ', model2%species_names
  WRITE (*,*) 'Number of parameters:', model2%nparameters
  WRITE (*,*) 'Parameter values:', model2%parameter_val
  WRITE (*,*) 'Stoichiometry:'
  DO i = 1,model2%nspecies
     DO j = 1, model2%nreactions
        WRITE (*,fmt = '(I2 4X)',advance='no') model2%stoichiometry(i,j)
     ENDDO
     WRITE (*,*) ''
  END DO
  model2%customprop=>prop
  error = 0.0d0
  ! Now test the accuracy of the read-from-file propensities
  DO i = 1, 50
     DO j = 1,50
        DO r = 1, model2%nreactions
           p1 = prop( (/i, j /),r)
           p2 = model2%propensity((/i, j/), r)
           error = error + ABS(p1-p2)
           WRITE (*, 201) i, j, r, p1, p2, ABS(p1-p2)
        ENDDO
     ENDDO
  ENDDO

  PRINT*,'Propensity testing completed for model 2. Accumulated error = ',error

CONTAINS
  DOUBLE PRECISION FUNCTION prop(state, reaction, parameters)
    ! propensity function for the toggle switch model
    IMPLICIT NONE
    INTEGER, INTENT(in) :: state(:), reaction
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: parameters(:)
    DOUBLE PRECISION c
    !------Generate individual values of propensity functions--
    !----------------------------------------------------------
    SELECT CASE (reaction)
    CASE (1)
       c=5000.0d0
       prop=c/(1.0d0+DBLE(state(2))**(2.5d0))
    CASE (2)
       c=1600.0d0
       prop=c/(1.0d0+DBLE(state(1))**(1.5d0))
    CASE (3)
       c=1.0d0*DBLE(state(1))
       prop=c
    CASE (4)
       c=1.0d0*DBLE(state(2))
       prop=c
    END SELECT
  END FUNCTION prop

END
