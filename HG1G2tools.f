! Fortran 2003 implementation of H,G1,G2 tools
!
! Can be compiled with, at least, gfortran verion 4.7 using
! options -std=f2003 -ffree-form
!
! Antti Penttilä, 2012
! Department of Physics, University of Helsinki
! see http://wiki.helsinki.fi/display/PSR/HG1G2+tools
!
MODULE HG1G2tools
  IMPLICIT NONE
  
  PRIVATE
  
  ! Error codes
  INTEGER, PARAMETER :: ECODE_FILE_NOT_FOUND = -1, ECODE_FILE_WRONG_FORMAT = -2, ECODE_NOT_IMPLEMENTED = -99, &
    ECODE_ALLOCATION = -3, ECODE_NOT_INITED = -4, ECODE_BELOW_RANGE = -5, ECODE_ABOVE_RANGE = -6, ECODE_OVERFLOW = -7, &
    ECODE_SINGULAR_MATRIX = -8, ECODE_CHOLESKY_FAIL = -9
  ! Double precision
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(15)
  ! Other parameters
  INTEGER, PARAMETER :: max_error_str = 512, func_type_constant = 1, func_type_linear = 2, &
    func_type_spline = 3, default_error_sim_n = 100000
  
  ! Degree (to radians) and radian (to degrees)
  REAL(dbl), PARAMETER :: degree = 0.0174532925199433_dbl, radian = 57.2957795130823_dbl, &
    pi = 3.14159265358979_dbl, DEFAULT_MAG_ERROR = 0.03_dbl, &
    coef_G1_small = 0.7527_dbl, const_G1_small = 0.06164_dbl, coef_G1_large = 0.9529_dbl, const_G1_large = 0.02162_dbl, &
    coef_G2_small = -0.9612_dbl, const_G2_small = 0.6270_dbl, coef_G2_large = -0.6125_dbl, const_G2_large = 0.5572_dbl
  
  ! Common error strings
  CHARACTER(LEN=*), PARAMETER :: ESTR_NOT_IMPLEMENTED = "This feature is not yet implemented", &
    ESTR_ALLOCATION = "Error in memory allocation", ESTR_NOT_INITED = "Base functions are not initialized", &
    ESTR_SINGULAR_MATRIX = "Datamatrix is singular", &
    ESTR_CHOLESKY_FAIL = "Cholesky decomposition of covariance matris failed, continuing with uncorrelated paramaters"
  
  ! Version number of the base functions
  INTEGER :: HG1G2version
  INTEGER :: statcode
  
  LOGICAL :: base_inited = .FALSE.

  TYPE error_type
    INTEGER :: code = 0
    CHARACTER(LEN=max_error_str) :: msg = ""
  END TYPE error_type
  
  TYPE spline_type
    INTEGER :: n_knots
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: knots
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: coef
  END TYPE spline_type
  
  TYPE generic_func_type
    INTEGER :: func_type
    REAL(dbl) :: low_limit, high_limit
    REAL(dbl) :: constant, coef
    TYPE(spline_type) :: spline  
  END TYPE generic_func_type
  
  TYPE base_func_type
    INTEGER :: no_of_funcs
    REAL(dbl) :: low_limit, high_limit
    TYPE(generic_func_type), DIMENSION(:), ALLOCATABLE :: funcs
  END TYPE base_func_type
  
  TYPE(error_type) :: last_error
  TYPE(base_func_type), DIMENSION(3) :: base
  
  INTERFACE form_base
    MODULE PROCEDURE form_base_default, form_base_from_file
  END INTERFACE form_base
  
  PUBLIC :: form_base, get_last_error, give_base_values, system_version, fit_HG1G2, fit_HG12, &
    simulate_errors, &
    dbl, degree, radian, max_error_str, DEFAULT_MAG_ERROR, &
    ECODE_FILE_NOT_FOUND, ECODE_FILE_WRONG_FORMAT, ECODE_NOT_IMPLEMENTED, ECODE_ALLOCATION, ECODE_NOT_INITED, &
    ECODE_BELOW_RANGE, ECODE_ABOVE_RANGE, ECODE_OVERFLOW, ECODE_SINGULAR_MATRIX, ECODE_CHOLESKY_FAIL
  
CONTAINS

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Public functions
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize base functions a1,a2,a3 with default system version 20101000
SUBROUTINE form_base_default(OP_STAT)
  IMPLICIT NONE
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  
  LOGICAL :: stat_present
  INTEGER :: n
  
  stat_present = PRESENT(OP_STAT)

  HG1G2version = 20101000

  ! base function a1
  base(1)%no_of_funcs = 2
  base(1)%low_limit = 0.0_dbl
  base(1)%high_limit = 150.0_dbl * degree
  ALLOCATE(base(1)%funcs(2), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(stat_present) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  ! first part of a1
  base(1)%funcs(1)%func_type = func_type_linear
  base(1)%funcs(1)%low_limit = 0.0_dbl
  base(1)%funcs(1)%high_limit = 7.5_dbl * degree
  base(1)%funcs(1)%constant = 1.0_dbl
  base(1)%funcs(1)%coef = -1.90985931710274_dbl
  ! second part of a1
  base(1)%funcs(2)%func_type = func_type_spline
  base(1)%funcs(2)%low_limit = 7.5_dbl * degree
  base(1)%funcs(2)%high_limit = 150.0_dbl * degree
  n = 6
  base(1)%funcs(2)%spline%n_knots = n
  ALLOCATE(base(1)%funcs(2)%spline%knots(n), base(1)%funcs(2)%spline%coef(n-1,4), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(stat_present) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  CALL spline_coefs(degree*(/ 7.5_dbl, 30.0_dbl, 60.0_dbl, 90.0_dbl, 120.0_dbl, 150.0_dbl /), &
    (/ 0.75_dbl,0.33486_dbl,0.134106_dbl,0.0511048_dbl,0.0214657_dbl,0.0036397_dbl /), &
    (/ -1.90986_dbl, -0.0913286_dbl /), 1, 2)

  ! base function a2
  base(2)%no_of_funcs = 2
  base(2)%low_limit = 0.0_dbl
  base(2)%high_limit = 150.0_dbl * degree
  ALLOCATE(base(2)%funcs(2), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(stat_present) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  ! first part of a2
  base(2)%funcs(1)%func_type = func_type_linear
  base(2)%funcs(1)%low_limit = 0.0_dbl
  base(2)%funcs(1)%high_limit = 7.5_dbl * degree
  base(2)%funcs(1)%constant = 1.0_dbl
  base(2)%funcs(1)%coef = -0.572957795130823_dbl
  ! second part of a2
  base(2)%funcs(2)%func_type = func_type_spline
  base(2)%funcs(2)%low_limit = 7.5_dbl * degree
  base(2)%funcs(2)%high_limit = 150.0_dbl * degree
  n = 6
  base(2)%funcs(2)%spline%n_knots = n
  ALLOCATE(base(2)%funcs(2)%spline%knots(n), base(2)%funcs(2)%spline%coef(n-1,4), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(stat_present) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  CALL spline_coefs(degree*(/ 7.5_dbl, 30.0_dbl, 60.0_dbl, 90.0_dbl, 120.0_dbl, 150.0_dbl /), &
    (/ 0.925_dbl,0.628842_dbl,0.317555_dbl,0.127164_dbl,0.0223739_dbl,0.000165057_dbl /), &
    (/ -0.572958_dbl, -8.6573138e-8_dbl /), 2, 2)

  ! base function a3
  base(3)%no_of_funcs = 2
  base(3)%low_limit = 0.0_dbl
  base(3)%high_limit = 150.0_dbl * degree
  ALLOCATE(base(3)%funcs(2), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(stat_present) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  ! second part of a3
  base(3)%funcs(2)%func_type = func_type_constant
  base(3)%funcs(2)%low_limit = 30.0_dbl * degree
  base(3)%funcs(2)%high_limit = 150.0_dbl * degree
  base(3)%funcs(2)%constant = 0.0_dbl
  ! first part of a3
  base(3)%funcs(1)%func_type = func_type_spline
  base(3)%funcs(1)%low_limit = 0.0_dbl
  base(3)%funcs(1)%high_limit = 30.0_dbl * degree
  n = 9
  base(3)%funcs(1)%spline%n_knots = n
  ALLOCATE(base(3)%funcs(1)%spline%knots(n), base(3)%funcs(1)%spline%coef(n-1,4), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(stat_present) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  CALL spline_coefs(degree*(/ 0.0_dbl, 0.3_dbl, 1.0_dbl, 2.0_dbl, 4.0_dbl, 8.0_dbl, 12.0_dbl, 20.0_dbl, 30.0_dbl /), &
    (/ 1._dbl,0.833812_dbl,0.577354_dbl,0.421448_dbl,0.231742_dbl,0.103482_dbl,0.0617335_dbl,0.016107_dbl,0.0_dbl /), &
    (/ -0.106301_dbl, 0.0_dbl /), 3, 1)

  base_inited = .TRUE.
  IF(PRESENT(OP_STAT)) OP_STAT = 0

END SUBROUTINE form_base_default

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Not implemented yet
SUBROUTINE form_base_from_file(fn, OP_STAT)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: fn
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT

  last_error%code = ECODE_NOT_IMPLEMENTED
  last_error%msg = ESTR_NOT_IMPLEMENTED
  IF(PRESENT(OP_STAT)) THEN
    OP_STAT = ECODE_NOT_IMPLEMENTED
    RETURN
  END IF
  STOP

END SUBROUTINE form_base_from_file

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fit H,G1,G2 system to observations.
FUNCTION fit_HG1G2(d, DEGR, WEIGHTED, RMS, COVMAT, OP_STAT) RESULT(params)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:,:), INTENT(IN) :: d
  LOGICAL, INTENT(IN), OPTIONAL :: DEGR
  REAL(dbl), INTENT(IN), OPTIONAL :: WEIGHTED
  REAL(dbl), INTENT(OUT), OPTIONAL :: RMS
  REAL(dbl), DIMENSION(3,3), INTENT(OUT), OPTIONAL :: COVMAT
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  REAL(dbl), DIMENSION(3) :: params
  
  INTEGER :: i, j, n
  REAL(dbl) :: x, y
  REAL(dbl), DIMENSION(2) :: resvar
  REAL(dbl), DIMENSION(3) :: as
  REAL(dbl), DIMENSION(SIZE(d,1)) :: xvec, yvec, wvec
  REAL(dbl), DIMENSION(SIZE(d,1),3) :: Amat
  REAL(dbl), DIMENSION(3,3) :: Cmat

  IF(.NOT. base_inited) THEN
    last_error%code = ECODE_NOT_INITED
    last_error%msg = ESTR_NOT_INITED
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = ECODE_NOT_INITED
      RETURN
    END IF
    STOP
  END IF
  
  IF(PRESENT(DEGR)) THEN
    IF(DEGR) THEN
      xvec(:) = degree * d(:,1)
    ELSE
      xvec(:) = d(:,1)
    END IF
  ELSE
    xvec(:) = degree * d(:,1)
  END IF

  yvec(:) = 10.0_dbl**(-0.4_dbl * d(:,2))
  
  IF(PRESENT(WEIGHTED)) THEN
    IF(SIZE(d,2) == 3) THEN
      ! If data includes errors in magnitude space, always use these.
      ! Change to flux space
      wvec(:) = 1.0_dbl/(yvec(:)*(10**(0.4_dbl*d(:,3))-1.0_dbl))
    ELSE IF(WEIGHTED < 0) THEN
      ! Default errors in magnitude space, sigma~DEFAULT_MAG_ERROR
      x = (10.0_dbl**(0.4_dbl*DEFAULT_MAG_ERROR)-1.0_dbl)
      wvec(:) = 1.0_dbl/(yvec(:)*x)
    ELSE
      ! Constant error in magnitude space, given by user
      x = (10.0_dbl**(0.4_dbl*WEIGHTED)-1.0_dbl)
      wvec(:) = 1.0_dbl/(yvec(:)*x)
    END IF
  END IF

  n = SIZE(d,1)
  DO i=1,n
    Amat(i,1) = get_value(xvec(i), 1, statcode)
    IF(statcode /= 0) THEN
      IF(PRESENT(OP_STAT)) THEN
        OP_STAT = statcode
        RETURN
      ELSE
        STOP
      END IF
    END IF
    Amat(i,2) = get_value(xvec(i), 2, statcode)
    IF(statcode /= 0) THEN
      IF(PRESENT(OP_STAT)) THEN
        OP_STAT = statcode
        RETURN
      ELSE
        STOP
      END IF
    END IF
    Amat(i,3) = get_value(xvec(i), 3, statcode)
    IF(statcode /= 0) THEN
      IF(PRESENT(OP_STAT)) THEN
        OP_STAT = statcode
        RETURN
      ELSE
        STOP
      END IF
    END IF      
  END DO

  IF(PRESENT(WEIGHTED)) THEN
    DO i=1,n
      Amat(i,:) = wvec(i) * Amat(i,:)
      yvec(i) = wvec(i) * yvec(i)
    END DO
  END IF
  
  CALL least_squares(Amat, yvec, as, Cmat)
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = statcode
      RETURN
    ELSE
      STOP
    END IF
  END IF
  
  IF(PRESENT(RMS) .OR. PRESENT(COVMAT)) THEN
    resvar(:) = residual_variance(Amat, yvec, as)
    IF(PRESENT(RMS)) RMS = SQRT(resvar(2))
    IF(PRESENT(COVMAT)) COVMAT(:,:) = Cmat(:,:)
  END IF

  params(:) = a1a2a3_to_HG1G2(as)
  
  IF(PRESENT(OP_STAT)) OP_STAT = 0

END FUNCTION fit_HG1G2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fit H,G12 system to observations.
FUNCTION fit_HG12(d, DEGR, WEIGHTED, RMS, COVMAT, OP_STAT) RESULT(params)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:,:), INTENT(IN) :: d
  LOGICAL, INTENT(IN), OPTIONAL :: DEGR
  REAL(dbl), INTENT(IN), OPTIONAL :: WEIGHTED
  REAL(dbl), INTENT(OUT), OPTIONAL :: RMS
  REAL(dbl), DIMENSION(2,2), INTENT(OUT), OPTIONAL :: COVMAT
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  REAL(dbl), DIMENSION(2) :: params
  
  INTEGER :: i, j, n
  REAL(dbl) :: x, y, b1, b2, b3
  REAL(dbl), DIMENSION(2) :: resvar_small, resvar_large, resvar
  REAL(dbl), DIMENSION(2) :: as_small, as_large
  REAL(dbl), DIMENSION(SIZE(d,1)) :: xvec, yvec, wvec
  REAL(dbl), DIMENSION(SIZE(d,1),2) :: Amat_small, Amat_large
  REAL(dbl), DIMENSION(2,2) :: Cmat_small, Cmat_large, Cmat

  IF(.NOT. base_inited) THEN
    last_error%code = ECODE_NOT_INITED
    last_error%msg = ESTR_NOT_INITED
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = ECODE_NOT_INITED
      RETURN
    END IF
    STOP
  END IF
  
  IF(PRESENT(DEGR)) THEN
    IF(DEGR) THEN
      xvec(:) = degree * d(:,1)
    ELSE
      xvec(:) = d(:,1)
    END IF
  ELSE
    xvec(:) = degree * d(:,1)
  END IF

  yvec(:) = 10.0_dbl**(-0.4_dbl * d(:,2))
  
  IF(PRESENT(WEIGHTED)) THEN
    IF(SIZE(d,2) == 3) THEN
      ! If data includes errors in magnitude space, always use these.
      ! Change to flux space
      wvec(:) = yvec(:)*(10**(0.4_dbl*d(:,3))-1.0_dbl)
    ELSE IF(WEIGHTED < 0) THEN
      ! Default errors in magnitude space, sigma~DEFAULT_MAG_ERROR
      x = (10.0_dbl**(0.4_dbl*DEFAULT_MAG_ERROR)-1.0_dbl)
      wvec(:) = 1.0_dbl/(yvec(:)*x)
    ELSE
      ! Constant error in magnitude space, given by user
      x = (10.0_dbl**(0.4_dbl*WEIGHTED)-1.0_dbl)
      wvec(:) = 1.0_dbl/(yvec(:)*x)
    END IF
  END IF

  n = SIZE(d,1)
  DO i=1,n
    b1 = get_value(xvec(i), 1, statcode)
    IF(statcode /= 0) THEN
      IF(PRESENT(OP_STAT)) THEN
        OP_STAT = statcode
        RETURN
      ELSE
        STOP
      END IF
    END IF
    b2 = get_value(xvec(i), 2, statcode)
    IF(statcode /= 0) THEN
      IF(PRESENT(OP_STAT)) THEN
        OP_STAT = statcode
        RETURN
      ELSE
        STOP
      END IF
    END IF
    b3 = get_value(xvec(i), 3, statcode)
    IF(statcode /= 0) THEN
      IF(PRESENT(OP_STAT)) THEN
        OP_STAT = statcode
        RETURN
      ELSE
        STOP
      END IF
    END IF
    Amat_small(i,1) = const_G1_small*b1 + const_G2_small*b2 + (1-const_G1_small-const_G2_small)*b3
    Amat_large(i,1) = const_G1_large*b1 + const_G2_large*b2 + (1-const_G1_large-const_G2_large)*b3
    Amat_small(i,2) = coef_G1_small*b1 + coef_G2_small*b2 - (coef_G1_small+coef_G2_small)*b3
    Amat_large(i,2) = coef_G1_large*b1 + coef_G2_large*b2 - (coef_G1_large+coef_G2_large)*b3
  END DO

  IF(PRESENT(WEIGHTED)) THEN
    DO i=1,n
      Amat_small(i,:) = wvec(i) * Amat_small(i,:)
      Amat_large(i,:) = wvec(i) * Amat_large(i,:)
      yvec(i) = wvec(i) * yvec(i)
    END DO
  END IF
  
  CALL least_squares2(Amat_small, yvec, as_small, Cmat_small)
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = statcode
      RETURN
    ELSE
      STOP
    END IF
  END IF
  CALL least_squares2(Amat_large, yvec, as_large, Cmat_large)
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = statcode
      RETURN
    ELSE
      STOP
    END IF
  END IF
  
  resvar_small(:) = residual_variance(Amat_small, yvec, as_small)
  resvar_large(:) = residual_variance(Amat_large, yvec, as_large)
  
  ! Select best fit
  IF(resvar_small(1) < resvar_large(1)) THEN
    params(:) = a1a2_to_HG12(as_small)
    resvar(:) = resvar_small(:)
    Cmat(:,:) = Cmat_small(:,:)
  ELSE
    params(:) = a1a2_to_HG12(as_large)
    resvar(:) = resvar_large(:)
    Cmat(:,:) = Cmat_large(:,:)
  END IF
  
  IF(PRESENT(RMS) .OR. PRESENT(COVMAT)) THEN
    IF(PRESENT(RMS)) RMS = SQRT(resvar(2))
    IF(PRESENT(COVMAT)) COVMAT(:,:) = Cmat(:,:)
  END IF
  
  IF(PRESENT(OP_STAT)) OP_STAT = 0

END FUNCTION fit_HG12

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Suits both for HG1G2 and HG12 -systems
FUNCTION simulate_errors(params, covmat, SIM_N, INIT_RANDOM, OP_STAT) RESULT(errmat)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:), INTENT(IN) :: params
  REAL(dbl), DIMENSION(:,:), INTENT(IN) :: covmat
  INTEGER, INTENT(IN), OPTIONAL :: SIM_N
  LOGICAL, INTENT(IN), OPTIONAL :: INIT_RANDOM
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  REAL(dbl), DIMENSION(4,8) :: errmat
  
  INTEGER :: sn, i, pern, k
  REAL(dbl), DIMENSION(SIZE(params,1)) :: as
  REAL(dbl), DIMENSION(3) :: vec
  REAL(dbl), DIMENSION(7) :: percentiles
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: simH, simG1, simG2, simk
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: sim_as
  
  ! HG1G2 (k=3) or HG12 (k=2)
  k = SIZE(params,1)
  
  IF(PRESENT(SIM_N)) THEN
    sn = SIM_N
  ELSE
    sn = default_error_sim_n
  END IF
  
  IF(k == 2) THEN
    as(:) = HG12_to_a1a2(params)
  ELSE
    as(:) = HG1G2_to_a1a2a3(params)
  END IF
  
  ALLOCATE(simH(sn), simG1(sn), simG2(sn), simk(sn), sim_as(sn,k), STAT=statcode)
  IF(statcode /= 0) THEN
    last_error%code = ECODE_ALLOCATION
    last_error%msg = ESTR_ALLOCATION
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = ECODE_ALLOCATION
      RETURN
    END IF
    STOP
  END IF
  
  IF(PRESENT(INIT_RANDOM)) THEN
    IF(INIT_RANDOM) CALL init_rng()
  ELSE
    CALL init_rng()
  END IF

  statcode = 0
  IF(k == 2) THEN
    CALL normal2_random(as, covmat, sim_as, statcode)
  ELSE
    CALL normal3_random(as, covmat, sim_as, statcode)
  END IF
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      ! do not stop, only raise code
      OP_STAT = statcode
    ELSE
      STOP
    END IF
  END IF

  ! Optimize, same loop with IF outside
  IF(k == 2) THEN
    DO i=1,sn
      vec(:) = a1a2_to_HG1G2(sim_as(i,:))
      simH(i) = vec(1)
      simG1(i) = vec(2)
      simG2(i) = vec(3)
      simk(i) = -(30.0_dbl*vec(2)+9.0_dbl*vec(3)) / (5.0_dbl*pi*(vec(2)+vec(3)))
    END DO
  ELSE
    DO i=1,sn
      vec(:) = a1a2a3_to_HG1G2(sim_as(i,:))
      simH(i) = vec(1)
      simG1(i) = vec(2)
      simG2(i) = vec(3)
      simk(i) = -(30.0_dbl*vec(2)+9.0_dbl*vec(3)) / (5.0_dbl*pi*(vec(2)+vec(3)))
    END DO
  END IF

  CALL Qsort(simH, 1, sn)
  CALL Qsort(simG1, 1, sn)
  CALL Qsort(simg2, 1, sn)
  CALL Qsort(simk, 1, sn)

  ! Means
  errmat(1,1) = SUM(simH)/sn
  errmat(2,1) = SUM(simG1)/sn
  errmat(3,1) = SUM(simG2)/sn
  errmat(4,1) = SUM(simk)/sn
  ! Median and other percentiles, corresponding to -3*sigma, -2*sigma, -sigma, sigma, 2*sigma and 3*sigma for
  ! normal distribution, i.e. 99.73%, 95.45% and 68.3% confidence intervals
  percentiles(:) = (/ 0.5_dbl, 0.0013499_dbl, 0.0227501_dbl, 0.158655_dbl, 0.841345_dbl, 0.97725_dbl, 0.99865_dbl /)
  DO i=2,8
    pern = MAX(1,NINT(percentiles(i-1)*sn))
    errmat(1,i) = simH(pern)
    errmat(2,i) = simG1(pern)
    errmat(3,i) = simG2(pern)
    errmat(4,i) = simk(pern)
  END DO
  
END FUNCTION simulate_errors

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION system_version(LOWER, UPPER) RESULT(ver)
  IMPLICIT NONE
  REAL(dbl), INTENT(OUT), OPTIONAL :: LOWER, UPPER
  INTEGER :: ver
  
  IF(.NOT. base_inited) THEN
    last_error%code = ECODE_NOT_INITED
    last_error%msg = ESTR_NOT_INITED
    ver = -1
    RETURN
  END IF
  
  IF(PRESENT(LOWER)) THEN
    LOWER = MAXVAL(base(:)%low_limit)
  END IF
  IF(PRESENT(UPPER)) THEN
    UPPER = MINVAL(base(:)%high_limit)
  END IF
  
  ver = HG1G2version
  
END FUNCTION system_version

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Public only for debuggin purposes, you can check 
! the base functions
SUBROUTINE give_base_values(x, b1, b2, b3, OP_STAT)
  IMPLICIT NONE
  REAL(dbl), INTENT(IN) :: x
  REAL(dbl), INTENT(OUT) :: b1, b2, b3
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  
  IF(.NOT. base_inited) THEN
    last_error%code = ECODE_NOT_INITED
    last_error%msg = ESTR_NOT_INITED
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = ECODE_NOT_INITED
      RETURN
    END IF
    STOP
  END IF
  
  b1 = get_value(x, 1, statcode)
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = statcode
      RETURN
    END IF
  END IF
  b2 = get_value(x, 2, statcode)
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = statcode
      RETURN
    END IF
  END IF
  b3 = get_value(x, 3, statcode)
  IF(statcode /= 0) THEN
    IF(PRESENT(OP_STAT)) THEN
      OP_STAT = statcode
      RETURN
    END IF
  END IF
  
  IF(PRESENT(OP_STAT)) OP_STAT = 0
  
END SUBROUTINE give_base_values

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_last_error(code, MSG)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: code
  CHARACTER(LEN=max_error_str), INTENT(OUT), OPTIONAL :: MSG
  
  code = last_error%code
  IF(PRESENT(MSG)) THEN
    WRITE(MSG, '(A)') last_error%msg
  END IF

END SUBROUTINE get_last_error

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private functions
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute residual variance in both flux (rv(1)) and magnitude (rv(2)) space
! Suits for both HG1G2 and HG12 fits
FUNCTION residual_variance(Amat, yvec, as) RESULT(rv)
  IMPLICIT NONE
  REAL(dbl) :: f
  REAL(dbl), DIMENSION(:,:), INTENT(IN) :: Amat
  REAL(dbl), DIMENSION(:), INTENT(IN) :: yvec
  REAL(dbl), DIMENSION(:), INTENT(IN) :: as
  REAL(dbl), DIMENSION(2) :: rv
  
  INTEGER :: i, j, n, k

  n = SIZE(yvec,1)
  k = SIZE(as,1)
 
  rv(:) = 0.0_dbl
  DO i=1,n
    f = 0.0_dbl
    DO j=1,k
      f = f + as(j)*Amat(i,j)
    END DO
    rv(1) = rv(1) + (yvec(i) - f)**2
    rv(2) = rv(2) + LOG10(yvec(i)/f)**2
  END DO
  rv(1) = rv(1)/(n-k)
  rv(2) = 6.25_dbl*rv(2)/(n-k)

END FUNCTION residual_variance

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Linear parameters a1,a2,a3 to H,G1,G2.
FUNCTION a1a2a3_to_HG1G2(as) RESULT(pvec)
  IMPLICIT NONE
  REAL(dbl), INTENT(IN), DIMENSION(3) :: as
  REAL(dbl), DIMENSION(3) :: pvec
  
  REAL(dbl) :: x
  
  x = as(1)+as(2)+as(3)
  pvec(1) = -2.5_dbl * LOG10(x)
  pvec(2) = as(1)/x
  pvec(3) = as(2)/x
  
END FUNCTION a1a2a3_to_HG1G2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Noninear parameters H,G1,G2 to  a1,a2,a3.
FUNCTION HG1G2_to_a1a2a3(params) RESULT(res)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(3), INTENT(IN) :: params
  REAL(dbl), DIMENSION(3) :: res
  
  REAL(dbl) :: x
  
  x = EXP(-0.9210340371976184_dbl*params(1))
  res(1) = x * params(2)
  res(2) = x * params(3)
  res(3) = x * (1.0_dbl-params(2)-params(3))
  
END FUNCTION HG1G2_to_a1a2a3

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Linear parameters a1,a2 to H,G12.
FUNCTION a1a2_to_HG12(as) RESULT(pvec)
  IMPLICIT NONE
  REAL(dbl), INTENT(IN), DIMENSION(2) :: as
  REAL(dbl), DIMENSION(2) :: pvec
  
  pvec(1) = -2.5_dbl * LOG10(as(1))
  pvec(2) = as(2)/as(1)
  
END FUNCTION a1a2_to_HG12

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Noninear parameters H,G12 to  a1,a2.
FUNCTION HG12_to_a1a2(params) RESULT(res)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(2), INTENT(IN) :: params
  REAL(dbl), DIMENSION(2) :: res
  
  res(1) = EXP(-0.9210340371976184_dbl*params(1))
  res(2) = res(1) * params(2)
  
END FUNCTION HG12_to_a1a2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Linear parameters a1,a2 to H,G1,G2.
FUNCTION a1a2_to_HG1G2(as) RESULT(pvec)
  IMPLICIT NONE
  REAL(dbl), INTENT(IN), DIMENSION(2) :: as
  REAL(dbl), DIMENSION(3) :: pvec
  
  pvec(1) = -2.5_dbl * LOG10(as(1))
  IF(as(2) < 0.2_dbl) THEN
    pvec(2) = coef_G1_small*as(2) + const_G1_small
    pvec(3) = coef_G2_small*as(2) + const_G2_small
  ELSE
    pvec(2) = coef_G1_large*as(2) + const_G1_large
    pvec(3) = coef_G2_large*as(2) + const_G2_large
  END IF

END FUNCTION a1a2_to_HG1G2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! beta = (X'X)^-1 X' y, X=Amat, y=yvec, Cmat=(X'X)^-1
SUBROUTINE least_squares(Amat, yvec, as, Cmat)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:,:), INTENT(IN) :: Amat
  REAL(dbl), DIMENSION(:), INTENT(IN) :: yvec
  REAL(dbl), DIMENSION(3), INTENT(OUT) :: as
  REAL(dbl), DIMENSION(3,3), INTENT(OUT) :: Cmat

  INTEGER :: i, j, k, n
  REAL(dbl) :: x
  REAL(dbl), DIMENSION(3,3) :: XtrX
  
  n = SIZE(yvec,1)
  
  ! Computing X'X
  DO i=1,3
    DO j=1,i
      XtrX(i,j) = 0.0_dbl
      DO k=1,n
        XtrX(i,j) = XtrX(i,j) + Amat(k,i)*Amat(k,j)
      END DO
    END DO
  END DO
  XtrX(1,2) = XtrX(2,1)
  XtrX(1,3) = XtrX(3,1)
  XtrX(2,3) = XtrX(3,2)
  
  ! Computing (X'X)^-1
  x = -2.0_dbl*XtrX(2,1)*XtrX(3,1)*XtrX(3,2) + XtrX(1,1)*XtrX(3,2)**2 + XtrX(2,1)**2*XtrX(3,3) + &
    XtrX(2,2)*(XtrX(3,1)**2 - XtrX(1,1)*XtrX(3,3))
  IF(x == 0.0_dbl) THEN
    statcode = ECODE_SINGULAR_MATRIX
    last_error%code =  ECODE_SINGULAR_MATRIX
    last_error%msg =  ESTR_SINGULAR_MATRIX
    RETURN
  END IF
  Cmat(1,1) = (XtrX(3,2)**2 - XtrX(2,2)*XtrX(3,3))/x
  Cmat(1,2) = (-(XtrX(3,1)*XtrX(3,2)) + XtrX(2,1)*XtrX(3,3))/x
  Cmat(2,1) = Cmat(1,2)
  Cmat(1,3) = (XtrX(2,2)*XtrX(3,1) - XtrX(2,1)*XtrX(3,2))/x
  Cmat(3,1) = Cmat(1,3)
  Cmat(2,2) = (XtrX(3,1)**2 - XtrX(1,1)*XtrX(3,3))/x
  Cmat(2,3) = (-(XtrX(2,1)*XtrX(3,1)) + XtrX(1,1)*XtrX(3,2))/x
  Cmat(3,2) = Cmat(2,3)
  Cmat(3,3) = (XtrX(2,1)**2 - XtrX(1,1)*XtrX(2,2))/x
  
  ! Computing (X'X)^-1 X' y
  DO k=1,3
    as(k) = 0.0_dbl
    DO i=1,n
      x = 0.0_dbl
      DO j=1,3
        x = x + Cmat(k,j)*Amat(i,j)
      END DO
      as(k) = as(k) + x*yvec(i)
    END DO
  END DO

END SUBROUTINE least_squares

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! beta = (X'X)^-1 X' y, X=Amat, y=yvec, Cmat=(X'X)^-1
SUBROUTINE least_squares2(Amat, yvec, as, Cmat)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:,:), INTENT(IN) :: Amat
  REAL(dbl), DIMENSION(:), INTENT(IN) :: yvec
  REAL(dbl), DIMENSION(2), INTENT(OUT) :: as
  REAL(dbl), DIMENSION(2,2), INTENT(OUT) :: Cmat

  INTEGER :: i, j, k, n
  REAL(dbl) :: x
  REAL(dbl), DIMENSION(2,2) :: XtrX
  
  n = SIZE(yvec,1)
  
  ! Computing X'X
  DO i=1,2
    DO j=1,i
      XtrX(i,j) = 0.0_dbl
      DO k=1,n
        XtrX(i,j) = XtrX(i,j) + Amat(k,i)*Amat(k,j)
      END DO
    END DO
  END DO
  XtrX(1,2) = XtrX(2,1)
  
  ! Computing (X'X)^-1
  x = XtrX(1,1)*XtrX(2,2)-XtrX(1,2)**2
  IF(x == 0.0_dbl) THEN
    statcode = ECODE_SINGULAR_MATRIX
    last_error%code =  ECODE_SINGULAR_MATRIX
    last_error%msg =  ESTR_SINGULAR_MATRIX
    RETURN
  END IF
  Cmat(1,1) = XtrX(2,2)/x
  Cmat(1,2) = -XtrX(1,2)/x
  Cmat(2,1) = Cmat(1,2)
  Cmat(2,2) = XtrX(1,1)/x
  
  ! Computing (X'X)^-1 X' y
  DO k=1,2
    as(k) = 0.0_dbl
    DO i=1,n
      x = 0.0_dbl
      DO j=1,2
        x = x + Cmat(k,j)*Amat(i,j)
      END DO
      as(k) = as(k) + x*yvec(i)
    END DO
  END DO

END SUBROUTINE least_squares2

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spline coefficients from system A.coefs=r, where A is tri-diagonal.
! Tri-diagonal matrix solver adapted from Numerical Recipies.
! Note, dim(xval,yval)=n, but dim(spline coefs)=n-1
SUBROUTINE spline_coefs(xval, yval, deriv, base_i, func_i)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:), INTENT(IN) :: xval, yval
  REAL(dbl), DIMENSION(2), INTENT(IN) :: deriv
  INTEGER, INTENT(IN) :: base_i, func_i
  
  INTEGER :: i, j, n
  REAL(dbl) :: bet
  REAL(dbl), DIMENSION(SIZE(xval,1)) :: a, b, c, d, r, gam, u
  
  n = SIZE(xval,1)
  
  ! Spline system matrix A coefficients
  a(:) = 1.0_dbl
  a(1) = 0.0_dbl
  a(n) = 0.0_dbl
  b(:) = 4.0_dbl
  b(1) = 1.0_dbl
  b(n) = 1.0_dbl
  c(:) = 1.0_dbl
  c(1) = 0.0_dbl
  c(n) = 0.0_dbl
  
  r(1) = deriv(1)*(xval(2)-xval(1))
  DO i=2,n-1
    r(i) = 3.0_dbl * (yval(i+1)-yval(i-1))
  END DO
  r(n) = deriv(2)*(xval(n)-xval(n-1))
  
  ! Tri-diagonal solver
  bet=b(1)
  u(1)=r(1)/bet
  DO i=2,n
    gam(i)=c(i-1)/bet
    bet=b(i)-a(i)*gam(i)
    u(i)=(r(i)-a(i)*u(i-1))/bet
  END DO
  
  DO i=n-1,1,-1
    u(i)=u(i)-gam(i+1)*u(i+1)
  END DO
  
  ! System solved, convert to spline coefficients (a,b,c,d), when
  ! spline is y(t) = a + b*t + c*t^2 + d*t^3
  base(base_i)%funcs(func_i)%spline%knots(:) = xval(:)
  base(base_i)%funcs(func_i)%spline%coef(:,1) = yval(1:n-1)
  base(base_i)%funcs(func_i)%spline%coef(:,2) = u(1:n-1)
  DO i=1,n-1
    base(base_i)%funcs(func_i)%spline%coef(i,3) = 3.0_dbl*(yval(i+1)-yval(i)) - 2.0_dbl*u(i) - u(i+1)
    base(base_i)%funcs(func_i)%spline%coef(i,4) = 2.0_dbl*(yval(i)-yval(i+1)) + u(i) + u(i+1)
  END DO

END SUBROUTINE spline_coefs

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get value from basis function i at x.
FUNCTION get_value(x, base_i, op_code) RESULT(y)
  IMPLICIT NONE
  REAL(dbl), INTENT(IN) :: x
  INTEGER, INTENT(IN) :: base_i
  INTEGER, INTENT(OUT) :: op_code
  REAL(dbl) :: y
  
  INTEGER :: i, j
  REAL(dbl) :: t
  
  IF(x < base(base_i)%low_limit) THEN
    last_error%code = ECODE_BELOW_RANGE
    WRITE(last_error%msg, *) "Value ", x, " is below the lower limit ", base(base_i)%low_limit, " of the base function."
    op_code = ECODE_BELOW_RANGE
    RETURN
  ELSE IF(x > base(base_i)%high_limit) THEN
    last_error%code = ECODE_ABOVE_RANGE
    WRITE(last_error%msg, *) "Value ", x, " is above the upper limit ", base(base_i)%high_limit, " of the base function."
    op_code = ECODE_ABOVE_RANGE
    RETURN
  END IF

  ! find right interval
  DO i=1,base(base_i)%no_of_funcs
    IF(base(base_i)%funcs(i)%low_limit <= x .AND. x <= base(base_i)%funcs(i)%high_limit) EXIT
  END DO

  SELECT CASE(base(base_i)%funcs(i)%func_type)
  CASE(func_type_constant)
    y = base(base_i)%funcs(i)%constant
  CASE(func_type_linear)
    y = base(base_i)%funcs(i)%constant + x*base(base_i)%funcs(i)%coef
  CASE(func_type_spline)
    ! Find right interval
    DO j=1,base(base_i)%funcs(i)%spline%n_knots-1
      IF(base(base_i)%funcs(i)%spline%knots(j) <= x .AND. x <= base(base_i)%funcs(i)%spline%knots(j+1)) EXIT
    END DO
    t = (x - base(base_i)%funcs(i)%spline%knots(j)) / &
      (base(base_i)%funcs(i)%spline%knots(j+1)-base(base_i)%funcs(i)%spline%knots(j))
    y = base(base_i)%funcs(i)%spline%coef(j,1) + base(base_i)%funcs(i)%spline%coef(j,2)*t + &
      base(base_i)%funcs(i)%spline%coef(j,3)*t**2 + base(base_i)%funcs(i)%spline%coef(j,4)*t**3
  END SELECT
  
  op_code = 0

END FUNCTION get_value

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Intianlize random number generator from system clock.
SUBROUTINE init_rng()
  IMPLICIT NONE
  INTEGER :: ssize
  INTEGER, DIMENSION(8) :: t
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL(dbl) :: x

  CALL RANDOM_SEED(SIZE=ssize)
  ALLOCATE(seed(ssize))
  CALL DATE_AND_TIME(VALUES=t)
  seed = 100*t(7) + t(8)/10
  CALL RANDOM_SEED(PUT=seed)
  CALL RANDOM_NUMBER(x) ! First number is not good

END SUBROUTINE init_rng

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Draw 3-dimensional normal random numbers.
SUBROUTINE normal3_random(as, covmat, sim_as, OP_STAT)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(3), INTENT(IN) :: as
  REAL(dbl), DIMENSION(3,3), INTENT(IN) :: covmat
  REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: sim_as
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  
  INTEGER :: n, tn, i
  REAL(dbl) :: u, v, w, ch11, ch12, ch13, ch22, ch23, ch33
  LOGICAL :: even, fallback
  
  n = SIZE(sim_as,1)
  
  ! Independent numbers
  even = (MOD(n,2) == 0)
  IF(even) THEN
    tn = n
  ELSE
    tn = n-1
  END IF
  DO i=1,tn
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(v)
    u = SQRT(-2.0_dbl*LOG(u))
    v = 2*pi*v
    sim_as(i,1) = u*COS(v)
    sim_as(i+1,1) = u*SIN(v)
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(v)
    u = SQRT(-2.0_dbl*LOG(u))
    v = 2*pi*v
    sim_as(i,2) = u*COS(v)
    sim_as(i+1,2) = u*SIN(v)
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(v)
    u = SQRT(-2.0_dbl*LOG(u))
    v = 2*pi*v
    sim_as(i,3) = u*COS(v)
    sim_as(i+1,3) = u*SIN(v)
  END DO
  IF(.NOT. even) THEN
    i = i+1
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(v)
    u = SQRT(-2.0_dbl*LOG(u))
    v = 2*pi*v
    sim_as(i,1) = u*COS(v)
    sim_as(i,2) = u*SIN(v)
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(v)
    u = SQRT(-2.0_dbl*LOG(u))
    v = 2*pi*v
    sim_as(i,3) = u*COS(v)
  END IF
  
  ! Cholesky decomposition of covariance matrix
  fallback = .TRUE.
  ch11 = SQRT(covmat(1,1))
  ch12 = covmat(2,1)/SQRT(covmat(1,1))
  ch13 = covmat(3,1)/SQRT(covmat(1,1))
  u = -(covmat(2,1)**2/covmat(1,1)) + covmat(2,2)
  ! Make sure that matrix is not too singular
  IF(u > 0.0_dbl) THEN
    ch22 = SQRT(u)
    u = covmat(1,1)*(-covmat(2,1)**2+covmat(1,1)*covmat(2,2))
    IF(u > 0.0_dbl) THEN
      ch23 = (-(covmat(2,1)*covmat(3,1)) + covmat(1,1)*covmat(3,2))/SQRT(u)
      u = ABS(covmat(2,1)**2 - covmat(1,1)*covmat(2,2))*covmat(1,1)
      IF(u /= 0.0_dbl) THEN
        v = -(covmat(3,1)**2/covmat(1,1)) - (covmat(2,1)*covmat(3,1) - covmat(1,1)*covmat(3,2))**2 / u + covmat(3,3)
        IF(v > 0.0_dbl) THEN
          ch33 = SQRT(v)
          fallback = .FALSE.
        END IF
      END IF
    END IF
  END IF
  
  IF(fallback) THEN
    ! Only independent random numbers
    u = SQRT(covmat(1,1))
    v = SQRT(covmat(2,2))
    w = SQRT(covmat(3,3))
    DO i=1,n
      sim_as(i,1) = as(1) + u*sim_as(i,1)
      sim_as(i,2) = as(2) + v*sim_as(i,2)
      sim_as(i,3) = as(3) + w*sim_as(i,3)
    END DO
    last_error%code = ECODE_CHOLESKY_FAIL
    last_error%msg = ESTR_CHOLESKY_FAIL
    IF(PRESENT(OP_STAT)) OP_STAT = ECODE_CHOLESKY_FAIL
    RETURN
  END IF
  ! Independent * cholesky = dependent
  DO i=1,n
    u = sim_as(i,1)
    v = sim_as(i,2)
    w = sim_as(i,3)
    sim_as(i,1) = as(1) + ch11*u
    sim_as(i,2) = as(2) + ch12*u+ch22*v
    sim_as(i,3) = as(3) + ch13*u+ch23*v+ch33*w
  END DO
  
  IF(PRESENT(OP_STAT)) OP_STAT = 0

END SUBROUTINE normal3_random

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Draw 2-dimensional normal random numbers.
SUBROUTINE normal2_random(as, covmat, sim_as, OP_STAT)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(2), INTENT(IN) :: as
  REAL(dbl), DIMENSION(2,2), INTENT(IN) :: covmat
  REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: sim_as
  INTEGER, INTENT(OUT), OPTIONAL :: OP_STAT
  
  INTEGER :: n, i
  REAL(dbl) :: u, v, ch11, ch12, ch22
  LOGICAL :: even, fallback
  
  n = SIZE(sim_as,1)
  
  ! Independent numbers
  DO i=1,n
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(v)
    u = SQRT(-2.0_dbl*LOG(u))
    v = 2*pi*v
    sim_as(i,1) = u*COS(v)
    sim_as(i,2) = u*SIN(v)
  END DO
  
  ! Cholesky decomposition of covariance matrix
  fallback = .TRUE.
  ch11 = SQRT(covmat(1,1))
  ch12 = covmat(1,2)/SQRT(covmat(1,1))
  u = -(covmat(1,2)**2/covmat(1,1)) + covmat(2,2)
  ! Make sure that matrix is not too singular
  IF(u > 0.0_dbl) THEN
    ch22 = SQRT(u)
    fallback = .FALSE.
  END IF
  
  IF(fallback) THEN
    ! Only independent random numbers
    u = SQRT(covmat(1,1))
    v = SQRT(covmat(2,2))
    DO i=1,n
      sim_as(i,1) = as(1) + u*sim_as(i,1)
      sim_as(i,2) = as(2) + v*sim_as(i,2)
    END DO
    last_error%code = ECODE_CHOLESKY_FAIL
    last_error%msg = ESTR_CHOLESKY_FAIL
    IF(PRESENT(OP_STAT)) OP_STAT = ECODE_CHOLESKY_FAIL
    RETURN
  END IF
  ! Independent * cholesky = dependent
  DO i=1,n
    u = sim_as(i,1)
    v = sim_as(i,2)
    sim_as(i,1) = as(1) + ch11*u
    sim_as(i,2) = as(2) + ch12*u+ch22*v
  END DO
  
  IF(PRESENT(OP_STAT)) OP_STAT = 0

END SUBROUTINE normal2_random

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recursive, in-place quicksort
RECURSIVE SUBROUTINE Qsort(vec, left, right)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: n, pivot

  n = right-left+1
  IF(MOD(n,2) == 0) THEN
    pivot = left-1 + n/2
  ELSE
    pivot = left-1 + (n+1)/2
  END IF

  pivot = partition(vec, pivot, left, right)

  IF(left < pivot-1) CALL Qsort(vec, left, pivot-1)
  IF(right > pivot+1) CALL Qsort(vec, pivot+1, right)

END SUBROUTINE Qsort

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Partition of real array according to pivot value
FUNCTION partition(vec, pivot, left, right) RESULT(new_pivot)
  IMPLICIT NONE
  REAL(dbl), DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: new_pivot

  INTEGER :: pivot, i
  REAL(dbl) :: rpv_value, rtemp

  rpv_value = vec(pivot)
  vec(pivot) = vec(right)
  vec(right) = rpv_value

  new_pivot = left
  DO i=left,right-1
    IF(vec(i) < rpv_value) THEN
      rtemp = vec(i)
      vec(i) = vec(new_pivot)
      vec(new_pivot) = rtemp
      new_pivot = new_pivot+1
    END IF
  END DO
  rtemp = vec(new_pivot)
  vec(new_pivot) = vec(right)
  vec(right) = rtemp

END FUNCTION partition

  
END MODULE HG1G2tools
