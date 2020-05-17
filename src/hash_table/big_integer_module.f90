MODULE big_integer_module

  !  Copyright (c) 1993-2002 Unicomp, Inc.
  !
  !  Developed at Unicomp, Inc.
  !
  !  Permission to use, copy, modify, and distribute this
  !  software is freely granted, provided that this notice
  !  is preserved.

  !  The module named BIG_INTEGERS defines a new data type
  !  BIG_INTEGER.  This data type represents nonnegative integers
  !  up to 10**n - 1, where n is the parameter (named constant)
  !  NR_OF_DECIMAL_DIGITS.  This value may be changed, but the
  !  module must then be recompiled (it is not dynamic).

  !  The following operations are implemented.
  !  b represents a big integer, c a character string,
  !  and i an ordinary integer.

  !  integer, parameter :: nr_of_decimal_digits

  !  big (i)
  !  big (c)
  !  int (b)
  !  int (c)
  !  char (b)
  !  char (i)

  !  b = i
  !  b = c
  !  i = b
  !  i = c
  !  c = b
  !  c = i

  !  b ? i, i ? b, and b ? b, where ? is
  !    +, -, *, /,
  !    <, <=, >, >=, ==, or /=

  !  b ** i

  !  modulo (b, i)  [result is integer]
  !  modulo (i, b)  [result is big_integer]
  !  modulo (b, b)  [result is big_integer]

  !  huge (b)
  !  sqrt (b)

  !  call print_big (b)
  !  call random_number (b, low, high)

  !  Many operations of the form b ? i, where i < base,
  !  are implemented to be efficient as special cases.

  IMPLICIT NONE

  PUBLIC :: big
  PUBLIC :: int
  PUBLIC :: char
  PUBLIC :: print_big
  PUBLIC :: ASSIGNMENT (=)
  PUBLIC :: OPERATOR (+)
  PUBLIC :: OPERATOR (-)
  PUBLIC :: OPERATOR (*)
  PUBLIC :: OPERATOR (/)
  PUBLIC :: OPERATOR (**)
  PUBLIC :: modulo
  PUBLIC :: huge
  PUBLIC :: sqrt
  PUBLIC :: random_number
  PUBLIC :: OPERATOR (==)
  PUBLIC :: OPERATOR (/=)
  PUBLIC :: OPERATOR (<=)
  PUBLIC :: OPERATOR (<)
  PUBLIC :: OPERATOR (>=)
  PUBLIC :: OPERATOR (>)

  PRIVATE :: big_gets_int, &
       big_int, &
       big_gets_char, &
       big_char, &
       int_gets_big, &
       int_big, &
       int_gets_char, &
       int_char, &
       char_gets_big, &
       char_big, &
       char_gets_int, &
       char_int

  PRIVATE :: big_plus_int, &
       int_plus_big, &
       big_plus_big, &
       big_minus_int, &
       int_minus_big, &
       big_minus_big, &
       big_times_int, &
       int_times_big, &
       big_times_big, &
       big_div_int, &
       int_div_big, &
       big_div_big, &
       big_power_int, &
       modulo_big_int, &
       modulo_int_big, &
       modulo_big_big

  PRIVATE :: big_eq_int, &
       int_eq_big, &
       big_eq_big, &
       big_ne_int, &
       int_ne_big, &
       big_ne_big, &
       big_le_int, &
       int_le_big, &
       big_le_big, &
       big_ge_int, &
       int_ge_big, &
       big_ge_big, &
       big_lt_int, &
       int_lt_big, &
       big_lt_big, &
       big_gt_int, &
       int_gt_big, &
       big_gt_big

  PRIVATE :: huge_big, &
       big_base_to_power, &
       print_big_base, &
       sqrt_big, &
       msd, &
       random_number_big

  INTRINSIC :: char
  INTRINSIC :: int
  INTRINSIC :: modulo
  INTRINSIC :: huge
  INTRINSIC :: sqrt
  INTRINSIC :: random_number
  INTRINSIC :: radix
  INTRINSIC :: digits

  ! This indicates the maximum number of decimal digits
  ! that a big integer may contain.

  INTEGER, PARAMETER, PUBLIC :: nr_of_decimal_digits = 150

  ! If the radix (returned by "radix(0)" of the integers on
  ! your system is not 2 change the following constant to
  ! the logarithm in the base 10 of the radix: log10(radix)

  REAL, PARAMETER, PRIVATE :: log_base_10_of_radix = 0.30103

  INTEGER, PARAMETER, PRIVATE :: &
       d = DIGITS (0) / 2, &
       r = RADIX (0), &
       base = r ** d, &
       nr_of_digits = nr_of_decimal_digits / (log_base_10_of_radix * d) + 1

  ! The base of the number system is r ** d,
  ! so that each "digit" is 0 to r**d - 1

  TYPE, PUBLIC :: big_integer
     PRIVATE
     INTEGER, DIMENSION (0 : nr_of_digits) :: digit
  END TYPE big_integer

  INTERFACE big
     MODULE PROCEDURE big_char, &
          big_int
  END INTERFACE

  INTERFACE int
     MODULE PROCEDURE int_char, &
          int_big
  END INTERFACE

  INTERFACE char
     MODULE PROCEDURE char_big, &
          char_int
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE big_gets_int, &
          big_gets_char, &
          int_gets_big, &
          int_gets_char, &
          char_gets_big, &
          char_gets_int
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE big_plus_int, &
          big_plus_big
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE big_minus_int, &
          int_minus_big, &
          big_minus_big
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE big_times_int, &
          int_times_big, &
          big_times_big
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE big_div_int, &
          int_div_big, &
          big_div_big
  END INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE big_power_int
  END INTERFACE

  INTERFACE modulo
     MODULE PROCEDURE modulo_big_int, &
          modulo_int_big, &
          modulo_big_big
  END INTERFACE

  INTERFACE OPERATOR (==)
     MODULE PROCEDURE big_eq_int, &
          int_eq_big, &
          big_eq_big
  END INTERFACE

  INTERFACE OPERATOR (/=)
     MODULE PROCEDURE big_ne_int, &
          int_ne_big, &
          big_ne_big
  END INTERFACE

  INTERFACE OPERATOR (<=)
     MODULE PROCEDURE big_le_int, &
          int_le_big, &
          big_le_big
  END INTERFACE

  INTERFACE OPERATOR (>=)
     MODULE PROCEDURE big_ge_int, &
          int_ge_big, &
          big_ge_big
  END INTERFACE

  INTERFACE OPERATOR (<)
     MODULE PROCEDURE big_lt_int, &
          int_lt_big, &
          big_lt_big
  END INTERFACE

  INTERFACE OPERATOR (>)
     MODULE PROCEDURE big_gt_int, &
          int_gt_big, &
          big_gt_big
  END INTERFACE

  INTERFACE huge
     MODULE PROCEDURE huge_big
  END INTERFACE

  INTERFACE sqrt
     MODULE PROCEDURE sqrt_big
  END INTERFACE

  INTERFACE random_number
     MODULE PROCEDURE random_number_big
  END INTERFACE

CONTAINS

  PURE FUNCTION big_char (c) RESULT (b)

    CHARACTER (len=*), INTENT (in) :: c
    TYPE (big_integer) :: b
    INTEGER :: temp_digit, n

    IF (LEN (c) > nr_of_decimal_digits) THEN
       b = HUGE (b)
       RETURN
    END IF
    b % digit = 0
    DO n = 1, LEN (c)
       temp_digit = INDEX ("0123456789", c (n:n)) - 1
       IF (temp_digit < 0) THEN
          b = HUGE (b)
       END IF
       b = b * 10 + temp_digit
    END DO

  END FUNCTION big_char

  PURE SUBROUTINE big_gets_char (b, c)

    TYPE (big_integer), INTENT (out) :: b
    CHARACTER (len=*), INTENT (in) :: c

    b = big_char (TRIM (c))

  END SUBROUTINE big_gets_char

  PURE FUNCTION big_int (i) RESULT (b)

    INTEGER, INTENT (in) :: i
    TYPE (big_integer) :: b
    INTEGER :: temp_i, n

    IF (i < 0) THEN
       b = HUGE (b)
    END IF

    b % digit = 0
    temp_i = i
    DO n = 0, nr_of_digits - 1
       IF (temp_i == 0) THEN
          RETURN
       END IF
       b % digit (n) = MODULO (temp_i, base)
       temp_i = temp_i / base
    END DO

    IF (temp_i /= 0) THEN
       b = HUGE (b)
    END IF

  END FUNCTION big_int

  PURE SUBROUTINE big_gets_int (b, i)

    TYPE (big_integer), INTENT (out) :: b
    INTEGER, INTENT (in) :: i

    b = big (i)

  END SUBROUTINE big_gets_int

  PURE FUNCTION int_char (c) RESULT (i)

    CHARACTER (len=*), INTENT (in) :: c
    INTEGER :: i

    i = INT (big (TRIM (c)))

  END FUNCTION int_char

  PURE SUBROUTINE int_gets_char (i, c)

    INTEGER, INTENT (out) :: i
    CHARACTER (len=*), INTENT (in) :: c

    i = int_char (TRIM (c))

  END SUBROUTINE int_gets_char

  PURE FUNCTION int_big (b) RESULT (i)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER :: i

    IF (msd (b) > 1) THEN
       i = HUGE (i)
    ELSE
       i = base * b % digit (1) + b % digit (0)
    END IF

  END FUNCTION int_big

  PURE SUBROUTINE int_gets_big (i, b)

    INTEGER, INTENT (out) :: i
    TYPE (big_integer), INTENT (in) :: b

    i = INT (b)

  END SUBROUTINE int_gets_big

  PURE FUNCTION char_big (b) RESULT (c)

    TYPE (big_integer), INTENT (in) :: b
    CHARACTER (len=nr_of_decimal_digits+9) :: c
    TYPE (big_integer) :: temp_big
    INTEGER :: n, remainder
    CHARACTER (len = *), PARAMETER :: digit_chars = "0123456789"

    temp_big = b
    c = REPEAT (" ", LEN(c)-1) // "0"
    DO n = LEN (c), 1, -1
       IF (temp_big == 0) THEN
          EXIT
       END IF
       remainder = MODULO (temp_big, 10) + 1
       temp_big = temp_big / 10
       c (n:n) = digit_chars (remainder:remainder)
    END DO

    c = ADJUSTL (c)

  END FUNCTION char_big

  PURE SUBROUTINE char_gets_big (c, b)

    TYPE (big_integer), INTENT (in) :: b
    CHARACTER (len=*), INTENT (out) :: c

    c = CHAR (b)

  END SUBROUTINE char_gets_big

  PURE FUNCTION char_int (i) RESULT (c)

    INTEGER, INTENT (in) :: i
    CHARACTER (len=nr_of_decimal_digits+9) :: c

    c = big (i)

  END FUNCTION char_int

  PURE SUBROUTINE char_gets_int (c, i)

    INTEGER, INTENT (in) :: i
    CHARACTER (len=*), INTENT (out) :: c

    c = big (i)

  END SUBROUTINE char_gets_int

  PURE FUNCTION msd (x) RESULT (msd_result)

    ! Find most significan digit of x

    TYPE (big_integer), INTENT (in) :: x
    INTEGER :: msd_result
    INTEGER :: n

    DO n = nr_of_digits, 0, -1
       IF (x % digit (n) /= 0) THEN
          msd_result = n
          RETURN
       END IF
    END DO

    msd_result = -1

  END FUNCTION msd

  PURE FUNCTION big_plus_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    TYPE (big_integer) :: bi
    INTEGER :: n, summ, carry

    IF (i < base) THEN
       carry = i
       DO n = 0, nr_of_digits - 1
          summ = b % digit (n) + carry
          bi % digit (n) = MODULO (summ, base)
          carry = summ / base
          IF (carry == 0) THEN
             bi % digit (n+1:) = b % digit (n+1:)
             RETURN
          END IF
       END DO
       IF (n==nr_of_digits) THEN
          bi = HUGE (bi)
       END IF
    ELSE
       bi = b + big (i)
    END IF

  END FUNCTION big_plus_int

  PURE FUNCTION int_plus_big (i, b) RESULT (bi)

    INTEGER, INTENT (in) :: i
    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: bi

    bi = b + i

  END FUNCTION int_plus_big

  PURE FUNCTION big_plus_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    TYPE (big_integer) :: bb
    INTEGER :: carry, temp_digit, n, m

    carry = 0
    m = MAX (msd (x), msd (y))
    DO n = 0, m
       temp_digit = &
            x % digit (n) + y % digit (n) + carry
       bb % digit (n) = MODULO (temp_digit, base)
       carry = temp_digit / base
    END DO

    bb % digit (m+1) = carry
    bb % digit (m+2:nr_of_digits) = 0
    IF (bb % digit (nr_of_digits) /= 0) THEN
       bb = HUGE (bb)
    END IF

  END FUNCTION big_plus_big

  PURE FUNCTION big_minus_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    TYPE (big_integer) :: bi
    INTEGER :: n, borrow, diff, msdb

    bi % digit = 0
    msdb = msd (b)
    IF (msdb<1 .AND. b % digit (0) < i) THEN
       RETURN
    END IF

    IF (i < base) THEN
       borrow = i
       DO n = 0, nr_of_digits - 1
          diff = b % digit (n) - borrow
          bi % digit (n) = MODULO (diff, base)
          borrow = (base - diff) / base
          IF (borrow == 0) THEN
             bi % digit (n+1:msdb) = b % digit (n+1:msdb)
             RETURN
          END IF
       END DO
    ELSE
       bi = b - big (i)
    END IF

  END FUNCTION big_minus_int

  PURE FUNCTION int_minus_big (i, b) RESULT (ib)

    INTEGER, INTENT (in) :: i
    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: ib

    ib = big (i) - b

  END FUNCTION int_minus_big

  PURE FUNCTION big_minus_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    TYPE (big_integer) :: bb
    TYPE (big_integer) :: temp_big
    INTEGER :: n

    temp_big = x
    DO n = 0, nr_of_digits - 1
       bb % digit (n) = temp_big % digit (n) - y % digit (n)
       IF (bb % digit (n) < 0) THEN
          bb % digit (n) = bb % digit (n) + base
          temp_big % digit (n + 1) = temp_big % digit (n + 1) - 1
       END IF
    END DO

    IF (temp_big % digit (nr_of_digits) < 0) THEN
       bb % digit = 0
    ELSE
       bb % digit (nr_of_digits) = 0
    END IF

  END FUNCTION big_minus_big

  PURE FUNCTION big_times_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    TYPE (big_integer) :: bi
    INTEGER :: ib, prod, carry

    IF (i < base) THEN
       bi % digit = 0
       carry = 0
       DO ib = 0, msd (b)
          prod = b % digit (ib) * i + carry
          bi % digit (ib) = MODULO (prod, base)
          carry = prod / base
       END DO
       IF (ib==nr_of_digits .AND. carry /= 0) THEN
          bi = HUGE (bi)
       ELSE
          bi % digit (ib) = carry
       END IF
    ELSE
       bi = b * big (i)
    END IF

  END FUNCTION big_times_int

  PURE FUNCTION int_times_big (i, b) RESULT (bi)

    INTEGER, INTENT (in) :: i
    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: bi

    bi = b * i

  END FUNCTION int_times_big

  PURE FUNCTION big_times_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    TYPE (big_integer) :: bb

    INTEGER :: ix, iy, ib, carry, prod

    bb % digit = 0

    DO ix = 0, msd (x)
       carry = 0
       ib = ix
       DO iy = 0, msd (y)
          prod = x % digit (ix) * y % digit (iy) + bb % digit (ib) + carry
          carry = prod / base
          bb % digit (ib) = MODULO (prod, base)
          IF (ib == nr_of_digits) THEN
             bb = HUGE (bb)
             RETURN
          END IF
          ib = ib + 1
       END DO
       bb % digit (ib) = bb % digit (ib) + carry
    END DO

  END FUNCTION big_times_big

  PURE FUNCTION big_base_to_power (n)  RESULT (b)

    INTEGER, INTENT (in) :: n
    TYPE (big_integer) :: b

    IF (n < 0) THEN
       b = 0
    ELSE IF (n >= nr_of_digits) THEN
       b = HUGE (b)
    ELSE
       b % digit = 0
       b % digit (n) = 1
    END IF

  END FUNCTION big_base_to_power

  PURE FUNCTION big_div_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    TYPE (big_integer) :: bi
    INTEGER :: n, temp_int, remainder

    IF (i == 0) THEN
       bi = HUGE (bi)
    ELSE IF (i < base) THEN
       bi % digit = 0
       remainder = 0
       DO n = msd(b), 0, -1
          temp_int = base * remainder + b % digit (n)
          bi % digit (n) = temp_int / i
          remainder = MODULO (temp_int, i)
       END DO
    ELSE
       bi = b / big (i)
    END IF

  END FUNCTION big_div_int

  PURE FUNCTION int_div_big (i, b) RESULT (ib)

    INTEGER, INTENT (in) :: i
    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: ib

    ib = big (i) / b

  END FUNCTION int_div_big

  PURE FUNCTION big_div_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    TYPE (big_integer) :: bb

    TYPE (big_integer) :: tx, ty

    INTEGER :: msdx, msdy, ix, iy
    INTEGER :: v1, v2, u0, u1, u2
    INTEGER :: dd, bi, car, bar, prd

    IF (y == 0) THEN
       bb = HUGE (bb)
       RETURN
    END IF

    msdx = msd(x)
    msdy = msd(y)

    IF (msdy == 0) THEN
       bb = x / y % digit (0)
       RETURN
    END IF

    bb % digit = 0

    IF (msdy < msdy) THEN
       RETURN
    END IF

    tx = x
    ty = y

    car = 0
    bar = 0
    prd = 0
    dd = base / (ty % digit (msdy) + 1)
    IF (dd /= 1) THEN
       DO ix = 0, msdx
          tx % digit (ix) = tx % digit (ix) * dd + car
          car = tx % digit (ix) / base
          tx % digit (ix) = tx % digit (ix) - base * car
       END DO
       tx % digit (msdx+1) = car
       car = 0
       DO iy = 0, msdy
          ty % digit (iy) = ty % digit (iy) * dd + car
          car = ty % digit (iy) / base
          ty % digit (iy) = ty % digit (iy) - base * car
       END DO
    END IF

    msdx = msdx + 1

    v1 = ty % digit (msdy)
    v2 = ty % digit (msdy-1)
    bb % digit = 0

    DO msdx = msdx, msdy + 1, -1

       u0 = tx % digit (msdx)
       u1 = tx % digit (msdx-1)
       u2 = tx % digit (msdx-2)

       IF (u0 == v1) THEN
          bi = base - 1
       ELSE
          bi = (u0*base + u1) / v1
       END IF

       DO
          IF (v2*bi <= (u0*base + u1 - bi*v1) * base + u2) THEN
             EXIT
          END IF
          bi = bi - 1
       END DO

       IF (bi > 0) THEN
          car = 0
          bar = 0
          ix = msdx - msdy - 1
          DO iy = 0, msdy
             prd = bi * ty % digit (iy) + car
             car = prd / base
             prd = prd - base * car
             tx % digit (ix) = tx % digit (ix) - (prd + bar)
             IF (tx % digit (ix) < 0) THEN
                bar = 1
                tx % digit (ix) = tx % digit (ix) + base
             ELSE
                bar = 0
             END IF
             ix = ix + 1
          END DO
          IF (tx % digit (msdx) < car + bar) THEN
             car = 0
             bi = bi -1
             ix = msdx - msdy - 1
             DO iy = 0, msdy
                tx % digit (ix) = tx % digit (ix) + ty % digit (iy) + car
                IF (tx % digit (ix) > base) THEN
                   car = 1
                   tx % digit (ix) = tx % digit (ix) - base
                ELSE
                   car = 0
                END IF
                ix = ix + 1
             END DO
          END IF
       END IF
       tx % digit (msdx) = 0
       bb % digit (1:nr_of_digits) = bb % digit (0:nr_of_digits-1)
       bb % digit (0) = bi
    END DO

  END FUNCTION big_div_big

  PURE FUNCTION modulo_big_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    INTEGER :: bi
    INTEGER :: n

    IF (i == 0) THEN
       bi = HUGE (bi)
    ELSE IF (i < base) THEN
       bi = 0
       DO n = msd(b), 0, -1
          bi = MODULO (base * bi + b % digit (n), i)
       END DO
    ELSE
       bi = MODULO (b, big (i))
    END IF

  END FUNCTION modulo_big_int

  PURE FUNCTION modulo_int_big (ii, b) RESULT (ib)

    INTEGER, INTENT (in) :: ii
    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: ib

    ib = MODULO (big (ii), b)

  END FUNCTION modulo_int_big

  PURE FUNCTION modulo_big_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    TYPE (big_integer) :: bb

    bb = x - x / y * y

  END FUNCTION modulo_big_big

  PURE FUNCTION big_eq_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) == i

  END FUNCTION big_eq_int

  PURE FUNCTION int_eq_big (i, b) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) == i

  END FUNCTION int_eq_big

  PURE FUNCTION big_eq_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    LOGICAL :: bb

    bb = ALL (x % digit == y % digit)

  END FUNCTION big_eq_big

  PURE FUNCTION big_ne_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) /= i

  END FUNCTION big_ne_int

  PURE FUNCTION int_ne_big (i, b) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) /= i

  END FUNCTION int_ne_big

  PURE FUNCTION big_ne_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    LOGICAL :: bb

    bb = ANY (x % digit /= y % digit)

  END FUNCTION big_ne_big

  PURE FUNCTION big_le_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) <= i

  END FUNCTION big_le_int

  PURE FUNCTION int_le_big (i, b) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = i <= INT (b)

  END FUNCTION int_le_big

  PURE FUNCTION big_le_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    LOGICAL :: bb
    INTEGER :: n

    bb = .TRUE.
    DO n = nr_of_digits, 0, -1
       IF (x % digit (n) /= y % digit (n)) THEN
          bb = (x % digit (n) < y % digit (n))
          EXIT
       END IF
    END DO

  END FUNCTION big_le_big

  PURE FUNCTION big_gt_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) > i

  END FUNCTION big_gt_int

  PURE FUNCTION int_gt_big (i, b) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = i > INT (b)

  END FUNCTION int_gt_big

  PURE FUNCTION big_gt_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    LOGICAL :: bb
    INTEGER :: n

    bb = .TRUE.
    DO n = nr_of_digits, 0, -1
       IF (x % digit (n) /= y % digit (n)) THEN
          bb = (x % digit (n) < y % digit (n))
          EXIT
       END IF
    END DO

    bb = .NOT. bb

  END FUNCTION big_gt_big

  PURE FUNCTION big_lt_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) < i

  END FUNCTION big_lt_int

  PURE FUNCTION int_lt_big (i, b) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = i < INT (b)

  END FUNCTION int_lt_big

  PURE FUNCTION big_lt_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    LOGICAL :: bb
    INTEGER :: n

    bb = .FALSE.
    DO n = nr_of_digits, 0, -1
       IF (x % digit (n) /= y % digit (n)) THEN
          bb = (x % digit (n) < y % digit (n))
          EXIT
       END IF
    END DO

  END FUNCTION big_lt_big

  PURE FUNCTION big_ge_int (b, i) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = INT (b) >= i

  END FUNCTION big_ge_int

  PURE FUNCTION int_ge_big (i, b) RESULT (bi)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    LOGICAL :: bi

    bi = i >= INT (b)

  END FUNCTION int_ge_big

  PURE FUNCTION big_ge_big (x, y) RESULT (bb)

    TYPE (big_integer), INTENT (in) :: x, y
    LOGICAL :: bb
    INTEGER :: n

    bb = .FALSE.
    DO n = nr_of_digits, 0, -1
       IF (x % digit (n) /= y % digit (n)) THEN
          bb = (x % digit (n) < y % digit (n))
          EXIT
       END IF
    END DO

    bb = .NOT. bb

  END FUNCTION big_ge_big

  PURE FUNCTION huge_big (b) RESULT (hb)

    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: hb

    hb % digit (0) = b % digit (0) ! to avoid diagnostic
    hb % digit = base - 1
    hb % digit (nr_of_digits) = 0

  END FUNCTION huge_big

  PURE FUNCTION sqrt_big (b) RESULT (sb)

    TYPE (big_integer), INTENT (in) :: b
    TYPE (big_integer) :: sb
    TYPE (big_integer) :: old_sqrt_big, new_sqrt_big
    INTEGER :: i, n

    n = -1
    DO i = nr_of_digits, 0, -1
       IF (b % digit (i) /= 0) THEN
          n = i
          EXIT
       END IF
    END DO

    IF (n == -1) THEN
       sb = 0
    ELSE IF (n == 0) THEN
       sb = INT (SQRT (REAL (b % digit (0))))
    ELSE
       old_sqrt_big = 0
       IF (MODULO (n, 2) == 0) THEN
          old_sqrt_big % digit (n / 2) = INT (SQRT (REAL (b % digit (n))))
       ELSE
          old_sqrt_big % digit ((n - 1) / 2) =  &
               INT (SQRT (REAL (base * b % digit (n) + b % digit (n-1))))
       END IF

       DO
          new_sqrt_big = (old_sqrt_big + b / old_sqrt_big) / 2
          IF (new_sqrt_big == old_sqrt_big .OR.  &
               new_sqrt_big == old_sqrt_big + 1 .OR.  &
               new_sqrt_big == 0) THEN
             EXIT
          ELSE
             old_sqrt_big = new_sqrt_big
          END IF
       END DO
       sb = old_sqrt_big
    END IF

  END FUNCTION sqrt_big

  RECURSIVE FUNCTION big_power_int (b, i)  &
       RESULT (big_power_int_result)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER, INTENT (in) :: i
    TYPE (big_integer) :: big_power_int_result
    TYPE (big_integer) :: temp_big

    IF (i <= 0) THEN
       big_power_int_result = 1
    ELSE
       temp_big = big_power_int (b, i / 2)
       IF (MODULO (i, 2) == 0) THEN
          big_power_int_result = temp_big * temp_big
       ELSE
          big_power_int_result = temp_big * temp_big * b
       END IF
    END IF

  END FUNCTION big_power_int

  SUBROUTINE print_big (b)

    TYPE (big_integer), INTENT (in) :: b

    WRITE (unit = *, fmt = "(a)", advance = "no") TRIM (CHAR (b))

  END SUBROUTINE print_big

  SUBROUTINE print_big_base (b)

    TYPE (big_integer), INTENT (in) :: b
    INTEGER :: n

    PRINT *, "base: ", base
    DO n = nr_of_digits, 1, -1
       IF (b % digit (n) /= 0) THEN
          EXIT
       END IF
    END DO
    PRINT "(10i9)", b % digit (n:0:-1)

  END SUBROUTINE print_big_base

  SUBROUTINE random_number_big (r, low, high)

    !  Generate by linear congruence x' = ax + c mod m
    !  where m is huge (b) + 1 = base ** nr_of_digits

    TYPE (big_integer), INTENT (out) :: r
    TYPE (big_integer), INTENT (in) :: low, high
    INTEGER :: n, i, carry, prod, summ
    TYPE (big_integer), SAVE :: x = big_integer ( (/ (1, i=0,nr_of_digits-1), 0 /) )
    TYPE (big_integer), PARAMETER :: h = big_integer ( (/ (base-1, i=0,nr_of_digits-1), 0 /) )
    INTEGER, PARAMETER :: a = 16907, c = 8191

    !  Multiply by a
    carry = 0
    DO n = 0, nr_of_digits - 1
       prod = x % digit (n) * a + carry
       x % digit (n) = MODULO (prod, base)
       carry = prod / base
    END DO

    !  Add c
    carry = c
    DO n = 0, nr_of_digits - 1
       summ = x % digit (n) + carry
       x % digit (n) = MODULO (summ, base)
       carry = summ / base
    END DO

    r = x / (h / (high -low + 1)) + low

  END SUBROUTINE random_number_big

END MODULE big_integer_module
