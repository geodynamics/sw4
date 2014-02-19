      subroutine fdump ( )

c*********************************************************************72
c
cc FDUMP produces a symbolic dump.
c
c  Discussion:
c
c    This routine is intended to be replaced by a locally written
c    version which produces a symbolic dump.  Failing this,
c    it should be replaced by a version which prints the
c    subprogram nesting list.
c
c    Normally, the dump information should be printed to all the
c    active error output units.  The number and value of these
c    units can be determined by calling XGETUA.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    None
c
      implicit none

      return
      end
      function i1mach ( i )

c*********************************************************************72
c
cc I1MACH returns integer machine dependent constants.
c
c  Discussion:
c
c    This routine can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    subroutine with one (input) argument, and can be called
c    as follows, for example
c
c      K = I1MACH(I)
c
c    where I=1,...,16.  The output value of K above is
c    determined by the input value of I.  The results for
c    various values of I are discussed below.
c
c    I/O unit numbers.
c
c    I1MACH( 1) = the standard input unit.
c    I1MACH( 2) = the standard output unit.
c    I1MACH( 3) = the standard punch unit.
c    I1MACH( 4) = the standard error message unit.
c
c    Words.
c
c    I1MACH( 5) = the number of bits per integer storage unit.
c    I1MACH( 6) = the number of characters per integer storage unit.
c
c    Integers.
c
c    Assume integers are represented in the S-digit, base-A form
c
c      sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
c
c    where 0 .LE. X(I) .LT. A for I=0,...,S-1.
c
c    I1MACH( 7) = A, the base.
c    I1MACH( 8) = S, the number of base-A digits.
c    I1MACH( 9) = A**S - 1, the largest magnitude.
c
c    Floating-Point Numbers.
c
c    Assume floating-point numbers are represented in the T-digit,
c    base-B form
c
c      sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c    where 0 <= X(I) < B for I=1,...,T,
c    0 < X(1), and EMIN <= E <= EMAX.
c
c    I1MACH(10) = B, the base.
c
c    Single-Precision
c
c    I1MACH(11) = T, the number of base-B digits.
c    I1MACH(12) = EMIN, the smallest exponent E.
c    I1MACH(13) = EMAX, the largest exponent E.
c
c    Double-Precision
c
c    I1MACH(14) = T, the number of base-B digits.
c    I1MACH(15) = EMIN, the smallest exponent E.
c    I1MACH(16) = EMAX, the largest exponent E.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, integer I1MACH, the value of the constant.
c
      implicit none

      integer i
      integer i1mach
      integer imach(16)
c
c  Machine constants for IEEE arithmetic machines.
c
      data imach( 1) /    5 /
      data imach( 2) /    6 /
      data imach( 3) /    7 /
      data imach( 4) /    6 /
      data imach( 5) /   32 /
      data imach( 6) /    4 /
      data imach( 7) /    2 /
      data imach( 8) /   31 /
      data imach( 9) / 2147483647 /
      data imach(10) /    2 /
      data imach(11) /   24 /
      data imach(12) / -125 /
      data imach(13) /  128 /
      data imach(14) /   53 /
      data imach(15) / -1021 /
      data imach(16) /  1024 /

      if ( i .lt. 1  .or.  16 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  I out of bounds.'
        stop
      end if

      i1mach = imach(i)

      return
      end
      function j4save ( which, value, set )

c*********************************************************************72
c
cc J4SAVE saves and recalls global variables.
c
c  Discussion:
c
c    The internal parameters are initialized to the following values:
c
c    #1 =  0, NERR, the index of the most recent error;
c    #2 =  0, KONTRL, error control flag (0 means only level 2 errors are fatal,
c             and get a printout, while lower level errors get no printout.)
c    #3 =  0, IUNIT, the main error output unit (0 means use standard output).
c    #4 = 10, MAXMES, the maximum number of times any message is printed.
c    #5 =  1, NUNIT, total number of error output units in use.
c    #6 = -1, second error output unit (-1 means not being used).
c    #7 = -1, third error output unit (-1 means not being used).
c    #8 = -1, fourth error output unit (-1 means not being used).
c    #9 = -1, fifth error output unit (-1 means not being used).
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer WHICH, the index of the desired item.
c    Legal values of WHICH are between 1 and 9.
c    *1, NERR, the index of the most recent error;
c    *2, KONTRL, the current error control flag;
c    *3, IUNIT, the main error output unit.  (0 means use the standard output);
c    *4, MAXMES, the maximum number of times any message is to be printed.
c    *5, NUNIT, the total number of error output units.
c    *6, the second unit for error messages, if being used;
c    *7, the third unit for error messages, if being used;
c    *8, the fourth unit for error messages, if being used;
c    *9, the fifth unit for error messages, if being used.
c
c    Input, integer VALUE, will be used to reset the value of the
c    WHICH-th internal parameter, if SET is TRUE.
c
c    Input, logical SET, is TRUE if the call is being made in order
c    to set the WHICH-th parameter to VALUE.  If SET is FALSE,
c    then no internal assignment is made; presumably the call is
c    being made to inquire about the current internal value of
c    the WHICH-th parameter.
c
c    Output, integer J4SAVE, the value of the WHICH-th parameter.
c    If the input value of SET was TRUE, then the value of the
c    WHICH-th parameter has just been changed.  In that case, the
c    value returned in J4SAVE is the old value of this parameter.
c
      implicit none

      integer j4save
      integer param(9)
      logical set
      integer value
      integer which

      save param

      data param /
     &  0, 2, 0, 10, 1, -1, -1, -1, -1 /

      if ( which .lt. 1 .or. 9 .lt. which ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'J4SAVE - Fatal error!'
        write ( *, '(a,i10)' ) 
     &    '  Illegal input value of WHICH = ', which
        stop
      end if

      j4save = param ( which )

      if ( set ) then
        param(which) = value
      end if

      return
      end
      function numxer ( nerr )

c*********************************************************************72
c
cc NUMXER returns the most recent error number.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Output, integer NERR, the most recent error number.
c
c    Output, integer NUMXER, the most recent error number.
c
      implicit none

      integer j4save
      integer nerr
      integer numxer
      logical set
      integer value
      integer which

      which = 1
      value = 0
      set = .false.

      nerr = j4save ( which, value, set )

      numxer = nerr

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine xerabt ( messg, nmessg )

c*********************************************************************72
c
cc XERABT aborts program execution.
c
c  Discussion:
c
c    This routine is called to abort execution of a running program,
c    indicated by the occurrence of a fatal error.
c
c    The error message associated with the fatal condition is provided
c    in the calling sequence.
c
c    This routine is used when the error message handlers XERROR and
c    XERRWV are employed.  The similar routine XERHLT is to be used 
c    when the more modern error message handler XERMSG is used.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*NMESSG MESSG, the error message.
c
c    Input, integer NMESSG, the number of characters in the error message.
c    If NMESSG <= 0, it is presumed that no error message has been supplied.
c
      implicit none

      character*(*) messg
      integer nmessg

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XERABT - Program abort!'

      if ( 0 < nmessg ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Associated error message:'
        write ( *, '(a)' ) messg(1:nmessg)

      end if

      stop
      end
      subroutine xerbla ( subrou, nerr )

c*********************************************************************72
c
cc XERBLA is an error handler for the Level 2 and Level 3 BLAS routines.
c
c  Discussion:
c
c    This routine is called by Level 2 and 3 BLAS routines if an input 
c    parameter is invalid.
c
c  Modified:
c
c    06 April 2007
c
c  Parameters:
c
c    Input, character*(*) SUBROU, the name of the routine which
c    called XERBLA.  The name will not be more than 6 characters.
c
c    Input, integer NERR, the error number, which here is used to
c    indicate the position of the invalid parameter in the 
c    parameter-list of the calling routine.
c
      implicit none

      integer level
      character*6 librar
      character*60 message
      integer nerr
      character*(*) subrou

      librar = 'SLATEC'
      write ( message, '(a,a,a,i2,a)' ) 'On entry to ', subrou,
     &  ' parameter number ', nerr, ' had an illegal value.'
      level = 1

      call xermsg ( librar, subrou, message, nerr, level )

      return
      end
      subroutine xerclr ( )

c*********************************************************************72
c
cc XERCLR resets the most recent error number to 0.
c
c  Discussion:
c
c    This routine simply resets the most recent error index to 0.
c
c    This may be necessary in order to determine that a certain
c    error has occurred again.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    None
c
      implicit none

      integer j4save
      integer junk
      logical set
      integer value
      integer which

      which = 1
      value = 0
      set = .true.

      junk = j4save ( which, value, set )

      return
      end
      subroutine xercnt ( librar, subrou, messg, nerr, level, kontrl )

c*********************************************************************72
c
cc XERCNT allows user control over the handling of errors.
c
c  Description:
c
c    This routine allows user control over handling of individual errors.
c
c    This routine is to be used when the error message routine XERMSG
c    is employed.  The similar routine XERCTL is to be used for the
c    older error message routines XERROR and XERRWV.
c
c    Just after each message is recorded, but before it is
c    processed any further (i.e., before it is printed or
c    a decision to abort is made), a call is made to XERCNT.
c
c    If the user has replaced this default, dummy version of XERCNT
c    with a customized routine, it can then be used to override the 
c    value of KONTROL used in processing this message by redefining its value.
c
c    KONTRL may be set to any value from -2 to 2.
c
c    The meanings for KONTRL are the same as in XSETF, except
c    that the value of KONTRL changes only for this message.
c
c    If KONTRL is set to a value outside the range from -2 to 2,
c    it will be moved back into that range.
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*(*) LIBRAR, the library or software package
c    from which the error message is coming.
c
c    Input, character*(*) SUBROU, the subroutine or function within
c    the library, from which the error message is coming.
c
c    Input, character*(*) MESSG, the error message.
c
c    Input, integer NERR, the error number.
c
c    Input, integer LEVEL, the error severity level.
c    * 2, this is an unconditionally fatal error.
c    * 1, this is a recoverable error.  It is normally non-fatal, unless
c         KONTRL has been reset by XSETF.
c    * 0, this is a warning message only.
c    *-1, this is a warning message which is to be printed at most once, 
c         regardless of how many times this call is executed.
c
c    Input/output, integer KONTRL.  This routine receives the current
c    value of KONTRL, and may reset it.  The change is effective only
c    for the current error message.  This allows the user to suppress
c    or force printing of certain messages, for instance.
c
      implicit none

      integer kontrl
      integer level
      character*(*) librar
      character*(*) messg
      integer nerr
      character*(*) subrou

      return
      end
      subroutine xerctl ( messg, nmessg, nerr, level, kontrl )

c*********************************************************************72
c
cc XERCTL allows the user control over individual errors.
c
c  Discussion:
c
c    This routine gives the user control over handling of individual errors.
c
c    This routine is to be used when the error message routines XERROR
c    and XERRWV are used.  The similar routine XERCNT is to be used for 
c    the newer error message routine XERMSG.
c
c    This routine is called just after each message has been recorded, 
c    but before it is processed any further; that is, before the
c    message is printed or a decision to abort is made.
c
c    If the user wishes to influence the behavior of the error package
c    with respect to certain errors, then this dummy version of the
c    routine should be replaced by a routine that carries out the
c    actions the user desires.
c
c    In particular, the user can override the value of KONTRL used 
c    in processing this message by redefining its value.
c
c    KONTRL may be set to any value from -2 to 2.
c    The meanings for KONTRL are the same as in XSETF, except
c    that the value of KONTRL changes only for this message.
c
c    If KONTRL is set to a value outside the range from -2 to 2,
c    it will be moved back into that range.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*NMESSG MESSG, the error message.
c
c    Input, integer NMESSG, the number of characters in the error message.
c    If NMESSG <= 0, it is presumed that no error message has been supplied.
c
c    Input, integer NERR, the error number.
c
c    Input, integer LEVEL, the error severity level.
c
c    Input/output, integer KONTRL.  This routine receives the current
c    value of KONTRL, and may reset it.  The change is effective only
c    for the current error message.  This allows the user to suppress
c    or force printing of certain messages, for instance.
c
      implicit none

      integer kontrl
      integer level
      character*(*) messg
      integer nerr
      integer nmessg

      return
      end
      subroutine xerdmp ( )

c*********************************************************************72
c
cc XERDMP prints the error tables and then clears them.
c
c  Discussion:
c
c    The error handling package keeps a record of the most recent
c    errors, and the number of times they occurred.  Calling this
c    routine will print out a summary table of these records.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    None
c
      implicit none

      integer count
      integer level
      character*(1) messg
      integer nerr
      integer nmessg
  
      messg = ' '
      nmessg = 0
      nerr = 0
      level = 0
      count = 0

      call xersav ( messg, nmessg, nerr, level, count )

      return
      end
      subroutine xerhlt ( messg )

c*********************************************************************72
c
cc XERHLT aborts program execution.
c
c  Discussion:
c
c    This routine aborts the execution of the program.
c
c    The error message causing the abort is given in the calling
c    sequence.
c
c    This routine is used when the error message handler XERMSG is
c    employed.  The similar routine XERABT is to be used when the
c    older error message handlers XERROR and XERRWV are used.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*(*) MESSG, the error message associated
c    with the halt in execution.
c
      implicit none

      character*(*) messg

      stop
      end
      subroutine xermax ( maxmes )

c*********************************************************************72
c
cc XERMAX sets the maximum number of times an error message is printed.
c
c  Discussion:
c
c    This routine sets the maximum number of times any error message
c    is to be printed.  That is, a non-fatal message associated with
c    a particular numbered error should not be be printed more than
c    MAXMES times.
c
c    Most error messages won't be printed at all if the error printout
c    suppression mode has been set.  That is the case if the variable
c    KONTRL has been set to zero.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer MAXMES, the maximum number of times that
c    any error message should be printed.
c
      implicit none

      integer j4save
      integer junk
      integer maxmes
      logical set
      integer value
      integer which

      which = 4
      value = maxmes
      set = .true.

      junk = j4save ( which, value, set )

      return
      end
      subroutine xermsg ( librar, subrou, messg, nerr, level )

c*********************************************************************72
c
cc XERMSG processes error messages.
c
c  Description:
c
c    This routine processes a diagnostic message in a manner determined by the
c    value of LEVEL and the current value of the library error control
c    flag, KONTRL.
c
c    See subroutine XSETF for details on KONTRL.
c
c  Modified:
c
c    08 April 2007
c
c  Author:
c
c    Kirby Fong
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*(*) LIBRAR, the name of the library from which the
c    error message was generated. 
c
c    Input, character*(*) SUBROU, the name of the subroutine or function
c    from which the error message was generated.
c
c    Input, character(*) MESSG, the text of the error or warning message.
c    In the example below, the message is a character constant that 
c    contains a generic message.
c
c      call xermsg ('SLATEC', 'MMPY',
c      'The order of the matrix exceeds the row dimension', 3, 1)
c
c    It is possible (and is sometimes desirable) to generate a
c    specific message--e.g., one that contains actual numeric
c    values.  Specific numeric values can be converted into
c    character strings using formatted WRITE statements into
c    character variables.  This is called standard Fortran
c    internal file I/O and is exemplified in the first three
c    lines of the following example.  You can also catenate
c    substrings of characters to construct the error message.
c    Here is an example showing the use of both writing to
c    an internal file and catenating character strings.
c
c      character*5 charn, charl
c      write (charn,'(i5)') n
c      write (charl,'(i5)') lda
c      call xermsg ('SLATEC', 'MMPY', 'The order'//charn//
c      ' of the matrix exceeds its row dimension of'// charl, 3, 1)
c
c    There are two subtleties worth mentioning.  One is that
c    the // for character catenation is used to construct the
c    error message so that no single character constant is
c    continued to the next line.  This avoids confusion as to
c    whether there are trailing blanks at the end of the line.
c    The second is that by catenating the parts of the message
c    as an actual argument rather than encoding the entire
c     message into one large character variable, we avoid
c    having to know how long the message will be in order to
c    declare an adequate length for that large character
c    variable.  XERMSG calls XERPRN to print the message using
c    multiple lines if necessary.  If the message is very long,
c    XERPRN will break it into pieces of 72 characters (as
c    requested by XERMSG) for printing on multiple lines.
c    Also, XERMSG asks XERPRN to prefix each line with ' *  '
c    so that the total line length could be 76 characters.
c    Note also that XERPRN scans the error message backwards
c    to ignore trailing blanks.  Another feature is that
c    the substring '$$' is treated as a new line sentinel
c    by XERPRN.  If you want to construct a multiline
c    message without having to count out multiples of 72
c    characters, just use '$$' as a separator.  '$$'
c    obviously must occur within 72 characters of the
c    start of each line to have its intended effect since
c    XERPRN is asked to wrap around at 72 characters in
c    addition to looking for '$$'.
c
c    Input, integer NERR, the error number, chosen by the library routine's
c    author.  It must be in the range -99 to 999 (three printable digits).  
c    Each distinct error should have its own error number.  These error 
c    numbers should be described in the machine readable documentation 
c    for the routine.  The error numbers need be unique only within each 
c    routine, so it is reasonable for each routine to start enumerating
c    errors from 1 and proceeding to the next integer.
c
c    Input, integer LEVEL, a value in the range 0 to 2 that indicates the
c    level (severity) of the error.  Their meanings are
c    * -1: A warning message.  This is used if it is not clear
c    that there really is an error, but the user's attention
c    may be needed.  An attempt is made to only print this
c    message once.
c    * 0: A warning message.  This is used if it is not clear
c    that there really is an error, but the user's attention
c    may be needed.
c    * 1: A recoverable error.  This is used even if the error is
c    so serious that the routine cannot return any useful
c    answer.  If the user has told the error package to
c    return after recoverable errors, then XERMSG will
c    return to the Library routine which can then return to
c    the user's routine.  The user may also permit the error
c    package to terminate the program upon encountering a
c    recoverable error.
c
c    * 2: A fatal error.  XERMSG will not return to its caller
c    after it receives a fatal error.  This level should
c    hardly ever be used; it is much better to allow the
c    user a chance to recover.  An example of one of the few
c    cases in which it is permissible to declare a level 2
c    error is a reverse communication Library routine that
c    is likely to be called repeatedly until it integrates
c    across some interval.  If there is a serious error in
c    the input such that another step cannot be taken and
c    the Library routine is called again without the input
c    error having been corrected by the caller, the Library
c    routine will probably be called forever with improper
c    input.  In this case, it is reasonable to declare the
c    error to be fatal.
c
      implicit none

      integer i
      integer j4save
      integer kdummy
      integer kount
      integer lerr
      integer level
      character*20  lfirst
      character*(*) librar
      integer lkntrl
      integer llevel
      integer ltemp
      integer maxmes
      character*(*) messg
      integer mkntrl
      integer nerr
      logical set
      logical skip
      character*(*) subrou
      character*72 temp
      integer value
      integer which
      character*8 xlibr
      character*8 xsubr

      which = 2
      value = 0
      set = .false.

      lkntrl = j4save ( which, value, set )

      which = 4
      value = 0
      set = .false.

      maxmes = j4save ( which, value, set )
c
c  LKNTRL is a local copy of the control flag KONTRL.
c
c  MAXMES is the maximum number of times any particular message
c  should be printed.
c
c  We print a fatal error message and terminate for an error in
c  calling XERMSG.  The error number should be positive,
c  and LEVEL should be between 0 and 2.
c
      if ( nerr .lt. -9999999 .or. 
     &     99999999 .lt. nerr .or. 
     &     nerr .eq. 0 .or.
     &     level .lt. -1 .or. 
     &     2 .lt. level ) then

        call xerprn ( ' ***', -1, 
     &    'Fatal error in...$$ ' //
     &    'XERMSG -- Invalid error number or level$$ ' //
     &    'Job abort due to fatal error.', 72 )

        call xersve ( ' ', ' ', ' ', 0, 0, 0, kdummy )
        call xerhlt ( ' ***XERMSG -- Invalid input' )
        return

      end if
c
c  Record the message.
c
      which = 1
      value = nerr
      set = .true.

      i = j4save ( which, value, set )

      call xersve ( librar, subrou, messg, 1, nerr, level, kount )
c
c  Handle print-once warning messages.
c
      if ( level .eq. -1 .and. 1 .lt. kount ) then
        return
      end if
c
c  Allow temporary user override of the control flag.
c
      xlibr  = librar
      xsubr  = subrou
      lfirst = messg
      lerr   = nerr
      llevel = level
      call xercnt ( xlibr, xsubr, lfirst, lerr, llevel, lkntrl )

      lkntrl = max ( -2, min ( 2, lkntrl ) )
      mkntrl = abs ( lkntrl )
c
c  Skip printing if the control flag value as reset in xercnt is
c  zero and the error is not fatal.
c
      skip = .false.

      if ( level .lt. 2 .and. lkntrl .eq. 0 ) then
        skip = .true.
      end if

      if ( level .eq. 0 .and. maxmes .lt. kount ) then
        skip = .true.
      end if

      if ( level .eq. 1 .and. 
     &     maxmes .lt. kount .and. 
     &     mkntrl .eq. 1 ) then
        skip = .true.
      end if

      if ( level .eq. 2 .and. max ( 1, maxmes ) .lt. kount ) then
        skip = .true.
      end if

      if ( .not. skip ) then
c
c  Announce the names of the library and subroutine by building a
c  message in character variable TEMP (not exceeding 66 characters)
c  and sending it out via XERPRN.  Print only if control flag
c  is not zero.
c
      if ( lkntrl .ne. 0 ) then
        temp(1:21) = 'Message from routine '
        i = min ( len ( subrou ), 16 )
        temp(22:21+i) = subrou(1:i)
        temp(22+i:33+i) = ' in library '
        ltemp = 33 + i
        i = min ( len ( librar ), 16)
        temp(ltemp+1:ltemp+i) = librar (1:i)
        temp(ltemp+i+1:ltemp+i+1) = '.'
        ltemp = ltemp + i + 1
        call xerprn ( ' ***', -1, temp(1:ltemp), 72 )
      end if
c
c  If LKNTRL is positive, print an introductory line before
c  printing the message.  The introductory line tells the choice
c  from each of the following three options.
c
c  1.  Level of the message
c
c  'Informative message'
c  'Potentially recoverable error'
c  'Fatal error'
c
c  2.  Whether control flag will allow program to continue
c
c  'Prog continues'
c  'Prog aborted'
c
c  3.  Whether or not a traceback was requested.  (The traceback
c  may not be implemented at some sites, so this only tells
c  what was requested, not what was delivered.)
c
c  'Traceback requested'
c  'Traceback not requested'
c
c  Notice that the line including four prefix characters will not
c  exceed 74 characters.
c  We skip the next block if the introductory line is not needed.
c
      if ( 0 .lt. lkntrl ) then
c
c  The first part of the message tells about the level.
c
        if ( level .le. 0 ) then
          temp(1:20) = 'Informative message,'
          ltemp = 20
        else if ( level .eq. 1 ) then
          temp(1:30) = 'Potentially recoverable error,'
          ltemp = 30
        else
          temp(1:12) = 'Fatal error,'
          ltemp = 12
        end if
c
c  Then whether the program will continue.
c
        if ( ( mkntrl .eq. 2 .and. 1 .le. level ) .or.
     &       ( mkntrl .eq. 1 .and. level .eq. 2 ) ) then
          temp(ltemp+1:ltemp+14) = ' Prog aborted,'
          ltemp = ltemp + 14
        else
          temp(ltemp+1:ltemp+16) = ' Prog continues,'
          ltemp = ltemp + 16
        end if
c
c  Finally tell whether there should be a traceback.
c
        if ( 0 .lt. lkntrl ) then
          temp(ltemp+1:ltemp+20) = ' Traceback requested'
          ltemp = ltemp + 20
        else
          temp(ltemp+1:ltemp+24) = ' Traceback not requested'
          ltemp = ltemp + 24
        end if

        call xerprn ( ' ***', -1, temp(1:ltemp), 72 )

      end if
c
c  Now send out the message.
c
      call xerprn ( ' *  ', -1, messg, 72 )
c
c  IF LKNTRL is positive, write the error number and request a
c  traceback.
c
      if ( 0 .lt. lkntrl ) then

        write ( temp, '(a,i8)' ) '  Error number = ', nerr

        call xerprn ( ' *  ', -1, temp, 72 )
        call fdump ( )

      end if
c
c  IF LKNTRL is not zero, print a blank line and an end of message.
c
      if ( lkntrl .ne. 0 ) then
        call xerprn ( ' *  ', -1, ' ', 72 )
        call xerprn ( ' ***', -1, 'End of message', 72 )
        call xerprn ( '    ',  0, ' ', 72 )
      end if
c
c  If the error is not fatal or the error is recoverable and the
c  control flag is set for recovery, then return.
c
      end if

      if ( level .le. 0 .or. 
     &  ( level .eq. 1 .and. mkntrl .le. 1 ) ) then
        return
      end if
c
c  The program will be stopped due to an unrecovered error or a
c  fatal error.  Print the reason for the abort and the error
c  summary if the control flag and the maximum error count permit.
c
      if ( 0 .lt. lkntrl .and. kount .lt. max ( 1, maxmes ) ) then

        if ( level .eq. 1 ) then
          call xerprn
     &         ( ' ***', -1, 'Job abort due to unrecovered error.', 72 )
        else
          call xerprn ( ' ***', -1, 'Job abort due to fatal error.',
     &      72 )
        end if

        call xersve ( ' ', ' ', ' ', -1, 0, 0, kdummy )
        call xerhlt ( ' ' )

      else

        call xerhlt ( messg )

      end if

      return
      end
      subroutine xerprn ( prefix, npref, messg, nwrap )

c*********************************************************************72
c
cc XERPRN prints error messages processed by XERMSG.
c
c  Description:
c
c  Discussion:
c
c    This routine is used by the error handling routine XERMSG.  A related
c    routine, XERPRT, is used by the older error handling routines 
c    XERROR and XERRWV.
c
c    This routine sends one or more lines to each of the (up to five)
c    logical units to which error messages are to be sent.  This routine
c    is called several times by XERMSG, sometimes with a single line to
c    print and sometimes with a (potentially very long) message that may
c    wrap around into multiple lines.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Kirby Fong
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*(*) PREFIX, a string to be put at the beginning of 
c    each line before the body of the message.  No more than 16 characters
c    of PREFIX will be used.
c
c    Input, integer NPREF, the number of characters to use from PREFIX. 
c    If it is negative, the intrinsic function LEN is used to determine 
c    its length.  If it is zero, PREFIX is not used.  If it exceeds 16 or if
c    LEN(PREFIX) exceeds 16, only the first 16 characters will be
c    used.  If NPREF is positive and the length of PREFIX is less
c    than NPREF, a copy of PREFIX extended with blanks to length
c    NPREF will be used.
c
c    Input, character*(*) MESSG, the error message.  If it is a long message, 
c    it will be broken into pieces for printing on multiple lines.  Each line
c    will start with the appropriate prefix and be followed by a piece of the 
c    message.  NWRAP is the number of characters per piece; that is, after 
c    each NWRAP characters, we break and start a new line.  In addition,
c    the characters '$$' embedded in MESSG are a sentinel for a new line.  
c    The counting of characters up to NWRAP starts over for each new line.  
c    The value of NWRAP typically used by XERMSG is 72 since many
c    older error messages in the SLATEC Library are laid out to rely on 
c    wrap-around every 72 characters.
c
c    Input, integer NWRAP, the maximum size piece into which to break MESSG 
c    for printing on multiple lines.  An embedded '$$' ends a line, and the 
c    count restarts at the following character.  If a line break does not occur
c    on a blank (it would split a word) that word is moved to the next line.  
c    Values of NWRAP less than 16 will be treated as 16.  Values of NWRAP 
c    greater than 132 will be treated as 132.  The actual line length will 
c    be NPREF + NWRAP after NPREF has been adjusted to fall between 0 and 16 
c    and NWRAP has been adjusted to fall between 16 and 132.
c
      implicit none

      character*148 cbuff
      integer i
      integer i1mach
      integer idelta
      integer iu(5)
      integer lenmsg
      integer lpiece
      integer lpref
      integer lwrap
      character*(*) messg
      integer n
      character*2 newlin
      parameter ( newlin = '$$' )
      integer nextc
      integer npref
      integer nunit
      integer nwrap
      character*(*) prefix

      call xgetua ( iu, nunit )
c
c  A zero value for a logical unit number means to use the standard
c  error message unit instead.  I1MACH(4) retrieves the standard
c  error message unit.
c
      n = i1mach(4)
      do i = 1, nunit
        if ( iu(i) .eq. 0 ) then
          iu(i) = n
        end if
      end do
c
c  LPREF is the length of the prefix.  The prefix is placed at the
c  beginning of CBUFF, the character buffer, and kept there during
c  the rest of this routine.
c
      if ( npref .lt. 0 ) then
        lpref = len ( prefix )
      else
        lpref = npref
      end if

      lpref = min ( 16, lpref )

      if ( lpref .ne. 0 ) then
        cbuff(1:lpref) = prefix
      end if
c
c  LWRAP is the maximum number of characters we want to take at one
c  time from MESSG to print on one line.
c
      lwrap = max ( 16, min ( 132, nwrap ) )
c
c  Set LENMSG to the length of MESSG, ignore any trailing blanks.
c
      lenmsg = len ( messg )
      n = lenmsg
      do i = 1, n
        if ( messg(lenmsg:lenmsg) .ne. ' ' ) then
          go to 30
        end if
        lenmsg = lenmsg - 1
      end do

   30 continue
c
c  If the message is all blanks, then print one blank line.
c
      if ( lenmsg .eq. 0 ) then
        cbuff(lpref+1:lpref+1) = ' '
        do i = 1, nunit
          write ( iu(i), '(a)' ) cbuff(1:lpref+1)
        end do
        return
      end if
c
c  Set NEXTC to the position in MESSG where the next substring
c  starts.  From this position we scan for the new line sentinel.
c  When NEXTC exceeds LENMSG, there is no more to print.
c  We loop back to label 50 until all pieces have been printed.
c
c  We look for the next occurrence of the new line sentinel.  The
c  INDEX intrinsic function returns zero if there is no occurrence
c  or if the length of the first argument is less than the length
c  of the second argument.
c
c  There are several cases which should be checked for in the
c  following order.  We are attempting to set LPIECE to the number
c  of characters that should be taken from MESSG starting at
c  position NEXTC.
c
c  * LPIECE == 0
c  The new line sentinel does not occur in the remainder of the 
c  character string.  LPIECE should be set to LWRAP or LENMSG+1-NEXTC,
c  whichever is less.
c
c  * LPIECE == 1
c  The new line sentinel starts at MESSG(NEXTC:NEXTC).  LPIECE is effectively 
c  zero, and we print nothing to avoid producing unnecessary blank lines.
c  This takes care of the situation where the library routine has a message of
c  exactly 72 characters followed by a new line sentinel followed by more
c  characters.  NEXTC should be incremented by 2.
c
c  * LWRAP + 1 < LPIECE
c  Reduce LPIECE to LWRAP.
c
c  * Otherwise
c  This last case means 2 <= LPIECE <= LWRAP+1.  Reset LPIECE = LPIECE-1.
c  Note that this properly handles the end case where LPIECE = LWRAP+1.  
c  That is, the sentinel falls exactly at the end of a line.
c
      nextc = 1

   50 continue

      lpiece = index ( messg(nextc:lenmsg), newlin )
c
c  There was no new line sentinel found.
c
      if ( lpiece .eq. 0 ) then

        idelta = 0
        lpiece = min ( lwrap, lenmsg + 1 - nextc )

        if ( lpiece .lt. lenmsg + 1 - nextc ) then
          do i = lpiece+1, 2, -1
            if ( messg(nextc+i-1:nextc+i-1) .eq. ' ' ) then
              lpiece = i - 1
              idelta = 1
              go to 54
            end if
          end do
        end if

   54   continue

        cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
        nextc = nextc + lpiece + idelta
c
c  We have a new line sentinel at MESSG(NEXTC:NEXTC+1).
c  Don't print a blank line.
c
      else if ( lpiece .eq. 1 ) then

        nextc = nextc + 2
        go to 50
c
c  LPIECE should be set down to LWRAP.
c
      else if ( lwrap + 1 .lt. lpiece ) then

        idelta = 0
        lpiece = lwrap

        do i = lpiece + 1, 2, -1
          if ( messg(nextc+i-1:nextc+i-1) .eq. ' ' ) then
            lpiece = i - 1
            idelta = 1
            go to 58
          end if
        end do

   58   continue

        cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
        nextc = nextc + lpiece + idelta
c
c  If we arrive here, it means 2 <= LPIECE <= LWRAP+1.
c  We should decrement LPIECE by one.
c
      else

        lpiece = lpiece - 1
        cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
        nextc = nextc + lpiece + 2
      end if
c
c  Print
c
      do i = 1, nunit
        write ( iu(i), '(a)' ) cbuff(1:lpref+lpiece)
      end do

      if ( nextc .le. lenmsg ) then
        go to 50
      end if

      return
      end
      subroutine xerprt ( messg, nmessg )

c*********************************************************************72
c
cc XERPRT prints a message on each file receiving error messages.
c
c  Discussion:
c
c    This routine is used by the error handling routines XERROR and
c    XERRWV.  A related routine, XERPRN, is used by the more modern
c    error handling routine XERMSG.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*(*) MESSG, the message to be printed.
c
c    Input, integer NMESSG, the length of the message.
c    This parameter is actually ignored; the length of MESSG is
c    determined by the LEN function.
c
      implicit none

      integer ichar
      integer iunit
      integer kunit
      integer last
      integer lun(5)
      character*(*) messg
      integer messg_len
      integer nmessg
      integer nunit
c
c  Obtain output unit numbers.
c
      call xgetua ( lun, nunit )

      messg_len = len ( messg )

      do kunit = 1, nunit

        iunit = lun(kunit)

        do ichar = 1, messg_len, 72

          last = min ( ichar + 71, messg_len )

          if ( iunit .eq. 0 ) then
            write ( *, '(a)' ) messg(ichar:last)
          else
            write ( iunit, '(a)' ) messg(ichar:last)
          end if

        end do

      end do

      return
      end
      subroutine xerror ( messg, nmessg, nerr, level )

c*********************************************************************72
c
cc XERROR processes a diagnostic error message.
c
c  Discussion:
c
c    This routine processes a diagnostic message, in a manner determined 
c    by the value of LEVEL and the current value of the library error 
c    control flag KONTRL.
c
c    See XSETF for details about KONTRL.
c
c  Example:
c
c    call xerror ( 'SMOOTH -- NUM was zero.', 23, 1, 2 )
c
c    call xerror ( 'INTEG  -- Less than full accuracy achieved.', 43, 2, 1 )
c
c    call xerror ( 'ROOTER -- Actual zero of f found before interval
c    fully collapsed.', 65, 3, 0 )
c
c    call xerror ( 'EXP    -- Underflows being set to zero.', 39, 1, -1 )
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*NMESSG MESSG, the error message.
c
c    Input, integer NMESSG, the number of characters in the error message.
c    If NMESSG <= 0, it is presumed that no error message has been supplied.
c
c    Input, integer NERR, the error number.  NERR should not be 0.
c
c    Input, integer LEVEL, the error severity level.
c    * 2, this is an unconditionally fatal error.
c    * 1, this is a recoverable error.  It is normally non-fatal, unless
c         KONTRL has been reset by XSETF.
c    * 0, this is a warning message only.
c    *-1, this is a warning message which is to be printed at most once, 
c         regardless of how many times this call is executed.
c
      implicit none

      integer nmessg

      integer i1
      integer i2
      integer level
      character*(*) messg
      integer nerr
      integer ni
      integer nr
      real r1
      real r2

      ni = 0
      i1 = 0
      i2 = 0
      nr = 0
      r1 = 0.0E+00
      r2 = 0.0E+00

      call xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2 )

      return
      end
      subroutine xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, 
     &  nr, r1, r2 )

c*********************************************************************72
c
cc XERRWV processes a diagnostic error message.
c
c  Discussion:
c
c    This routine processes a diagnostic message, in a manner determined
c    by the value of LEVEL and the current value of the library error
c    control flag KONTRL.
c
c    See XSETF for details about KONTRL.
c
c    In addition, up to two integer values and two real values may be 
c    printed along with the message.
c
c  Example:
c
c    call xerrwv ( 'SMOOTH -- NUM (=I1) was zero.', 
c    29, 1, 2, 1, num, 0, 0, 0.0, 0.0 )
c
c    call xerrwv (
c    'QUADXY -- Requested error (R1) less than minimum (R2)',
c    54, 77, 1, 0, 0, 0, 2, errreq, errmin )
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*NMESSG MESSG, the error message.
c
c    Input, integer NMESSG, the number of characters in the error message.
c    If NMESSG <= 0, it is presumed that no error message has been supplied.
c
c    Input, integer NERR, the error number.  NERR should not be 0.
c
c    Input, integer LEVEL, the error severity level.
c    * 2, this is an unconditionally fatal error.
c    * 1, this is a recoverable error.  It is normally non-fatal, unless
c         KONTRL has been reset by XSETF.
c    * 0, this is a warning message only.
c    *-1, this is a warning message which is to be printed at most once, 
c         regardless of how many times this call is executed.
c
c    Input, integer NI, the number of integers to print, 0, 1 or 2.
c
c    Input, integer N1, N2, the two integers that may be printed.
c
c    Input, integer NR, the number of single precision real values
c    to print, 0, 1 or 2.
c
c    Input, real R1, R2, the two single precision real values that
c    may be printed.
c
      implicit none

      character*37 form
      integer i
      integer i1
      integer i1mach
      integer i2
      integer ifatal
      integer isizef
      integer isizei
      integer iunit
      integer j4save
      integer junk
      integer kdummy
      integer kount
      integer kunit
      integer level
      integer lerr
      character*20 lfirst
      integer lifatal
      integer lkntrl
      integer llevel
      integer lmessg
      integer lun(5)
      integer maxmes
      character*(*) messg
      integer mkntrl
      integer nerr
      integer ni
      integer nmessg
      integer nr
      integer nunit
      real r1
      real r2
      logical set
      logical skip
      integer value
      integer which
c
c  Get flags.
c
      which = 2
      value = 0
      set = .false.

      lkntrl = j4save ( which, value, set )

      which = 4
      value = 0
      set = .false.

      maxmes = j4save ( which, value, set )
c
c  Check for valid input.
c
      if ( nmessg .le. 0 .or. 
     &     nerr .eq. 0 .or.
     &     level .lt. -1 .or.
     &     2 .lt. level ) then

        if ( 0 .lt. lkntrl ) then
          call xerprt ( 'Fatal error in...', 17 )
        end if

        call xerprt ( 'XERROR -- Invalid input', 23 )

        if ( 0 .lt. lkntrl ) then

          call fdump ( )

          call xerprt ( 'Job abort due to fatal error.', 29 )

          call xersav ( ' ', 0, 0, 0, kdummy )

        end if

        call xerabt ( 'XERROR -- invalid input', 23 )

        return

      end if
c
c  Record message.
c
      which = 1
      value = nerr
      set = .true.

      junk = j4save ( which, value, set )

      call xersav ( messg, nmessg, nerr, level, kount )
c
c  Let the user override.
c
      lmessg = nmessg
      lerr = nerr
      llevel = level

      call xerctl ( messg, lmessg, lerr, llevel, lkntrl )
c
c  Reset to original values.
c
      lmessg = nmessg
      lerr = nerr
      llevel = level

      lkntrl = max ( -2, min ( 2, lkntrl ) )
      mkntrl = abs ( lkntrl )
c
c  Decide whether to print message.
c
      skip = 
     &  ( llevel .lt. 2 .and. 
     &    lkntrl .eq .0 ) .or.
     &  ( llevel .eq. -1  .and.
     &    kount  .gt. min ( 1, maxmes ) ) .or.
     &  ( llevel .eq.  0  .and.
     &    kount  .gt. maxmes ) .or.
     &  ( llevel .eq.  1  .and. 
     &    kount  .gt. maxmes .and.
     &    mkntrl .eq. 1 ) .or.
     &  ( llevel .eq.  2  .and.
     &    kount  .gt. max ( 1, maxmes ) )

      if ( .not. skip ) then
c
c  Introduction.
c
        if ( 0 .lt. lkntrl ) then

          if ( llevel .eq. (-1) ) then

            call xerprt ( ' ', 1 )
            call xerprt ( 'One time warning message...', 27 )

          else if ( llevel .eq. 0 ) then

            call xerprt ( ' ', 1 )
            call xerprt ( 'Warning in...', 13 )

          else if ( llevel .eq. 1 ) then

            call xerprt ( ' ', 1 )
            call xerprt ( 'Recoverable error in...', 23 )

          else if ( llevel .eq. 2 ) then

            call xerprt ( ' ', 1 )
            call xerprt ( 'Fatal error in...', 17 )

          end if

        end if
c
c  Message.
c
        call xerprt ( messg, lmessg )
        call xgetua ( lun, nunit )

        isizei = 1 + int ( log10 ( real ( i1mach(9) ) ) )
        isizef = 1 + int ( log10 ( real ( i1mach(10) )**i1mach(11) ) )

        do kunit = 1, nunit

          iunit = lun(kunit)

          do i = 1, min ( ni, 2 )
            write ( form, 21 ) i, isizei
   21       format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
            if ( iunit .eq. 0 ) then
              if ( i .eq. 1 ) then
                write ( *, form ) i1
              else if ( i .eq. 2 ) then
                write ( *, form ) i2
              end if
            else
              if ( i .eq. 1 ) then
                write (iunit,form) i1
              else if ( i .eq. 2 ) then
                write (iunit,form) i2
              end if
            end if
          end do

          do i = 1, min ( nr, 2 )
            write ( form, 23 ) i, isizef+10, isizef
   23       format ('(11x,21hin above message, r',i1,'=,e',
     1      i2,'.',i2,')')
            if ( iunit .eq. 0 ) then
              if ( i .eq. 1 ) then
                write ( *, form ) r1
              else if ( i .eq. 2 ) then
                write ( *, form ) r2
              end if
            else
              if ( i .eq. 1 ) then
                write ( iunit, form ) r1
              else if ( i .eq. 2 ) then
                write ( iunit, form ) r2
              end if
            end if
          end do
c
c  Print the error number.
c
          if ( 0 .lt. lkntrl ) then

            if ( iunit .eq. 0 )then
              write ( *, '(a,i10)' ) '  Error number = ', lerr
            else
              write ( iunit, '(a,i10)') '  Error number = ', lerr
            end if

          end if

        end do
c
c  Trace-back.
c
        if ( 0 .lt. lkntrl ) then
          call fdump ( )
        end if

      end if

      if ( llevel .eq. 2 .or.
     &   ( llevel .eq. 1 .and. mkntrl .eq. 2 ) ) then
        ifatal = 1
      else
        ifatal = 0
      end if
c
c  Quit here if message is not fatal.
c
      if ( ifatal .le. 0 ) then
        return
      end if
c
c  Print reason for abort and error summary.
c
      if ( 0 .lt. lkntrl .and. kount .le. max ( 1, maxmes ) ) then

        if ( llevel .eq. 1 ) then
          call xerprt  ( 'Job abort due to unrecovered error.', 35 )
        else if ( llevel .eq. 2 ) then
          call xerprt ( 'Job abort due to fatal error.', 29 )
        end if

        call xersav ( ' ', -1, 0, 0, kdummy )

      end if
c
c  Abort.
c
      if ( llevel .eq. 2 .and.
     &  max ( 1, maxmes ) .lt. kount ) then
        lmessg = 0
      end if

      call xerabt ( messg, lmessg )

      return
      end
      subroutine xersav ( messg, nmessg, nerr, level, count )

c*********************************************************************72
c
cc XERSAV records that a particular error has occurred.
c
c  Discussion:
c
c    This routine maintains a table of up to 10 error messages,
c    and records the number of times each error occurred.  
c
c    This same routine can be "persuaded" to print out the
c    contents of the error table, by the trick of passing a
c    value of NMESSG that is less than or equal to zero.
c
c    If NMESSG is exactly zero, then once they have been printed,
c    the error tables are also cleared, that is, reset as though
c    no errors had occurred.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*NMESSG MESSG, the error message.
c
c    Input, integer NMESSG, the number of characters in the error message.
c    If NMESSG = 0, the tables will be dumped and cleared.
c    If NMESSG < 0, the tables will be dumped, but not cleared.
c
c    Input, integer NERR, the error number.  NERR should not be 0.
c
c    Input, integer LEVEL, the error severity level.
c    * 2, this is an unconditionally fatal error.
c    * 1, this is a recoverable error.  It is normally non-fatal, unless
c         KONTRL has been reset by XSETF.
c    * 0, this is a warning message only.
c    *-1, this is a warning message which is to be printed at most once, 
c         regardless of how many times this call is executed.
c
c    Output, integer COUNT, the number of times this message has
c    been seen, or zero if the table has overflowed and does not contain 
c    this message specifically.
c    When NMESSG = 0, COUNT will not be altered.
c
      implicit none

      integer count
      integer i
      integer i1mach
      integer ii
      integer iunit
      integer kount(10)
      integer kountx
      integer kunit
      integer level
      integer levtab(10)
      integer lun(5)
      character*20 mes
      character*(*) messg
      character*20 mestab(10)
      integer nerr
      integer nertab(10)
      integer nmessg
      integer nunit

      save kount
      save kountx
      save levtab
      save mestab
      save nertab

      data kount /0,0,0,0,0,0,0,0,0,0/
      data kountx / 0 /
      data levtab /0,0,0,0,0,0,0,0,0,0/
      data mestab /' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/
      data nertab /0,0,0,0,0,0,0,0,0,0/
c
c  Dump the table.
c
      if ( nmessg .le. 0 ) then

        if ( kount(1) .eq. 0 ) then
          return
        end if
c
c  Print to each unit.
c
        call xgetua ( lun, nunit )

        do kunit = 1, nunit

          iunit = lun(kunit)
          if ( iunit .eq. 0 ) then
            iunit = i1mach(4)
          end if
c
c  Print table header.
c
          write ( iunit, '(a)' ) ' '
          write ( iunit, '(a)' )
     &  '         Error message summary'
          write ( iunit, '(a)' )
     &  'Message start             NERR     Level     Count'
c
c  Print body of table.
c
          do i = 1, 10

            if ( kount(i) .eq. 0 ) then
              go to 30
            end if

            write ( iunit,'(a20,3i10)' ) 
     &        mestab(i), nertab(i), levtab(i), kount(i)

          end do

   30     continue
c
c  Print number of other errors.
c
          if ( kountx .ne. 0 ) then

            write ( iunit, '(a)' ) ' '
            write ( iunit, '(a,i10)' ) 
     &        'Other errors not individually tabulated = ', kountx

          end if

          write ( iunit, '(a)' ) ' '

        end do

        if ( nmessg .lt. 0 ) then
          return
        end if
c
c  Clear the error tables.
c
        do i = 1, 10
          kount(i) = 0
        end do

        kountx = 0
c
c  Process a message.
c
c  Search for this message, or else an empty slot for this message,
c  or else determine that the error table is full.
c
      else

        mes(1:20) = messg(1:20)

        do i = 1, 10

          ii = i
c
c  Empty slot found for new message.
c
          if ( kount(i) .eq. 0 ) then
            mestab(ii) = mes
            nertab(ii) = nerr
            levtab(ii) = level
            kount(ii)  = 1
            count = 1
            return
          end if
c
c  Message found in table.
c
          if ( mes .eq. mestab(i) .and.
     &         nerr .eq. nertab(i) .and.
     &         level .eq. levtab(i) ) then
            kount(ii) = kount(ii) + 1
            count = kount(ii)
            return
          end if

        end do
c
c  Table is full.
c
        kountx = kountx + 1
        count = 1

      end if

      return
      end
      subroutine xersve ( librar, subrou, messg, kflag, nerr, level,
     &  icount )

c*********************************************************************72
c
cc XERSVE records that an error has occurred.
c
c  Discussion:
c
c    This routine is used by the error handling routines associated
c    with XERMSG.  It is a revised version of the routine XERSAV, which
c    was used with the older pair of error handling routines XERROR
c    and XERRWV.
c
c  Modified:
c
c    06 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character*(*) LIBRAR, the name of the library or software package
c    from which the error message comes.
c
c    Input, character*(*) SUBROU, the name of the subroutine or function
c    from which the error message comes.
c
c    Input, character*(*) MESSG, the error message.
c
c    Input, integer KFLAG, indicates the action to be performed.
c    0 < KFLAG, the message in MESSG is saved.
c    KFLAG=0 the tables will be dumped and cleared.
c    KFLAG < 0, the tables will be dumped and not cleared.
c
c    Input, integer NERR, the error number.
c
c    Input, integer LEVEL, the error severity.
c
c    Output, integer ICOUNT, the number of times this message has been seen,
c    or zero if the table has overflowed and does not contain this message 
c    specifically.  When KFLAG=0, ICOUNT will not be altered from its
c    input value.
c
      implicit none

      integer lentab
      parameter ( lentab = 10 )

      integer i
      integer i1mach
      integer icount
      integer iunit
      integer kflag
      integer kount(lentab)
      integer kountx
      integer kunit
      integer level
      integer levtab(lentab)
      character*8 lib
      character*(*) librar
      character*8 libtab(lentab)
      integer lun(5)
      character*20 mes
      character*(*) messg
      character*20 mestab(lentab)
      integer nerr
      integer nertab(lentab)
      integer nmsg
      integer nunit
      character*8 sub
      character*(*) subrou
      character*8 subtab(lentab)

      save kount
      save kountx
      save levtab
      save libtab
      save mestab
      save nertab
      save nmsg
      save subtab

      data kountx /0/
      data nmsg /0/

      if ( kflag .le. 0 ) then
c
c  Dump the table.
c
        if ( nmsg .eq. 0 ) then
          return
        end if
c
c  Print to each unit.
c
        call xgetua ( lun, nunit )

        do kunit = 1, nunit

          iunit = lun(kunit)
          if ( iunit .eq. 0 ) then
            iunit = i1mach(4)
          end if
c
c  Print the table header.
c
          write ( iunit, '(a)' ) ' '
          write ( iunit, '(a)' ) '          Error message summary'
          write ( iunit, '(a,a)' )
     &      'Library    Subroutine Message start             NERR',
     &      '     Level     Count'

c
c  Print body of table.
c
          do i = 1, nmsg
            write ( iunit, '(a,3x,a,3x,a,3i10)' ) 
     &        libtab(i), subtab(i), mestab(i), nertab(i), levtab(i), 
     &        kount(i)
          end do
c
c  Print the number of other errors.
c
          if ( kountx .ne. 0 ) then
            write ( iunit, '(a)' ) ' '
            write ( iunit, '(a,i10)' )
     &        'Other errors not individually tabulated = ', kountx
          end if

          write ( iunit, '(1x)' )

        end do
c
c  Clear the error tables.
c
        if ( kflag .eq. 0 ) then
          nmsg = 0
          kountx = 0
        end if

      else
c
c  Process a message.
c
c  Search for this message, or else an empty slot for this message,
c  or else determine that the error table is full.
c
        lib = librar
        sub = subrou
        mes = messg

        do i = 1, nmsg

          if ( 
     &      lib .eq. libtab(i) .and. 
     &      sub .eq. subtab(i) .and.
     &      mes .eq. mestab(i) .and. 
     &      nerr .eq. nertab(i) .and.
     &      level .eq. levtab(i) ) then
            kount(i) = kount(i) + 1
            icount = kount(i)
            return
          end if

        end do
c
c  Empty slot found for new message.
c
        if ( nmsg .lt. lentab ) then

          nmsg = nmsg + 1
          libtab(i) = lib
          subtab(i) = sub
          mestab(i) = mes
          nertab(i) = nerr
          levtab(i) = level
          kount(i)  = 1
          icount    = 1
c
c  Table is full.
c
        else

          kountx = kountx + 1
          icount = 0

        end if

      end if

      return
      end
      subroutine xgetf ( kontrl )

c*********************************************************************72
c
cc XGETF returns the error control flag.
c
c  Discussion:
c
c    This routine returns the current value of KONTRL, the error 
c    control flag.
c
c    The amount of output printed for a given error is determined
c    by LEVEL, the level of severity, and KONTRL, which controls
c    how much output the user wants to see.
c
c    The following table shows how each message is treated,
c    depending on the values of KONTRL and LEVEL.
c
c    If KONTRL is zero or negative, no information other than the
c    message itself (including numeric values, if any) will be
c    printed.  If KONTRL is positive, introductory messages,
c    trace-backs, and so on, will be printed in addition to the message.
c
c            KONTRL   0            -1/+1          -2/+2
c        LEVEL
c          2        fatal          fatal          fatal
c          1     not printed      printed         fatal
c          0     not printed      printed        printed
c         -1     not printed      printed once   printed once
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Output, integer KONTRL, the current error control flag.
c
      implicit none

      integer j4save
      integer kontrl
      logical set
      integer value
      integer which

      which = 2
      value = 0
      set = .false.

      kontrl = j4save ( which, value, set )

      return
      end
      subroutine xgetua ( iunit, nunit )

c*********************************************************************72
c
cc XGETUA reports the unit numbers of the files receiving error messages.
c
c  Discussion:
c
c    These unit numbers may have been set by a call to XSETUN,
c    or a call to XSETUA, or may be default values.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c  Parameters:
c
c    Output, integer IUNIT(5), an array into which the routine will
c    store the values of the NUNIT units to which the error messages
c    are being sent.  The value of NUNIT is never more than 5, so
c    using an array of dimension 5 will be sufficient.
c
c    Output, integer NUNIT, the number of units to which the
c    error messages are being sent.  NUNIT will be in the
c    range from 1 to 5.
c
      implicit none

      integer i
      integer iunit(5)
      integer j4save
      integer nunit
      logical set
      integer value
      integer which

      which = 5
      value = 0
      set = .false.

      nunit = j4save ( which, value, set)

      if ( nunit .lt. 1 ) then
        return
      end if

      which = 3
      value = 0
      set = .false.

      iunit(1) = j4save ( which, value, set )

      do i = 2, nunit

        which = i + 4
        value = 0
        set = .false.

        iunit(i) = j4save ( which, value, set )

      end do

      return
      end
      subroutine xgetun ( iunit )

c*********************************************************************72
c
cc XGETUN reports the first unit number of the files receiving error messages.
c
c  Discussion:
c
c    This routine returns the unit number associated with the first or
c    main file to which error messages are being sent.  
c
c    To find out if more than one file is being used for error output, 
c    one must use the XGETUA routine.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Output, integer IUNIT, the logical unit number of the first unit
c    to which error messages are being sent.
c
      implicit none

      integer iunit
      integer j4save
      logical set
      integer value
      integer which

      which = 3
      value = 0
      set = .false.

      iunit = j4save ( which, value, set )

      return
      end
      subroutine xsetf ( kontrl )

c*********************************************************************72
c
cc XSETF sets the error control flag.
c
c  Discussion:
c
c    This routine sets the error control flag value to KONTRL.
c
c    The following table shows how each message is treated,
c    depending on the values of KONTRL and LEVEL.  See XERROR
c    for description of LEVEL.
c
c    If KONTRL is zero or negative, no information other than the
c    message itself (including numeric values, if any) will be
c    printed.  If KONTRL is positive, introductory messages,
c    trace-backs, and so on, will be printed in addition to the message.
c
c            KONTRL   0            -1/+1          -2/+2
c        LEVEL
c          2        fatal          fatal          fatal
c          1     not printed      printed         fatal
c          0     not printed      printed        printed
c         -1     not printed      printed once   printed once
c
c  Modified:
c
c    08 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer KONTRL, the new value to be assigned the error 
c    control flag.
c
      implicit none

      integer j4save
      integer junk
      integer kontrl
      integer level
      integer nerr
      logical set
      integer value
      integer which
      character*8 xern1

      if ( kontrl .lt. -2 .or. 2 .lt. kontrl ) then

         write ( xern1, '(i8)' ) kontrl
         nerr = 1
         level = 2

         call xermsg ( 'XERROR', 'XSETUA',
     &      'Invalid value of KONTRL = ' // xern1, nerr, level )

        return

      end if

      which = 2
      value = kontrl
      set = .true.

      junk = j4save ( which, value, set )

      return
      end
      subroutine xsetua ( iunita, nunit )

c*********************************************************************72
c
cc XSETUA sets the unit numbers of the files receiving error messages.
c
c  Discussion:
c
c    This routine may be called to declare a list of up to five
c    logical units, each of which is to receive a copy of
c    each error message processed by this package.
c
c    The existence of multiple error output units makes it possible for
c    the user to ensure simultaneous printing of each error message to, 
c    say, a main output file, an interactive terminal, and other files 
c    such as graphics communication files.
c
c  Modified:
c
c    08 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer IUNIT(NUNIT), unit numbers to which the error messages
c    should be printed.  Normally these numbers should all be different
c    but duplicates are not prohibited.
c
c    Input, integer NUNIT, the number of unit numbers provided in IUNIT.
c    1 <= N <= 5.
c
      implicit none

      integer nunit

      integer i
      integer iunita(nunit)
      integer j4save
      integer junk
      integer level
      character*38 messg
      integer nerr
      logical set
      integer value
      integer which
      character*8 xern1

      if ( nunit .lt. 1 .or. nunit .gt. 5 ) then

         write ( xern1, '(i8)' ) nunit
         nerr = 1
         level = 2

         call xermsg ( 'XERROR', 'XSETUA',
     &      'Invalid number of units, NUNIT = ' // xern1, nerr, level )

         return

      end if
c
c  Set the main error output unit.
c
      which = 3
      value = iunita(1)
      set = .true.

      junk = j4save ( which, value, set )
c
c  If 1 < NUNIT, set auxilliary output units.
c
      do i = 2, nunit

        which = i + 4
        value = iunita(i)
        set = .true.

        junk = j4save ( which, value, set )

      end do

      which = 5
      value = nunit
      set = .true.

      junk = j4save ( which, value, set )

      return
      end
      subroutine xsetun ( iunit )

c*********************************************************************72
c
cc XSETUN sets the unit number of the main error output file.
c
c  Discussion:
c
c    This routine sets the unit number associated with the main error
c    output file.  If auxilliary error output units were defined, 
c    this routine suppresses them, as well.
c
c    Note that this error package initializes this unit number to the standard
c    output file, a reasonable choice.
c
c    Common choices for the unit number to be associated with the main
c    error output file might be 1, 6 or 9.  FORTRAN generally requires 
c    that such unit numbers be between 1 and 99.  Depending on the 
c    compiler, units -1 or 0 might also be legal choices.  It may 
c    be the case that some unit numbers cannot be chosen, because 
c    they are reserved by the compiler for other purposes.  In 
c    particular, some compilers reserve unit 5 for input.
c
c    Copies of the error output may also be sent to up to four auxilliary 
c    units, which can be defined by calling XSETUA.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer IUNIT, the unit number of the main error file.
c
      implicit none

      integer iunit
      integer j4save
      integer junk
      integer nunit
      logical set
      integer value
      integer which
c
c  Set the main error output unit.
c
      which = 3
      value = iunit
      set = .true.

      junk = j4save ( which, value, set )
c
c  Suppress all the other error output units.
c
      nunit = 1

      which = 5
      value = nunit
      set = .true.

      junk = j4save ( which, value, set )

      return
      end
