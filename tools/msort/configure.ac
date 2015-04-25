AC_INIT(msort, 1.1-1)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE

AC_CANONICAL_HOST
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_STDC_HEADERS

is_static=0
case ${prefix} in
  NONE);;
  *) is_static=1
     LDFLAGS="-L${prefix}/lib"
     CPPFLAGS="-I${prefix}/include";;
esac

# set CFLAGS and LDFLAGS

CFLAGS="-g -O2 -Wall -W -m64"
LDFLAGS="-m64 -lm"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
