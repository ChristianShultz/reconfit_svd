dnl Balint Joo, 13/12/2002
dnl George T. Fleming, 03/03/2003

dnl PAC_ADAT_LINK_CXX_FUNC(
dnl   ADAT_CXXFLAGS,
dnl   ADAT_LDFLAGS,
dnl   ADAT_LIBS,
dnl   ADAT_VARS,
dnl   ADAT_FUNC,
dnl   [action if working],
dnl   [action if not working]
dnl )
dnl
dnl  ADAT_CXXFLAGS for the necessary includes paths (-I)
dnl  ADAT_LDFLAGS  for the necessary library search paths (-L)
dnl  ADAT_LIBS     for the libraries (-l<lib> etc)
dnl  ADAT_VARS     for the declaration of variables needed
dnl                   to call ADAT_FUNC code fragment
dnl  ADAT_FUNC     for the body of a QDP++ function call or even general code
dnl                 fragment on which to run a compile/link test.
dnl                 If ADAT_VARS and ADAT_FUNC are empty, a basic test
dnl                 of compiling and linking a ADAT program is run.
dnl
AC_DEFUN(
  PAC_ADAT_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl - so input variables can be CXXFLAGS, LDFLAGS or LIBS
    pac_ADAT_CXXFLAGS="$1"
    pac_ADAT_LDFLAGS="$2"
    pac_ADAT_LIBS="$3"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_ADAT_CXXFLAGS"
    LDFLAGS="$LDFLAGS $pac_ADAT_LDFLAGS"
    LIBS="$LIBS $pac_ADAT_LIBS"
    AC_TRY_LINK(
      [
        #include <libxml/xmlmemory.h>
	#include <libxml/parser.h>
      ], [
        int argc ; char **argv ;
        xmlDocPtr doc;
	char *docname="foo";	
	doc = xmlParseFile(docname);
        $4 ;
        $5 ;
      ],
      [pac_libxml2_working=yes],
      [pac_libxml2_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_libxml2_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)
