#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/src \
               @CXXFLAGS@ @SEMBLE_CXXFLAGS@ @JACKFITTER_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/src \
               @LDFLAGS@ @SEMBLE_LDFLAGS@ @JACKFITTER_LDFLAGS@ @ADAT_LDFLAGS@ 

LDADD = -lreconfit @SEMBLE_LIBS@ @JACKFITTER_LIBS@ @ADAT_LIBS@  @ITPP_LIBS@ @LIBS@

#
# Library
#
lib_LIBRARIES = libreconfit.a

libreconfit_a_SOURCES = \
    correlator_reader_factory.cc \
    correlator_util.cc \
    redstar_reader.cc \
    redstarSUN_reader.cc \
    colorvec_db_reader.cc \
    semble_fit_ini_xml.cc

#    semble_load_correlators.cc

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
#noinst_PROGRAMS = spectrum

bin_PROGRAMS = semble_vfit
noinst_PROGRAMS = 
check_PROGRAMS =  effective_mass_onecorr recon_util proj_op_util dispersion_fit 
EXTRA_PROGRAMS = 


#
# The program and its dependencies
#

semble_vfit_SOURCES = semble_vfit.cc

effective_mass_onecorr_SOURCES = effective_mass_onecorr.cc 

recon_util_SOURCES = recon_util.cc 

proj_op_util_SOURCES = proj_op_util.cc

dispersion_fit_SOURCES = dispersion_fit.cc multi_ensemble_fitter.cc

#
# The following rule will always try to rebuild libreconfit.a when you 
# try to compile anything else
#
DEPENDENCIES = ${lib_LIBRARIES}
${noinst_PROGRAMS}: ${DEPENDENCIES}
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}
