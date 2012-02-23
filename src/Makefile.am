#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/src \
               @CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ 

AM_LDFLAGS  = -L@top_builddir@/src \
               @LDFLAGS@ @ADAT_LDFLAGS@

LDADD = @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
#noinst_PROGRAMS = spectrum

bin_PROGRAMS = semble_vfit timeslice_covariance fit_prin_corr fit_multiple_princorrs fit_corr_exp fit_corr_cosh fit_corr_exp_and_const fit_corr_cosh_and_const dispersion_fit generate_ensem fit_Z_const_and_exp
noinst_PROGRAMS = 
check_PROGRAMS =
EXTRA_PROGRAMS = 


#
# The program and its dependencies
#

semble_vfit_SOURCES = semble_vfit.cc semble_load_correlators.cc semble_fit_ini_xml.cc semble_file_management.cc semble_histogram.cc fit_correlators.cc fit_forms.cc jackknife_fitter.cc ensem_data.cc plot.cc linear_algebra.cc semble_linear_algebra.cc

timeslice_covariance_SOURCES = timeslice_covariance.cc linear_algebra.cc ensem_data.cc

fit_prin_corr_SOURCES = fit_prin_corr.cc fit_correlators.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc

fit_multiple_princorrs_SOURCES = fit_multiple_princorrs.cc multi_princorr_fit.cc fit_correlators.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc

fit_corr_exp_SOURCES = fit_corr_exp.cc fit_correlators.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc  plot.cc fit_forms.cc

fit_corr_cosh_SOURCES = fit_corr_cosh.cc fit_correlators.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc

fit_corr_exp_and_const_SOURCES = fit_corr_exp_and_const.cc fit_correlators.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc

fit_corr_cosh_and_const_SOURCES = fit_corr_cosh_and_const.cc fit_correlators.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc

#fitZ_experiment_SOURCES = fitZ_experiment.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc

dispersion_fit_SOURCES = dispersion_fit.cc multi_ensemble_fitter.cc jackknife_linear_fitter.cc jackknife_fitter.cc ensem_data.cc plot.cc linear_algebra.cc semble_linear_algebra.cc

generate_ensem_SOURCES = generate_ensem.cc

fit_Z_const_and_exp_SOURCES = fit_Z_const_and_exp.cc jackknife_fitter.cc ensem_data.cc linear_algebra.cc plot.cc fit_forms.cc


#
# The following rule will always try to rebuild libadat.a when you 
# try to compile anything else
#
####DEPENDENCIES = build_lib ${lib_LIBRARIES}
DEPENDENCIES = ${lib_LIBRARIES}
${noinst_PROGRAMS}: ${DEPENDENCIES}
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
