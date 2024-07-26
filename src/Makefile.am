bin_PROGRAMS = trf
trf_SOURCES = trf.c trfrun.c trfclean.c tr30dat.c indexlist.c
trf_DEPENDENCIES = tr30dat.h trf.c trfclean.h trfrun.h indexlist.h
trf_CFLAGS = -fno-align-functions -fno-align-loops -Wall -Wextra -pedantic -O2

install-exec-hook:
	$(LN_S) $(DESTDIR)$(bindir)/trf$(EXEEXT) \
	$(DESTDIR)$(bindir)/trf@PACKAGE_VERSION@.@TARGET_NAME@.exe

uninstall-hook:
	rm $(DESTDIR)$(bindir)/trf@PACKAGE_VERSION@.@TARGET_NAME@.exe
