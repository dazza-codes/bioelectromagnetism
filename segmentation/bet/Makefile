include $(FSLCONFDIR)/default.mk

PROJNAME = bet

LIBS = -lss_$$FDT -lm

XFILES = bet_8UI bet_8SI bet_16UI bet_16SI bet_32UI bet_32SI bet_32R

RUNAVWS = bet

RUNTCLS = Bet

SCRIPTS = betfunc betpremask

all:	bet_8UI

bet_8UI: bet.c
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R ; do $(DATATYPES) ; \
        $(CC) $(CFLAGS) -DFDT="$$STR" -o bet_$$FDT bet.c $(LDFLAGS) $(LIBS) ; \
	done
