# Makefile for MinGW on Windows
CC  = c:/JTSDK/Qt55/Tools/mingw492_32/bin/gcc
FC = c:/JTSDK/Qt55/Tools/mingw492_32/bin/gfortran
CXX = c:/JTSDK/Qt55/Tools/mingw492_32/bin/g++

FFLAGS = -O2 -fbounds-check -Wall -Wno-conversion -Wno-precision-loss \
       -fno-second-underscore -fbacktrace \
       -ffpe-summary=none
#       -ffpe-trap=invalid,zero,overflow,underflow,denormal
CFLAGS = -I. -fbounds-check

# Default rules
%.o: %.c
	${CC} ${CFLAGS} -c $<
%.o: %.f
	${FC} ${FFLAGS} -c $<
%.o: %.F
	${FC} ${FFLAGS} -c $<
%.o: %.f90
	${FC} ${FFLAGS} -c $<
%.o: %.F90
	${FC} ${FFLAGS} -c $<

all:    packjt.o wavhdr.o libm65.a m65.exe mapsim.exe libastro.a

OBJS1 = trimlist.o display.o getdphi.o pctile.o ccf65.o \
	decode1a.o sort.o filbig.o fil6521.o afc65b.o \
	twkfreq.o decode65b.o indexx.o ssort.o fchisq.o setup65.o \
	extract.o deep65.o ccf2.o demod64a.o chkhist.o graycode.o \
	interleave63.o encode65.o igray.o set.o shell.o qra64_subs.o \
	grid2k.o getpfx2.o q65b.o qra64c.o twkfreq_xy.o qra64zap.o \
	deg2grid.o getpfx1.o k2grid.o ftrsd2.o graycode65.o \
	wrapkarn.o nchar.o init_rs.o encode_rs.o decode_rs.o \
	four2a.o rfile3a.o grid2deg.o pfxdump.o dpol.o sync64.o spec64.o \
	astro.o tm2.o sun.o moondop.o coord.o tmoonsub.o \
	geocentric.o moon2.o toxyz.o dot.o dcoord.o f77_wisdom.o \
	gen65.o chkmsg.o ptt.o astrosub.o astro0.o recvpkt.o symspec.o \
	iqcal.o iqfix.o timf2.o s3avg.o packjt.o badmsg.o fmtmsg.o \
	qracodes.o qra64.o smo.o averms.o lorentzian.o npfwht.o pdmath.o \
	qra13_64_64_irr_e.o fchisq0.o gen_q65_wave.o wavhdr.o

packjt.o: packjt.f90
	$(FC) -c packjt.f90

libm65.a: $(OBJS1)
	ar cr libm65.a $(OBJS1) 
	ranlib libm65.a

qra64_subs.o: ./qra64/qra64_subs.c
	gcc -c -O2 -o qra64_subs.o ./qra64/qra64_subs.c

qracodes.o: qracodes/qracodes.c
	gcc -c -O2 -o qracodes.o qracodes/qracodes.c

qra64.o: qra64/qra64.c
	gcc -c -O2 -o qra64.o qra64/qra64.c

qra13_64_64_irr_e.o: qracodes/qra13_64_64_irr_e.c
	gcc -c -O2 -o qra13_64_64_irr_e.o  \
	qracodes/qra13_64_64_irr_e.c

npfwht.o: qracodes/npfwht.c
	gcc -c -O2 -o npfwht.o qracodes/npfwht.c

pdmath.o: qracodes/pdmath.c
	gcc -c -O2 -o pdmath.o qracodes/pdmath.c

OBJS3 = m65.o m65a.o map65a.o symspec.o decode0.o ftninit.o ftnquit.o \
	timer.o ipcomm.o sec_midn.o usleep.o
# cutil.o 

LIBS3 = -L'C:/JTSDK/Qt55/5.5/mingw492_32/lib' -lQt5Core

m65.exe: $(OBJS3) libm65.a
	$(CXX) -o m65.exe $(OBJS3) $(LIBS3) libm65.a ../libfftw3f_win.a \
	-lgfortran
#	cp m65.exe ../../map65_install
	cp m65.exe /c/MAP65_2.9

OBJS7 = astrosub.o astro0.o astro.o tm2.o grid2deg.o sun.o moondop.o \
	coord.o dot.o moon2.o tmoonsub.o toxyz.o geocentric.o \
	dcoord.o

libastro.a: $(OBJS7)
	ar cr libastro.a $(OBJS7)
	ranlib libastro.a

OBJS6 = mapsim.o cgen65.o noisegen.o gran.o
mapsim.exe: $(OBJS6) libm65.a
	$(FC) -o mapsim.exe $(OBJS6) libm65.a

INCPATH = -I'C:/JTSDK/Qt55/5.5/mingw492_32/include/QtCore' \
	-I'C:/JTSDK/Qt55/5.5/mingw492_32/include' 

ipcomm.o: ipcomm.cpp
	$(CXX) -c $(INCPATH) ipcomm.cpp

#m65a.o:	m65a.f90
#	$(FC) -c -fno-second-underscore -cpp m65a.f90

sec_midn.o: sec_midn.f90
	$(FC) -c -fno-second-underscore sec_midn.f90

#symspec.o: ../symspec.f90
#	$(FC) -c $(FFLAGS) -o symspec.o ../symspec.f90

OBJS4 = tastro.o astro0.o libm65.a
tastro.exe: $(OBJS4)
	$(FC) $(FFLAGS) -o tastro.exe $(OBJS4) libm65.a 

OBJS5 = t1.o timer.o libm65.a
t1.exe: $(OBJS5)
	$(FC) $(FFLAGS) -o t1.exe $(OBJS5) libm65.a 

#astro0.o: ../astro0.f90
#	$(FC) -c $(FFLAGS) -o astro0.o ../astro0.f90

init_rs.o: init_rs.c
	$(CC) -c -DBIGSYM=1 -o init_rs.o init_rs.c

encode_rs.o: encode_rs.c
	$(CC) -c -DBIGSYM=1 -o encode_rs.o encode_rs.c

decode_rs.o: decode_rs.c
	$(CC) -c -DBIGSYM=1 -o decode_rs.o decode_rs.c

.PHONY : clean

clean:
	rm -f *.o *.mod libm65.a m65.exe
