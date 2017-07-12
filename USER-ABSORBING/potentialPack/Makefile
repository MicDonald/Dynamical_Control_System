.DEFAULT_GOAL := LIB

CC_PATH = /opt/gcc-5.2.0
CC = $(CC_PATH)/bin/g++
CFLAGS = -std=c++14 -Wall -Wextra
OBJDIR = ../build
LIBDIR = ..
LIBFILE = libpotentialPack.a

CPPFILES = Math_Triangular.cpp \
	Potential.cpp \
	PotentialPair.cpp \
	PotentialAngle.cpp \
	PotentialBondAngle.cpp \
	PotentialDihedral.cpp \
	PotentialManybody.cpp \
	Potential_HarmonicBond.cpp \
	Potential_LennardJones.cpp \
	Potential_Buckingham.cpp \
	Potential_Morse.cpp \
	Potential_HarmonicAngle.cpp \
	Potential_HarmonicCosineAngle.cpp \
	Potential_CosineAngle.cpp \
	Potential_QuarticAngle.cpp \
	Potential_ScreenedHarmonic.cpp \
	Potential_Compass.cpp \
	Potential_HarmonicDihedral.cpp \
	Potential_CosineDihedral.cpp \
	Potential_HarmonicCosineDihedral.cpp \
	Potential_StillingerWeber.cpp 
#	Potential_Tersoff.cpp 

OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(CPPFILES))


$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o
	rm -f $(OBJDIR)/console

LIB: $(OBJDIR) $(LIBDIR) $(OBJS)
	ar rvs $(LIBDIR)/$(LIBFILE) $(OBJS)


test: $(OBJDIR)/main.o $(OBJS)
	$(CC) $(OBJDIR)/main.o $(OBJS) -Wl,--rpath=$(CC_PATH)/lib64 -o $(OBJDIR)/console

