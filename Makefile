# get the type of OS currently running
OS=$(shell uname -s)
PWD=$(shell pwd)

INC         = -Isrc
LIBS        = -L./lib -lQuartic
DEFS        =
STATIC_EXT  = .a
DYNAMIC_EXT = .so
AR          = ar rcs
LDCONFIG    = sudo ldconfig

WARN=-Wall -Wno-sign-compare
#-Weverything -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command 

# default values

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  LIBS     = -static -L./lib -lQuartic
  CXXFLAGS = -std=c++11 $(WARN) -O3 -fPIC
  AR       = ar rcs
  LDCONFIG = sudo ldconfig
endif

# check if the OS string contains 'MINGW'
ifneq (,$(findstring MINGW, $(OS)))
  LIBS     = -static -L./lib -lQuartic
  CXXFLAGS = -std=c++11 $(WARN) -O3
  AR       = ar rcs
  LDCONFIG = sudo ldconfig
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN        = -Wall -Weverything -Wno-sign-compare -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command 
  LIBS        = -L./lib -lQuartic
  CXXFLAGS    = $(WARN) -O3 -fPIC
  AR          = libtool -static -o
  LDCONFIG    =
  DYNAMIC_EXT = .dylib
endif

LIB_QUARTIC = libQuartic

SRCS = \
src/PolynomialRoots-1-Quadratic.cc \
src/PolynomialRoots-2-Cubic.cc \
src/PolynomialRoots-3-Quartic.cc \
src/PolynomialRoots-Jenkins-Traub.cc \
src/PolynomialRoots-Utils.cc

OBJS  = $(SRCS:.cc=.o)
DEPS  = src/PolynomialRoots-Utils.hh src/PolynomialRoots.hh
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Quartic

all: bin

travis: bin

bin: lib
	@$(MKDIR) bin
	$(CXX) $(INC) $(CXXFLAGS) -o bin/check_1_quadratic test/check_1_quadratic.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/check_2_cubic     test/check_2_cubic.cc     $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/check_3_quartic   test/check_3_quartic.cc   $(LIBS)

lib: lib/$(LIB_QUARTIC)$(STATIC_EXT) lib/$(LIB_QUARTIC)$(DYNAMIC_EXT)

include_local:
	@rm -rf lib/include
	$(MKDIR) lib
	$(MKDIR) lib/include
	@cp -f src/*.hh lib/include

src/%.o: src/%.cc $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@ 

src/%.o: src/%.c $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/libQuartic.a: $(OBJS) include_local
	@$(MKDIR) lib
	$(AR) lib/libQuartic.a $(OBJS) 

lib/libQuartic.dylib: $(OBJS) include_local
	@$(MKDIR) lib
	$(CXX) -shared -o lib/libQuartic.dylib $(OBJS) 

lib/libQuartic.so: $(OBJS) include_local
	@$(MKDIR) lib
	$(CXX) -shared -o lib/libQuartic.so $(OBJS) 

install: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include
	cp src/*.hh                $(PREFIX)/include
	cp lib/$(LIB_QUARTIC).*   $(PREFIX)/lib
	@$(LDCONFIG) $(PREFIX)/lib

install_as_framework: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/*.hh           $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_QUARTIC) $(PREFIX)/lib

run:
	./bin/check_1_quadratic
	./bin/check_2_cubic
	./bin/check_3_quartic

doc:
	doxygen

clean:
	rm -f lib/include/* lib/libQuartic.* src/*.o
	rm -rf bin
