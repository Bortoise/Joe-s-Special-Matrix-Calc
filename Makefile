CFLAGS :=  -std=c++17 #-fPIC # do not include -fPIC if static library
LDFLAGS := -lginac -lcln -lgmp -static-libstdc++ # this is more portable, but idk if it will work for all installs

UNAME := $(shell uname)
ifneq ($(UNAME),Darwin)
	LDFLAGS += -Wl,-Bstatic
endif

SOURCES=$(src/)
HEADERS=$(src/headers/)

#Makefile Targets
out/%.o: src/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

out/library.a: out/utils.o out/lie_algebra.o out/lin_alg.o
# $(CXX) $(CFLAGS) -c -shared -o $@ $^ $(LDFLAGS)
# mv $@ out/$@
	ar rvs $@ $^ 
# We need to change this to get a shared library if matlab wants it
#mv $@ out/$@

out/test.o: src/test.cpp out/library.a #out/library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
#mv $@ out/$@

.PHONY: all test clear
all: out/utils.o out/lin_alg.o out/lie_algebra.o out/library.a
test: out/test.o
clear:
	rm out/*
#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .a and remove -shared flag