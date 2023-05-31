CFLAGS := -fPIC -std=c++17# do not include -fPIC if static library
LDFLAGS := -Wl,-Bstatic -lginac -lcln -lgmp -static-libstdc++ # this is more portable, but idk if it will work for all installs
SOURCES=$(src/)
HEADERS=$(src/headers/)

#Makefile Targets
out/%.o: src/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

library.a: out/utils.o out/lie_algebra.o out/lin_alg.o
# $(CXX) $(CFLAGS) -c -shared -o $@ $^ $(LDFLAGS)
# mv $@ out/$@
	ar rcs $@ $^ # We need to change this to get a shared library if matlab wants it

test.o: src/test.cpp library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	mv $@ out/$@

.PHONY: all
all: out/utils.o out/lin_alg.o out/lie_algebra.o library.a

#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .a and remove -shared flag