#CFLAGS := -I/usr/local/include/ginac/ -fPIC # do not include -fPIC if static library
LDFLAGS := -Wl,-Bstatic -lginac -lcln -lgmp -static-libstdc++ # this is more portable, but idk if it will work for all installs
SOURCES=$(src/)
HEADERS=$(src/headers/)

#Makefile Targets
library.so: src/utils.cpp
	$(CXX) $(CFLAGS) -c -shared -o $@ $^ $(LDFLAGS)
	mv $@ out/$@

test.o: src/test.cpp out/library.so
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	mv $@ out/$@


#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .o and remove -shared flag