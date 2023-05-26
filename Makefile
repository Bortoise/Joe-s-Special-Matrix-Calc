#CFLAGS := -I/usr/local/include/ginac/ -fPIC # do not include -fPIC if static library
LDFLAGS := -Wl,-Bstatic -lginac -lcln -lgmp -static-libstdc++ # this is more portable, but idk if it will work for all installs

LieAlgebra.so: LieAlgebra.cpp
	$(CXX) $(CFLAGS) -c -shared -o $@ $^ $(LDFLAGS)

test.o: test.cpp LieAlgebra.so
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)


#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .o and remove -shared flag