#CFLAGS := -I/usr/local/include/ginac/ -fPIC # do not include -fPIC if static library
LDFLAGS := /usr/local/lib/libginac.a /usr/local/lib/libcln.a -Wl,-Bstatic -lgmp #replace the path with the path on your machine
#LDFLAGS := -Wl,-Bstatic -lgmp -lginac -lcln # this is more portable, but idk if it will work for all installs

LieAlgebra.so: LieAlgebra.cpp
	$(CXX) $(CFLAGS) -c -shared -o $@ $^ $(LDFLAGS)

test.o: test.cpp LieAlgebra.so
	$(CXX) $(CFLAGS) -o $@ $^ -Wl,-rpath,LieAlgebra.so


#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .o and remove -shared flag