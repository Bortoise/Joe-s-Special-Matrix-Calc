CFLAGS := -I/usr/local/include/ginac/ -fPIC # do not include -fPIC if static library
LDFLAGS := /usr/local/lib/libginac.a /usr/local/lib/libcln.a -Wl,-Bstatic -lgmp

LieAlgebra.so: LieAlgebra.cpp
	$(CXX) $(CFLAGS) -c -shared -o $@ $^ $(LDFLAGS)

#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .o and remove -shared flag