SOURCES=$(src/)
HEADERS=$(src/headers/)

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CFLAGS :=  -std=c++17 -arch x86_64
	LDFLAGS := -Wl,-v -lginac -lcln -lgmp -static-libstdc++
else
	CFLAGS :=  -std=c++17
	LDFLAGS := -Wl,-Bstatic -lginac -lcln -lgmp -static-libstdc++ # this is more portable, but idk if it will work for all installs
endif

#Makefile Targets
out/%.o: src/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

out/library.a: out/utils.o out/lin_alg.o out/lie_algebra.o
#	$(CXX) $(CFLAGS) -shared -o $@ $^ # $(LDFLAGS)
	ar rcs $@ $^

out/test.o: src/test.cpp out/library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

out/example.o: src/example.cpp out/library.so #out/library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
#mv $@ out/$@

.PHONY: all test clear
all: out/utils.o out/lin_alg.o out/lie_algebra.o out/library.so
test: out/test.o
example: out/example.o
clear:
	rm out/*
#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .a and remove -shared flag
