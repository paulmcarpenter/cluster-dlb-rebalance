BINARIES=fakeruntime

all:${BINARIES}

%:%.cpp
	g++ --std=c++11 -lpthread -o $@ $<

clean:
	rm -rf ${BINARIES}
