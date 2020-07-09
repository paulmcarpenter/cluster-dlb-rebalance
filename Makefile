BINARIES=fakeruntime

all:${BINARIES}

%:%.cpp
	g++ -lpthread -o $@ $<

clean:
	rm -rf ${BINARIES}
