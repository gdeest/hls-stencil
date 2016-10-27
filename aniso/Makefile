APPSOURCES = \
src/main.cpp \
src/inner_tile.cpp \
src/perf_counter.cpp \

EXECUTABLE = jacobi2d.elf

PLATFORM = -sds-pf zc706
# PLATFORM = -sds-pf zybo
VERBOSE  = 
OLEVEL   = -O3

# SDSFLAGS = \
# ${PLATFORM}


SDSFLAGS = \
${PLATFORM} \
-sds-hw inner_tiles inner_tile.cpp -sds-end\
${VERBOSE}

INC=-I$(XILINX_SDX)/include -I$(XILINX_SDX)/Vivado_HLS/include/

CPP    = sds++ $(SDSFLAGS)
CFLAGS = -Wall ${OLEVEL} -c
CFLAGS += -MMD -MP -MF"$(@:%.o=%.d)"
LFLAGS = ${OLEVEL} -poll-mode 1

CXX    = g++

OBJECTS := $(APPSOURCES:.cpp=.o)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all

all: test ${EXECUTABLE}

test.exe: ${APPSOURCES}
	${CXX} -o test.exe ${INC} ${APPSOURCES} -lm -lpthread

test: test.exe
	./test.exe 512 512

${EXECUTABLE}: ${OBJECTS} ${INCSOURCES}
	${CPP} ${LFLAGS} ${OBJECTS} -o $@

-include ${DEPS}

%.o: %.cpp
	${CPP} ${CFLAGS} $< -o $@

clean:
	${RM} ${EXECUTABLE} ${OBJECTS} *.d

ultraclean: clean
	${RM} ${EXECUTABLE}.bit ${EXECUTABLE}
	${RM} -rf _sds sd_card
