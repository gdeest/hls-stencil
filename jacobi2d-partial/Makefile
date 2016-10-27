APPSOURCES = \
src/main.cpp \
src/inner_tile.cpp \
src/perf_counter.cpp

EXECUTABLE = jacobi2d.elf

PLATFORM = -sds-pf zc706
# PLATFORM = -sds-pf zybo
VERBOSE  = -verbose
OLEVEL   = -O3

# SDSFLAGS = \
# ${PLATFORM} \
# ${VERBOSE}


SDSFLAGS = \
${PLATFORM} \
-sds-hw inner_tiles inner_tile.cpp -sds-end\
${VERBOSE}

INCDIR=/run/media/huginn/b85bf520-705c-4507-b6d8-262e11ce3f5b/soft/SDSoC/2016.2/Vivado_HLS/2016.2/include

CPP    = sds++ $(SDSFLAGS)
CFLAGS = -Wall ${OLEVEL} -c
CFLAGS += -MMD -MP -MF"$(@:%.o=%.d)"
LFLAGS = ${OLEVEL} -poll-mode 1

CXX    = g++

OBJECTS := $(APPSOURCES:.cpp=.o)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all

all: ${EXECUTABLE}

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
