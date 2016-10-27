APPSOURCES = \
src/main.cpp \
src/inner_tile.cpp \
src/perf_counter.cpp \
src/outer_tile.cpp \
src/utils.cpp

INCSOURCES = \
src/compute.cpp \
src/read_inputs.cpp \
src/read_coords.cpp \
src/write_outputs.cpp

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

all: tests ${EXECUTABLE}

validate_sw: ${APPSOURCES}
	${CXX} -o validate_sw -I${INCDIR} -DSW -DVALIDATE ${APPSOURCES} -lm -lpthread
validate: ${APPSOURCES}
	${CXX} -o validate -I${INCDIR} -DSW -DVALIDATE -DSIMULATE ${APPSOURCES} -lm -lpthread
simulate_debug: ${APPSOURCES}
	${CXX} -o simulate_debug -I${INCDIR} -DSIMULATE -DDEBUG ${APPSOURCES} -lm -lpthread
simulate: ${APPSOURCES}
	${CXX} -o simulate -I${INCDIR} -DSIMULATE ${APPSOURCES} -lm -lpthread

tests: validate_sw validate simulate
	./validate_sw 500 500 500 > validate_sw.log
	./validate 500 500 500 > validate.log
	# ./simulate_debug 50 50 200 > simulate_debug.log
	./simulate 500 500 500 > simulate.log

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
