CC = gcc
CFLAGS = -O3 -Wall -mavx2 -mfma -flto -g -Wno-stringop-truncation

NVCC = nvcc
CUDA_ARCH ?= sm_89
NVCC_FLAGS = -arch=$(CUDA_ARCH) -lineinfo

SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

# FTYPE=FLOAT uses fp32
FTYPE ?= DOUBLE
# Deprecated variable, needs to be removed
VEC ?= EXPL
# TARGET=CUDA compiles with CUDA support
TARGET ?=
# TIMEIT=ON logs TIMEIT() calls
TIMEIT ?=

DEFINE = -D$(FTYPE) -D$(VEC)

ifeq ($(TIMEIT), ON)
	DEFINE += -DTIMEIT_ON
endif

INCLUDE = -I$(INC_DIR) -I$(SRC_DIR)
LIBS = -lm

SOLVER_OBJS = solver.o momentum.o pressure.o output.o thread-array.o
UNIT_TEST_OBJS = unit-test.o momentum-test.o pressure-test.o
CONVERGENCE_TEST_OBJS = $(SOLVER_OBJS) convergence-test.o

CUDA_OBJS = momentum.o

solver: mkdir-build $(BUILD_DIR)/solver
tests: mkdir-build $(BUILD_DIR)/convergence-test

mkdir-build:
	mkdir -p $(BUILD_DIR)/objs
ifeq ($(TARGET), CUDA)
	mkdir -p $(BUILD_DIR)/objs/cuda
endif

LINKER = $(CC)

SOLVER_OBJS_PATH = $(addprefix $(BUILD_DIR)/objs/, $(SOLVER_OBJS) main.o)
CONVERGENCE_TEST_OBJS_PATH = $(addprefix $(BUILD_DIR)/objs/, $(CONVERGENCE_TEST_OBJS))

ifeq ($(TARGET), CUDA)
	LINKER = $(NVCC)
	DEFINE += -DTARGET_CUDA
	CUDA_OBJS_PATH = $(addprefix $(BUILD_DIR)/objs/cuda/, $(CUDA_OBJS))
	SOLVER_OBJS_PATH += $(CUDA_OBJS_PATH)
	CONVERGENCE_TEST_OBJS_PATH += $(CUDA_OBJS_PATH)
endif

$(BUILD_DIR)/solver: $(SOLVER_OBJS_PATH)
	$(LINKER) $^ $(LIBS) -o $@

$(BUILD_DIR)/unit-test: $(addprefix $(BUILD_DIR)/objs/, $(UNIT_TEST_OBJS))
	$(CC) $^ $(LIBS) -o $@

$(BUILD_DIR)/convergence-test: $(CONVERGENCE_TEST_OBJS_PATH)
	$(LINKER) $^ $(LIBS) -o $@

$(BUILD_DIR)/convergence-pressure-test: $(BUILD_DIR)/objs/convergence-pressure-test.o
	$(CC) $^ $(LIBS) -o $@

$(BUILD_DIR)/objs/%.o: $(SRC_DIR)/%.c
	$(CC) -c $^ $(CFLAGS) $(INCLUDE) $(DEFINE) -o $@

$(BUILD_DIR)/objs/%.o: tests/%.c
	$(CC) -c $^ $(CFLAGS) $(INCLUDE) $(DEFINE) -o $@

$(BUILD_DIR)/objs/cuda/%.o: $(SRC_DIR)/cuda/%.cu
	$(NVCC) -c $^ $(INCLUDE) $(DEFINE) $(NVCC_FLAGS) -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
