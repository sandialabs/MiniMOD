# Compiler
CC = mpicc

# Compiler flags
CFLAGS = -Wall -O2 -lm -g

# Source files
SRCS = all_to_all_rma.c broadcast_rma.c allgather_rma.c earlycoll.c

# Header files
HDRS = earlycoll.h

# Object directory
OBJDIR = build

# Results directory
RESULTSDIR = results

# Object files
OBJS = $(SRCS:%.c=$(OBJDIR)/%.o)

# Executable names
DRIVER_EXEC = earlycoll
BENCHMARK_EXEC = benchmark
BENCH_ALLGATHER_EXEC = benchmark_allgather
TEST_ALLTOALL_EXEC = test_alltoall 
TEST_ALLGATHER_EXEC = test_allgather 
TEST_BROADCAST_EXEC = test_broadcast

# Check for persistent collectives support
PERSISTENT_COLL_SUPPORT = $(shell ./build_scripts/check_persistent_collectives.sh 2>/dev/null)

# Conditionally define PERSISTENT_COLL based on the script output
ifeq ($(PERSISTENT_COLL_SUPPORT),true)
    CFLAGS += -DPERSISTENT_COLL
endif

# Default target
all: $(DRIVER_EXEC) $(TEST_ALLTOALL_EXEC) $(TEST_ALLGATHER_EXEC) $(TEST_BROADCAST_EXEC) $(BENCHMARK_EXEC) $(BENCH_ALLGATHER_EXEC)

# Rule to build the driver executable
$(DRIVER_EXEC): $(OBJS) $(OBJDIR)/main.o
	$(CC) $(CFLAGS) -o $@ $^

# Rule to build the benchmark executable
$(BENCHMARK_EXEC): $(OBJS) $(OBJDIR)/benchmark2.o
	$(CC) $(CFLAGS) -o $@ $^

# Rule to build the benchmark executable
$(BENCH_ALLGATHER_EXEC): $(OBJS) $(OBJDIR)/benchmark_allgather.o
	$(CC) $(CFLAGS) -o $@ $^

# Rule to build the test_alltoall executable
$(TEST_ALLTOALL_EXEC): $(OBJS) $(OBJDIR)/test_alltoall.o
	$(CC) $(CFLAGS) -o $@ $^

# Rule to build the test_allgather executable
$(TEST_ALLGATHER_EXEC): $(OBJS) $(OBJDIR)/test_allgather.o
	$(CC) $(CFLAGS) -o $@ $^

# Rule to build the test_broadcast executable
$(TEST_BROADCAST_EXEC): $(OBJS) $(OBJDIR)/test_broadcast.o
	$(CC) $(CFLAGS) -o $@ $^

# Rule to build object files
$(OBJDIR)/%.o: %.c $(HDRS)
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target
clean:
	rm -rf $(OBJDIR) $(DRIVER_EXEC) $(BENCHMARK_EXEC) $(TEST_ALLTOALL_EXEC) $(TEST_ALLGATHER_EXEC) $(TEST_BROADCAST_EXEC) $(RESULTSDIR)

# Phony targets
.PHONY: all clean

