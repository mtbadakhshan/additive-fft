CC=    gcc

CFLAGS= -O3 -fopenmp -funroll-loops -mavx2 -mpclmul -std=c23 -Wextra -Wall -I .

CANTOR_DIR = cantor
CANTOR_SRC = $(wildcard $(CANTOR_DIR)/*.c)
CANTOR_OBJ = $(CANTOR_SRC:.c=.o)

DYADIC_DIR = dyadic
DYADIC_SRC = $(filter-out $(DYADIC_DIR)/tune_dyadic.c,$(wildcard $(DYADIC_DIR)/*.c))
DYADIC_OBJ = $(DYADIC_SRC:.c=.o)

BITPOLYMUL_DIR = bitpolymul
BITPOLYMUL_SRC = bc.c bitpolymul.c encode.c gf2128_tower_iso.c butterfly_net.c gf16.c ska.c gf2256_tower_iso.c gf2128_cantor_iso.c btfy.c trunc_btfy_tab.c
BITPOLYMUL_OBJ = $(patsubst %.c, $(BITPOLYMUL_DIR)/%.o, $(BITPOLYMUL_SRC))

UTILS_DIR = utils
UTILS_SRC = $(wildcard $(UTILS_DIR)/*.c)
UTILS_OBJ = $(UTILS_SRC:.c=.o)

TARGET  = main.out

.PHONY: all clean tune

all: $(TARGET)
main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o
# Dyadic-only tuner (no LCH): needs dyadic_fft + Cantor iso table it links.
DYADIC_TUNE_DEPS = $(DYADIC_DIR)/dyadic_fft.o $(BITPOLYMUL_DIR)/gf2128_cantor_iso.o
tune.out: $(DYADIC_DIR)/tune_dyadic.c $(DYADIC_TUNE_DEPS)
	$(CC) $(CFLAGS) -o $@ $^ -lm
# gf2x-style: measure B on this machine, write dyadic_tune_params.h, rebuild
tune: tune.out
	./tune.out
	$(MAKE) $(DYADIC_DIR)/dyadic_fft.o $(TARGET)
$(CANTOR_DIR)/%.o: $(CANTOR_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@	
$(DYADIC_DIR)/%.o: $(DYADIC_DIR)/%.c $(DYADIC_DIR)/dyadic_tune_params.h
	$(CC) $(CFLAGS) -c $< -o $@	
$(BITPOLYMUL_DIR)/%.o: $(BITPOLYMUL_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(UTILS_DIR)/%.o: $(UTILS_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@	
$(TARGET): main.o $(CANTOR_OBJ) $(DYADIC_OBJ) $(BITPOLYMUL_OBJ) $(UTILS_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f main.o $(CANTOR_DIR)/*.o $(DYADIC_DIR)/*.o $(BITPOLYMUL_DIR)/*.o $(UTILS_DIR)/*.o $(TARGET) tune.out
