CC=    gcc

CFLAGS= -O3 -funroll-loops -mavx2 -mpclmul -std=c23 -Wextra -Wall -I .

UTILS_DIR = utils
UTILS_SRC = $(wildcard $(UTILS_DIR)/*.c)
UTILS_OBJ = $(UTILS_SRC:.c=.o)

BITPOLYMUL_DIR = bitpolymul
BITPOLYMUL_SRC = bc.c bitpolymul.c encode.c gf2128_tower_iso.c butterfly_net.c gf16.c ska.c gf2256_tower_iso.c gf2128_cantor_iso.c btfy.c trunc_btfy_tab.c
BITPOLYMUL_OBJ = $(patsubst %.c, $(BITPOLYMUL_DIR)/%.o, $(BITPOLYMUL_SRC))

TARGET  = main

.PHONY: all clean

all: $(TARGET)
main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o
$(BITPOLYMUL_DIR)/%.o: $(BITPOLYMUL_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(UTILS_DIR)/%.o: $(UTILS_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@	
$(TARGET): main.o $(BITPOLYMUL_OBJ) $(UTILS_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f main.o $(BITPOLYMUL_DIR)/*.o $(UTILS_DIR)/*.o $(TARGET)