# -------------- #
# -- Makefile -- #
# -------------- #

# -- Lile list ----------

FILE_NRC    = nrutil.c 
FILE_MOTION = sigmadelta.c motion.c motion_test.c
FILE_MORPHO = morpho_max.c morpho_min.c morpho_ouverture.c morpho_test.c
FILE_MAIN   = main.c
FILE_SWP    = swp.c swp_test.c
FILE = $(FILE_MAIN) $(FILE_NRC) $(FILE_SWP) $(FILE_MORPHO) $(FILE_MOTION)

# -- paths ----------
SRC_PATH = src
OBJ_PATH = obj
EXE_PATH = exe
INC_PATH = include

# -- Macros ----------
CC = gcc

# -- Flags ----------
C_OPTIMISATION_FLAGS = -std=c99 -O3 -fstrict-aliasing -Wno-comment
C_OPTIMISATION_FLAGS = -std=c99 -O0 -fstrict-aliasing -Wno-comment
C_ARCH_FLAGS = -mtune=native -march=native
C_ARCH_FLAGS =
C_INC_FLAGS = -I$(INC_PATH)

CFLAGS  = $(C_OPTIMISATION_FLAGS) $(C_ARCH_FLAGS) $(C_INC_FLAGS) $(LIB_INC_PATH)
LDFLAGS = $(C_OPTIMISATION_FLAGS) $(C_ARCH_FLAGS) $(C_INC_FLAGS) $(LIB_LIB_PATH)

# -- Final product ----------
PRODUCT   = motion_acc.exe

# -- src and obj List ----------
SRC = $(addprefix ${SRC_PATH}/, $(FILE))
OBJ = $(addprefix ${OBJ_PATH}/, $(addsuffix .o, $(basename $(FILE))))

# -- Base rules ----------
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
   
#-----Main rule ----------
$(EXE_PATH)/$(PRODUCT): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(INC) $(LIB_LIB_PATH) -lm

# -- Other stuff ----------
clean:
	rm -f $(OBJ)
	rm -f $(PRODUCT)
	#rm -f ${LIB_PATH}/${PRODUCT}

