# ---------------------------------------------------------------------------- #
#                                                                              #
# LICENSE NOTICE.                                                              #
#                                                                              #
#     LICENSE: GPL-3.0                                                         #
#                                                                              #
#     IMPORTANT: THIS IS FREE SOFTWARE WITHOUT ANY WARRANTY. THE USER IS FREE  #
#                TO MODIFY AND REDISTRIBUTE THIS SOFTWARE UNDER THE TERMS OF   #
#                THE LICENSE LISTED ABOVE PUBLISHED BY THE FREE SOFTWARE       #
#                FOUNDATION. THE PUBLISHER, SIMON EUCHNER, IS NOT RESPONSIBLE  #
#                FOR ANY NEGATIVE EFFECTS THIS SOFTWARE MAY CAUSE.             #
#                                                                              #
# ---------------------------------------------------------------------------- #
#                                                                              #
# Makefile for the EIGS C-libary                                               #
#                                                                              #
# ---------------------------------------------------------------------------- #


### Variables

# Compiler and Linker (GNU C compiler/linker)
CC = gcc -fPIC
LD = gcc -shared
FLAGS = -Wall -Wextra -pedantic -std=c99 -fPIC
OLVL = -O1

# Paths
SRC = ./src.d
OBJ = ./obj.d
LIB = ./lib.d

# Target (library *eigs*)
TARGET = ${LIB}/libeigs.so

# Define file names
F1  = eigs
F2  = zgeigsf
F3  = zgeigsa
F4  = dgeigsa
F5  = zheigsa
F6  = dseigsa
F7  = dgeigsf

# Specify object code files
OBJ_FILENAMES = ${F1}.o ${F2}.o ${F3}.o ${F4}.o ${F5}.o ${F6}.o ${F7}.o
OBJ_FILES = ${foreach file, ${OBJ_FILENAMES}, ${OBJ}/${file}}


### Print info
all: ${TARGET}
	@echo -e "\nBuilding process complete\n"


### Link
${TARGET}: ${OBJ_FILES}
	${LD} ${FLAGS} ${OLVL} -o ${LIB}/libeigs.so ${wildcard ${OBJ}/*.o}


### Compile

# eigs.c
${OBJ}/${F1}.o: ${SRC}/${F1}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F1}.o -c ${SRC}/${F1}.c

# zgeigsf.c
${OBJ}/${F2}.o: ${SRC}/${F2}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F2}.o -c ${SRC}/${F2}.c

# zgeigsa.c
${OBJ}/${F3}.o: ${SRC}/${F3}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F3}.o -c ${SRC}/${F3}.c

# dgeigsa.c
${OBJ}/${F4}.o: ${SRC}/${F4}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F4}.o -c ${SRC}/${F4}.c

# zheigsa.c
${OBJ}/${F5}.o: ${SRC}/${F5}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F5}.o -c ${SRC}/${F5}.c

# dseigsa.c
${OBJ}/${F6}.o: ${SRC}/${F6}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F6}.o -c ${SRC}/${F6}.c

# dgeigsf.c
${OBJ}/${F7}.o: ${SRC}/${F7}.c
	${CC} ${FLAGS} ${OLVL} -o ${OBJ}/${F7}.o -c ${SRC}/${F7}.c


### Cleanup

clean:
	rm ${OBJ}/*.o
	rm ${LIB}/libeigs.so

.PHONY: clean
