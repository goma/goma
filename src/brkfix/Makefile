#/************************************************************************ *
#* Goma - Multiphysics finite element software                             *
#* Sandia National Laboratories                                            *
#*                                                                         *
#* Copyright (c) 2014 Sandia Corporation.                                  *
#*                                                                         *
#* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
#* the U.S. Government retains certain rights in this software.            *
#*                                                                         *
#* This software is distributed under the GNU General Public License.      *
#\************************************************************************/

# Makefile for brk/fix, a tool to add nodes to an EXODUS FE database

OBJ_DIR := ${addprefix ../, $(BRKFIX_OBJ_DIR)}

RM = rm -f

CC = gcc

ACCESS = /projects/seacas/linux_rhel6/current

_INC = -I../$(INCLUDE_DIR)/brkfix -I$(ACCESS)/inc

INCLUDES = $(_INC) -I.

CDEBUGFLAGS = -g -Wall 

COPTIMIZEFLAGS = -O -Wall

DEFINES = -DHAVE_CONFIG_H=1 -DBUILD_LIB
# DEFINES = -DHAVE_CONFIG_H=1 -DPOTATO_CHIP_DIET

CFLAGS = $(CDEBUGFLAGS) $(COPTIMIZEFLAGS) $(CCOPTIONS) $(DEFINES) $(INCLUDES)

_LIB = -L$(ACCESS)/lib

LIB_CHACO = -lchaco

LIB_NETCDF = -lnetcdf -lhdf5_hl -lhdf5

LIB_EXODUS = -lexoIIv2c

LIB_SYS = -lm -lz -ldl

LIBS = $(_LIB) $(LIB_EXODUS) $(LIB_NETCDF) $(LIB_CHACO) $(LIB_SYS)

FOBJS = fix.o rd_exo.o wr_exo.o rd_dpi.o wr_dpi.o utils.o \
	aalloc.o eh.o bbb.o exo_utils.o

OBJS = brk.o brk_exo_file.o fix_exo_file.o bbb.o rd_exo.o wr_exo.o rd_in.o wr_graph_file.o mk_dm.o \
	wr_coords.o ppi.o sam_perea.o wr_dpi.o \
	utils.o aalloc.o eh.o exo_utils.o exo_conn.o emuck.o \
	rd_dpi.o

FOBJS := ${addprefix $(OBJ_DIR)/, $(FOBJS)}
OBJS := ${addprefix $(OBJ_DIR)/, $(OBJS)}

LIBRARY = libbrkfix.a

all: brk fix $(LIBRARY)

fix: $(FOBJS)
	$(RM) $@;\
	$(CC) $(COPTIMIZEFLAGS) -o $@ $(FOBJS) $(LIBS)

brk: $(OBJS)
	$(RM) $@;\
	$(PURE) $(CC) $(COPTIMIZEFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: library

library: $(LIBRARY)

$(OBJ_DIR)/%.o: %.c
	@mkdir -p $(OBJ_DIR)
	$(CC) -c $(CFLAGS) -o $@ $^

$(LIBRARY): $(FOBJS) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?; \
	$(RANLIB) $@

clean:
	-rm -f *.o *~ tmp* core junk* *.BAK *%

realclean:
	-rm -f *.fem *.pcon *.psplot *.o *~ tmp* junk* *.BAK *%
