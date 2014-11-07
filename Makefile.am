##  Include the m4 directory for m4 macros; required when using autoreconf
ACLOCAL_AMFLAGS = -I m4

##  Compiler
CC = $(MPICC)

##  Source flags
bin_PROGRAMS = plsa
plsa_SOURCES = \
  plsa-defn.h \
  em-steps.h \
  em-steps.c \
  input.h \
  output.h \
  wmalloc.h \
  input.c \
  main.c \
  output.c \
  wmalloc.c \
  run.h \
  run.c \
  parameters.h \
  parameters.c \
  debug.h \
  debug.c \
  comm.h \
  comm.c
plsa_LDADD = -lm

dist_doc_DATA = AUTHORS COPYING ChangeLog INSTALL NEWS README VERSION \
                docs/manual.pdf

