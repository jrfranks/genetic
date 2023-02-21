# +
# Makefile for building genetic.o
#
# Copyright (c) 2022 Sveltesoft. All rights reserved.
#
# This software is the confidential and proprietary information of SvelteSoft.
# You shall not disclose such confidential information and shall use it only in
# accordance with the terms of the license agreement you entered into with
# SvelteSoft.
#-

CC = gcc
CFLAGS = -c -Wall

all: genetic.o

genetic.o: genetic.c genetic.h
	$(CC) $(CFLAGS) genetic.c

clean:
	rm -f *.o
