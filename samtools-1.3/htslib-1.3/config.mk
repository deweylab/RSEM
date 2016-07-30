#  Optional configure Makefile overrides for htslib.
#
#    Copyright (C) 2015 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# This is config.mk.  Generated from config.mk.in by configure.
#
# If you use configure, this file overrides variables and augments rules
# in the Makefile to reflect your configuration choices.  If you don't run
# configure, the main Makefile contains suitable conservative defaults.

prefix       = /usr/local
exec_prefix  = ${prefix}
bindir       = ${exec_prefix}/bin
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib
libexecdir   = ${exec_prefix}/libexec
datarootdir  = ${prefix}/share
mandir       = ${datarootdir}/man

CC     = gcc
RANLIB = ranlib

CPPFLAGS = 
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = 

# Lowercase here indicates these are "local" to config.mk
plugin_OBJS =
noplugin_LDFLAGS =
noplugin_LIBS =

# ifeq/.../endif, +=, and target-specific variables are GNU Make-specific.
# If you don't have GNU Make, comment out this conditional and note that
# to enable iRODS you will need to implement the following elsewhere.
ifeq "iRODS-disabled" "iRODS-enabled"

IRODS_HOME ?= /disabled

IRODS_CPPFLAGS = \
	-I$(IRODS_HOME)/lib/api/include \
	-I$(IRODS_HOME)/lib/core/include \
	-I$(IRODS_HOME)/lib/md5/include \
	-I$(IRODS_HOME)/lib/sha1/include \
	-I$(IRODS_HOME)/server/core/include \
	-I$(IRODS_HOME)/server/drivers/include \
	-I$(IRODS_HOME)/server/icat/include

IRODS_LDFLAGS = -L$(IRODS_HOME)/lib/core/obj
IRODS_LIBS = -lRodsAPIs -lgssapi_krb5

plugin_OBJS += hfile_irods.o

hfile_irods.o hfile_irods.pico: CPPFLAGS += $(IRODS_CPPFLAGS)
hfile_irods.so hfile_irods.bundle: LDFLAGS += $(IRODS_LDFLAGS)
hfile_irods.so hfile_irods.bundle: LIBS += $(IRODS_LIBS)

noplugin_LDFLAGS += $(IRODS_LDFLAGS)
noplugin_LIBS += $(IRODS_LIBS)

endif

ifeq "libcurl-disabled" "libcurl-enabled"

LIBCURL_LIBS = -lcurl 

plugin_OBJS += hfile_libcurl.o

hfile_libcurl.so hfile_libcurl.bundle: LIBS += $(LIBCURL_LIBS)

noplugin_LIBS += $(LIBCURL_LIBS)

endif

ifeq "plugins-no" "plugins-yes"

plugindir  = $(libexecdir)/htslib
pluginpath = $(libexecdir)/htslib

LIBHTS_OBJS += plugin.o
PLUGIN_OBJS += $(plugin_OBJS)

plugin.o plugin.pico: CPPFLAGS += -DPLUGINPATH=\"$(pluginpath)\"

else

LIBHTS_OBJS += $(plugin_OBJS)
LDFLAGS += $(noplugin_LDFLAGS)
LIBS += $(noplugin_LIBS)

endif
