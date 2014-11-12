#!/bin/bash
./configure \
	--enable-debug \
	--enable-atmo \
	--enable-rwgt \
	--enable-vle-extension \
	--with-optimiz-level=O0 \
	--with-log4cpp-inc=$LOG4CPP_INC \
	--with-log4cpp-lib=$LOG4CPP_LIB >& log.config
