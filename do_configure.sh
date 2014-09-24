#!/bin/bash
./configure --enable-debug --enable-test --enable-numi --enable-t2k --enable-atmo --enable-rwgt --enable-vle-extension --enable-validation-tools --with-optimiz-level=O0 --with-log4cpp-inc=$LOG4CPP_INC --with-log4cpp-lib=$LOG4CPP_LIB >& log.config
