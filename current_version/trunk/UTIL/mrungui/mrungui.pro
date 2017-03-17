#!/bin/ksh
#------------------------------------------------------------------------------#
# Current revisions:
# -----------------
# 
# 
# Former revisions:
# -----------------
# $Id: mrungui.pro 1805 2016-04-05 16:32:34Z maronga $
#
# 1804 2016-04-05 16:30:18Z maronga
# Update for use with qt5
#
# 793 2011-12-12 15:15:24Z maronga
# Initial revision
#
# Description:
# ------------
# Project file for mrungui (qmake instructions)
#------------------------------------------------------------------------------#

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = mrungui.x
DESTDIR = ../../SCRIPTS/
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui \
    help.ui \
    about.ui
