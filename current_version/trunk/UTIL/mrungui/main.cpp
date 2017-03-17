//------------------------------------------------------------------------------//
// This file is part of PALM.
//
// PALM is free software: you can redistribute it and/or modify it under the terms
// of the GNU General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// PALM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// PALM. If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 1997-2012  Leibniz University Hannover
//--------------------------------------------------------------------------------//
//
// Current revisions:
// -----------------
//
//
// Former revisions:
// -----------------
// $Id: main.cpp 1805 2016-04-05 16:32:34Z maronga $
//
// 1804 2016-04-05 16:30:18Z maronga 
// Update for use with qt5
//
// 1046 2012-11-09 14:38:45Z maronga
// code put under GPL (PALM 3.9)
//
// 793 2011-12-12 15:15:24Z maronga
// Initial revision
//
// Description:
// ------------
// Graphical user-interface (GUI) for the mrun script. All routines are placed
// in mainwindow.cpp. There are three windows: mainwindow; and two childs:
// about window and a help window (files ending with .ui). The initial code was
// written using Qt Creator 2.3.0.
//----------------------------------------------------------------------------//

#include <QApplication>
#include "mainwindow.h"
#include <QProcess>

int main(int argc, char *argv[])
{

    QApplication app(argc, argv);
//  app.setWindowIcon(QIcon("./mrungui.png"));

    MainWindow w;
    w.show();

    return app.exec();
}


