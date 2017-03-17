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
// $Id: mainwindow.cpp 1805 2016-04-05 16:32:34Z maronga $
//
// 1804 2016-04-05 16:30:18Z maronga 
// Removed parameter file check
// Update for use with qt5
//
// 1723 2015-11-16 15:25:51Z boeske
// Added checkbox for cyclic fill flag in mrungui
//
// 1611 2015-07-07 12:23:22Z maronga
// Added routine start_watchdog.
// 
// 1046 2012-11-09 14:38:45Z maronga
// code put under GPL (PALM 3.9)
//
// mainwindow.cpp 920 2012-06-05 09:56:53Z maronga
// Added -Z option (disable combine_plot_fields)
//
// 818 2012-02-08 16:11:23Z maronga
// New: flag parameter -z implemented (used for skipping parameter file check)
//
// 811 2012-01-31 09:45:54Z maronga
// Bugfix: waitForFinished time exceeded limits, adjusted to 1h waiting time
//
// 809 2012-01-30 13:32:58Z marong
// Bugfix: waiting for 20 minutes does not suffice for local runs
//
// 793 2011-12-12 15:15:24Z maronga
// Initial revision
//
// Description:
// ------------
// All subroutines for mrungui
//----------------------------------------------------------------------------//

#include <QApplication>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "ui_about.h"
#include "ui_help.h"
#include <QString>
#include <QFileDialog>
#include <QDateTime>
#include "stdlib.h"
#include <QProcess>
#include <QTextStream>


// Define username on host as global variable
QString username = getenv("USER");
// Define variable branch (location of mrungui, e.g.
// "/home/user/palm/current_version"
QString branch;


//****************************************************************************//
//  Initialization of the main window
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

//  Empty the mrun command line (hereafter mrunline)
    ui->commandline->setText("");

//  Get path of the program and set default file
    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);


    QFile file(branch+"/.mrun.gui.default");

//  Read default settings
    if ( file.exists() == true && file.size() > 10)
    {

       QString mrunline;
       if ( file.open(QIODevice::ReadOnly | QIODevice::Text ) )
       {

//        File opened successfully
          QTextStream in(&file);

          QString line = in.readLine();
          while (!line.isNull())
          {
             mrunline = line;
             line = in.readLine();
          }

          file.close();
       }

       mrunline = mrunline.right(mrunline.length() - 17);
       ui->commandline->setText(mrunline);

       setup_gui(mrunline);

    }

//  Load jobs into the recent jobs list
    recent_jobs(10);
}


//****************************************************************************//
//  Routine for listing jobs in the recent jobs list
int MainWindow::recent_jobs(int noj)
{
    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);

    QFile file(branch+"/.mrun.history");
    QString listitem, timestamp;

//  Read mrun history and select the last used jobs
    if ( file.exists() == true && file.size() > 10)
    {

//     Open history
       QStringList history, tmphistory;
       if ( file.open(QIODevice::ReadOnly | QIODevice::Text ) )
       {

//        file opened successfully
          QTextStream in(&file);

          QString line = in.readLine();
          while (!line.isNull())
          {
             history.append(line);
             line = in.readLine();
          }

          file.close();
       }

       int j = 0;

       ui->list_jobname->clear();

//     Read history entries and append to recent job list
       for (int i=history.count(); i>=1; i--)
       {
           timestamp = history[i-1].left(16);
           listitem = history[i-1].right(history[i-1].length() - 17);
           listitem = listitem.split("-d ", QString::SkipEmptyParts)[1];
           listitem = listitem.split(" -", QString::SkipEmptyParts)[0];
           listitem = listitem.replace(" ","");

           QList<QListWidgetItem *> matchitems = \
                   ui->list_jobname->findItems(listitem, Qt::MatchExactly);

           if ( matchitems.count() == 0 )
           {
              ui->list_jobname->addItem(listitem);
              tmphistory.append(listitem+" ("+timestamp+")");
              j++;
           }
           if ( j == noj )
           {
               break;
           }
       }

//     Send to list
       ui->list_jobname->clear();
       for (int i=tmphistory.count(); i>=1; i--)
       {
           ui->list_jobname->addItem(tmphistory[i-1]);
       }

    }
    return 0;
}

//****************************************************************************//
//  Exit program
MainWindow::~MainWindow()
{

    delete ui;
}


//****************************************************************************//
//  Start the mrun command via xterm
int MainWindow::startmrun()
{

    QString xtermoutput;
    QString mrunline_save;
    QString history_line;
    QString mrunline = ui->commandline->text();
    QString userline = ui->line_user->text();

//  Check for empty line
    mrunline_save = mrunline;
    if (userline != "")
    {
        mrunline = mrunline + " " + userline;
    }
    history_line = mrunline;

//  Disable the main window
    ui->group_job->setEnabled(false);
    ui->group_execution->setEnabled(false);
    ui->group_runcontrol->setEnabled(false);
    ui->group_advanced->setEnabled(false);
    ui->check_advanced->setEnabled(false);
    ui->check_verbose->setEnabled(false);
    ui->button_start->setEnabled(false);
    ui->button_start->setText("wait...");

    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);

    ui->commandline->setText("Executing MRUN in xterm...");

//  Wait until all commands have been executed (ugly)
    for(int i=0; i<20; i++)
       qApp->processEvents();

//  Start xterm as QProcess
    QProcess mrun;
    mrun.setProcessChannelMode(QProcess::MergedChannels);
    mrun.setWorkingDirectory(branch);

    mrunline = " -title \"Executing MRUN...\" -geometry \"100x55+970+0\" -e \""\
          +mrunline.replace("\"","\'")\
          +";echo -n '--> Press Enter to continue...';read yesno\"";

    mrun.start("xterm"+mrunline);

    if(!mrun.waitForStarted())
        return 0;

//    while(mrun.waitForReadyRead())
//        xtermoutput = xtermoutput+mrun.readAllStandardOutput();

//  Wait until mrun has finished or wait for 200 minutes
    mrun.waitForFinished(3600000);


//  Jobs has been submitted or aborted. Continuing...
//  Save the mrun command to history file
    QString filename = branch+"/.mrun.history";

    QDateTime time = QDateTime::currentDateTime();
    QString tmptime = time.toString("yyyy/MM/dd hh:mm");

    QFile file(filename);
    file.open(QIODevice::Append | QIODevice::Text);
    QTextStream out(&file);
    out << tmptime + " " + history_line + "\n";
    file.close();

//  Enable main window again
    ui->group_job->setEnabled(true);
    ui->group_execution->setEnabled(true);
    ui->group_runcontrol->setEnabled(true);
    if ( ui->check_advanced->isChecked() == true)
    {
       ui->group_advanced->setEnabled(true);
    }
    ui->check_advanced->setEnabled(true);
    ui->check_verbose->setEnabled(true);
    ui->button_start->setEnabled(true);
    ui->action_save->setEnabled(true);
    ui->button_start->setText("MRUN Start");
    ui->commandline->setText(mrunline_save);

//  Reload recent jobs
    recent_jobs(10);

    return 0;
}


//****************************************************************************//
//  Disable/Enable advanced settings dialog
int MainWindow::enable_advanced()
{
    bool status;

    status = ui->group_advanced->isEnabled();
    if (status == true)
    {
       ui->group_advanced->setEnabled(false);
    }
    else
    {
       ui->group_advanced->setEnabled(true);
    }

    return 0;
}


//****************************************************************************//
//  Disable/enable dialog for coupled runs
int MainWindow::enable_coupled()
{
    QString status;

    status = ui->drop_job->currentText();
    if (status == "Coupled restart")
    {
       ui->label_coupled1->setEnabled(true);
       ui->label_coupled2->setEnabled(true);
       ui->label_coupled3->setEnabled(true);
       ui->line_PE_atmos->setEnabled(true);
       ui->line_PE_ocean->setEnabled(true);

       QString pe_total = ui->line_pe->text();
       QString pe_atmos = ui->line_PE_atmos->text();
       QString pe_ocean = ui->line_PE_ocean->text();
       if (pe_total != "")
       {
           int PE_split = pe_total.toInt()/2;
           ui->line_PE_atmos->setText(QString::number(PE_split));
           ui->line_PE_ocean->setText(QString::number(PE_split));
       }
    }
    else
    {
        ui->label_coupled1->setEnabled(false);
        ui->label_coupled2->setEnabled(false);
        ui->label_coupled3->setEnabled(false);
        ui->line_PE_atmos->setEnabled(false);
        ui->line_PE_ocean->setEnabled(false);
    }

    return 0;
}


//****************************************************************************//
//  Choose job from dialog
int MainWindow::choosejob()
{
    QString filename;
    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);

    QString jobdir = branch+"/JOBS";


//  Pick a job from dialog
    filename = QFileDialog::getExistingDirectory(this, \
                            tr("Choose a job directory"), branch+"/JOBS", \
               QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks | \
               QFileDialog::DontUseNativeDialog | QFileDialog::ReadOnly | \
               QFileDialog::DontConfirmOverwrite);





//  If a file was selected, load it into mainwindow
    if ( filename != "")
    {
        filename = filename.right(filename.length() - jobdir.length() - 1);

        ui->line_jobname->setText(filename);
        ui->list_jobname->clearSelection();

//      Change mrunline accordingly
        change_commandline("d","");
        change_commandline("r","");
        return 0;
    }
    else
    {
        return 1;
    }
}


//****************************************************************************//
//  Choose job from the recent job list and load all settings
int MainWindow::choosejob_list()
{
    QString filename, timestamp, jobname;

//  Get selected item from list
    filename = ui->list_jobname->currentItem()->text();
    int itemint = ui->list_jobname->currentRow();

//  Reload list
    ui->setupUi(this);
    recent_jobs(10);

//  Set selected item to jobname
    ui->list_jobname->item(itemint)->setSelected(true);

    timestamp = filename.right(17).left(16);
    jobname = filename.left(filename.length() - 19);

    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);
    QFile file(branch+"/.mrun.history");
    QString listitem;

//  Load history
    if ( file.exists() == true && file.size() > 10)
    {

       QStringList history;
       if ( file.open(QIODevice::ReadOnly | QIODevice::Text ) )
       {

//        file opened successfully
          QTextStream in(&file);

          QString line = in.readLine();
          while (!line.isNull())
          {
             history.append(line);
             line = in.readLine();
          }

          file.close();
       }

       for (int i=history.count(); i>=1; i--)
       {
           listitem = history[i-1].right(history[i-1].length() - 17);
           listitem = listitem.split("-d ", QString::SkipEmptyParts)[1];
           listitem = listitem.split(" -", QString::SkipEmptyParts)[0];
           listitem = listitem.replace(" ","");


//         Select command line with correct timestamp and jobname
           if (history[i-1].left(16) == timestamp && listitem == jobname)
           {
              QString mrunline = history[i-1];
              mrunline = mrunline.right(mrunline.length() - 17);
              ui->commandline->setText(mrunline);

              setup_gui(mrunline);
           }
       }
    }

     return 0;
}


//****************************************************************************//
//  Change run identifer (initial, restart, coupled...)
int MainWindow::change_rc_list()
{


   QString drop_job = ui->drop_job->currentText();

   change_commandline("r","");

// Enable PE distribution for atmosphere/ocean
   if ( drop_job == "Coupled restart")
   {
       QString drop_atmos = ui->line_PE_atmos->text();
       QString drop_ocean = ui->line_PE_ocean->text();

       change_commandline("Y",drop_atmos+" "+drop_ocean);
   }

// Check of ocean runs
   else
   {
      delete_commandline("Y");
      if (drop_job == "Precursor run (Ocean)")
      {
          activate_flag("y");
      }
      else
      {
          deactivate_flag("y");
      }
   }


    return 0;
}


//****************************************************************************//
//  Routine for processing any changes in the mainwindow settings
int MainWindow::change_commandline(QString id,QString fwstring)
{

//  First get the mrunline
    QString newmrunline;
    bool    initialize=false;

    QString mrunline = ui->commandline->text();

    QStringList splitline = mrunline.split(" -"+id+"", QString::SkipEmptyParts);
    if ( splitline.count() == 1)
    {
        splitline.append(" ");
    }
    else if ( splitline.count() == 0 )
    {
       splitline.append("mrun");
       splitline.append(" ");
       initialize = true;
    }

    QStringList param = splitline[1].split("-");

//  Change in parameter "d" (jobname)
    if (id == "d")
    {
       QString filename = ui->line_jobname->text();

       param[0] = filename.replace(" ","");

       if ( initialize == true && ui->group_runcontrol->isEnabled() == true )
       {
           ui->group_runcontrol->setEnabled(true);
           ui->group_execution->setEnabled(true);
           ui->drop_job->setEnabled(true);
           ui->check_advanced->setEnabled(true);
           ui->button_start->setEnabled(true);
           ui->action_save->setEnabled(true);
       }
       else if ( param[0] == "")
       {
           ui->group_runcontrol->setEnabled(false);
           ui->group_execution->setEnabled(false);
           ui->drop_job->setEnabled(false);
           ui->button_start->setEnabled(false);
           ui->action_save->setEnabled(false);
           ui->group_advanced->setEnabled(false);
           ui->check_advanced->setEnabled(false);
           delete_commandline("d");
           change_commandline("r","remove");
           ui->label_usercode->setText("");
           return 1;
       }
       else
       {

//         Check if user code is available
           branch = QCoreApplication::applicationDirPath();
           branch = branch.left(branch.length() - 14);
           QDir usercode = branch+"/USER_CODE/"+param[0];
           if (usercode.exists() == true)
           {
               ui->label_usercode->setText("User code found.");
           }
           else
           {
               ui->label_usercode->setText("Warning: no user code found!");
           }

//         Check if _pdf file is available, otherwise notice user
           if (ui->check_restarts->checkState() == 2)
           {
              QString jobname = ui->line_jobname->text();
              QFile restartfile(branch+"/JOBS/"+jobname+"/INPUT/"+jobname+"_p3df");
              if (restartfile.exists() == true)
              {
                 ui->label_restart->setText("");
              }
              else
              {
                 ui->label_restart->setText("Warning: No p3df file found!");
              }
           }

           ui->group_runcontrol->setEnabled(true);
           ui->group_execution->setEnabled(true);
           ui->drop_job->setEnabled(true);
           ui->check_advanced->setEnabled(true);

           if ( ui->check_advanced->isChecked() == true )
           {
              ui->group_advanced->setEnabled(true);
              if ( ui->line_i->text() != "" || ui->line_o->text() != "")
              {
                 change_commandline("r","remove");
                 ui->group_runcontrol->setEnabled(false);
              }
           }
       }

    }

//  Change in parameter "r" (run control list)
    else if (id == "r")
    {
        int status_ts = ui->check_ts->checkState();
        int status_pr = ui->check_pr->checkState();
        int status_xy = ui->check_xy->checkState();
        int status_xz = ui->check_xz->checkState();
        int status_yz = ui->check_yz->checkState();
        int status_3d = ui->check_3d->checkState();
        int status_ma = ui->check_ma->checkState();
        int status_sp = ui->check_sp->checkState();
        int status_pts = ui->check_pts->checkState();
        int status_prt = ui->check_prt->checkState();

        QString drop_job = ui->drop_job->currentText();
        QString rc_flag = "#";

        if (drop_job == "Initial run")
        {
            rc_flag = "#";
        }
        else if (drop_job == "Restart run")
        {
            rc_flag = "f";
        }
        else if (drop_job == "Precursor run (Atmosphere)")
        {
            rc_flag = "#";
        }
        else if (drop_job == "Precursor run (Ocean)")
        {
           rc_flag = "o#";
        }
        else if (drop_job == "Coupled restart")
        {
           rc_flag = "f";
        }

        param[0] = "\"d3"+rc_flag;

        if (status_ts == 2)
        {
            param[0] = param[0]+" ts"+rc_flag;
        }
        if (status_pr == 2)
        {
           param[0] = param[0]+" pr"+rc_flag;
        }
        if (status_xy == 2)
        {
           param[0] = param[0]+" xy"+rc_flag;
        }
        if (status_xz == 2)
        {
           param[0] = param[0]+" xz"+rc_flag;
        }
        if (status_yz == 2)
        {
           param[0] = param[0]+" yz"+rc_flag;
        }
        if (status_3d == 2)
        {
           param[0] = param[0]+" 3d"+rc_flag;
        }
        if (status_ma == 2)
        {
           param[0] = param[0]+" ma"+rc_flag;
        }
        if (status_sp == 2)
        {
           param[0] = param[0]+" sp"+rc_flag;
        }
        if (status_prt == 2)
        {
           param[0] = param[0]+" prt"+rc_flag;
        }
        if (status_pts == 2)
        {
           param[0] = param[0]+" pts"+rc_flag;
        }

        if (drop_job == "Coupled restart")
        {
            rc_flag = "of";
            param[0] = param[0]+" d3"+rc_flag;

            if (status_ts == 2)
            {
                param[0] = param[0]+" ts"+rc_flag;
            }
            if (status_pr == 2)
            {
               param[0] = param[0]+" pr"+rc_flag;
            }
            if (status_xy == 2)
            {
               param[0] = param[0]+" xy"+rc_flag;
            }
            if (status_xz == 2)
            {
               param[0] = param[0]+" xz"+rc_flag;
            }
            if (status_yz == 2)
            {
               param[0] = param[0]+" yz"+rc_flag;
            }
            if (status_3d == 2)
            {
               param[0] = param[0]+" 3d"+rc_flag;
            }
            if (status_ma == 2)
            {
               param[0] = param[0]+" ma"+rc_flag;
            }
            if (status_sp == 2)
            {
               param[0] = param[0]+" sp"+rc_flag;
            }
            if (status_prt == 2)
            {
               param[0] = param[0]+" prt"+rc_flag;
            }
            if (status_pts == 2)
            {
               param[0] = param[0]+" pts"+rc_flag;
            }
        }

        int status_restarts = ui->check_restarts->checkState();



        if (status_restarts == 2)
        {
            param[0]=param[0]+" restart";

//          Check if _pdf file is available, otherwise notice user
            if (ui->check_restarts->checkState() == 2)
            {
               QString jobname = ui->line_jobname->text();
               branch = QCoreApplication::applicationDirPath();
               branch = branch.left(branch.length() - 14);
               QFile restartfile(branch+"/JOBS/"+jobname+"/INPUT/"+jobname+ \
                                 "_p3df");
               if (restartfile.exists() == true)
               {
                  ui->label_restart->setText("");
               }
               else
               {
                  ui->label_restart->setText("Warning: No p3df file found!");
               }
            }

        }
        else  {
          ui->label_restart->setText("");
        }

        
        
        int status_cycfill = ui->check_cycfill->checkState();

        if (status_cycfill == 2)
        {
            param[0]=param[0]+" cycfill";
        }


        
        param[0]=param[0]+"\" ";

        if ( fwstring == "remove")
        {
            delete_commandline(id);
            return 1;
        }
        else
        {
           ui->button_start->setEnabled(true);
           ui->action_save->setEnabled(true);
        }
    }
//  Change in any other parameter
    else
    {
        if ( fwstring != "")
        {
           param[0] = "\""+fwstring+"\"";
        }
        else
        {
            delete_commandline(id);
            return 1;
        }

    }
    param[0] = param[0] + " ";

//  Join the new mrunline
    splitline[1]= param.join("-");
    newmrunline = splitline.join(" -"+id+" ");

//  Print the new mrunline to mainwindow
    newmrunline.replace("  "," ");
    ui->commandline->setText(newmrunline);

    return 0;
}


//****************************************************************************//
//  Get all signals from mainwindow and perform change via change_commandline()
int MainWindow::change_lineinput()
{
    QString tmptext;

    if ( sender() == ui->line_host )
    {
        tmptext = ui->line_host->text();
        change_commandline("h",tmptext);
    }
    else if ( sender() == ui->line_jobname)
    {
        tmptext = ui->line_jobname->text();
        change_commandline("d",tmptext);
    }
    else if ( sender() == ui->line_q)
    {
        tmptext = ui->line_q->text();
        change_commandline("q",tmptext);
    }
    else if ( sender() == ui->line_account)
    {
        tmptext = ui->line_account->text();
        change_commandline("u",tmptext);
    }
    else if ( sender() ==  ui->line_pe)
    {
        tmptext = ui->line_pe->text();
        change_commandline("X",tmptext);
    }
    else if ( sender() == ui->line_tpn)
    {
        tmptext = ui->line_tpn->text();
        change_commandline("T",tmptext);
    }
    else if ( sender() == ui->line_branch)
    {
        tmptext = ui->line_branch->text();
        change_commandline("K",tmptext);
    }
    else if ( sender() == ui->line_time)
    {
        tmptext = ui->line_time->text();
        change_commandline("t",tmptext);
    }
    else if ( sender() == ui->line_M)
    {
        tmptext = ui->line_M->text();
        change_commandline("M",tmptext);
    }
    else if ( sender() == ui->line_m)
    {
        tmptext = ui->line_m->text();
        change_commandline("m",tmptext);
    }
    else if ( sender() == ui->line_a)
    {
        tmptext = ui->line_a->text();
        change_commandline("a",tmptext);
    }
    else if ( sender() == ui->line_D)
    {
        tmptext = ui->line_D->text();
        change_commandline("D",tmptext);
    }
    else if ( sender() == ui->line_c)
    {
        tmptext = ui->line_c->text();
        if ( tmptext == ".mrun.config")
        {
            tmptext = "";
        }
        change_commandline("c",tmptext);
    }
    else if ( sender() == ui->line_p)
    {
        tmptext = ui->line_p->text();
        change_commandline("p",tmptext);
    }
    else if ( sender() == ui->line_s)
    {
        tmptext = ui->line_s->text();
        change_commandline("s",tmptext);
    }
    else if ( sender() == ui->line_i)
    {
        tmptext = ui->line_i->text();
        if ( tmptext != "")
        {
            change_commandline("r","remove");
            ui->group_runcontrol->setEnabled(false);

            ui->button_start->setEnabled(true);
            ui->action_save->setEnabled(true);
        }
        else if (ui->line_o->text() == "" )
        {
           ui->group_runcontrol->setEnabled(true);
           change_commandline("r","");
        }

        change_commandline("i",tmptext);
    }
    else if ( sender() == ui->line_o)
    {
        tmptext = ui->line_o->text();
        if ( tmptext != "")
        {
            change_commandline("r","remove");
            ui->button_start->setEnabled(true);
            ui->action_save->setEnabled(true);
            ui->group_runcontrol->setEnabled(false);
        }
        else if (ui->line_i->text() == "" )
        {
           ui->group_runcontrol->setEnabled(true);
           change_commandline("r","");
        }

        change_commandline("o",tmptext);

    }
    else if ( sender() == ui->line_w)
    {
        tmptext = ui->line_w->text();
        change_commandline("w",tmptext);
    }
    else if ( sender() == ui->line_PE_atmos || sender() == ui->line_PE_ocean)
    {
        tmptext = ui->line_PE_atmos->text() + " "  + ui->line_PE_ocean->text();
        change_commandline("Y",tmptext);
        int PE_total = ui->line_PE_atmos->text().toInt() + \
                        ui->line_PE_ocean->text().toInt();
        ui->line_pe->setText(QString::number(PE_total));
        change_commandline("X",QString::number(PE_total));
    }

    else if ( sender() == ui->combo_n)
    {
        tmptext = ui->combo_n->currentText();
        if ( tmptext == "default")
        {
            tmptext = "";
        }
        change_commandline("n",tmptext);
    }
    else if ( sender() == ui->line_user)
    {
       return 0;
    }

    return 0;
}


//****************************************************************************//
//  Delete a parameter from mrunline
int MainWindow::delete_commandline(QString id)
{

//  Read mrunline
    QString newmrunline;
    QString mrunline = ui->commandline->text();
    QStringList splitline = mrunline.split(" -"+id, QString::SkipEmptyParts);
    if ( splitline.count() == 1)
    {
       return 0;
    }
    else
    {
       QStringList param = splitline[1].split("-");
       param[0] = "";
       splitline[1]= param.join(" -");
       newmrunline = splitline.join("");

//     Print new mrunline to screen
       ui->commandline->setText(newmrunline);

       return 0;
    }
}


//****************************************************************************//
//  Check routine for all flag parameters
int MainWindow::check_flags()
{

    int status = ui->check_delete_tmp_files->checkState();

    if (status == 2)
    {
        activate_flag("B");
    }
    else
    {
        deactivate_flag("B");
    }

    status = ui->check_verbose->checkState();

    if (status == 2)
    {
        activate_flag("v");
    }
    else
    {
        deactivate_flag("v");
    }

    status = ui->check_A->checkState();

    if (status == 2)
    {
        activate_flag("A");
    }
    else
    {
        deactivate_flag("A");
    }

    status = ui->check_b->checkState();

    if (status == 2)
    {
        activate_flag("b");
    }
    else
    {
        deactivate_flag("b");
    }

    status = ui->check_F->checkState();

    if (status == 2)
    {
        activate_flag("F");
    }
    else
    {
        deactivate_flag("F");
    }
    status = ui->check_I->checkState();

    if (status == 2)
    {
        activate_flag("I");
    }
    else
    {
        deactivate_flag("I");
    }
    status = ui->check_k->checkState();

    if (status == 2)
    {
        activate_flag("k");
    }
    else
    {
        deactivate_flag("k");
    }
    status = ui->check_O->checkState();

    if (status == 2)
    {
        activate_flag("O");
    }
    else
    {
        deactivate_flag("O");
    }
    status = ui->check_S->checkState();

    if (status == 2)
    {
        activate_flag("S");
    }
    else
    {
        deactivate_flag("S");
    }
    status = ui->check_x->checkState();

    if (status == 2)
    {
        activate_flag("x");
    }
    else
    {
        deactivate_flag("x");
    }
    
    status = ui->check_Z->checkState();

    if (status == 2)
    {
        activate_flag("Z");
    }
    else
    {
        deactivate_flag("Z");
    }
    return 0;

}


//****************************************************************************//
//  Activate flag
int MainWindow::activate_flag(QString id)
{
    QString newmrunline;

    QString mrunline = ui->commandline->text();

    QStringList splitline = mrunline.split(" -"+id, QString::SkipEmptyParts);
    if ( splitline.count() == 1)
    {
        splitline.append("");
        newmrunline = splitline.join(" -"+id);

 //     print new commandline to screen
        newmrunline.replace("  "," ");
        ui->commandline->setText(newmrunline);

       return 0;
    }
    else
    {
       return 0;
    }
}


//****************************************************************************//
//  Deactivate flag
int MainWindow::deactivate_flag(QString id)
{
    QString newmrunline;

    QString mrunline = ui->commandline->text();

    QStringList splitline = mrunline.split(" -"+id, QString::SkipEmptyParts);

    newmrunline = splitline.join("");
//  print new commandline to screen
    newmrunline.replace("  "," ");
    ui->commandline->setText(newmrunline);

    return 0;
}


//****************************************************************************//
//  Mainwindow reset
int MainWindow::reset_window()
{
    ui->setupUi(this);
    recent_jobs(10);
    return 0;
}


//****************************************************************************//
//  Save current setting as default
int MainWindow::save_default()
{

//  Get mrunline
    QString mrunline_save;
    QString mrunline = ui->commandline->text();
    QString userline = ui->line_user->text();

    mrunline_save = mrunline;
    if (userline != "")
    {
        mrunline = mrunline + " " + userline;
    }


//  Define filename
    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);

    QString filename = branch+"/.mrun.gui.default";

//  Prepare data output
    QDateTime time = QDateTime::currentDateTime();
    QString tmptime = time.toString("yyyy/MM/dd hh:mm");

//  Write to file
    QFile file(filename);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    if ( mrunline == "" || mrunline == " ")
    {
        out << "";
    }
    else
    {
       out << tmptime + " " + mrunline;
    }

    file.close();

    return 0;
}


//****************************************************************************//
//  Save current settings to a file of choice
int MainWindow::save_to_file()
{

//  Get mrunline
    QString mrunline_save;
    QString mrunline = ui->commandline->text();
    QString userline = ui->line_user->text();

    mrunline_save = mrunline;
    if (userline != "")
    {
        mrunline = mrunline + " " + userline;
    }

//  Define a filename
    QString filename = QFileDialog::getSaveFileName(this, tr("Save to file"), \
                                    "", tr("Save files (*.sav)"));
    QString extension = filename.right(4);

    if ( extension != ".sav" )
    {
        filename = filename + ".sav";
    }

//  Prepare data output
    QDateTime time = QDateTime::currentDateTime();
    QString tmptime = time.toString("yyyy/MM/dd hh:mm");

//  Write to file
    QFile file(filename);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out << tmptime + " " + mrunline;
    file.close();

    return 0;
}

//****************************************************************************//
//  Start watchdog (palm_wd)
int MainWindow::start_watchdog()
{
    return system("nohup palm_wd >> /dev/null 2>&1 &");

}


//****************************************************************************//
//  Open job from file (previously saved by the user)
int MainWindow::open_from_file()
{

//  Select filename and open it
    QString filename = QFileDialog::getOpenFileName(this, tr("Open File"), \
                                    "", tr("Save files (*.sav)"));

    QFile file(filename);
    QString mrunline;
    if ( filename != "")
    {
       if ( file.open(QIODevice::ReadOnly | QIODevice::Text ) )
       {

//        File opened successfully
          QTextStream in(&file);

          QString line = in.readLine();
          while (!line.isNull())
          {
             mrunline = line;
             line = in.readLine();
          }
          file.close();
       }

//     In case a mrunline was found, load it to mainwindow
       if ( mrunline != "")
       {

          mrunline = mrunline.right(mrunline.length() - 17);
          ui->commandline->setText(mrunline);

          setup_gui(mrunline);
       }
    }
    return 0;
}


//****************************************************************************//
//  Open the last submitted job
int MainWindow::open_last()
{
    branch = QCoreApplication::applicationDirPath();
    branch = branch.left(branch.length() - 14);

    QFile file(branch+"/.mrun.history");

//  Load mrun history
    QString mrunline;
    if ( file.open(QIODevice::ReadOnly | QIODevice::Text ) )
    {

//     file opened successfully
       QTextStream in(&file);

       QString line = in.readLine();
       while (!line.isNull())
       {
          mrunline = line;
          line = in.readLine();
       }

       file.close();
    }

//  Select the last submitted job and print it to mainwindow
    if ( mrunline != "" )
    {
       mrunline = mrunline.right(mrunline.length() - 17);
       ui->commandline->setText(mrunline);

       setup_gui(mrunline);
    }
    return 0;
}

//****************************************************************************//
//  Display help window
int MainWindow::help()
{

    QDialog *child_help = new QDialog;
    Ui::child_help ui;
    ui.setupUi(child_help);
    child_help->show();


    return 0;
}

//****************************************************************************//
//  Display about window
int MainWindow::about_gui()
{

    QDialog *about = new QDialog;
    Ui::about ui;
    ui.setupUi(about);
    about->show();

    return 0;
}

//****************************************************************************//
//  Setup mainwindow according to a given mrunline
int MainWindow::setup_gui(QString mrunline)
{

//  Some initial settings
    QString user = "";

    bool coupled_run = false;
    bool ocean_run   = false;
    bool advanced    = false;
    bool nojob       = false;
    bool rc_manual   = false;

//  Split parameters in mrunline
    QStringList splitline = mrunline.split(" -", QString::SkipEmptyParts);

    if ( splitline[0] != "mrun")
    {
        return 1;

    }
    else
    {
        ui->group_job->setEnabled(true);

//      Loop for all parameters in mrunline
        for(int i=splitline.count()-1; i>=1; i--)
        {

//          Determine parameter
            QStringList splitparameter = splitline[i].split(" ", \
                                                      QString::SkipEmptyParts);

            QString parameter = splitparameter[0];
            splitparameter.removeFirst();
            QString options = splitparameter.join(" ");
            options = options.replace("\"","");

//          Check for suitable switch
            if ( parameter == "d")
            {
               if ( options != "")
               {
                  ui->line_jobname->setText(options);
                  nojob = false;
               }
               else
               {
                  nojob = true;
               }
            }
            else if ( parameter == "h")
            {
               ui->line_host->setText(options);
            }
            else if ( parameter == "q")
            {
               ui->line_q->setText(options);
            }
            else if ( parameter == "K")
            {
               ui->line_branch->setText(options);
            }
            else if ( parameter == "u")
            {
               ui->line_account->setText(options);
            }
            else if ( parameter == "X")
            {
               ui->line_pe->setText(options);
            }
            else if ( parameter == "T")
            {
               ui->line_tpn->setText(options);
            }
            else if ( parameter == "t")
            {
               ui->line_time->setText(options);
            }
            else if ( parameter == "B")
            {
               ui->check_delete_tmp_files->setChecked(true);
            }
            else if ( parameter == "v")
            {
               ui->check_verbose->setChecked(true);
            }
            else if ( parameter == "A")
            {
               ui->check_A->setChecked(true);
            }
            else if ( parameter == "b")
            {
               ui->check_b->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "F")
            {
               ui->check_F->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "I")
            {
               ui->check_I->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "k")
            {
               ui->check_k->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "O")
            {
               ui->check_O->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "S")
            {
               ui->check_S->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "x")
            {
               ui->check_x->setChecked(true);
               advanced = true;
            }
            else if ( parameter == "Z")
            {
               ui->check_Z->setChecked(true);
               advanced = true;
            }  
            else if ( parameter == "m")
            {
               ui->line_m->setText(options);
               advanced = true;
            }
            else if ( parameter == "M")
            {
               ui->line_M->setText(options);
               advanced = true;
            }
            else if ( parameter == "a")
            {
               ui->line_a->setText(options);
               advanced = true;
            }
            else if ( parameter == "D")
            {
               ui->line_D->setText(options);
               advanced = true;
            }
            else if ( parameter == "c")
            {
               ui->line_c->setText(options);
               advanced = true;
            }
            else if ( parameter == "p")
            {
               ui->line_p->setText(options);
               advanced = true;
            }
            else if ( parameter == "s")
            {
               ui->line_s->setText(options);
               advanced = true;
            }
            else if ( parameter == "i")
            {
               ui->line_i->setText(options);
               advanced = true;
               rc_manual = true;
            }
            else if ( parameter == "o")
            {
               ui->line_o->setText(options);
               advanced = true;
               rc_manual = true;
            }
            else if ( parameter == "w")
            {
               ui->line_w->setText(options);
               advanced = true;
            }

//          Determine settings for coupled restart runs
            else if ( parameter == "Y")
            {
                QStringList optionssplit = options.split(" ", \
                                                   QString::SkipEmptyParts);

                ui->line_PE_atmos->setEnabled(true);
                ui->line_PE_ocean->setEnabled(true);
                ui->label_coupled1->setEnabled(true);
                ui->label_coupled2->setEnabled(true);
                ui->label_coupled3->setEnabled(true);
                ui->label_coupling->setEnabled(true);

                if (optionssplit.count() == 2)
                {
                   ui->line_PE_atmos->setText(optionssplit[0]);
                   ui->line_PE_ocean->setText(optionssplit[1]);
                }
                else
                {
                   ui->line_PE_atmos->setText("");
                   ui->line_PE_ocean->setText("");
                }
                coupled_run = true;
            }
            else if ( parameter == "n")
            {
                if ( options == "shared")
                {
                    ui->combo_n->setCurrentIndex(1);
                }
                else if ( options == "non_shared")
                {
                    ui->combo_n->setCurrentIndex(2);
                }
                else
                {
                    ui->combo_n->setCurrentIndex(0);
                }
            }
            else if ( parameter == "y")
            {
               ui->drop_job->setCurrentIndex(3);
            }

//          Determine settings for the run control list
            else if ( parameter == "r")
            {
               QStringList optionssplit = options.split(" ", \
                                          QString::SkipEmptyParts);

               QString options_2;
               QString options_all;
               for (int j=0; j<optionssplit.count(); j++)
               {
                   options_all = optionssplit[j];
                   options_2 = optionssplit[j].left(2);
                   if (options_2 == "ts")
                   {
                       ui->check_ts->setChecked(true);
                   }
                   if (options_2 == "pr" && options_all.left(3) != "prt")
                   {
                       ui->check_pr->setChecked(true);
                   }
                   if (options_2 == "xy")
                   {
                       ui->check_xy->setChecked(true);
                   }
                   if (options_2 == "xz")
                   {
                       ui->check_xz->setChecked(true);
                   }
                   if (options_2 == "yz")
                   {
                       ui->check_yz->setChecked(true);
                   }
                   if (options_2 == "3d")
                   {
                       ui->check_3d->setChecked(true);
                   }
                   if (options_2 == "ma")
                   {
                       ui->check_ma->setChecked(true);
                   }
                   if (options_2 == "sp")
                   {
                       ui->check_sp->setChecked(true);
                   }
                   if (options_all.left(3) == "prt")
                   {
                       ui->check_prt->setChecked(true);
                   }
                   if (options_all.left(3) == "pts")
                   {
                       ui->check_pts->setChecked(true);
                   }
                   if (options_2 == "d3")
                   {
                       if (options_all.left(3).right(1) == "#")
                       {
                          ui->drop_job->setCurrentIndex(0);
                       }
                       else if (options_all.left(3).right(1) == "f")
                       {
                          ui->drop_job->setCurrentIndex(1);
                       }
                       else if (options_all.left(3).right(1) == "o")
                       {
                          ocean_run = true;
                       }
                   }
                   if (options_all == "restart")
                   {
                      ui->check_restarts->setChecked(true);
//                    Check if _pdf file is available, otherwise notice user
                      QString jobname = ui->line_jobname->text();
                      branch = QCoreApplication::applicationDirPath();
                      branch = branch.left(branch.length() - 14);

                      QFile restartfile(branch+"/JOBS/"+jobname+"/INPUT/"+ \
                                        jobname+"_p3df");
                      if (restartfile.exists() == true)
                      {
                         ui->label_restart->setText("");
                      }
                      else
                      {
                         ui->label_restart->setText("Warning: No p3df file \
                                                    found!");
                      }
                   }

               }

            }

//          All unknown parameters are set as extra user parameters
            else
            {
                user = user + "-" + parameter + " \"" + options + "\" ";
                splitline.removeAt(i);
            }
        }

//      Change drop box state in case of ocean precursor or coupled restart runs
        if ( ocean_run == true )
        {
           if ( coupled_run == true )
           {
             ui->drop_job->setCurrentIndex(4);
           }
           else
           {
              ui->drop_job->setCurrentIndex(3);
           }
        }

        if ( user != "")
        {
           ui->line_user->setText(user);
        }

//      Join mrunline and post it to mainwindow
        mrunline = splitline.join(" -");
        ui->commandline->setText(mrunline);

//      Check if advanced settings are used and enable advanced dialog
        if ( advanced == true )
        {
           ui->check_advanced->setChecked(true);
        }

//      Disable mainwindow if no job was found, otherwise enable
        if ( nojob == true )
        {
           ui->group_execution->setEnabled(false);
           ui->group_runcontrol->setEnabled(false);
           ui->button_start->setEnabled(false);
           ui->action_save->setEnabled(false);
           ui->drop_job->setEnabled(false);
           ui->group_advanced->setEnabled(false);
           ui->check_advanced->setEnabled(false);
        }
        else
        {
           ui->group_execution->setEnabled(true);
           ui->group_runcontrol->setEnabled(true);
           ui->button_start->setEnabled(true);
           ui->action_save->setEnabled(true);
           ui->check_advanced->setEnabled(true);
           ui->drop_job->setEnabled(true);
           if ( advanced == true )
           {
              ui->group_advanced->setEnabled(true);
           }
           else
           {
              ui->group_advanced->setEnabled(false);
           }
        }

//      Disable run control box, if manual settings of -i and -o are used
        if ( rc_manual == true )
        {
           ui->group_runcontrol->setEnabled(false);
           change_commandline("r","remove");
           ui->button_start->setEnabled(true);
           ui->action_save->setEnabled(true);
        }

    }

    return 0;
}
