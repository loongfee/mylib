#pragma once
#ifndef FILE_UTIL_HPP
#define FILE_UTIL_HPP

#ifdef USE_QT
#include <QFileInfo>
#include <QDir>
#endif

#include <direct.h>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <wtypes.h>
#include <atltypes.h>
#include <list>
using namespace std;

namespace mylib{
#ifdef USE_QT
/************************************************************************/
/* Qt                                                                     */
/************************************************************************/
bool QFindFile(const QString& path, const QStringList& filter, QStringList &vFiles);

bool QFindFile(const QString& path, const QStringList& filter, list<QString> &vFiles);

bool QCopyFolder(QString pathFrom, QString pathTo, QStringList filters = QStringList("*.*"), bool bOverWrite = true);
#endif

char* getExePath();
}; // end of namespace mylib
#endif