#pragma once
#ifndef FUNC_H
#define FUNC_H

//////// QT
#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QStandardItemModel>
#include <QGridLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QCheckBox>
#include <QToolButton>
#include <QStringList>
#include <QString>
#include <QMouseEvent>
#include <QPoint>
#include <QSettings>
#include <QMessageBox>
#include <QProgressBar>
#include <vector>
#include <QProgressDialog>
#include <QFileDialog>
#include <QThread> 
#include <QPoint> 
#include <QSettings>
#include <QDateEdit>
#include <QThread> 
#include <QFileInfo>
#include <QDir>
#include <QPointer>
#include "qfileinfo.h"
#include "qstring.h"
#include "qstringlist.h"
#include "qiodevice.h"
#include "qfile.h"
#include "qdir.h"
#include "qtextstream.h"
#include "qbytearray.h"

/////// ossim
#include <rspf/base/rspfDirectory.h>
#include <rspf/base/rspfDirectoryTree.h>
#include <rspf/base/rspfString.h>
#include "rspf/base/rspfFilename.h"
#include "rspf/base/rspfString.h"
#include "rspf/base/rspferrorcodes.h"
#include "rspf/imaging/rspfImageHandlerRegistry.h"
#include "rspf/imaging/rspfImageHandler.h"
#include "rspf/imaging/rspfImageFileWriter.h"
#include "rspf/imaging/rspfImageWriterFactoryRegistry.h"
#include <rspf/base/rspfDirectory.h>
#include <rspf/base/rspfThreeParamDatum.h>
#include <rspf/imaging/rspfFilterResampler.h>
#include "rspf/projection/rspfProjection.h"
#include <rspf/projection/rspfRpcProjection.h> 
#include "rspf/projection/rspfUtmProjection.h"
#include "rspf/projection/rspfTransMercatorProjection.h"
#include "rspf/projection/rspfProjectionFactoryRegistry.h"
#include "rspf/imaging/rspfImageRenderer.h"
#include "rspf/init/rspfInit.h"
#include "rspf/projection/rspfIkonosRpcModel.h"
#include "rspf/projection/rspfquickbirdrpcmodel.h"
#include "rspf/projection/rspfLandSatModel.h"
#include "rspf/support_data/rspfFfL5.h"
#include <rspf/support_data/rspfSpotDimapSupportData.h>
#include <rspf/projection/rspfSpot5Model.h> 
#include <rspf/projection/rspfProjection.h>
#include <rspf/projection/rspfMapProjectionFactory.h>
#include <rspf/projection/rspfProjectionFactoryRegistry.h>
#include <rspf/base/rspfStdOutProgress.h>
#include "rspf/base/rspfGpt.h"
#include "rspf/base/rspfDpt.h"
#include <rspf/imaging/rspfBandSelector.h>
#include <rspf/imaging/rspfPolyCutter.h>


//////gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include <direct.h>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <wtypes.h>
#include <atltypes.h>


bool leap_year(int iyear);
void parseDate(QString strDay, QDate& dt);

bool FindFile(const QString& path, const QStringList& filter, list<QString> &vFiles);
bool FindFile(const QString& path, const QStringList& filter, QStringList &vFiles);
QString QBeforeLast(const QString& str, QChar ch);
QString QAfterFirst(const QString& str, QChar ch);
QString QAfterLast(const QString& str, QChar ch);
QString QBeforeFirst(const QString& str, QChar ch);

bool CopyFolder(QString pathFrom, QString pathTo, QStringList filters = QStringList("*.*"), bool bOverWrite = true);

bool bandMerge(std::vector<string>& fileList,string createtiffile);

static QString getLastErrorMsg() {
	LPWSTR bufPtr = NULL;
	DWORD err = GetLastError();
	FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL, err, 0, (LPWSTR)&bufPtr, 0, NULL);
	const QString result = 
		(bufPtr) ? QString::fromUtf16((const ushort*)bufPtr).trimmed() :
		QString("Unknown Error %1").arg(err);
	LocalFree(bufPtr);
	return result;
}

class CopierWorker : public QThread { // only to be used by the Copier object
	BOOL m_stop;
	QString m_src, m_dst;
	QPointer<QObject> m_object;
	static DWORD CALLBACK copyProgress(
		LARGE_INTEGER totalSize, LARGE_INTEGER totalTransferred,
		LARGE_INTEGER streamSize, LARGE_INTEGER streamTransferred,
		DWORD streamNo, DWORD callbackReason, HANDLE src, HANDLE dst,
		LPVOID data)
	{
		Q_UNUSED(streamSize) Q_UNUSED(streamTransferred)
			Q_UNUSED(streamNo) Q_UNUSED(callbackReason)
			Q_UNUSED(src) Q_UNUSED(dst)
			QObject * object = static_cast<QObject*>(data);
		const QString text = QString("Transferred %1 of %2 bytes").
			arg(totalTransferred.QuadPart).arg(totalSize.QuadPart);
		QMetaObject::invokeMethod(object, "newStatus", Qt::QueuedConnection,
			Q_ARG(QString, text));
		return PROGRESS_CONTINUE;
	}
	void run() {
		m_stop = FALSE;
		BOOL rc = CopyFileExW((LPCWSTR)m_src.utf16(), (LPCWSTR)m_dst.utf16(),
			&copyProgress, m_object, &m_stop, 0);
		if (!rc) {
			QMetaObject::invokeMethod(m_object, "newStatus", Qt::QueuedConnection,
				Q_ARG(QString, getLastErrorMsg()));
		}
	}
	CopierWorker(const QString & src, const QString & dst, QObject * obj) :
		m_src(src), m_dst(dst), m_object(obj) {}
	void stop() { m_stop = TRUE; }
	friend class Copier;
};

class Copier : public QObject {
	Q_OBJECT
		QPointer<CopierWorker> m_worker;
public:
	Copier(const QString & src, const QString & dst) : m_worker(new CopierWorker(src, dst, this)) {
		connect(m_worker, SIGNAL(finished()), SIGNAL(finished()));
		connect(m_worker, SIGNAL(finished()), m_worker, SLOT(deleteLater()));
	}
	~Copier() {
		if (!m_worker) return;
		m_worker->stop();
		if (!m_worker->isRunning()) delete m_worker;
	}
	Q_SIGNAL void newStatus(const QString &);
	Q_SIGNAL void finished();
	Q_SLOT void stop() { m_worker->stop(); }
	void copy() { m_worker->start(); }
};

#endif