/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef OPENGLSCENE_H
#define OPENGLSCENE_H

#include <GL/glew.h>
#include <QGraphicsScene>
#include <QGraphicsProxyWidget>
#include <QPlainTextEdit>
#include <wrap/gui/trackball.h>
#include <wrap/gl/trimesh.h>
#include <vcg/math/histogram.h>

#include "qanimationwrapper.h"
#include "qmeshsequencefieldswrapper.h"

#ifndef QT_NO_CONCURRENT
    #include <QFutureWatcher>
#endif

class OpenGLScene : public QGraphicsScene
{
    Q_OBJECT
public:
    typedef enum{
        AVERAGE,
        MAXIMUM
    } RestFieldsKinds;

    explicit OpenGLScene(QObject *parent = 0);

    /// static object to send messages
    static const QString msg;
    
signals:
    void noOpenGLPainter(const QString &msg);
    void sendLogMessage(const QString &text);
//    void controlWidgetClosed(bool closed);
    void meshesOpening(bool opening);
    void meshesOpened(bool opened);
    void projectOpening(bool opening);
    void projectOpened(bool opened);
    void projectStoring(bool storing);
    void projectStored(bool stored);
    void settingsStoring(bool storing);
    void settingsStored(bool stored);
    void settingsOpening(bool opening);
    void settingsOpened(bool opened);
    void meshSelect(int i);
    void meshMaxIndex(int m);
    void persistentHardConstraintThresholdMaxIndex(int m);
    void computingFields(bool computing);
    void computedFields(bool computed);
    void computingRestFields(bool computing);
    void computedRestFields(bool computed);
    void smoothingFields(bool smoothing);
    void smoothedFields(bool smoothed);
    void restPoseSelect(int i);
    void restPoseToMeshSelect(int i);
    void restMeshFieldsKindChanged(int kind);
    void restMeshFieldsStretchTypeChanged(int stretchT);
    void updatedSingularValue2LowerBound(float l);
    void updatedSingularValue1UpperBound(float u);
    void updatedCurrentStretchMeasureLowerBound(float l);
    void updatedCurrentStretchMeasureUpperBound(float u);
    void anisotropyThresholdChanged(double a);
    void averageThresholdChanged(double a);
    void smootherAlphaChanged(double a);
    void stretchThresholdChanged(double t);
    void persistentHardConstraintThresholdChanged(int t);
    void updatingPersistentHardConstraintsVector(bool updating);
    void ffieldsStoring(bool storing);
    void ffieldsStored(bool stored);
    
public slots:
//    void showControlWidget();
//    void emitControlWidgetClosed();
    void showTrackBall(bool show);
    void setMeshType(int type);
    void openMeshes();
    void openProject();
    void storeProject();
    void storeSettings();
    void openSettings();
    void closeMeshes();
    void selectMesh(int i);
    void selectRestPose(int i);
    void selectRestPoseToMesh(int i);
    void selectColorThreshold(double t);
    void setAverageThreshold(double t);
    void setSmootherAlpha(double alpha);
    void setStretchThreshold(double t);
    void setAnisotropyThreshold(double t);
    void setPersistentHardConstraintThreshold(int t);
    void setLowerBound(float l);
    void setUpperBound(float u);
    void showRestFieldsSmoothed(bool show);
    void showQuadSingularities(bool show);
    void showQuadColorGold(bool gold);
    void storeFFields();
    void printCurrentScene();
    
protected:
    void drawBackground(QPainter *painter, const QRectF &rect); // draw OPENGL here
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);
    void wheelEvent(QGraphicsSceneWheelEvent *event);
    //    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event);

private:
    // widgets
    QGraphicsProxyWidget *controlWidget;
    QPlainTextEdit logPlainTextEdit;

    // for logs
    static OpenGLScene * scene;
    mutable QMutex mutex;
    void textToLog(const QString text);
    static void toLog(const string text);

    // objects of the scene
    QAnimationWrapper animation;
    QMeshSequenceFieldsWrapper fields;
    std::vector< std::pair<TMesh::CoordType,TMesh::CoordType> > smoothedFieldsVector;
    std::list<int> smoothedFieldsSingularitiesList;
    std::vector<unsigned long> persistentHardConstraintsVector;
    TMesh *currentTriMesh;
//    PMesh *currentPolyMesh;

    // objects and functions for drawing the scene
    vcg::Trackball track;
    vcg::GlTrimesh<TMesh> glTriWrap;
    bool viewTrackBall;
    bool viewFields;
    bool viewFaceColors;
    bool viewCurvature;
    bool viewMeshesCurvature;
    bool viewPersistentCurvature;
    bool viewRestFieldsSmoothed;
    bool viewFieldsAsOnlyPassingSoftConstraint;
    bool viewQuadSingularities;
    vcg::Color4b quadColor;
//    vector<vcg::Histogramf> histograms;
//    vector<unsigned long> maxCountsOfHistograms;
//    void updateHistograms();
    void setFacesColor(TMesh *mesh, unsigned long index);

    // state variables
    QString folder;
    QString projectName;
    vector<string> fnames;
    MeshTypes currentMeshType;
    RestFieldsKinds restFieldsKind;
    QMeshSequenceFieldsWrapper::StretchType stretchType;
    bool meshesLoaded;
    bool fieldsComputing;
    bool fieldsComputed;
    bool restFieldsComputing;
    bool restFieldsComputed;
    bool restFieldsSmoothing;
    bool restFieldsSmoothed;
    unsigned long currentMeshIndex;
    unsigned long restPoseTriMeshIndex;
    unsigned long importedRestPoseTriMeshIndex;
    unsigned long currentRestPoseToTriMeshIndex;
    unsigned long persistentHardConstraintThreshold;
    bool loadAll;
    float colorThreshold;
    double averageThreshold;
    double smootherAlpha;
    double stretchThreshold;
    double anisotropyThreshold;

    // loading and computing functions
    void computeMeshSequenceFields();
    void computeFieldsOnRestPoseMesh();
    void computeSmoothFields();
    void persistentHardConstraintThresholding();
    void saveFField(const char * ffieldFilename);

#ifndef QT_NO_CONCURRENT
    // functions for multi-threading
    static void loadMeshesParallel();
    static void loadProjectParallel(fstream::pos_type pos);
    static void storeProjectParallel();
    static void computeTriMeshesFieldsParallel();
    static void computeRestFieldsParallel();
    static void smoothFieldsParallel();
    static void persistentHardConstraintThresholdingParallel();

    // object for multi-threading
    QFutureWatcher<void> loadFutureWatcher;
    QFutureWatcher<void> computeFutureWatcher;
    QFutureWatcher<void> computeRestFutureWatcher;
    QFutureWatcher<void> smoothFutureWatcher;
    QFutureWatcher<void> persistentConstraintsFutureWatcher;
    QFutureWatcher<void> openFutureWatcher;
    QFutureWatcher<void> storeFutureWatcher;
#endif

public slots:
    // slot for loading and computing functions
    void importedMeshes();
    void computeTriMeshesFields();
    void computedTriMeshesFields();
    void setRestMeshFieldsKind(int kind);
    void setRestMeshFieldsStretchType(int stretchT);
    void computeRestMeshFields();
    void computedRestMeshFields();
    void smoothFields();
    void fieldsSmoothed();
    void importedProject();
    void exportedProject();
    void storedFFields();
    void thresholdPersistentHardContraints();
    void persistentHardConstraintThresholded();

    // slots for changing the mode of drawing the scene
    void selectDrawMode(int controlIndex);
    void selectFieldsDrawMode(int drawIndex);
    void showMeshesCurvature(bool showCurvature);
    void showFieldsAsOnlyPassingSoftConstraint(bool show);
//    void showFaceColors(bool showFaceColors);
};

#endif // OPENGLSCENE_H
