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

#include "openglscene.h"
#include "spinbox.h"
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QImage>
#include <QToolBox>
#include <QPlainTextEdit>
#include <QGroupBox>
#include <QFormLayout>
#include <QStyleFactory>
#include <QComboBox>
#include <QButtonGroup>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QSlider>
#include <QLabel>
#include <QPainter>
#include <QPaintEngine>
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QTime>
#include <QtAlgorithms>
#ifndef Q_NO_CONCURRENT
    #include <QtConcurrentRun>
#endif
//#include <qopengl.h>
#include <vcg/math/matrix33.h>
#include <wrap/qt/trackball.h>
#include <wrap/io_trimesh/import.h>
#include <vcg/complex/algorithms/bitquad_support.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <wrap/gl/addons.h>
#include "compute.h"
#include "SmootherWrapper.h"

/************************ static stuff **********************************************************************************/

/// initialize static const private member
const QString OpenGLScene::msg = OpenGLScene::tr("OpenGLScene needs OpenGL (>=2) as paint engine.");

/// initialize static private reference to 'this' OpenGLScene object (it will be changed when instantiating OpenGLScene)
OpenGLScene * OpenGLScene::scene = 0;


/**
 * @brief OpenGLScene::toLog appends text to logPlainTextEdit.
 * @param text The text to append.
 */
void OpenGLScene::toLog(const string text) {
    if (scene)
        scene->textToLog(QString::fromStdString(text));
}

#ifndef Q_NO_CONCURRENT
/**
 * @brief OpenGLScene::loadMeshesParallel
 */
void OpenGLScene::loadMeshesParallel() {
    if (scene)
        scene->animation.load(scene->fnames);
}

/**
 * @brief computedTriMeshesFieldsParallel
 */
void OpenGLScene::computeTriMeshesFieldsParallel() {
    if (scene) {
        // compute fields
        scene->computeMeshSequenceFields();

        // update histograms
//        scene->updateHistograms();
    }
}

/**
 * @brief OpenGLScene::computeRestFieldsParallel
 */
void OpenGLScene::computeRestFieldsParallel() {
    if (scene) {
        // compute rest fields
        scene->computeFieldsOnRestPoseMesh();

        // update histograms
//        scene->updateHistograms();
    }
}

/**
 * @brief OpenGLScene::smoothFieldsParallel
 */
void OpenGLScene::smoothFieldsParallel() {
    if (scene) {
        // smooth fields
        scene->computeSmoothFields();
    }
}

/**
 * @brief OpenGLScene::loadFieldsParallel
 * @param pos
 */
void OpenGLScene::loadProjectParallel(fstream::pos_type pos) {
    if (scene) {
        // if all meshes must be kept in memory
        if (scene->loadAll)
            // load all
            scene->animation.load(scene->fnames);
        else
            // else load only names
            scene->animation.setNames(scene->fnames);

        // load fields
        scene->fields.loadFields(scene->projectName.toStdString(), pos);

//        scene->fieldsComputed = true;
        // update histograms
//        scene->updateHistograms();
    }
}

/**
 * @brief OpenGLScene::storeProjectParallel
 */
void OpenGLScene::storeProjectParallel() {
    if (scene)
        scene->fields.storeFields(scene->projectName.toStdString());
}

/**
 * @brief OpenGLScene::persistentHardConstraintThresholdingParallel
 */
void OpenGLScene::persistentHardConstraintThresholdingParallel() {
    if (scene)
        scene->persistentHardConstraintThresholding();
}
#endif

/**
 * @brief setColorScaled scales color such that for min you have BLUE, and for max you have RED.
 * @param c
 * @param min
 * @param max
 * @param v
 */
template < typename ColorType >
void setColorScaled(vcg::Color4<ColorType> &c, float min, float max, float v, float th) {
    // check min and max
    if (max == min){
        c = vcg::Color4<ColorType>(vcg::Color4<ColorType>::Blue);
        return;
    }
    if (max < min) {
        setColorScaled(c, max, min, v, th);
        return;
    }

    // check threshold
    assert(th >= min && th <= max);

    // check value
    assert(v >= min && v <= max);

    // generate components
    float R = v <= th ? 0.0 : (v - th) / (max - th);
    float G = v < th ? (v - min) / (th - min) : v > th ? (max - v) / (max - th) : 1.0;
    float B = v < th ? (th - v) / (th - min) : 0.0;
    float A = 1.0;

    // set color
    c.Import(vcg::Color4f(R, G, B, A));
}




/******************** non-static stuff ***********************************************************************************/

/**
 * @brief OpenGLScene::OpenGLScene Contructor make up a scene to bee shown into a QGraphicsView.
 * @param parent The parent object.
 */
OpenGLScene::OpenGLScene(QObject *parent) :
    QGraphicsScene(parent)
{
    currentMeshType = NO_MESHTYPE;
    currentTriMesh = 0;
//    currentPolyMesh = 0;
    currentMeshIndex = 0;
    restPoseTriMeshIndex = 0;
    importedRestPoseTriMeshIndex = 0;
    currentRestPoseToTriMeshIndex = 0;
    restFieldsKind = MAXIMUM;
    stretchType = QMeshSequenceFieldsWrapper::RATIO;
    folder = QDir::homePath();
    meshesLoaded = false;
    loadAll = false;
    fieldsComputing = false;
    fieldsComputed = false;
    restFieldsComputing = false;
    restFieldsComputed = false;
    restFieldsSmoothed = false;
    viewTrackBall = true;
    viewFields = false;
    viewCurvature = false;
    viewMeshesCurvature = false;
    viewPersistentCurvature = false;
    viewRestFieldsSmoothed = false;
    viewFaceColors = false;
    viewFieldsAsOnlyPassingSoftConstraint = false;
    viewQuadSingularities = false;
    colorThreshold = 0.5;
    smootherAlpha = 0.5;
    stretchThreshold = 0.1;
    anisotropyThreshold = 0.8;
    scene = this;
    animation.setLogFunction(toLog);
    fields.setLogFunction(toLog);
    fields.setAbsSingularValue2LowerBound(0.1);
    fields.setAbsSingularValue1UpperBound(10);
    quadColor = vcg::Color4b(184,184,184,255);

#ifndef Q_NO_CONCURRENT
    connect(&loadFutureWatcher, SIGNAL(finished()), this, SLOT(importedMeshes()));
    connect(&computeFutureWatcher, SIGNAL(finished()), this, SLOT(computedTriMeshesFields()));
    connect(&computeRestFutureWatcher, SIGNAL(finished()), this, SLOT(computedRestMeshFields()));
    connect(&smoothFutureWatcher, SIGNAL(finished()), this, SLOT(fieldsSmoothed()));
    connect(&openFutureWatcher, SIGNAL(finished()), this, SLOT(importedProject()));
    connect(&storeFutureWatcher, SIGNAL(finished()), this, SLOT(exportedProject()));
    connect(&persistentConstraintsFutureWatcher, SIGNAL(finished()), this, SLOT(persistentHardConstraintThresholded()));
#endif

//    QToolBox *controlToolBox = new QToolBox;

//    QVBoxLayout *boxLayout = new QVBoxLayout;

//    QRadioButton *triMeshTypeRadioButton = new QRadioButton(tr("triangles"));
//    QRadioButton *polyMeshTypeRadioButton = new QRadioButton(tr("polygons"));
//    triMeshTypeRadioButton->setChecked(true);
//    QGroupBox *meshTypeGroupBox = new QGroupBox(tr("Mesh type"));
//    QButtonGroup *meshTypeButtonGroup = new QButtonGroup(meshTypeGroupBox);
//    meshTypeButtonGroup->addButton(triMeshTypeRadioButton, 1);
//    meshTypeButtonGroup->addButton(polyMeshTypeRadioButton, 2);
//    connect(meshTypeButtonGroup, SIGNAL(buttonPressed(int)), this, SLOT(setMeshType(int)));
//    connect(this, SIGNAL(meshesOpened(bool)), triMeshTypeRadioButton, SLOT(setEnabled(bool)));
//    connect(this, SIGNAL(meshesOpened(bool)), polyMeshTypeRadioButton, SLOT(setEnabled(bool)));
//    triMeshTypeRadioButton->setEnabled(false);
//    polyMeshTypeRadioButton->setEnabled(false);
//    QVBoxLayout *meshTypeGroupBoxLayout = new QVBoxLayout;
//    meshTypeGroupBoxLayout->addWidget(triMeshTypeRadioButton);
//    meshTypeGroupBoxLayout->addWidget(polyMeshTypeRadioButton);
//    meshTypeGroupBox->setLayout(meshTypeGroupBoxLayout);
//    boxLayout->addWidget(meshTypeGroupBox);

//    SpinBox *meshSelectionSpinBox = new SpinBox;
//    meshSelectionSpinBox->setSingleStep(1);
//    meshSelectionSpinBox->setMinimum(0);
//    meshSelectionSpinBox->setMaximum(0);
//    QVBoxLayout *meshSelectionGroupBoxLayout = new QVBoxLayout;
//    meshSelectionGroupBoxLayout->addWidget(meshSelectionSpinBox);
//    QGroupBox *meshSelectionGroupBox = new QGroupBox(tr("Mesh selection"));
//    meshSelectionGroupBox->setLayout(meshSelectionGroupBoxLayout);
//    boxLayout->addWidget(meshSelectionGroupBox);
//    meshSelectionSpinBox->setEnabled(false);
//    connect(this, SIGNAL(meshesOpened(bool)), meshSelectionSpinBox, SLOT(setEnabled(bool)));
//    connect(this, SIGNAL(meshSelect(int)), meshSelectionSpinBox, SLOT(setValue(int)));
//    connect(meshSelectionSpinBox, SIGNAL(valueChanged(int)), this, SLOT(selectMesh(int)));
//    connect(this, SIGNAL(meshMaxIndex(int)), meshSelectionSpinBox, SLOT(maximumChanged(int)));

//    QGroupBox *meshFeaturesGroupBox = new QGroupBox;
//    meshFeaturesGroupBox->setLayout(boxLayout);
//    controlToolBox->addItem(meshFeaturesGroupBox, tr("Mesh features"));

//    SpinBox *restPoseSelectionSpinBox = new SpinBox;
//    restPoseSelectionSpinBox->setSingleStep(1);
//    restPoseSelectionSpinBox->setMinimum(0);
//    restPoseSelectionSpinBox->setMaximum(0);
//    restPoseSelectionSpinBox->setEnabled(false);
//    connect(this, SIGNAL(meshesOpened(bool)), restPoseSelectionSpinBox, SLOT(setEnabled(bool)));
//    connect(this, SIGNAL(restPoseSelect(int)), restPoseSelectionSpinBox, SLOT(setValue(int)));
//    connect(restPoseSelectionSpinBox, SIGNAL(valueChanged(int)), this, SLOT(selectRestPose(int)));
//    connect(this, SIGNAL(meshMaxIndex(int)), restPoseSelectionSpinBox, SLOT(maximumChanged(int)));
//    connect(this, SIGNAL(computingFields(bool)), restPoseSelectionSpinBox, SLOT(setDisabled(bool)));
//    connect(this, SIGNAL(computedFields(bool)), restPoseSelectionSpinBox, SLOT(setDisabled(bool)));
//    QVBoxLayout *restPoseSelectionGroupBoxLayout = new QVBoxLayout;
//    restPoseSelectionGroupBoxLayout->addWidget(restPoseSelectionSpinBox);
//    QGroupBox *restPoseSelectionGroupBox = new QGroupBox(tr("Rest pose selection"));
//    restPoseSelectionGroupBox->setLayout(restPoseSelectionGroupBoxLayout);
//    QPushButton *computePushButton = new QPushButton(tr("Compute fields"));
//    computePushButton->setEnabled(false);
//    connect(this, SIGNAL(meshesOpened(bool)), computePushButton, SLOT(setEnabled(bool)));
//    connect(this, SIGNAL(computingFields(bool)), computePushButton, SLOT(setDisabled(bool)));
//    connect(this, SIGNAL(computedFields(bool)), computePushButton, SLOT(setDisabled(bool)));
//    connect(computePushButton, SIGNAL(clicked()), this, SLOT(computeTriMeshesFields()));
//    SpinBox *restPoseToMeshSelectionSpinBox = new SpinBox;
//    restPoseToMeshSelectionSpinBox->setSingleStep(1);
//    restPoseToMeshSelectionSpinBox->setMinimum(0);
//    restPoseToMeshSelectionSpinBox->setMaximum(0);
//    restPoseToMeshSelectionSpinBox->setEnabled(false);
//    connect(this, SIGNAL(meshesOpened(bool)), restPoseToMeshSelectionSpinBox, SLOT(setEnabled(bool)));
//    connect(this, SIGNAL(restPoseToMeshSelect(int)), restPoseToMeshSelectionSpinBox, SLOT(setValue(int)));
//    connect(restPoseToMeshSelectionSpinBox, SIGNAL(valueChanged(int)), this, SLOT(selectRestPoseToMesh(int)));
//    connect(this, SIGNAL(meshMaxIndex(int)), restPoseToMeshSelectionSpinBox, SLOT(maximumChanged(int)));
//    QVBoxLayout *restPoseToMeshSelectionGroupBoxLayout = new QVBoxLayout;
//    restPoseToMeshSelectionGroupBoxLayout->addWidget(restPoseToMeshSelectionSpinBox);
//    QGroupBox *restPoseToMeshSelectionGroupBox = new QGroupBox(tr("Reference mesh selection"));
//    restPoseToMeshSelectionGroupBox->setLayout(restPoseToMeshSelectionGroupBoxLayout);
//    QSlider *colorThresholdSlider = new QSlider(Qt::Horizontal);
//    colorThresholdSlider->setMinimum(0);
//    colorThresholdSlider->setMaximum(10);
//    colorThresholdSlider->setSingleStep(1);
//    colorThresholdSlider->setSliderPosition(5);
//    connect(colorThresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(selectColorThreshold(int)));
//    connect(this, SIGNAL(computedFields(bool)), colorThresholdSlider, SLOT(setEnabled(bool)));
//    QVBoxLayout *colorThresholdSelectionLayout = new QVBoxLayout;
//    colorThresholdSelectionLayout->addWidget(colorThresholdSlider);
//    QGroupBox *colorThresholdSelectionGroupBox = new QGroupBox(tr("Color Threshold"));
//    colorThresholdSelectionGroupBox->setLayout(colorThresholdSelectionLayout);
//    QVBoxLayout *computeGroupBoxLayout = new QVBoxLayout;
//    computeGroupBoxLayout->addWidget(restPoseSelectionGroupBox);
//    computeGroupBoxLayout->addWidget(computePushButton);
//    computeGroupBoxLayout->addWidget(restPoseToMeshSelectionGroupBox);
//    computeGroupBoxLayout->addWidget(colorThresholdSelectionGroupBox);
//    QGroupBox *computeGroupBox = new QGroupBox;
//    computeGroupBox->setLayout(computeGroupBoxLayout);
//    controlToolBox->addItem(computeGroupBox, tr("Compute Fields"));

//    logPlainTextEdit.setReadOnly(true);
//    logPlainTextEdit.ensureCursorVisible();
//    controlToolBox->addItem(&logPlainTextEdit, tr("Log"));

//    connect(controlToolBox, SIGNAL(currentChanged(int)), this, SLOT(selectDrawMode(int)));

//    controlWidget = addWidget(controlToolBox, Qt::Tool);
//    controlWidget->setWindowTitle(tr("Control"));
//    controlWidget->setOpacity(0.8);
//    controlWidget->setCacheMode(QGraphicsItem::DeviceCoordinateCache);
//    connect(controlWidget, SIGNAL(visibleChanged()), this, SLOT(emitControlWidgetClosed()));
//    showControlWidget();
//    QRectF geo = controlWidget->geometry();
//    geo.setWidth(geo.width() + 10);
//    geo.setHeight(geo.height() + 100);
//    controlWidget->setGeometry(geo.toRect());
}

/**
 * @brief OpenGLScene::drawBackground Draw on background OpenGL commands.
 * @param painter Painter the engine of which must be of type OpenGL.
 * @param rect The exposed rectangle. Unused.
 *
 * @overload
 */
void OpenGLScene::drawBackground(QPainter *painter, const QRectF &rect) {
    QGraphicsScene::drawBackground(painter, rect);

    // verify opengl compatibility
    if (painter->paintEngine()->type() != QPaintEngine::OpenGL
            && painter->paintEngine()->type() != QPaintEngine::OpenGL2
            && painter->paintEngine()->type() != QPaintEngine::Raster)
        emit noOpenGLPainter(msg);

    glewInit();
    glClearColor(255, 255, 255, 255);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
//    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.3);
//    glEnable(GL_LINE_SMOOTH);
//    glEnable(GL_BLEND);
//    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
//    glHint(GL_FOG_HINT, GL_NICEST);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, width()/height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center = vcg::Point3f(0, 0, 0);
    track.radius = 1;
    track.GetView();
    track.Apply();

    if (meshesLoaded) {
        switch (currentMeshType) {
        case TRI_MESHTYPE:
            if (currentTriMesh->VN() > 0) {
                glPushMatrix();
                vcg::glScale(2.0f/currentTriMesh->bbox.Diag());
                glTranslate(-currentTriMesh->bbox.Center());
                glColor(vcg::Color4b(198,226,255,255));
                if (viewFields && viewFaceColors && (restFieldsComputed || currentRestPoseToTriMeshIndex != restPoseTriMeshIndex))
                    glTriWrap.Draw(GLW::DMSmooth,GLW::CMPerFace,GLW::TMNone);
                else
                    glTriWrap.Draw(GLW::DMFlatWire,GLW::CMNone,GLW::TMNone);

                // render curvature
                if ((viewFields && !viewFaceColors && viewCurvature) || viewMeshesCurvature) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glDisable(GL_LIGHTING);
                    vcg::Point3f c, v;

        //                // draw points
        //                glPointSize(8);
        //                glBegin(GL_POINTS);
        //                for (int i = 0; i < currentTriMesh->VN(); i++) {
        //                    float tau = currentTriMesh->vert[i].Q();

        //                    if (tau >= anisotropyThreshold) {
        //                        // set color
        //                        vcg::Color4b color4b;
        //                        setColorScaled(color4b, 0.0, 1.0, tau, 0.1);
        //                        glColor(color4b);

        //                        glVertex(currentTriMesh->vert[i].P());
        //                    }
        //                }
        //                glEnd();

                    // temporary variables
                    float tau;
                    bool curvAll;

                    // for each face
                    for (unsigned long i = 0; i < (unsigned long)currentTriMesh->FN(); i++) {
                        // get values
                        tau = currentTriMesh->face[i].Q();

                        if (viewFields && viewPersistentCurvature)
                            curvAll = persistentHardConstraintsVector[i] >= persistentHardConstraintThreshold;
                        else
                            curvAll = tau >= anisotropyThreshold;

                        if (curvAll) {
                            // pick the center of the face
                            c = currentTriMesh->face[i].P(0);
                            for (int j = 1; j < 3; j++)
                                c = c + currentTriMesh->face[i].P(j);
                            c = c / 3.0;

                            // set color
                            vcg::Color4b color4b;
                            setColorScaled(color4b, 0.0, 1.0, tau, 0.5);
                            glColor(color4b);

                            // set line width [1..9] based on anisotropy
                            int lineWidth = 2 + 6 * tau;
                            glLineWidth(1);
//                            glEnable(GL_LINE_SMOOTH);
//                            glEnable(GL_BLEND);
//                            glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
//                            glHint(GL_FOG_HINT, GL_NICEST);

                            glBegin(GL_LINES);
                            // construct end point of line
                            v = currentTriMesh->face[i].PD1();
                            v = v * (currentTriMesh->bbox.Dim()[currentTriMesh->bbox.MinDim()] / 600);
                            v = v * lineWidth;
                            glVertex(c);
                            glVertex(c+v);
                            glVertex(c);
                            glVertex(c-v);
                            // construct end point of line
                            v = currentTriMesh->face[i].PD2();
                            v = v * (currentTriMesh->bbox.Dim()[currentTriMesh->bbox.MinDim()] / 600);
                            v = v * lineWidth;
                            glVertex(c);
                            glVertex(c+v);
                            glVertex(c);
                            glVertex(c-v);
                            glEnd();
                        }
                    }
                    glPopAttrib();
                }

                // render fields
                if (fieldsComputed && currentMeshIndex == restPoseTriMeshIndex
                        && viewFields&& !viewFaceColors&& !viewCurvature
                        && (restFieldsComputed || currentRestPoseToTriMeshIndex != restPoseTriMeshIndex)) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glDisable(GL_LIGHTING);
        //                glLineWidth(2);
        //                glBegin(GL_LINES);
        //                glDepthRange(0,0.997);
                    vcg::Point3f c, v;

                    // get min and max
                    float min = fields.getAbsStretchMeasureLowerBound(stretchType);
                    float max = fields.getAbsStretchMeasureUpperBound(stretchType);
        //                float min = fields.getMinAbsStretchMeasureAt(currentRestPoseToTriMeshIndex, stretchType);
        //                float max = fields.getMaxAbsStretchMeasureAt(currentRestPoseToTriMeshIndex, stretchType);
        //                std::cout << min << " " << max << " " << colorThreshold << std::endl;

                    // for each face
                    for (int i = 0; i < currentTriMesh->FN(); i++) {
                        // get stretch as weight
                        float w = fields.getClampedAbsStretchMeasureAt(currentRestPoseToTriMeshIndex, (unsigned long)i, stretchType);
                        // normalize
                        w = (w - min) / (max - min > 0.0 ? max - min : 1.0);

                        // persistent?
                        bool persistent = persistentHardConstraintsVector[i] >= persistentHardConstraintThreshold;

                        // set soft constraint
                        if (viewFieldsAsOnlyPassingSoftConstraint)
                            if (w < stretchThreshold || smootherAlpha == 0.0 || persistent)
                                continue;

                        // pick the center of the face
                        c = currentTriMesh->face[i].P(0);
                        for (int j = 1; j < 3; j++)
                            c = c + currentTriMesh->face[i].P(j);
                        c = c / 3.0;

                        // get values
                        float s = fields.getClampedAbsStretchMeasureAt(currentRestPoseToTriMeshIndex, i, stretchType);

                        // set color
                        vcg::Color4b color4b;
        //                    color4b.SetColorRamp(min, max, s);
        //                    color4b.SetColorRamp(0, maxCountsOfHistograms[currentRestPoseToTriMeshIndex],
        //                                         histograms[currentRestPoseToTriMeshIndex].BinCount(s));
                        setColorScaled(color4b, min, max, std::fabs(s), colorThreshold);
                        glColor(color4b);

                        // set line width [1..9] based on stretch
                        float range = max - min;
                        int lineWidth = 1 + 3 * (std::fabs(s) - min) / (range > 0.0 ? range : 1);
                        glLineWidth(1);
                        glEnable(GL_LINE_SMOOTH);
                        glEnable(GL_BLEND);
                        glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
                        glHint(GL_FOG_HINT, GL_NICEST);

                        glBegin(GL_LINES);
                        // construct end point of line
                        if (restFieldsSmoothed && viewRestFieldsSmoothed && currentRestPoseToTriMeshIndex == restPoseTriMeshIndex)
                            v = smoothedFieldsVector[i].first;
                        else
                            fields.getFirstFieldAt(currentRestPoseToTriMeshIndex, i, v);
                        v = v * (currentTriMesh->bbox.Dim()[currentTriMesh->bbox.MinDim()] / 400);
                        v = v * lineWidth;
                        glVertex(c);
                        glVertex(c+v);
                        glVertex(c);
                        glVertex(c-v);
                        // construct end point of line
                        if (restFieldsSmoothed && viewRestFieldsSmoothed && currentRestPoseToTriMeshIndex == restPoseTriMeshIndex)
                            v = smoothedFieldsVector[i].second;
                        else
                            fields.getSecondFieldAt(currentRestPoseToTriMeshIndex, i, v);
                        v = v * (currentTriMesh->bbox.Dim()[currentTriMesh->bbox.MinDim()] / 400);
                        v = v * lineWidth;
                        glVertex(c);
                        glVertex(c+v);
                        glVertex(c);
                        glVertex(c-v);
                        glEnd();
                    }
        //                glEnd();
                    glPopAttrib();

                    // draw singularities
                    if (viewQuadSingularities && restFieldsSmoothed && viewRestFieldsSmoothed && currentRestPoseToTriMeshIndex == restPoseTriMeshIndex) {
                        glPushAttrib(GL_ALL_ATTRIB_BITS);
//                        glDisable(GL_LIGHTING);
//                        glPointSize(8);
                        vcg::glColor(vcg::Color4b(255,0,0,255));
//                        glBegin(GL_POINTS);
                            std::list<int>::iterator it = smoothedFieldsSingularitiesList.begin();
                            for (; it != smoothedFieldsSingularitiesList.end(); it++) {
//                                glNormal(currentTriMesh->vert[*it].N());
//                                glVertex(currentTriMesh->vert[*it].P());
                                vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(currentTriMesh->vert[*it].P(),
                                        currentTriMesh->bbox.Diag()/200);
                            }
//                        glEnd();
                        glPopAttrib();
                    }
                }
                glPopMatrix();
            }
            break;

        case POLY_MESHTYPE:
//            if (currentPolyMesh->VN() > 0) {
            if (currentTriMesh->VN() > 0) {
                glPushAttrib(GL_ALL_ATTRIB_BITS);
                glEnable(GL_LIGHTING);
                glEnable(GL_LIGHT0);
                glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glEnable(GL_NORMALIZE);
                glEnable(GL_COLOR_MATERIAL);
                glDisable(GL_TEXTURE_2D);
                glColor(quadColor);
//                glFrontFace(GL_CW);
                glDepthRange(0,0.998);
                glPushMatrix();
//                vcg::glScale(2.0f/currentPolyMesh->bbox.Diag());
//                glTranslate(-currentPolyMesh->bbox.Center());
                vcg::glScale(2.0f/currentTriMesh->bbox.Diag());
                glTranslate(-currentTriMesh->bbox.Center());

                glTriWrap.SetHint(vcg::GlTrimesh<TMesh>::HNIsPolygonal);
                glTriWrap.Draw(GLW::DMFlatWire,GLW::CMNone,GLW::TMNone);

                if (viewQuadSingularities) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
//                    glDisable(GL_LIGHTING);
//                    glPointSize(8);
//                    glBegin(GL_POINTS);
                    for (int i = 0; i < currentTriMesh->VN(); i++) {
                        if (vcg::tri::BitQuad<TMesh>::GetValency(&currentTriMesh->vert[i]) == 3) {
                            vcg::glColor(vcg::Color4b(255,0,0,255));
                            vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(currentTriMesh->vert[i].P(),
                                                                         currentTriMesh->bbox.Diag()/200);
                        }
                        if (vcg::tri::BitQuad<TMesh>::GetValency(&currentTriMesh->vert[i]) == 5) {
                            vcg::glColor(vcg::Color4b(0,0,255,255));
                            vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(currentTriMesh->vert[i].P(),
                                                                         currentTriMesh->bbox.Diag()/200);
                        }
                    }
//                    glEnd();
                    glPopAttrib();
                }

//                for (int i = 0; i < currentPolyMesh->FN(); i++)
//                {
//                    vcg::Point3f n(0,0,0);
//                    for (int j = 0; j < currentPolyMesh->face[i].VN(); j++) {
//                        PMesh::VertexPointer v = currentPolyMesh->face[i].V(j);
//                        n += v->N();
//                    }
//                    n /= currentPolyMesh->face[i].VN();

//                    glBegin(GL_POLYGON);
//                    for (int j = 0; j < currentPolyMesh->face[i].VN(); j++)
//                    {
//                        PMesh::VertexPointer v = currentPolyMesh->face[i].V(j);
//                        glNormal(n);
//                        glVertex(v->P());
//                    }
//                    glEnd();
//                }

//                glDepthRange(0,0.997);
//                glDisable(GL_LIGHTING);
//        //            glEnable(GL_COLOR_MATERIAL);
//                glColor3d(0,0,0);
//                for (int i = 0; i < currentPolyMesh->FN(); i++)
//                {
//                    glBegin(GL_LINE_LOOP);
//                    for (int j = 0; j < currentPolyMesh->face[i].VN(); j++)
//                    {
//                        PMesh::VertexPointer v = currentPolyMesh->face[i].V(j);
//                        glVertex(v->P());
//                    }
//                    glEnd();
//                }

                glPopMatrix();
                glPopAttrib();
            }
            break;

        default:
            break;
        }
    }

    if (viewTrackBall)
        track.DrawPostApply();
}

///**
// * @brief OpenGLScene::showControlWidget Show the control widget.
// */
//void OpenGLScene::showControlWidget() {
//    controlWidget->setPos(10 - controlWidget->boundingRect().x(), 10 - controlWidget->boundingRect().y());
//    controlWidget->setVisible(true);
//}

///**
// * @brief OpenGLScene::emitControlWidgetDestroyed Emit a signal that the control widget has been closed.
// */
//void OpenGLScene::emitControlWidgetClosed() {
//    emit controlWidgetClosed(true);
//}

/**
 * @brief OpenGLScene::showTrackBall
 * @param show
 */
void OpenGLScene::showTrackBall(bool show) {
    // set flag
    viewTrackBall = show;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::printCurrentScene
 */
void OpenGLScene::printCurrentScene() {
    // read pixels from buffer
    unsigned char *pixelsInv = new unsigned char[4*(int)width()*(int)height()];
    glReadPixels(0, 0, (GLint)width(), (GLint)height(), GL_BGRA, GL_UNSIGNED_BYTE, pixelsInv);

    // create an image
    QImage image(pixelsInv, width(), height(), QImage::Format_ARGB32);

    // swap upside down
    image = image.mirrored(false, true);

    // save the image on file
    image.save(QString::fromStdString(fnames[currentMeshIndex]) + ".png", "png");

    delete[] pixelsInv;
}

/**
 * @brief OpenGLScene::showQuadColorGold
 * @param gold
 */
void OpenGLScene::showQuadColorGold(bool gold) {
    // change color
    if (gold)
        quadColor = vcg::Color4b(255,193,37,255);
    else
        quadColor = vcg::Color4b(184,184,184,255);

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::selectDrawMode
 * @param controlIndex
 */
void OpenGLScene::selectDrawMode(int controlIndex) {
    // switch between control toolbox items
    switch (controlIndex) {
    case 0:     // first tab - mesh features
        viewFields = false;
        break;

    case 1:
        viewFields = true;
        if (currentMeshIndex != restPoseTriMeshIndex)
            selectRestPose(restPoseTriMeshIndex);
        break;

    default:
        break;
    }

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::selectFieldsDrawMode
 * @param drawIndex
 */
void OpenGLScene::selectFieldsDrawMode(int drawIndex) {
    // switch between draw mode indexes
    switch (drawIndex) {
    case 1:     // view fields
        viewFaceColors = false;
        viewCurvature = false;
        break;

    case 2:     // view face colors
        viewFaceColors = true;
        viewCurvature = false;
        break;

    case 3:     // view curvature
        viewFaceColors = false;
        viewCurvature = true;
        break;

    default:
        break;
    }

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::showMeshesCurvature
 * @param showCurvature
 */
void OpenGLScene::showMeshesCurvature(bool showCurvature) {
    // set flag
    viewMeshesCurvature = showCurvature;

    // update the scene
    update();
}

///**
// * @brief OpenGLScene::showFaceColors
// * @param showFaceColors
// */
//void OpenGLScene::showFaceColors(bool showFaceColors) {
//    // set flag
//    viewFaceColors = showFaceColors;

//    // update the scene
//    update();
//}

/**
 * @brief OpenGLScene::textToLog Append text to logPlainTextEdit.
 * @param text The text to append.
 */
void OpenGLScene::textToLog(const QString text) {
    QMutexLocker locker(&mutex);
//    QTextCursor cursor = logPlainTextEdit.textCursor();
//    int prevPos = cursor.position();
//    logPlainTextEdit.appendPlainText(text);
//    cursor.setPosition(prevPos + text.length());
//    while(!cursor.atEnd())
//        cursor.setPosition(cursor.position() + 1);
//    logPlainTextEdit.setTextCursor(cursor);

    emit sendLogMessage(text);
}

/**
 * @brief OpenGLScene::openMeshes Open mesh(es) the user select.
 */
void OpenGLScene::openMeshes() {
    // get file names
    QStringList filenames = QFileDialog::getOpenFileNames(0, tr("Load mesh(es)..."), folder, tr("Meshes (*.ply *.obj *.off *.stl *.vmi)"));
    unsigned long numOfMeshes = filenames.size();


    // if at least one is selected
    if (numOfMeshes > 0) {
        // update state
        closeMeshes();

        // extract std strings
        fnames.clear();
        fnames.resize(numOfMeshes);
        for (unsigned long i = 0; i < numOfMeshes; i++)
            fnames[i] = filenames[i].toStdString();

        // update UI
        emit meshesOpening(true);

        // ask if all meshes must be kept in memory
        QMessageBox::StandardButton responce = QMessageBox::question(0, tr("Loading method"), tr("Keep all meshes on memory?"),
                                                                     QMessageBox::StandardButtons(QMessageBox::No | QMessageBox::Yes),
                                                                     QMessageBox::Yes);
        if (responce == QMessageBox::Yes)
            loadAll = true;
        else
            loadAll = false;

#ifndef Q_NO_CONCURRENT
        // cancel and wait loading threads
        if (loadFutureWatcher.isRunning()) {
            loadFutureWatcher.cancel();
            loadFutureWatcher.waitForFinished();
        }
        // if all meshes must be kept in memory
        if (loadAll)
            // load in another thread
            loadFutureWatcher.setFuture(QtConcurrent::run(loadMeshesParallel));
        else {
            // else load only names
            animation.setNames(fnames);

            // terminate updating current state
            importedMeshes();
        }
#else
        // load animation
        animation.load(fnames);

        // terminate loading by updating current state
        importedMeshes();
#endif
    }
}

/**
 * @brief OpenGLScene::importedMeshes End of opening the meshes. Update current state.
 */
void OpenGLScene::importedMeshes() {
    // set current mesh
    if (currentMeshType == NO_MESHTYPE)
        currentMeshType = TRI_MESHTYPE;
    selectMesh(0);

    // update status
    meshesLoaded = true;

    // reset trackball
//    track.Reset();
    track.ButtonUp(QT2VCG(Qt::NoButton, Qt::ControlModifier));
    track.ButtonUp(QT2VCG(Qt::NoButton, Qt::ShiftModifier));
    track.ButtonUp(QT2VCG(Qt::NoButton, Qt::AltModifier));

    // log
    textToLog(QString("Meshes opened.\n"));

    // update folder
    QFileInfo fileInfo(QString::fromStdString(fnames.front()));
    folder = QDir::toNativeSeparators(fileInfo.path() + "/");

    // update UI
    emit meshMaxIndex(animation.getNumberOfMeshes() - 1);
    emit persistentHardConstraintThresholdMaxIndex(animation.getNumberOfMeshes() + 1);
    emit meshSelect(0);
    emit restPoseSelect(0);
    emit restPoseToMeshSelect(0);
    emit meshesOpening(false);
    emit meshesOpened(true);
    emit projectStored(false);

    // update scene
    update();
}

/**
 * @brief OpenGLScene::closeMeshes Delete all meshes.
 */
void OpenGLScene::closeMeshes() {
    if (animation.getNumberOfMeshes() == 0)
        return;

    // if not all meshes are kept in memory
    if (!loadAll && animation.getNumberOfMeshes() > 0) {
        // unlock current mesh (one for trimesh and one for polymesh)
        animation.unlockAt(currentMeshIndex);
        animation.unlockAt(currentMeshIndex);
    }

    // unload animation
    animation.unload();

    // clear fields
    fields.clear();

    // unset indexes and pointers
    currentMeshIndex = 0;
    restPoseTriMeshIndex = 0;
    currentRestPoseToTriMeshIndex = 0;
    currentTriMesh = 0;
//    currentPolyMesh = 0;

    // unset opengl wrapper
    glTriWrap.m = 0;
    glTriWrap.Update();

    // reset trackball
//    track.Reset();

    // reset state
    meshesLoaded = false;
    fieldsComputing = false;
    fieldsComputed = false;
    restFieldsComputing = false;
    restFieldsComputed = false;
    restFieldsSmoothed = false;

    // clear names
    fnames.clear();

    // clear persistent hard constraint vector
    persistentHardConstraintsVector.clear();

    // log
    textToLog(QString("Meshes closed.\n"));

    // update UI
    emit meshMaxIndex(0);
    emit persistentHardConstraintThresholdMaxIndex(1);
    emit meshSelect(0);
    emit restPoseSelect(0);
    emit restPoseToMeshSelect(0);
    emit meshMaxIndex(0);
    emit projectOpened(false);
    emit projectStored(false);
    emit computedFields(false);
    emit computedRestFields(false);
    emit ffieldsStored(true);
    emit meshesOpened(false);

    // update histograms
//    updateHistograms();

    // update scene
    update();
}

/**
 * @brief OpenGLScene::setMeshType Set the current mesh type.
 * @param type Type index. Must be either 0 or 1.
 */
void OpenGLScene::setMeshType(int type) {
    // update state
    switch(type) {
    case 1:
        currentMeshType = TRI_MESHTYPE;
        break;

    case 2:
        currentMeshType = POLY_MESHTYPE;
        break;

    default:
        break;
    }

    // update the scene
    update();
}


/**
 * @brief OpenGLScene::selectMesh Select one mesh from vectors.
 * @param i Index of mesh.
 */
void OpenGLScene::selectMesh(int i) {
    // of course it must be inside of range
    if ((unsigned long)i < animation.getNumberOfMeshes()) {
        // if not all meshes are kept in memory
        if (!loadAll) {
            // unlock current mesh (one for trimesh and one for polymesh)
            animation.unlockAt(currentMeshIndex);
            animation.unlockAt(currentMeshIndex);

            // when not computing, current mesh (if equal to rest) can be unloaded
            if ((!fieldsComputing && !restFieldsComputing) || currentMeshIndex != restPoseTriMeshIndex)
                animation.unloadMeshAt(currentMeshIndex);

            if (restFieldsSmoothed && currentMeshIndex == restPoseTriMeshIndex) {
                restFieldsSmoothed = false;
                emit smoothedFields(false);
            }
        }

        // set current mesh index
        currentMeshIndex = i;

        // set index of rest mesh
        if (!fieldsComputing && !fieldsComputed)
            restPoseTriMeshIndex = i;

        // (load and) get current tri mesh
        currentTriMesh = animation.getTriMeshAt(i);

        // set color to faces
        if (fieldsComputed)
            if (currentMeshIndex != currentRestPoseToTriMeshIndex || restFieldsComputed)
                setFacesColor(currentTriMesh, currentRestPoseToTriMeshIndex);

        // (already loaded by getTriMeshAt(i)) get current poly mesh
//        currentPolyMesh = animation.getPolyMeshAt(i);

        // set opengl wrapper for tri mesh
        glTriWrap.m = currentTriMesh;
        glTriWrap.Update();

        // update scene
        update();

        // update UI
        if (!fieldsComputing && !fieldsComputed)
            emit restPoseSelect(i);
    }
}

/**
 * @brief OpenGLScene::selectRestPose Select rest pose mesh.
 * @param i Index of mesh.
 */
void OpenGLScene::selectRestPose(int i) {
    // of course it must be inside of range
    if ((unsigned long)i < animation.getNumberOfMeshes()) {
        // if not all meshes are kept in memory
        if (!loadAll) {
            // unlock current mesh (one for trimesh and one for polymesh)
            animation.unlockAt(currentMeshIndex);
            animation.unlockAt(currentMeshIndex);

            // when not computing, current mesh (if equal to rest) can be unloaded
            if ((!fieldsComputing && !restFieldsComputing) || currentMeshIndex != restPoseTriMeshIndex)
                animation.unloadMeshAt(currentMeshIndex);
        }

        // set current mesh index
        currentMeshIndex = i;

        // set index of rest mesh
        restPoseTriMeshIndex = i;

        // (load and) get current tri mesh
        currentTriMesh = animation.getTriMeshAt(i);

        // (already loaded by getTriMeshAt(i)) get current poly mesh
//        currentPolyMesh = animation.getPolyMeshAt(i);

        // set opengl wrapper for tri mesh
        glTriWrap.m = currentTriMesh;
        glTriWrap.Update();

        // update scene
        update();

        // update UI
        emit meshSelect(i);
    }
}

/**
 * @brief OpenGLScene::selectRestPoseToMesh Select current mesh the rest pose to map to.
 * @param i Index of mesh.
 */
void OpenGLScene::selectRestPoseToMesh(int i) {
    // of course it must be inside of range
    if ((unsigned long)i < animation.getNumberOfMeshes()) {
        // set current reference mesh
        currentRestPoseToTriMeshIndex = i;
    }

    // set color to faces
    if (fieldsComputed)
        if (currentMeshIndex != currentRestPoseToTriMeshIndex || restFieldsComputed)
            setFacesColor(currentTriMesh, currentRestPoseToTriMeshIndex);

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setLowerBound
 * @param l
 */
void OpenGLScene::setLowerBound(float l) {
    // update color threshold
    fields.setAbsSingularValue2LowerBound(l);

    // update UI
    emit updatedCurrentStretchMeasureLowerBound(fields.getAbsStretchMeasureLowerBound(stretchType));

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setUpperBound
 * @param u
 */
void OpenGLScene::setUpperBound(float u) {
    // update color threshold
    fields.setAbsSingularValue1UpperBound(u);

    // update UI
    emit updatedCurrentStretchMeasureUpperBound(fields.getAbsStretchMeasureUpperBound(stretchType));

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::selectColorThreshold
 * @param t
 */
void OpenGLScene::selectColorThreshold(double t) {
    // check parameter
    assert(t >= fields.getAbsStretchMeasureLowerBound(stretchType)
           && t <= fields.getAbsStretchMeasureUpperBound(stretchType));

    // update color threshold
    colorThreshold = t;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setAverageThreshold
 * @param t
 */
void OpenGLScene::setAverageThreshold(double t) {
    // set threshold
    if (t >= fields.getAbsStretchMeasureLowerBound(stretchType)
            && t <= fields.getAbsStretchMeasureUpperBound(stretchType))
        averageThreshold = t;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setSmootherAlpha
 * @param alpha
 */
void OpenGLScene::setSmootherAlpha(double alpha) {
    // check parameter
    assert(alpha>= 0.0 && alpha <= 1.0);

    // set smoother alpha
    smootherAlpha = alpha;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setStretchThreshold
 * @param t
 */
void OpenGLScene::setStretchThreshold(double t) {
    // check parameter
    assert(t>= 0.0 && t <= 1.0);

    // set smoother alpha
    stretchThreshold = t;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setAnisotropyThreshold
 * @param t
 */
void OpenGLScene::setAnisotropyThreshold(double t) {
    // set anisotropy threshold
    if (t < 0.0)
        anisotropyThreshold = 0.0;
    else if (t > 1.0)
        anisotropyThreshold = 1.0;
    else
        anisotropyThreshold = t;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setPersistentHardConstraintThreshold
 * @param t
 */
void OpenGLScene::setPersistentHardConstraintThreshold(int t) {
    // check parameter
    if (t >= 0)
        persistentHardConstraintThreshold = (unsigned long) t;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setFacesColor
 * @param mesh
 * @param index
 */
void OpenGLScene::setFacesColor(TMesh *mesh, unsigned long index) {
    if ((unsigned long)mesh->FN() != fields.getNumberOfFieldsAt(index))
        return;

    // set color to mesh faces
    // get min and max
    float min = fields.getAbsStretchMeasureLowerBound(stretchType);
    float max = fields.getAbsStretchMeasureUpperBound(stretchType);
    float s;

    // for each face
    for (int i = 0; i < mesh->FN(); i++) {
        // get values
        s = fields.getClampedAbsStretchMeasureAt(index, i, stretchType);

        // set color
        vcg::Color4b color4b;
        setColorScaled(color4b, min, max, std::fabs(s), colorThreshold);
        mesh->face[i].C() = color4b;
    }
}

/**
 * @brief OpenGLScene::showRestFieldsSmoothed
 * @param show
 */
void OpenGLScene::showRestFieldsSmoothed(bool show) {
    // set flag
    viewRestFieldsSmoothed = show;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::showFieldsAsOnlyPassingSoftConstraint
 * @param show
 */
void OpenGLScene::showFieldsAsOnlyPassingSoftConstraint(bool show) {
    // set flag
    viewFieldsAsOnlyPassingSoftConstraint = show;

    // update the scene
    update();
}

void OpenGLScene::showQuadSingularities(bool show) {
    // set flag
    viewQuadSingularities = show;

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::computeTriMeshesFields
 */
void OpenGLScene::computeTriMeshesFields() {
    // check if already computed
    if (fieldsComputing)
        return;

    // update state
    fieldsComputing = true;

    // update UI
    emit computingFields(true);

#ifndef Q_NO_CONCURRENT
    // cancel and wait loading threads
    if (computeFutureWatcher.isRunning()) {
        computeFutureWatcher.cancel();
        computeFutureWatcher.waitForFinished();
    }

    // compute in another thread
    computeFutureWatcher.setFuture(QtConcurrent::run(computeTriMeshesFieldsParallel));
#else
    // compute fields in this thread
    computeMeshSequenceFields();

    // update histograms
//    updateHistograms();

    // update current state
    computedTriMeshesFields();
#endif
}

/**
 * @brief OpenGLScene::computeMeshSequenceFields
 */
void OpenGLScene::computeMeshSequenceFields() {
    // get number of meshes
    unsigned long numOfMeshes = animation.getNumberOfMeshes();

    // check if some meshes are present
    if (numOfMeshes == 0) {
        textToLog(tr("Error. No meshes to compute.\n"));
        return;
    }

    // start time
    QTime time;
    time.start();

    // reset mesh fields
    fields.resizeNumOfMeshes(numOfMeshes);

    // state variables
    bool wasRestMeshLoaded;
    bool wasMeshLoaded;

    // check if rest mesh is already loaded
    wasRestMeshLoaded = animation.isMeshLoadedAt(restPoseTriMeshIndex);

    // get rest mesh
    TMesh * restMesh = animation.getTriMeshAt(restPoseTriMeshIndex);

    // check rest mesh
    if (!restMesh) {
        textToLog(tr("Error. Rest mesh unavailable.\n"));
        return;
    }

    // resize persistent hard constraints vector
    persistentHardConstraintsVector.resize(restMesh->FN());

    // reference mesh variable
    TMesh * referenceMesh = 0;

    // for each mesh
    for (unsigned long i = 0; i < numOfMeshes; i++) {
        // don't compute for same mesh
        if (i == restPoseTriMeshIndex)
            continue;

        // check if i-th mesh is already loaded
        wasMeshLoaded = animation.isMeshLoadedAt(i);

        // get i-th mesh
        referenceMesh = animation.getTriMeshAt(i);

        // compute fields
        fields.computeMeshFieldsAt(i, *restMesh, *referenceMesh);

        // unlock i-th mesh
        animation.unlockAt(i);

        // unload if necessary
        if (!loadAll && !wasMeshLoaded)
            animation.unloadMeshAt(i);

        // unset reference mesh
        referenceMesh = 0;
    }

    // unlock rest mesh
    animation.unlockAt(restPoseTriMeshIndex);

    // unset rest mesh
    restMesh = 0;

    // unload rest mesh if necessary
    if (!loadAll && !wasRestMeshLoaded)
        animation.unloadMeshAt(restPoseTriMeshIndex);

    // get elapsed time
    QTime elapsedTime(0, 0, 0, 0);
    elapsedTime = elapsedTime.addMSecs(time.elapsed());
    // log
    textToLog(tr("\nFields computation time: ") + elapsedTime.toString("HH:mm:ss.zzz") + "\n");

    // update status
//    fieldsComputed = true;
}

/**
 * @brief OpenGLScene::computedTriMeshesFields When finished to compute the fields, update the scene and the UI.
 */
void OpenGLScene::computedTriMeshesFields() {
    // log
    textToLog(tr("Computed triangle meshes fields.\n"));

    // update state
    fieldsComputing = false;
    fieldsComputed = true;

    // update UI
    emit computingFields(false);
    emit computedFields(true);

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setRestMeshFieldsKind
 * @param kind
 */
void OpenGLScene::setRestMeshFieldsKind(int kind) {
    switch(kind) {
    case 1:
        restFieldsKind = MAXIMUM;
        break;

    case 2:
        restFieldsKind = AVERAGE;
        break;

    default:
        break;
    }

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::setRestMeshFieldsStretchType
 * @param stretchT
 */
void OpenGLScene::setRestMeshFieldsStretchType(int stretchT) {
    switch(stretchT) {
    case 1:
        stretchType = QMeshSequenceFieldsWrapper::RATIO;
        break;

    case 2:
        stretchType = QMeshSequenceFieldsWrapper::MAXWARP;
        break;

    case 3:
        stretchType = QMeshSequenceFieldsWrapper::SINGVAL1;
        break;

    default:
        break;
    }

    // update UI
    emit updatedCurrentStretchMeasureLowerBound(fields.getAbsStretchMeasureLowerBound(stretchType));
    emit updatedCurrentStretchMeasureUpperBound(fields.getAbsStretchMeasureUpperBound(stretchType));

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::computeRestMeshFields
 */
void OpenGLScene::computeRestMeshFields() {
    // check if fields are already computed
    if (!fieldsComputed) {
        textToLog(tr("Error. Fields not computed, impossible to compute average.\n"));
        return;
    }

    // check if already computed
    if (restFieldsComputing)
        return;

    // update state
    restFieldsComputing = true;
    restFieldsComputed = false;
    restFieldsSmoothed = false;

    // update UI
    emit computingRestFields(true);

    // update the scene
    update();

#ifndef Q_NO_CONCURRENT
    // cancel and wait loading threads
    if (computeRestFutureWatcher.isRunning()) {
        computeRestFutureWatcher.cancel();
        computeRestFutureWatcher.waitForFinished();
    }

    // compute in another thread
    computeRestFutureWatcher.setFuture(QtConcurrent::run(computeRestFieldsParallel));
#else
    // compute rest fields
    computeFieldsOnRestPoseMesh();

    // update histograms
//    updateHistograms();

    // update current state
    computeFieldsOnRestPoseMesh();
#endif
}

/**
 * @brief OpenGLScene::computeFieldsOnRestPoseMesh
 */
void OpenGLScene::computeFieldsOnRestPoseMesh() {
    // log
    textToLog(tr("Computing rest pose mesh fields ... "));

    // start time
    QTime time;
    time.start();

    // check if rest mesh is already loaded
    bool wasRestMeshLoaded = animation.isMeshLoadedAt(restPoseTriMeshIndex);

    // get rest mesh
    TMesh * restMesh = animation.getTriMeshAt(restPoseTriMeshIndex);

    // compute average fields on rest pose mesh
    switch(restFieldsKind) {
    case AVERAGE:
        fields.computeAverageFieldsAt(restPoseTriMeshIndex, *restMesh, stretchType, averageThreshold);
        break;

    case MAXIMUM:
        fields.computeMaximalFieldsAt(restPoseTriMeshIndex, *restMesh, stretchType);
        break;

    default:
        assert(false);
        break;
    }

    // set color to faces
    if (loadAll || wasRestMeshLoaded)
        setFacesColor(restMesh, restPoseTriMeshIndex);

    // unlock rest mesh
    animation.unlockAt(restPoseTriMeshIndex);

    // unset rest mesh
    restMesh = 0;

    // unload rest mesh if necessary
    if (!loadAll && !wasRestMeshLoaded)
        animation.unloadMeshAt(restPoseTriMeshIndex);

    // get elapsed time
    QTime elapsedTime(0, 0, 0, 0);
    elapsedTime = elapsedTime.addMSecs(time.elapsed());
    // log
    textToLog(tr("Rest pose mesh fields computation time: ") + elapsedTime.toString("HH:mm:ss.zzz") + "\n");
}


/**
 * @brief OpenGLScene::computedRestFields
 */
void OpenGLScene::computedRestMeshFields() {
    // log
    textToLog(tr("Computed rest mesh fields.\n"));

    // update state
    restFieldsComputing = false;
    restFieldsComputed = true;

    // update UI
    emit computingRestFields(false);
    emit computedRestFields(true);

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::smoothFields
 */
void OpenGLScene::smoothFields() {
    // check if already computed
    if (!restFieldsComputed)
        return;

    // update state
    restFieldsSmoothing = true;
    restFieldsSmoothed = false;

    // update UI
    emit smoothingFields(true);
    emit smoothedFields(false);

    // update the scene
    update();

#ifndef Q_NO_CONCURRENT
    // cancel and wait loading threads
    if (smoothFutureWatcher.isRunning()) {
        smoothFutureWatcher.cancel();
        smoothFutureWatcher.waitForFinished();
    }

    // compute in another thread
    smoothFutureWatcher.setFuture(QtConcurrent::run(smoothFieldsParallel));
#else
    // compute rest fields
    computeSmoothFields();

    // update current state
    fieldsSmoothed();
#endif
}

/**
 * @brief OpenGLScene::thresholdPersistentHardContraints
 */
void OpenGLScene::thresholdPersistentHardContraints() {
    // update UI
    emit updatingPersistentHardConstraintsVector(true);

#ifndef Q_NO_CONCURRENT
    // cancel and wait loading threads
    if (persistentConstraintsFutureWatcher.isRunning()) {
        persistentConstraintsFutureWatcher.cancel();
        persistentConstraintsFutureWatcher.waitForFinished();
    }

    // compute in another thread
    persistentConstraintsFutureWatcher.setFuture(QtConcurrent::run(persistentHardConstraintThresholdingParallel));
#else
    // threshold
    persistentHardConstraintThresholding();

    // update current state
    persistentHardConstraintThresholded();
#endif
}

/**
 * @brief OpenGLScene::persistentHardConstraintThresholding
 * @param v
 */
void OpenGLScene::persistentHardConstraintThresholding() {
    // temporary variables
    TMesh *mesh;
    bool wasLoaded;

    // reset persistent hard constraints vector
    std::fill(persistentHardConstraintsVector.begin(), persistentHardConstraintsVector.end(), 0);

    // for each mesh
    for (unsigned long i = 0; i < (unsigned long)animation.getNumberOfMeshes(); i++) {
        // check if i-th mesh has been already loaded
        wasLoaded = animation.isMeshLoadedAt(i);

        // get mesh
        mesh = animation.getTriMeshAt(i);

        // for each face
        for (unsigned long j = 0; j < (unsigned long)mesh->FN(); j++)
            // thresholding anisotropy
            if (mesh->face[j].Q() >= anisotropyThreshold)
                persistentHardConstraintsVector[j]++;

        // unlock mesh
        animation.unlockAt(i);

        // unset mesh
        mesh = 0;

        // unload mesh if necessary
        if (!loadAll && !wasLoaded)
            animation.unloadMeshAt(i);
    }

    // check if rest pose mesh has been already loaded
    wasLoaded = animation.isMeshLoadedAt(restPoseTriMeshIndex);

    // get mesh
    mesh = animation.getTriMeshAt(restPoseTriMeshIndex);

    // for each face
    for (unsigned long j = 0; j < (unsigned long)mesh->FN(); j++)
        // if selected, do not pass the thresholding
        if (mesh->face[j].IsS())
            persistentHardConstraintsVector[j] = 0;

    // unlock mesh
    animation.unlockAt(restPoseTriMeshIndex);

    // unset mesh
    mesh = 0;

    // unload mesh if necessary
    if (!loadAll && !wasLoaded)
        animation.unloadMeshAt(restPoseTriMeshIndex);
}

/**
 * @brief OpenGLScene::persistentHardConstraintThresholded
 */
void OpenGLScene::persistentHardConstraintThresholded() {
    // update status
    viewPersistentCurvature = true;

    // update UI
    emit updatingPersistentHardConstraintsVector(false);

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::smoothFieldsOnRestPoseMesh
 */
void OpenGLScene::computeSmoothFields() {
    // log
    textToLog(tr("Smoothing rest pose mesh fields ... "));

    // start time
    QTime time;
    time.start();

    // check if rest mesh is already loaded
    bool wasRestMeshLoaded = animation.isMeshLoadedAt(restPoseTriMeshIndex);

    // get rest mesh
    TMesh * restMesh = animation.getTriMeshAt(restPoseTriMeshIndex);

    // read hard constraints vertices
    vcg::tri::UpdateFlags<TMesh>::VertexClearS(*restMesh);
    std::string hardstr;
    animation.getMeshNameAt(restPoseTriMeshIndex, hardstr);
    hardstr += ".fl";
    ifstream in(hardstr.c_str());
    if (in.is_open()) {
        int nEdges, v1, v2, tmp;
        in >> nEdges;
        assert(in.good());
        for (int i = 0; i < nEdges; i++) {
            in >> tmp;  // ignore "2"
            assert(in.good());
            in >> v1;   // read v1
            assert(in.good());
            in >> v2;   // read v2
            assert(in.good());
            restMesh->vert[v1].SetS();
        }
        in.close();
    }

    // persistent hard constraints thresholding
    persistentHardConstraintThresholding();

    // create smoother
    SmootherWrapper<TMesh> smoother;

    // initialize it
    smoother.initializeSmoother(*restMesh);

    // set alpha
    smoother.setSoftAlpha(smootherAlpha);

    // get min and max
    float min = fields.getAbsStretchMeasureLowerBound(stretchType);
    float max = fields.getAbsStretchMeasureUpperBound(stretchType);

    // set constraints
    TMesh::CoordType f;
    float w;
//    float tau;
    bool persistent;
    for (int i = 0; i < restMesh->FN(); i++) {
        // if on border do not set its field as constraint
        if (restMesh->face[i].IsB(TFace::BORDER012) ||
                &restMesh->face[i] == restMesh->face[i].FFp(0)
                || &restMesh->face[i] == restMesh->face[i].FFp(1)
                || &restMesh->face[i] == restMesh->face[i].FFp(2)) {
//            std::cout << "Border at " << i << std::endl;
            continue;
        }

        // get field as soft constraint
        fields.getFirstFieldAt(restPoseTriMeshIndex, (unsigned long)i, f);

        // get stretch as weight
        w = fields.getClampedAbsStretchMeasureAt(restPoseTriMeshIndex, (unsigned long)i, stretchType);
        // normalize
        w = (w - min) / (max - min > 0.0 ? max - min : 1.0);

        // persistent?
        persistent = !restMesh->face[i].IsS() && persistentHardConstraintsVector[i] >= persistentHardConstraintThreshold;
        // manual persistent hard constraint:
        for (int j = 0; j < 3; j++) {
            if (restMesh->face[i].V(j)->IsS() && restMesh->face[i].V((j+1)%3)->IsS()) {
                restMesh->face[i].PD1() = (restMesh->face[i].P((j+1)%3) - restMesh->face[i].P(j)).Normalize();
                restMesh->face[i].Q() = 1.0;    // max anisotropy
                persistentHardConstraintsVector[i] = animation.getNumberOfMeshes();
                persistent = true;
                break;
            }
        }

        // get face curvatures to compute anisotropy
//        tau = QAnimationWrapper::computeAnisotropyTau(restMesh->face[i].K1(), restMesh->face[i].K2());
//        tau = restMesh->face[i].Q();

        // set soft constraint
        if (w >= stretchThreshold && smootherAlpha > 0.0 && !persistent)
            smoother.setConstraintSoft(i, w, f);

        // set main curvature dir as hard constraint if anisotropy is greater than threshold
        if (persistent)
            smoother.setConstraintHard(i, restMesh->face[i].PD1());
    }

    // solve
    smoother.solve();

    // get smoothed fields
    smoother.getFFieldPerFace(smoothedFieldsVector);

    // find singularities
    smoother.findCones();

    // get singuarities
    smoother.getSingularitiesIndexPerVertexList(smoothedFieldsSingularitiesList,
                                                std::numeric_limits<float>::epsilon());

    // unlock rest mesh
    animation.unlockAt(restPoseTriMeshIndex);

    // unset rest mesh
    restMesh = 0;

    // unload rest mesh if necessary
    if (!loadAll && !wasRestMeshLoaded)
        animation.unloadMeshAt(restPoseTriMeshIndex);

    // get elapsed time
    QTime elapsedTime(0, 0, 0, 0);
    elapsedTime = elapsedTime.addMSecs(time.elapsed());
    // log
    textToLog(tr("done.\n"));
    textToLog(tr("Rest pose mesh fields smoothing time: ") + elapsedTime.toString("HH:mm:ss.zzz") + "\n");
}

/**
 * @brief OpenGLScene::fieldsSmoothed
 */
void OpenGLScene::fieldsSmoothed() {
    // update state
    restFieldsSmoothing = false;
    restFieldsSmoothed = true;

    // update UI
    emit smoothingFields(false);
    emit smoothedFields(true);

    // update the scene
    update();
}

/**
 * @brief OpenGLScene::storeFFields
 */
void OpenGLScene::storeFFields() {
    // check if computed
    if (!fieldsComputed)
        return;

    // get file name
    QString ffieldName = QFileDialog::getSaveFileName(0, tr("Choose file to store FFIELDS."), folder, tr("FFIELDS file (*.ffield)"));

    // if a file was selected
    if (!ffieldName.isEmpty()) {
        // update UI
        emit ffieldsStoring(true);

        // save file
        saveFField(ffieldName.toStdString().c_str());

        // terminate
        storedFFields();

        // log
        textToLog(QString("File '") + ffieldName + QString("' saved.\n"));
    }
}

/**
 * @brief OpenGLScene::saveFField
 * @param ffieldFilename
 */
void OpenGLScene::saveFField(const char *ffieldFilename) {
    // create and open file
    ofstream output(ffieldFilename);

    // check if file is open
    if (!output.is_open()) {
        // log
        textToLog(QString("Error. Not able to open FFIELD file '") + QString(ffieldFilename) + QString("' for writing.\n"));

        // update UI
        emit projectStoring(false);

        return;
    }

    // store data
    if (viewRestFieldsSmoothed) {
        output << "# Generated with NRosy field generator. (panozzo@inf.ethz.ch)" << endl;
        output << "target_frame" << endl;
        output << smoothedFieldsVector.size() << endl;
        output << "k1	 k2	 k1v_x	 k1v_y	 k1v_z	 k2v_x	 k2v_y	 k2v_z" << endl;

        for(unsigned long i = 0; i < smoothedFieldsVector.size(); ++i) {
          output << "1 1 ";
          for(int j = 0; j < 3; ++j)
              output << smoothedFieldsVector[i].first[j] << " ";
          for(int j = 0; j < 3; ++j) {
              output << smoothedFieldsVector[i].second[j];
              if (j != 3)
                  output << " ";
          }
          output << endl;
        }
    } else {
        // store data
        output << "# Generated with NRosy field generator. (panozzo@inf.ethz.ch)" << endl;
        output << "target_frame" << endl;
        output << fields.getNumberOfFieldsAt(currentRestPoseToTriMeshIndex) << endl;
        output << "k1	 k2	 k1v_x	 k1v_y	 k1v_z	 k2v_x	 k2v_y	 k2v_z" << endl;

        TMesh::CoordType f;

        for(unsigned long i = 0; i < fields.getNumberOfFieldsAt(currentRestPoseToTriMeshIndex); ++i) {
          output << "1 1 ";
          fields.getFirstFieldAt(currentRestPoseToTriMeshIndex, i, f);
          for(int j = 0; j < 3; ++j)
              output << f[j] << " ";
          fields.getSecondFieldAt(currentRestPoseToTriMeshIndex, i, f);
          for(int j = 0; j < 3; ++j) {
              output << f[j];
              if (j != 3)
                  output << " ";
          }
          output << endl;
        }
    }

    // close stream
    output.close();
}

/**
 * @brief OpenGLScene::storedFFields
 */
void OpenGLScene::storedFFields() {
    // update UI
    emit ffieldsStoring(false);
    emit ffieldsStored(true);
}

///**
// * @brief OpenGLScene::updateHistograms
// */
//void OpenGLScene::updateHistograms() {
//    // clear histograms
//    histograms.clear();
//    maxCountsOfHistograms.clear();

//    // check if computed
//    if (!fieldsComputed)
//        return;

//    // get number of meshes
//    unsigned long numOfMeshes = fields.getNumberOfMeshes();

//    // resize vector of histograms
//    histograms.resize(numOfMeshes);
//    maxCountsOfHistograms.resize(numOfMeshes);

//    // for each mesh
//    for (unsigned long i = 0; i < numOfMeshes; i++) {
//        // set range of current histogram
//        histograms[i].SetRange(fields.getMinAbsStretchMeasureAt(i, QMeshSequenceFieldsWrapper::SINGVAL1),
//                               fields.getMaxAbsStretchMeasureAt(i, QMeshSequenceFieldsWrapper::SINGVAL1),
//                               256);

//        // for each face
//        for (unsigned long j = 0; j < fields.getNumberOfFieldsAt(i); j++)
//            // add singular value to histogram
//            histograms[i].Add(fields.getStretchMeasureAt(i, j, QMeshSequenceFieldsWrapper::SINGVAL1));

//        // update max count
//        maxCountsOfHistograms[i] = histograms[i].MaxCount();
//    }
//}

/**
 * @brief OpenGLScene::storeProject
 */
void OpenGLScene::storeProject() {
    // check if computed
    if (!fieldsComputed)
        return;

    // get file name
    projectName = QFileDialog::getSaveFileName(0, tr("Choose file to store current project data."), folder, tr("Project file (*.a2q)"));

    // if a file was selected
    if (!projectName.isEmpty()) {
        // update UI
        emit projectStoring(true);

        // create and open file
        ofstream output(projectName.toStdString().c_str());

        // check if file is open
        if (!output.is_open()) {
            // log
            textToLog(QString("Error. Not able to open project file '") + projectName + QString("' for writing.\n"));

            // update UI
            emit projectStoring(false);

            return;
        }

        // first write number of meshes
        output << "Number of meshes: " << fnames.size() << std::endl;

        // leave one row
        output << endl;

        // write folder
        output << "Folder: " << folder.toStdString() << endl;

        // write text
        output << "Meshes: " << endl;

        // write mesh filenames
        for (unsigned long i = 0; i < fnames.size(); i++)
            output << fnames[i].substr(folder.length()) << endl;

        // leave one raw
        output << endl;

        // write rest pose mesh index
        output << "Rest pose mesh index: " << restPoseTriMeshIndex << endl;

        // leave one row
        output << endl;

        // write rest mesh fields kind
        output << "Computed rest pose mesh fields kind: ";
        switch(restFieldsKind) {
        case MAXIMUM:
            output << "1 (maximum)";
            break;

        case AVERAGE:
            output << "2 (average)";
            break;

        default:
            break;
        }
        output << endl;

        // leave one row
        output << endl;

        // write rest mesh fields stretch type
        output << "Computed rest pose mesh fields stretch type: ";
        switch(stretchType) {
        case QMeshSequenceFieldsWrapper::RATIO:
            output << "1 (|s1 / s2| - 1)";
            break;

        case QMeshSequenceFieldsWrapper::MAXWARP:
            output << "2 (max(|s1|,|1/s2|))";
            break;

        case QMeshSequenceFieldsWrapper::SINGVAL1:
            output << "3 (|s1|)";
            break;

        case QMeshSequenceFieldsWrapper::SINGVAL2:
            output << "4 (|s2|)";
            break;

        default:
            break;
        }
        output << endl;

        // leave three rows
        output << endl << endl << endl;

        // close file
        output.close();

#ifndef Q_NO_CONCURRENT
        // cancel and wait loading threads
        if (storeFutureWatcher.isRunning()) {
            storeFutureWatcher.cancel();
            storeFutureWatcher.waitForFinished();
        }

        // compute in another thread
        storeFutureWatcher.setFuture(QtConcurrent::run(storeProjectParallel));
#else
        // store fields
        fields.storeFields(projectName.toStdString());

        // update current state
        exportedProject();
#endif
    }
}

/**
 * @brief OpenGLScene::exportedProject
 */
void OpenGLScene::exportedProject() {
    // log
    textToLog(QString("Project saved in file '") + projectName + QString("'.\n"));

    // update UI
    emit projectStoring(false);
    emit projectStored(true);
}

/**
 * @brief OpenGLScene::openProject
 */
void OpenGLScene::openProject() {
    // get file name
    projectName = QFileDialog::getOpenFileName(0, tr("Open file of project data."), folder, tr("Project file (*.a2q)"));

    // if a file was selected
    if (!projectName.isEmpty()) {
        // update UI
        emit projectOpening(true);

        // open file
        ifstream input(projectName.toStdString().c_str());

        // check if file is open
        if (!input.is_open()) {
            // log
            textToLog(QString("Error. Not able to open project file '") + projectName + QString("' for reading.\n"));

            // update UI
            emit projectOpening(false);

            return;
        }

        // log
        textToLog(QString("Opening project file '") + projectName + QString("' for reading...\n"));

        // close all
        closeMeshes();

        // temporary variables
        string tmp;
        char buffer[4096];

        // ignore text
        tmp = "Number of meshes: ";
        input.ignore(tmp.length());

        // read number of meshes
        unsigned long num;
        input >> num;
        assert(input.good());

        // reset list of names
        fnames.clear();
        fnames.resize(num);

        // ignore two new line char
        input.ignore(2);

        // ignore text
        tmp = "Folder: ";
        input.ignore(tmp.length());

        // read folder
        input.getline(buffer, 4096);
        QFileInfo fileInfo(buffer);
        folder = QDir::toNativeSeparators(fileInfo.path() + "/");

        // ignore text
        tmp = "Meshes: ";
        input.ignore(tmp.length());

        // ignore new line char
        input.ignore();

        // read lines as mesh filenames
        for (unsigned long i = 0; i < num; i++) {
            // read line
            input.getline(buffer, 4096);

            // construct string
            fnames[i] = folder.toStdString() + buffer;
        }

        // ignore one new line char
        input.ignore();

        // ignore text
        tmp = "Rest pose mesh index: ";
        input.ignore(tmp.length());

        // read rest pose index
        input >> importedRestPoseTriMeshIndex;
        assert(input.good());

        // ignore end line char
        input.ignore();

        // ignore one row
        input.ignore();

        // read rest mesh fields kind
        tmp = "Computed rest pose mesh fields kind: ";
        input.ignore(tmp.length());
        int kind;
        input >> kind;
        assert(input.good());
        switch(kind) {
        case 1:
            restFieldsKind = MAXIMUM;
            break;

        case 2:
            restFieldsKind = AVERAGE;
            break;

        default:
            break;
        }
        // update UI
        emit restMeshFieldsKindChanged(kind);

        // ignore to end of line
        input.getline(buffer, 4096);

        // ignore one row
        input.ignore();

        // read rest mesh fields stretch type
        tmp = "Computed rest pose mesh fields stretch type: ";
        input.ignore(tmp.length());
        int type;
        input >> type;
        assert(input.good());
        switch(type) {
        case 1:
            stretchType = QMeshSequenceFieldsWrapper::RATIO;
            break;

        case 2:
            stretchType = QMeshSequenceFieldsWrapper::MAXWARP;
            break;

        case 3:
            stretchType = QMeshSequenceFieldsWrapper::SINGVAL1;
            break;

        case 4:
            stretchType = QMeshSequenceFieldsWrapper::SINGVAL2;
            break;

        default:
            break;
        }
        // update UI
        emit restMeshFieldsStretchTypeChanged(type);

        // ignore to end of line
        input.getline(buffer, 4096);

        // ignore three new line chars
        input.ignore();
        input.ignore();
        input.ignore();

        // save current position on file
        fstream::pos_type pos = input.tellg();

        // close file
        input.close();

        // open meshes:

        // update UI
        emit meshesOpening(true);

        // ask if all meshes must be kept in memory
        QMessageBox::StandardButton responce = QMessageBox::question(0, tr("Loading method"), tr("Keep all meshes on memory?"),
                                                                     QMessageBox::StandardButtons(QMessageBox::No | QMessageBox::Yes),
                                                                     QMessageBox::Yes);
        if (responce == QMessageBox::Yes)
            loadAll = true;
        else
            loadAll = false;

        // log
        textToLog(QString("Loading fields...\n"));

#ifndef Q_NO_CONCURRENT
        // cancel and wait loading threads
        if (openFutureWatcher.isRunning()) {
            openFutureWatcher.cancel();
            openFutureWatcher.waitForFinished();
        }

        // compute in another thread
        openFutureWatcher.setFuture(QtConcurrent::run(loadProjectParallel, pos));
#else
        // load animation
        if (loadAll)
            animation.load(fnames);
        else
            animation.setNames(fnames);

        // store fields
        fields.loadFields(projectName.toStdString(), pos);

//        fieldsComputed = true;
        // update histograms
//        updateHistograms();

        // update current state
        importedProject();
#endif
    }
}

/**
 * @brief OpenGLScene::importedProject
 */
void OpenGLScene::importedProject() {
    // check if both have same size
    if (animation.getNumberOfMeshes() != fields.getNumberOfMeshes()) {
        // close all
        closeMeshes();

        // log
        textToLog(QString("Project from file '") + projectName + QString("' is corrupted. Error: different number of meshes.\n"));

        return;
    }

    // terminate loading meshes by updating current state
    importedMeshes();

    // set rest pose
    selectRestPose(importedRestPoseTriMeshIndex);

    // set reference pose
    selectRestPoseToMesh(importedRestPoseTriMeshIndex);

    // terminate loading fields by updating current state
    computedTriMeshesFields();

    // terminate loading rest pose mesh fields by updating current state
    computedRestMeshFields();

    // log
    textToLog(QString("Project loaded from file '") + projectName + QString("'.\n"));

    // update UI
    emit projectOpening(false);
    emit projectOpened(true);
    emit updatedSingularValue2LowerBound(fields.getAbsStretchMeasureLowerBound(QMeshSequenceFieldsWrapper::SINGVAL2));
    emit updatedSingularValue1UpperBound(fields.getAbsStretchMeasureUpperBound(QMeshSequenceFieldsWrapper::SINGVAL1));
    emit updatedCurrentStretchMeasureLowerBound(fields.getAbsStretchMeasureLowerBound(stretchType));
    emit updatedCurrentStretchMeasureUpperBound(fields.getAbsStretchMeasureUpperBound(stretchType));

    // update scene
    update();
}

/**
 * @brief OpenGLScene::storeSettings
 */
void OpenGLScene::storeSettings() {
    // get file name
    QString settingsName = QFileDialog::getSaveFileName(0, tr("Store file of settings."), folder, tr("Settings file (*.txt)"));

    // if a file was selected
    if (!settingsName.isEmpty()) {
        // update UI
        emit settingsStoring(true);

        // open file
        ofstream output(settingsName.toStdString().c_str());

        // check if file is open
        if (!output.is_open()) {
            // log
            textToLog(QString("Error. Not able to store settings file '") + settingsName + QString("' for writing.\n"));

            // update UI
            emit settingsStoring(false);

            return;
        }

        // log
        textToLog(QString("Opening settings file '") + settingsName + QString("' for writing..."));

        // write anisotropy threshold
        output << "Anisotropy threshold: " << anisotropyThreshold << endl;

        // write bounds
        output << "Singular values lower bound: " << fields.getAbsStretchMeasureLowerBound(QMeshSequenceFieldsWrapper::SINGVAL2) << endl;
        output << "Singular values upper bound: " << fields.getAbsStretchMeasureUpperBound(QMeshSequenceFieldsWrapper::SINGVAL1) << endl;

        // write average threshold
        output << "Average threshold: " << averageThreshold << endl;

        // write rest mesh fields kind
        output << "Computed rest pose mesh fields kind: ";
        switch(restFieldsKind) {
        case MAXIMUM:
            output << "1 (maximum)";
            break;

        case AVERAGE:
            output << "2 (average)";
            break;

        default:
            break;
        }
        output << endl;

        // leave one row
        output << endl;

        // write rest mesh fields stretch type
        output << "Computed rest pose mesh fields stretch type: ";
        switch(stretchType) {
        case QMeshSequenceFieldsWrapper::RATIO:
            output << "1 (|s1 / s2| - 1)";
            break;

        case QMeshSequenceFieldsWrapper::MAXWARP:
            output << "2 (max(|s1|,|1/s2|))";
            break;

        case QMeshSequenceFieldsWrapper::SINGVAL1:
            output << "3 (|s1|)";
            break;

        case QMeshSequenceFieldsWrapper::SINGVAL2:
            output << "4 (|s2|)";
            break;

        default:
            break;
        }
        output << endl;

        // write soft contraints alpha
        output << "Soft constraints alpha: " << smootherAlpha << endl;

        // write soft constraints minimum threshold
        output << "Soft contraints minimum threshold: " << stretchThreshold << endl;

        // write persistent hard constraints threshold
        output << "Persistent hard constraints threshold: " << persistentHardConstraintThreshold << endl;

        // close the file
        output.close();

        // update UI
        emit settingsStoring(false);
        emit settingsStored(true);

        // log
        textToLog(QString("done.\n"));
    }
}

/**
 * @brief OpenGLScene::openSettings
 */
void OpenGLScene::openSettings() {
    // get file name
    QString settingsName = QFileDialog::getOpenFileName(0, tr("Open file of settings data."), folder, tr("Settings file (*.txt)"));

    // if a file was selected
    if (!settingsName.isEmpty()) {
        // update UI
        emit settingsOpening(true);

        // open file
        ifstream input(settingsName.toStdString().c_str());

        // check if file is open
        if (!input.is_open()) {
            // log
            textToLog(QString("Error. Not able to open settings file '") + settingsName + QString("' for reading.\n"));

            // update UI
            emit settingsOpening(false);

            return;
        }

        // log
        textToLog(QString("Opening settings file '") + settingsName + QString("' for reading..."));

        // temporary variables
        string tmp;
        char buffer[4096];

        // ignore text
        tmp = "Anisotropy threshold: ";
        input.ignore(tmp.length());

        // read anisotropy threshold
        double num;
        input >> num;
        assert(input.good());
        anisotropyThreshold = num;

        // update UI
        emit anisotropyThresholdChanged(anisotropyThreshold);

        // ignore new line char
        input.ignore();

        // ignore text
        tmp = "Singular values lower bound: ";
        input.ignore(tmp.length());

        // read lower bound
        input >> num;
        assert(input.good());
        fields.setAbsSingularValue2LowerBound(num);

        // update UI
        emit updatedSingularValue2LowerBound(num);

        // ignore new line char
        input.ignore();

        // ignore text
        tmp = "Singular values upper bound: ";
        input.ignore(tmp.length());

        // read upper bound
        input >> num;
        assert(input.good());
        fields.setAbsSingularValue1UpperBound(num);

        // update UI
        emit updatedSingularValue1UpperBound(num);

        // ignore new line char
        input.ignore();

        // ignore text
        tmp = "Average threshold: ";
        input.ignore(tmp.length());

        // read lower bound
        input >> num;
        assert(input.good());
        averageThreshold = num;

        // update UI
        emit averageThresholdChanged(averageThreshold);

        // ignore new line char
        input.ignore();

        // read rest mesh fields kind
        tmp = "Computed rest pose mesh fields kind: ";
        input.ignore(tmp.length());
        int kind;
        input >> kind;
        assert(input.good());
        switch(kind) {
        case 1:
            restFieldsKind = MAXIMUM;
            break;

        case 2:
            restFieldsKind = AVERAGE;
            break;

        default:
            break;
        }
        // update UI
        emit restMeshFieldsKindChanged(kind);

        // ignore to end of line
        input.getline(buffer, 4096);

        // ignore one row
        input.ignore();

        // read rest mesh fields stretch type
        tmp = "Computed rest pose mesh fields stretch type: ";
        input.ignore(tmp.length());
        int type;
        input >> type;
        assert(input.good());
        switch(type) {
        case 1:
            stretchType = QMeshSequenceFieldsWrapper::RATIO;
            break;

        case 2:
            stretchType = QMeshSequenceFieldsWrapper::MAXWARP;
            break;

        case 3:
            stretchType = QMeshSequenceFieldsWrapper::SINGVAL1;
            break;

        case 4:
            stretchType = QMeshSequenceFieldsWrapper::SINGVAL2;
            break;

        default:
            break;
        }
        // update UI
        emit restMeshFieldsStretchTypeChanged(type);
        emit updatedCurrentStretchMeasureLowerBound(fields.getAbsStretchMeasureLowerBound(stretchType));
        emit updatedCurrentStretchMeasureUpperBound(fields.getAbsStretchMeasureUpperBound(stretchType));

        // ignore to end of line
        input.getline(buffer, 4096);

        // ignore text
        tmp = "Soft constraints alpha: ";
        input.ignore(tmp.length());

        // read soft constraints alpha
        input >> num;
        assert(input.good());
        smootherAlpha = num;

        // update UI
        emit smootherAlphaChanged(smootherAlpha);

        // ignore new line char
        input.ignore();

        // ignore text
        tmp = "Soft contraints minimum threshold: ";
        input.ignore(tmp.length());

        // read soft constraints minimum threshold
        input >> num;
        assert(input.good());
        stretchThreshold = num;

        // update UI
        emit stretchThresholdChanged(stretchThreshold);

        // ignore new line char
        input.ignore();

        // ignore text
        tmp = "Persistent hard constraints threshold: ";
        input.ignore(tmp.length());

        // read persistent hard constraints threshold
        unsigned long pers;
        input >> pers;
        assert(input.good());
        persistentHardConstraintThreshold = pers;

        // update UI
        emit persistentHardConstraintThresholdChanged(persistentHardConstraintThreshold);

        // update UI
        emit settingsOpening(false);
        emit settingsOpened(true);

        // log
        textToLog(QString("done.\n"));
    }
}



/********************************** mouse event handling stuff ***********************************************************/

void OpenGLScene::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    QGraphicsScene::mousePressEvent(event);
    if (event->isAccepted())
        return;

    track.MouseDown (event->scenePos().x(), height() - event->scenePos().y(), QT2VCG(event->button(), event->modifiers()));

    event->accept();
    update();
}

void OpenGLScene::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
    QGraphicsScene::mouseMoveEvent(event);
    if (event->isAccepted())
        return;

    track.MouseMove(event->scenePos().x(), height() - event->scenePos().y());

    event->accept();
    update();
}

void OpenGLScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    QGraphicsScene::mouseReleaseEvent(event);
    if (event->isAccepted())
        return;

    track.MouseUp(event->scenePos().x(), height() - event->scenePos().y(), QT2VCG(event->button(), event->modifiers()));

    event->accept();
    update();
}

void OpenGLScene::wheelEvent(QGraphicsSceneWheelEvent *event) {
    QGraphicsScene::wheelEvent(event);
    if (event->isAccepted())
        return;

    const int WHEEL_STEP = 120;
    track.MouseWheel(event->delta()/(float)WHEEL_STEP, QTWheel2VCG(event->modifiers()));

    event->accept();
    update();
}

void OpenGLScene::keyPressEvent(QKeyEvent *event) {
    QGraphicsScene::keyPressEvent(event);
    if (event->isAccepted())
        return;

    if (event->key() == Qt::Key_Control)
        track.ButtonDown(QT2VCG(Qt::NoButton, Qt::ControlModifier));
    if (event->key() == Qt::Key_Shift)
        track.ButtonDown(QT2VCG(Qt::NoButton, Qt::ShiftModifier));
    if (event->key() == Qt::Key_Alt)
        track.ButtonDown(QT2VCG(Qt::NoButton, Qt::AltModifier));

    event->accept();
    update();
}

void OpenGLScene::keyReleaseEvent(QKeyEvent *event) {
    QGraphicsScene::keyReleaseEvent(event);
    if (event->isAccepted())
        return;

    if (event->key() == Qt::Key_Control)
        track.ButtonUp(QT2VCG(Qt::NoButton, Qt::ControlModifier));
    if (event->key() == Qt::Key_Shift)
        track.ButtonUp(QT2VCG(Qt::NoButton, Qt::ShiftModifier));
    if (event->key() == Qt::Key_Alt)
        track.ButtonUp(QT2VCG(Qt::NoButton, Qt::AltModifier));

    event->accept();
    update();
}
