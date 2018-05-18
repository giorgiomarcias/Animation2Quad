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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <GL/glew.h>
#include <QGLWidget>
#include <QMessageBox>
#include <QGroupBox>
#include <QFormLayout>
#include <QStyleFactory>
#include <QComboBox>
#include <QButtonGroup>
#include <QPushButton>
#include <QRadioButton>
#include <QSlider>
#include <QDoubleValidator>
#include <QLineEdit>
#include <QLabel>
#include "graphicsview.h"
#include "spinbox.h"
#include "doublespinbox.h"

/**
 * @brief MainWindow::MainWindow
 * @param parent
 */
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    toolWindow = 0;
    ui->setupUi(this);
    scene = new OpenGLScene(this);
    connect(scene, SIGNAL(noOpenGLPainter(const QString&)), this, SLOT(showErrorMessage(const QString&)));
    connect(ui->actionTrackBall, SIGNAL(toggled(bool)), scene, SLOT(showTrackBall(bool)));
    connect(ui->actionControls, SIGNAL(triggered()), this, SLOT(showToolWindow()));
    connect(ui->actionOpen, SIGNAL(triggered()), scene, SLOT(openMeshes()));
    connect(ui->actionClose, SIGNAL(triggered()), scene, SLOT(closeMeshes()));
    connect(ui->actionOpen_Project, SIGNAL(triggered()), scene, SLOT(openProject()));
    connect(ui->actionSave_Project, SIGNAL(triggered()), scene, SLOT(storeProject()));
    connect(ui->actionSave_FFIELD, SIGNAL(triggered()), scene, SLOT(storeFFields()));
    connect(ui->actionOpen_settings, SIGNAL(triggered()), scene, SLOT(openSettings()));
    connect(ui->actionSave_settings, SIGNAL(triggered()), scene, SLOT(storeSettings()));
    connect(ui->actionPrint_image, SIGNAL(triggered()), scene, SLOT(printCurrentScene()));
    connect(scene, SIGNAL(meshesOpening(bool)), ui->actionOpen, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(meshesOpened(bool)), ui->actionClose, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(projectOpening(bool)), ui->actionOpen_Project, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(projectStoring(bool)), ui->actionSave_Project, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(projectStored(bool)), ui->actionSave_Project, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(computedFields(bool)), ui->actionSave_FFIELD, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(ffieldsStoring(bool)), ui->actionSave_FFIELD, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(ffieldsStored(bool)), ui->actionSave_FFIELD, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(meshesOpened(bool)), ui->actionOpen_settings, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(settingsOpening(bool)), ui->actionOpen_settings, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(settingsOpened(bool)), ui->actionOpen_settings, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(smoothedFields(bool)), ui->actionSave_settings, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(settingsStoring(bool)), ui->actionSave_settings, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(settingsStored(bool)), ui->actionSave_settings, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(computedRestFields(bool)), ui->actionSave_Project, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(meshesOpened(bool)), ui->actionPrint_image, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(sendLogMessage(QString)), this, SLOT(textToLog(QString)));

    view = new GraphicsView(scene);
    view->setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));
    view->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
    ui->horizontalLayout->addWidget(view);
    connect(view, SIGNAL(printSizes(QString)), ui->statusBar, SLOT(showMessage(QString)));

    createToolWindow();
    showMaximized();
}

/**
 * @brief MainWindow::showErrorMessage
 * @param msg
 */
void MainWindow::showErrorMessage(const QString &msg) {
    QMessageBox::critical(this, tr("Error"), msg);
    exit(1);
}

/**
 * @brief MainWindow::~MainWindow
 */
MainWindow::~MainWindow()
{
    delete ui;
}

/**
 * @brief MainWindow::on_actionShow_FullScreen_triggered
 */
void MainWindow::on_actionShow_FullScreen_triggered()
{
    if (isFullScreen()) {
        showNormal();
        ui->actionShow_FullScreen->setChecked(false);
    } else {
        showFullScreen();
        ui->actionShow_FullScreen->setChecked(true);
    }
}

/**
 * @brief MainWindow::createToolWindow
 */
void MainWindow::createToolWindow() {
    // generate window
    toolWindow = new QToolBox(this, Qt::Tool);

    QVBoxLayout *boxLayout = new QVBoxLayout;

    QRadioButton *triMeshTypeRadioButton = new QRadioButton(tr("triangles"));
    QRadioButton *polyMeshTypeRadioButton = new QRadioButton(tr("polygons"));
    triMeshTypeRadioButton->setChecked(true);
    QGroupBox *meshTypeGroupBox = new QGroupBox(tr("Mesh type"));
    QButtonGroup *meshTypeButtonGroup = new QButtonGroup(meshTypeGroupBox);
    meshTypeButtonGroup->addButton(triMeshTypeRadioButton, 1);
    meshTypeButtonGroup->addButton(polyMeshTypeRadioButton, 2);
    connect(meshTypeButtonGroup, SIGNAL(buttonPressed(int)), scene, SLOT(setMeshType(int)));
//    connect(scene, SIGNAL(meshesOpened(bool)), triMeshTypeRadioButton, SLOT(setEnabled(bool)));
//    connect(scene, SIGNAL(meshesOpened(bool)), polyMeshTypeRadioButton, SLOT(setEnabled(bool)));
//    triMeshTypeRadioButton->setEnabled(false);
//    polyMeshTypeRadioButton->setEnabled(false);
    QVBoxLayout *meshTypeGroupBoxLayout = new QVBoxLayout;
    meshTypeGroupBoxLayout->addWidget(triMeshTypeRadioButton);
    meshTypeGroupBoxLayout->addWidget(polyMeshTypeRadioButton);
    meshTypeGroupBox->setLayout(meshTypeGroupBoxLayout);
    boxLayout->addWidget(meshTypeGroupBox);
    meshTypeGroupBox->setEnabled(false);
    connect(scene, SIGNAL(meshesOpened(bool)), meshTypeGroupBox, SLOT(setEnabled(bool)));

    SpinBox *meshSelectionSpinBox = new SpinBox;
    meshSelectionSpinBox->setSingleStep(1);
    meshSelectionSpinBox->setMinimum(0);
    meshSelectionSpinBox->setMaximum(0);
    QCheckBox *showMeshesCurvatureCheckBox = new QCheckBox(tr("Show curvature"));
    showMeshesCurvatureCheckBox->setChecked(false);
    connect(showMeshesCurvatureCheckBox, SIGNAL(toggled(bool)), scene, SLOT(showMeshesCurvature(bool)));
    connect(triMeshTypeRadioButton, SIGNAL(toggled(bool)), showMeshesCurvatureCheckBox, SLOT(setEnabled(bool)));
    QLabel *anisotropyThresholdMinLabel = new QLabel(tr("0"));
    QDoubleSpinBox *anisotropyThresholdDoubleSpinBox = new QDoubleSpinBox;
    anisotropyThresholdDoubleSpinBox->setMinimum(0.0);
    anisotropyThresholdDoubleSpinBox->setMaximum(1.0);
    anisotropyThresholdDoubleSpinBox->setDecimals(2);
    anisotropyThresholdDoubleSpinBox->setSingleStep(0.01);
    anisotropyThresholdDoubleSpinBox->setValue(0.99);
    scene->setAnisotropyThreshold(0.99);
    connect(anisotropyThresholdDoubleSpinBox, SIGNAL(valueChanged(double)), scene, SLOT(setAnisotropyThreshold(double)));
    connect(scene, SIGNAL(anisotropyThresholdChanged(double)), anisotropyThresholdDoubleSpinBox, SLOT(setValue(double)));
    QLabel *anisotropyThresholdMaxLabel = new QLabel(tr("1"));
    QHBoxLayout *anisotropyThresholdLayout = new QHBoxLayout;
    anisotropyThresholdLayout->addWidget(anisotropyThresholdMinLabel);
    anisotropyThresholdLayout->addWidget(anisotropyThresholdDoubleSpinBox);
    anisotropyThresholdLayout->addWidget(anisotropyThresholdMaxLabel);
    QGroupBox *anisotropyThresholdGroupBox = new QGroupBox(tr("Anisotropy threshold"));
    anisotropyThresholdGroupBox->setLayout(anisotropyThresholdLayout);
    anisotropyThresholdGroupBox->setEnabled(false);
    connect(showMeshesCurvatureCheckBox, SIGNAL(toggled(bool)), anisotropyThresholdGroupBox, SLOT(setEnabled(bool)));
    connect(triMeshTypeRadioButton, SIGNAL(toggled(bool)), anisotropyThresholdGroupBox, SLOT(setEnabled(bool)));
    QCheckBox *showQuadColorGoldCheckBox = new QCheckBox(tr("Show gold"));
    showQuadColorGoldCheckBox->setChecked(false);
    showQuadColorGoldCheckBox->setEnabled(false);
    connect(showQuadColorGoldCheckBox, SIGNAL(toggled(bool)), scene, SLOT(showQuadColorGold(bool)));
    connect(polyMeshTypeRadioButton, SIGNAL(toggled(bool)), showQuadColorGoldCheckBox, SLOT(setEnabled(bool)));
    QCheckBox *showQuadSingularitiesCheckBox = new QCheckBox(tr("Show quad singularities"));
    showQuadSingularitiesCheckBox->setChecked(false);
    showQuadSingularitiesCheckBox->setEnabled(false);
    connect(showQuadSingularitiesCheckBox, SIGNAL(toggled(bool)), scene, SLOT(showQuadSingularities(bool)));
    connect(polyMeshTypeRadioButton, SIGNAL(toggled(bool)), showQuadSingularitiesCheckBox, SLOT(setEnabled(bool)));
    QVBoxLayout *meshSelectionGroupBoxLayout = new QVBoxLayout;
    meshSelectionGroupBoxLayout->addWidget(meshSelectionSpinBox);
    meshSelectionGroupBoxLayout->addWidget(showMeshesCurvatureCheckBox);
    meshSelectionGroupBoxLayout->addWidget(anisotropyThresholdGroupBox);
    meshSelectionGroupBoxLayout->addWidget(showQuadColorGoldCheckBox);
    meshSelectionGroupBoxLayout->addWidget(showQuadSingularitiesCheckBox);
    QGroupBox *meshSelectionGroupBox = new QGroupBox(tr("Mesh selection"));
    meshSelectionGroupBox->setLayout(meshSelectionGroupBoxLayout);
    boxLayout->addWidget(meshSelectionGroupBox);
//    meshSelectionSpinBox->setEnabled(false);
//    connect(scene, SIGNAL(meshesOpened(bool)), meshSelectionSpinBox, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(meshSelect(int)), meshSelectionSpinBox, SLOT(setValue(int)));
    connect(meshSelectionSpinBox, SIGNAL(valueChanged(int)), scene, SLOT(selectMesh(int)));
    connect(scene, SIGNAL(meshMaxIndex(int)), meshSelectionSpinBox, SLOT(maximumChanged(int)));
    meshSelectionGroupBox->setEnabled(false);
    connect(scene, SIGNAL(meshesOpened(bool)), meshSelectionGroupBox, SLOT(setEnabled(bool)));

    QGroupBox *meshFeaturesGroupBox = new QGroupBox;
    meshFeaturesGroupBox->setLayout(boxLayout);
    toolWindow->addItem(meshFeaturesGroupBox, tr("Mesh features"));

    SpinBox *restPoseSelectionSpinBox = new SpinBox;
    restPoseSelectionSpinBox->setSingleStep(1);
    restPoseSelectionSpinBox->setMinimum(0);
    restPoseSelectionSpinBox->setMaximum(0);
//    restPoseSelectionSpinBox->setEnabled(false);
//    connect(scene, SIGNAL(meshesOpened(bool)), restPoseSelectionSpinBox, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(restPoseSelect(int)), restPoseSelectionSpinBox, SLOT(setValue(int)));
    connect(restPoseSelectionSpinBox, SIGNAL(valueChanged(int)), scene, SLOT(selectRestPose(int)));
    connect(scene, SIGNAL(meshMaxIndex(int)), restPoseSelectionSpinBox, SLOT(maximumChanged(int)));
//    connect(scene, SIGNAL(computingFields(bool)), restPoseSelectionSpinBox, SLOT(setDisabled(bool)));
//    connect(scene, SIGNAL(computedFields(bool)), restPoseSelectionSpinBox, SLOT(setDisabled(bool)));
    QVBoxLayout *restPoseSelectionGroupBoxLayout = new QVBoxLayout;
    restPoseSelectionGroupBoxLayout->addWidget(restPoseSelectionSpinBox);
    QGroupBox *restPoseSelectionGroupBox = new QGroupBox(tr("Select Rest Pose:"));
    restPoseSelectionGroupBox->setLayout(restPoseSelectionGroupBoxLayout);

    QPushButton *computePushButton = new QPushButton(tr("Compute fields"));
//    computePushButton->setEnabled(false);
//    connect(scene, SIGNAL(meshesOpened(bool)), computePushButton, SLOT(setEnabled(bool)));
//    connect(scene, SIGNAL(computingFields(bool)), computePushButton, SLOT(setDisabled(bool)));
//    connect(scene, SIGNAL(computedFields(bool)), computePushButton, SLOT(setDisabled(bool)));
    connect(computePushButton, SIGNAL(clicked()), scene, SLOT(computeTriMeshesFields()));
    QVBoxLayout *computeFieldsLayout = new QVBoxLayout;
    computeFieldsLayout->addWidget(restPoseSelectionGroupBox);
    computeFieldsLayout->addWidget(computePushButton);
    QGroupBox *computeFieldsGroupBox = new QGroupBox(tr("Sequence Fields"));
    computeFieldsGroupBox->setLayout(computeFieldsLayout);
    computeFieldsGroupBox->setEnabled(false);
    connect(scene, SIGNAL(computingFields(bool)), computeFieldsGroupBox, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(computedFields(bool)), computeFieldsGroupBox, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(meshesOpened(bool)), computeFieldsGroupBox, SLOT(setEnabled(bool)));

    SpinBox *restPoseToMeshSelectionSpinBox = new SpinBox;
    restPoseToMeshSelectionSpinBox->setSingleStep(1);
    restPoseToMeshSelectionSpinBox->setMinimum(0);
    restPoseToMeshSelectionSpinBox->setMaximum(0);
//    restPoseToMeshSelectionSpinBox->setEnabled(false);
//    connect(scene, SIGNAL(computedFields(bool)), restPoseToMeshSelectionSpinBox, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(restPoseToMeshSelect(int)), restPoseToMeshSelectionSpinBox, SLOT(setValue(int)));
    connect(restPoseToMeshSelectionSpinBox, SIGNAL(valueChanged(int)), scene, SLOT(selectRestPoseToMesh(int)));
    connect(scene, SIGNAL(meshMaxIndex(int)), restPoseToMeshSelectionSpinBox, SLOT(maximumChanged(int)));
//    QCheckBox *showFaceColorCheckBox = new QCheckBox(tr("Show face color"));
//    connect(showFaceColorCheckBox, SIGNAL(toggled(bool)), scene, SLOT(showFaceColors(bool)));
    QVBoxLayout *restPoseToMeshSelectionGroupBoxLayout = new QVBoxLayout;
    restPoseToMeshSelectionGroupBoxLayout->addWidget(restPoseToMeshSelectionSpinBox);
//    restPoseToMeshSelectionGroupBoxLayout->addWidget(showFaceColorCheckBox);
    QGroupBox *restPoseToMeshSelectionGroupBox = new QGroupBox(tr("Select Reference Pose:"));
    restPoseToMeshSelectionGroupBox->setLayout(restPoseToMeshSelectionGroupBoxLayout);

    QGroupBox *fieldsAppearenceGroupBox = new QGroupBox(tr("Fields Appearence"));
    QRadioButton *fieldsAppearenceFieldsRadioButton = new QRadioButton(tr("Show fields"));
    fieldsAppearenceFieldsRadioButton->setChecked(true);
    QRadioButton *fieldsAppearenceFaceRadioButton = new QRadioButton(tr("Show face colors"));
    QRadioButton *fieldsAppearenceCurvatureRadioButton = new QRadioButton(tr("Show curvature"));
    QButtonGroup *fieldsAppearenceButtonGroup = new QButtonGroup(fieldsAppearenceGroupBox);
    fieldsAppearenceButtonGroup->addButton(fieldsAppearenceFieldsRadioButton, 1);
    fieldsAppearenceButtonGroup->addButton(fieldsAppearenceFaceRadioButton, 2);
    fieldsAppearenceButtonGroup->addButton(fieldsAppearenceCurvatureRadioButton, 3);
    connect(fieldsAppearenceButtonGroup, SIGNAL(buttonClicked(int)), scene, SLOT(selectFieldsDrawMode(int)));
    QVBoxLayout *fieldsAppearenceLayout = new QVBoxLayout;
    fieldsAppearenceLayout->addWidget(fieldsAppearenceFieldsRadioButton);
    fieldsAppearenceLayout->addWidget(fieldsAppearenceFaceRadioButton);
    fieldsAppearenceLayout->addWidget(fieldsAppearenceCurvatureRadioButton);
    fieldsAppearenceGroupBox->setLayout(fieldsAppearenceLayout);

    QDoubleValidator *validator = new QDoubleValidator(0.0, numeric_limits<double>::max(), 20);
    validator->setNotation(QDoubleValidator::StandardNotation);
    colorLineEditValidator = new QDoubleValidator(0.0, numeric_limits<double>::max(), 20);
    validator->setNotation(QDoubleValidator::StandardNotation);
    lowerBoundEdit = new QLineEdit;
    lowerBoundEdit->setValidator(validator);
    lowerBoundEdit->setText("0.5");
    lowerBoundEdit->setFixedWidth(50);
    connect(lowerBoundEdit, SIGNAL(editingFinished()), this, SLOT(setLowerBound()));
    QLabel *lowerBoundEditLabel = new QLabel("Lower:");
    upperBoundEdit = new QLineEdit;
    upperBoundEdit->setValidator(validator);
    upperBoundEdit->setText("2.0");
    upperBoundEdit->setFixedWidth(50);
    connect(upperBoundEdit, SIGNAL(editingFinished()), this, SLOT(setUpperBound()));
    QLabel *upperBoundEditLabel = new QLabel("Upper:");
    colorThresholdEdit = new QDoubleSpinBox;
    colorThresholdEdit->setSingleStep(0.01);
    connect(colorThresholdEdit, SIGNAL(valueChanged(double)), scene, SLOT(selectColorThreshold(double)));
    QLabel *colorThresholdEditLabel = new QLabel("Color th:");
    currentStretchMeasureLowerBoundLabel = new QLabel;
    QLabel *currentStretchMeasureLowerBoundLabelHeader = new QLabel("Current Lower:");
    currentStretchMeasureUpperBoundLabel = new QLabel;
    QLabel *currentStretchMeasureUpperBoundLabelHeader = new QLabel("Current Upper:");
    QFormLayout *thresholdsFormLayout = new QFormLayout;
    thresholdsFormLayout->addRow(lowerBoundEditLabel, lowerBoundEdit);
    thresholdsFormLayout->addRow(upperBoundEditLabel, upperBoundEdit);
    thresholdsFormLayout->addRow(currentStretchMeasureLowerBoundLabelHeader, currentStretchMeasureLowerBoundLabel);
    thresholdsFormLayout->addRow(currentStretchMeasureUpperBoundLabelHeader, currentStretchMeasureUpperBoundLabel);
    thresholdsFormLayout->addRow(colorThresholdEditLabel, colorThresholdEdit);
    QGroupBox *thresholdsSelectionGroupBox = new QGroupBox(tr("Bounds and Thresholds"));
    thresholdsSelectionGroupBox->setLayout(thresholdsFormLayout);
    QVBoxLayout *referencePoseLayout = new QVBoxLayout;
    referencePoseLayout->addWidget(restPoseToMeshSelectionGroupBox);
    referencePoseLayout->addWidget(fieldsAppearenceGroupBox);
    referencePoseLayout->addWidget(thresholdsSelectionGroupBox);
    QGroupBox *referencePoseGroupBox = new QGroupBox(tr("Fields Appearence"));
    referencePoseGroupBox->setLayout(referencePoseLayout);
    referencePoseGroupBox->setEnabled(false);
    connect(scene, SIGNAL(computedFields(bool)), referencePoseGroupBox, SLOT(setEnabled(bool)));

    restFieldsKindMaxRadioButton = new QRadioButton(tr("Maximum fields"));
//    connect(scene, SIGNAL(computedFields(bool)), restFieldsKindMaxRadioButton, SLOT(setEnabled(bool)));
//    restFieldsKindMaxRadioButton->setEnabled(false);
//    restFieldsKindMaxRadioButton->setChecked(true);
    restFieldsKindAvgRadioButton = new QRadioButton(tr("Average fields"));
//    connect(scene, SIGNAL(computedFields(bool)), restFieldsKindAvgRadioButton, SLOT(setEnabled(bool)));
//    restFieldsKindAvgRadioButton->setEnabled(false);
    restFieldsKindAvgRadioButton->setChecked(true);
    QGroupBox *restFieldsKindGroupBox = new QGroupBox(tr("Select rest fields kind:"));
    QButtonGroup *restFieldsKindButtonGroup = new QButtonGroup(restFieldsKindGroupBox);
    restFieldsKindButtonGroup->addButton(restFieldsKindMaxRadioButton, 1);
    restFieldsKindButtonGroup->addButton(restFieldsKindAvgRadioButton, 2);
    connect(restFieldsKindButtonGroup, SIGNAL(buttonClicked(int)), scene, SLOT(setRestMeshFieldsKind(int)));
    connect(scene, SIGNAL(restMeshFieldsKindChanged(int)), this, SLOT(selectRestFieldKind(int)));
    scene->setRestMeshFieldsKind(2);
    averageThresholdMinLabel = new QLabel(tr("0"));
    averageThresholdEdit = new QDoubleSpinBox;
    averageThresholdEdit->setSingleStep(0.01);
    averageThresholdEdit->setValue(0.0);
    scene->setAverageThreshold(0.0);
    connect(averageThresholdEdit, SIGNAL(valueChanged(double)), scene, SLOT(setAverageThreshold(double)));
    connect(scene, SIGNAL(averageThresholdChanged(double)), averageThresholdEdit, SLOT(setValue(double)));
    averageThresholdMaxLabel = new QLabel(tr("1"));
    QHBoxLayout *averageThresholdLayout = new QHBoxLayout;
    averageThresholdLayout->addWidget(averageThresholdMinLabel);
    averageThresholdLayout->addWidget(averageThresholdEdit);
    averageThresholdLayout->addWidget(averageThresholdMaxLabel);
    QGroupBox *averageThresholdGroupBox = new QGroupBox(tr("Average threshold"));
    averageThresholdGroupBox->setLayout(averageThresholdLayout);
    connect(restFieldsKindAvgRadioButton, SIGNAL(toggled(bool)), averageThresholdGroupBox, SLOT(setEnabled(bool)));
    QVBoxLayout *restFieldsKindLayout = new QVBoxLayout;
    restFieldsKindLayout->addWidget(restFieldsKindMaxRadioButton);
    restFieldsKindLayout->addWidget(restFieldsKindAvgRadioButton);
    restFieldsKindLayout->addWidget(averageThresholdGroupBox);
    restFieldsKindGroupBox->setLayout(restFieldsKindLayout);

    connect(scene, SIGNAL(updatedSingularValue2LowerBound(float)), this, SLOT(updateSingulaValue2LowerBound(float)));
    connect(scene, SIGNAL(updatedSingularValue1UpperBound(float)), this, SLOT(updateSingulaValue1UpperBound(float)));
    connect(scene, SIGNAL(updatedCurrentStretchMeasureLowerBound(float)), this, SLOT(updateCurrentStretchMeasureLowerBoundLabel(float)));
    connect(scene, SIGNAL(updatedCurrentStretchMeasureUpperBound(float)), this, SLOT(updateCurrentStretchMeasureUpperBoundLabel(float)));
    scene->setLowerBound(0.5);
    scene->setUpperBound(2.0);
    scene->selectColorThreshold(1.0);
    colorThresholdEdit->setValue(1.0);

    stretchTypeRatioRadioButton = new QRadioButton(tr("|s1 / s2| - 1"));
//    connect(scene, SIGNAL(computedFields(bool)), stretchTypeRatioRadioButton, SLOT(setEnabled(bool)));
//    stretchTypeRatioRadioButton->setEnabled(false);
    stretchTypeRatioRadioButton->setChecked(true);
    stretchTypeMaxWarpRadioButton = new QRadioButton(tr("max(|s1|, |1/s2|)"));
//    connect(scene, SIGNAL(computedFields(bool)), stretchTypeMaxWarpRadioButton, SLOT(setEnabled(bool)));
//    stretchTypeMaxWarpRadioButton->setEnabled(false);
    stretchTypeSing1RadioButton = new QRadioButton(tr("|s1|"));
//    connect(scene, SIGNAL(computedFields(bool)), stretchTypeSing1RadioButton, SLOT(setEnabled(bool)));
//    stretchTypeSing1RadioButton->setEnabled(false);
    QGroupBox *stretchTypeGroupBox = new QGroupBox(tr("Select stretch type:"));
    QButtonGroup *stretchTypeButtonGroup = new QButtonGroup(stretchTypeGroupBox);
    stretchTypeButtonGroup->addButton(stretchTypeRatioRadioButton, 1);
    stretchTypeButtonGroup->addButton(stretchTypeMaxWarpRadioButton, 2);
    stretchTypeButtonGroup->addButton(stretchTypeSing1RadioButton, 3);
    connect(stretchTypeButtonGroup, SIGNAL(buttonClicked(int)), scene, SLOT(setRestMeshFieldsStretchType(int)));
    connect(scene, SIGNAL(restMeshFieldsStretchTypeChanged(int)), this, SLOT(selectRestFieldStretchType(int)));
    QVBoxLayout *stretchTypeLayout = new QVBoxLayout;
    stretchTypeLayout->addWidget(stretchTypeRatioRadioButton);
    stretchTypeLayout->addWidget(stretchTypeMaxWarpRadioButton);
    stretchTypeLayout->addWidget(stretchTypeSing1RadioButton);
    stretchTypeGroupBox->setLayout(stretchTypeLayout);

    QPushButton *restFieldsPushButton = new QPushButton(tr("Compute Rest Fields"));
//    connect(scene, SIGNAL(computedFields(bool)), restFieldsPushButton, SLOT(setEnabled(bool)));
//    connect(scene, SIGNAL(computingRestFields(bool)), restFieldsPushButton, SLOT(setDisabled(bool)));
//    connect(scene, SIGNAL(computedRestFields(bool)), restFieldsPushButton, SLOT(setEnabled(bool)));
//    restFieldsPushButton->setEnabled(false);
    connect(restFieldsPushButton, SIGNAL(clicked()), scene, SLOT(computeRestMeshFields()));
    connect(scene, SIGNAL(computingRestFields(bool)), restFieldsPushButton, SLOT(setDisabled(bool)));
//    connect(scene, SIGNAL(computedRestFields(bool)), restFieldsPushButton, SLOT(setEnabled(bool)));

    QVBoxLayout *restFieldsLayout = new QVBoxLayout;
    restFieldsLayout->addWidget(restFieldsKindGroupBox);
    restFieldsLayout->addWidget(stretchTypeGroupBox);
    restFieldsLayout->addWidget(restFieldsPushButton);
    QGroupBox *restFieldsGroupBox = new QGroupBox(tr("Rest Mesh Fields"));
    restFieldsGroupBox->setLayout(restFieldsLayout);
    restFieldsGroupBox->setEnabled(false);
    connect(scene, SIGNAL(computedFields(bool)), restFieldsGroupBox, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(computingRestFields(bool)), restFieldsGroupBox, SLOT(setDisabled(bool)));

    viewSmoothCheckBox = new QCheckBox(tr("Show smooth"));
    connect(viewSmoothCheckBox, SIGNAL(toggled(bool)), scene, SLOT(showRestFieldsSmoothed(bool)));
    connect(scene, SIGNAL(computedRestFields(bool)), this, SLOT(toSmoothFields(bool)));

    QDoubleSpinBox *smootherAlphaDoubleSpinBox = new QDoubleSpinBox;
    smootherAlphaDoubleSpinBox->setMinimum(0.0);
    smootherAlphaDoubleSpinBox->setMaximum(1.0);
    smootherAlphaDoubleSpinBox->setDecimals(2);
    smootherAlphaDoubleSpinBox->setSingleStep(0.01);
    smootherAlphaDoubleSpinBox->setValue(0.3);
    scene->setSmootherAlpha(0.3);
    connect(smootherAlphaDoubleSpinBox, SIGNAL(valueChanged(double)), scene, SLOT(setSmootherAlpha(double)));
    connect(scene, SIGNAL(smootherAlphaChanged(double)), smootherAlphaDoubleSpinBox, SLOT(setValue(double)));
    QDoubleSpinBox *stretchThresholdDoubleSpinBox = new QDoubleSpinBox;
    stretchThresholdDoubleSpinBox->setMinimum(0.0);
    stretchThresholdDoubleSpinBox->setMaximum(1.0);
    stretchThresholdDoubleSpinBox->setDecimals(2);
    stretchThresholdDoubleSpinBox->setSingleStep(0.01);
    stretchThresholdDoubleSpinBox->setValue(0.01);
    scene->setStretchThreshold(0.01);
    connect(stretchThresholdDoubleSpinBox, SIGNAL(valueChanged(double)), scene, SLOT(setStretchThreshold(double)));
    connect(scene, SIGNAL(stretchThresholdChanged(double)), stretchThresholdDoubleSpinBox, SLOT(setValue(double)));
    SpinBox *persistentHardConstraintsThresholdSpinBox = new SpinBox;
    persistentHardConstraintsThresholdSpinBox->setMinimum(1);
    persistentHardConstraintsThresholdSpinBox->setSingleStep(1);
    persistentHardConstraintsThresholdSpinBox->setMaximum(1);
    connect(persistentHardConstraintsThresholdSpinBox, SIGNAL(valueChanged(int)), scene, SLOT(setPersistentHardConstraintThreshold(int)));
    connect(scene, SIGNAL(persistentHardConstraintThresholdChanged(int)), persistentHardConstraintsThresholdSpinBox, SLOT(setValue(int)));
    connect(scene, SIGNAL(persistentHardConstraintThresholdMaxIndex(int)), persistentHardConstraintsThresholdSpinBox, SLOT(maximumChanged(int)));
    connect(scene, SIGNAL(persistentHardConstraintThresholdMaxIndex(int)), persistentHardConstraintsThresholdSpinBox, SLOT(setValue(int)));
    QPushButton *updatePersistentPushButton = new QPushButton(tr("Update"));
    connect(updatePersistentPushButton, SIGNAL(clicked()), scene, SLOT(thresholdPersistentHardContraints()));
    connect(scene, SIGNAL(updatingPersistentHardConstraintsVector(bool)), updatePersistentPushButton, SLOT(setDisabled(bool)));
    QCheckBox *showFieldsAsOnlyPassingCheckBox = new QCheckBox(tr("Show"));
    connect(showFieldsAsOnlyPassingCheckBox, SIGNAL(toggled(bool)), scene, SLOT(showFieldsAsOnlyPassingSoftConstraint(bool)));
    connect(scene, SIGNAL(updatingPersistentHardConstraintsVector(bool)), showFieldsAsOnlyPassingCheckBox, SLOT(setDisabled(bool)));
    QFormLayout *smootherParametersFormLayout = new QFormLayout;
    smootherParametersFormLayout->addRow(tr("Soft Alpha:"), smootherAlphaDoubleSpinBox);
    smootherParametersFormLayout->addRow(tr("Soft th.:"), stretchThresholdDoubleSpinBox);
    smootherParametersFormLayout->addRow(tr("Hard th.:"), persistentHardConstraintsThresholdSpinBox);
    smootherParametersFormLayout->addRow(tr("Update persistents:"), updatePersistentPushButton);
    smootherParametersFormLayout->addRow(tr("Only soft cons.:"), showFieldsAsOnlyPassingCheckBox);
    QGroupBox *smootherParametersGroupBox = new QGroupBox(tr("Smoother parameters"));
    smootherParametersGroupBox->setLayout(smootherParametersFormLayout);
    smootherParametersGroupBox->setEnabled(false);
    connect(viewSmoothCheckBox, SIGNAL(toggled(bool)), smootherParametersGroupBox, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(smoothingFields(bool)), smootherParametersGroupBox, SLOT(setDisabled(bool)));

    QPushButton *smoothPushButton = new QPushButton(tr("Smooth"));
    smoothPushButton->setEnabled(false);
    connect(smoothPushButton, SIGNAL(clicked()), scene, SLOT(smoothFields()));
    connect(viewSmoothCheckBox, SIGNAL(toggled(bool)), smoothPushButton, SLOT(setEnabled(bool)));
    connect(scene, SIGNAL(smoothingFields(bool)), smoothPushButton, SLOT(setDisabled(bool)));
    connect(scene, SIGNAL(updatingPersistentHardConstraintsVector(bool)), smoothPushButton, SLOT(setDisabled(bool)));


    QVBoxLayout *viewSmoothLayout = new QVBoxLayout;
    viewSmoothLayout->addWidget(viewSmoothCheckBox);
    viewSmoothLayout->addWidget(smootherParametersGroupBox);
    viewSmoothLayout->addWidget(smoothPushButton);
    QGroupBox *viewSmoothGroupBox = new QGroupBox(tr("Smooth Fields"));
    viewSmoothGroupBox->setLayout(viewSmoothLayout);
    viewSmoothGroupBox->setEnabled(false);
    connect(scene, SIGNAL(computingRestFields(bool)), viewSmoothGroupBox, SLOT(setDisabled(bool)));

    QVBoxLayout *computeGroupBoxLayout = new QVBoxLayout;
    computeGroupBoxLayout->addWidget(computeFieldsGroupBox);
    computeGroupBoxLayout->addWidget(referencePoseGroupBox);
    computeGroupBoxLayout->addWidget(restFieldsGroupBox);
    computeGroupBoxLayout->addWidget(viewSmoothGroupBox);
    QGroupBox *computeGroupBox = new QGroupBox;
    computeGroupBox->setLayout(computeGroupBoxLayout);
    toolWindow->addItem(computeGroupBox, tr("Fields Extraction"));

    logPlainTextEdit.setReadOnly(true);
    logPlainTextEdit.ensureCursorVisible();
    toolWindow->addItem(&logPlainTextEdit, tr("Log"));

    connect(toolWindow, SIGNAL(currentChanged(int)), scene, SLOT(selectDrawMode(int)));

    toolWindow->setWindowTitle(tr("Manager"));
    toolWindow->setWindowOpacity(0.95);
    toolWindow->setContentsMargins(4, 4, 4, 4);
//    connect(controlWidget, SIGNAL(visibleChanged()), this, SLOT(emitControlWidgetClosed()));

    showToolWindow();
}

/**
 * @brief MainWindow::showToolWindow
 */
void MainWindow::showToolWindow() {
    QPoint pos = ui->centralWidget->mapToGlobal(view->pos());
//    pos += QPoint(0, ui->mainToolBar->height());
    toolWindow->move(pos + QPoint(30, 30));
    toolWindow->setVisible(true);
}

/**
 * @brief MainWindow::textToLog
 * @param text
 */
void MainWindow::textToLog(const QString &text) {
    QTextCursor cursor = logPlainTextEdit.textCursor();
    int prevPos = cursor.position();
    logPlainTextEdit.appendPlainText(text);
    cursor.setPosition(prevPos + text.length());
    while(!cursor.atEnd())
        cursor.setPosition(cursor.position() + 1);
    logPlainTextEdit.setTextCursor(cursor);
}

/**
 * @brief MainWindow::selectRestFieldKind
 * @param b
 */
void MainWindow::selectRestFieldKind(int b) {
    switch(b) {
    case 1:
        restFieldsKindMaxRadioButton->setChecked(true);
        break;

    case 2:
        restFieldsKindAvgRadioButton->setChecked(true);
        break;

    default:
        break;
    }
}

/**
 * @brief MainWindow::selectRestFieldStretchType
 * @param s
 */
void MainWindow::selectRestFieldStretchType(int s) {
    switch(s) {
    case 1:
        stretchTypeRatioRadioButton->setChecked(true);
        break;

    case 2:
        stretchTypeMaxWarpRadioButton->setChecked(true);
        break;

    case 3:
        stretchTypeSing1RadioButton->setChecked(true);
        break;

    default:
        break;
    }
}

/**
 * @brief MainWindow::setLowerBound
 */
void MainWindow::setLowerBound() {
    bool ok;
    float l = lowerBoundEdit->text().toFloat(&ok);
    assert(ok);
    ostringstream stream;
    stream << l;
    QString str = QString::fromStdString(stream.str());
    if (l > upperBoundEdit->text().toFloat(&ok)) {
        upperBoundEdit->setText(str);
        scene->setUpperBound(l);
    }
    scene->setLowerBound(l);
}

/**
 * @brief MainWindow::setUpperBound
 */
void MainWindow::setUpperBound() {
    bool ok;
    float u = upperBoundEdit->text().toFloat(&ok);
    assert(ok);
    ostringstream stream;
    stream << u;
    QString str = QString::fromStdString(stream.str());
    if (u < lowerBoundEdit->text().toFloat(&ok)) {
        lowerBoundEdit->setText(str);
        scene->setLowerBound(u);
    }
    scene->setUpperBound(u);
}

/**
 * @brief MainWindow::setColorThreshold
 */
void MainWindow::setColorThreshold() {
    bool ok;
    float t = colorThresholdEdit->text().toFloat(&ok);
    assert(ok);
    scene->selectColorThreshold(t);
}

/**
 * @brief MainWindow::updateSingulaValue2LowerBound
 * @param l
 */
void MainWindow::updateSingulaValue2LowerBound(float l) {
    ostringstream stream;
    stream << l;
    QString str = QString::fromStdString(stream.str());
    lowerBoundEdit->setText(str);
}

/**
 * @brief MainWindow::updateSingulaValue1UpperBound
 * @param u
 */
void MainWindow::updateSingulaValue1UpperBound(float u) {
    ostringstream stream;
    stream << u;
    QString str = QString::fromStdString(stream.str());
    upperBoundEdit->setText(str);
}

/**
 * @brief MainWindow::updateCurrentStretchMeasureLowerBoundLabel
 * @param l
 */
void MainWindow::updateCurrentStretchMeasureLowerBoundLabel(float l) {
    ostringstream stream;
    stream << l;
    QString str = QString::fromStdString(stream.str());
    currentStretchMeasureLowerBoundLabel->setText(str);
    averageThresholdMinLabel->setText(str);

    colorThresholdEdit->setMinimum(l);
    averageThresholdEdit->setMinimum(l);

    float t = colorThresholdEdit->value();
    if (t < l) {
        colorThresholdEdit->setValue(l);
        scene->selectColorThreshold(l);
    }

    t = averageThresholdEdit->value();
    if (t < l) {
        averageThresholdEdit->setValue(l);
        scene->setAverageThreshold(l);
    }
}

/**
 * @brief MainWindow::updateCurrentStretchMeasureUpperBoundLabel
 * @param u
 */
void MainWindow::updateCurrentStretchMeasureUpperBoundLabel(float u) {
    ostringstream stream;
    stream << u;
    QString str = QString::fromStdString(stream.str());
    currentStretchMeasureUpperBoundLabel->setText(str);
    averageThresholdMaxLabel->setText(str);

    colorThresholdEdit->setMaximum(u);
    averageThresholdEdit->setMaximum(u);

    float t = colorThresholdEdit->value();
    if (t > u) {
        colorThresholdEdit->setValue(u);
        scene->selectColorThreshold(u);
    }

    t = averageThresholdEdit->value();
    if (t > u) {
        averageThresholdEdit->setValue(u);
        scene->setAverageThreshold(u);
    }
}

/**
 * @brief MainWindow::toSmoothFields
 * @param toSmooth
 */
void MainWindow::toSmoothFields(bool toSmooth) {
    viewSmoothCheckBox->setChecked(!toSmooth);
}
