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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QToolBox>
#include <QPlainTextEdit>
#include <QLineEdit>
#include <QLabel>
#include <QDoubleValidator>
#include <QDoubleSpinBox>
#include <QRadioButton>
#include <QCheckBox>
#include "openglscene.h"
#include "graphicsview.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void showErrorMessage(const QString &msg);
    void textToLog(const QString &text);

private slots:
    void showToolWindow();
    void on_actionShow_FullScreen_triggered();
    void selectRestFieldKind(int b);
    void selectRestFieldStretchType(int s);
    void setLowerBound();
    void setUpperBound();
    void updateSingulaValue2LowerBound(float l);
    void updateSingulaValue1UpperBound(float u);
    void updateCurrentStretchMeasureLowerBoundLabel(float l);
    void updateCurrentStretchMeasureUpperBoundLabel(float u);
    void setColorThreshold();
    void toSmoothFields(bool toSmooth);

private:
    Ui::MainWindow *ui;
    OpenGLScene *scene;
    GraphicsView *view;
    QToolBox * toolWindow;
    QPlainTextEdit logPlainTextEdit;
    QRadioButton *restFieldsKindMaxRadioButton;
    QRadioButton *restFieldsKindAvgRadioButton;
    QRadioButton *stretchTypeRatioRadioButton;
    QRadioButton *stretchTypeMaxWarpRadioButton;
    QRadioButton *stretchTypeSing1RadioButton;
    QLineEdit *lowerBoundEdit;
    QLineEdit *upperBoundEdit;
    QLabel *currentStretchMeasureLowerBoundLabel;
    QLabel *currentStretchMeasureUpperBoundLabel;
    QDoubleValidator *colorLineEditValidator;
    QDoubleSpinBox *colorThresholdEdit;
    QDoubleSpinBox *averageThresholdEdit;
    QLabel *averageThresholdMinLabel;
    QLabel *averageThresholdMaxLabel;
    QCheckBox *viewSmoothCheckBox;
    void createToolWindow();
};

#endif // MAINWINDOW_H
