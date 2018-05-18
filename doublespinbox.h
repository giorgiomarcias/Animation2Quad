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

#ifndef DOUBLESPINBOX_H
#define DOUBLESPINBOX_H

#include <QDoubleSpinBox>

class DoubleSpinBox : public QDoubleSpinBox
{
    Q_OBJECT
public:
    /**
     * @brief DoubleSpinBox
     * @param parent
     */
    explicit DoubleSpinBox(QWidget *parent = 0) : QDoubleSpinBox(parent) {

    }
    
signals:
    
public slots:
    /**
     * @brief minimumChanged
     * @param min
     */
    void minimumChanged(double min) {
        setMinimum(min);
    }

    /**
     * @brief maximumChanged
     * @param max
     */
    void maximumChanged(double max) {
        setMaximum(max);
    }
};

#endif // DOUBLESPINBOX_H
