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

#ifndef SPINBOX_H
#define SPINBOX_H

#include <QSpinBox>

class SpinBox : public QSpinBox
{
    Q_OBJECT
public:
    /**
     * @brief SpinBox
     * @param parent
     */
    explicit SpinBox(QWidget *parent = 0) : QSpinBox(parent) {

    }
    
signals:
    
public slots:
    /**
     * @brief maximumChanged
     * @param m
     */
    void maximumChanged(int m)  {
        setMaximum(m);
    }

    /**
     * @brief minimumChanged
     * @param m
     */
    void minimumChanged(int m) {
        setMinimum(m);
    }
    
};

#endif // SPINBOX_H
