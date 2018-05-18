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

#include "graphicsview.h"
#include <QResizeEvent>
#include <sstream>

/**
 * @brief GraphicsView::GraphicsView Constructor simply call QGraphicsView constructor.
 * @param parent Parent widget.
 */
GraphicsView::GraphicsView(QWidget *parent) :
    QGraphicsView(parent)
{
}

/**
 * @brief GraphicsView::GraphicsView Contructor simply call QGraphicsView contructor.
 * @param scene The QGraphicsScene to add.
 * @param parent Parent widget.
 */
GraphicsView::GraphicsView(QGraphicsScene *scene, QWidget *parent) :
    QGraphicsView(scene, parent)
{
}

/**
 * @brief GraphicsView::resizeEvent Changes size of this view and that of its scene.
 * @param event The resize event.
 */
void GraphicsView::resizeEvent(QResizeEvent *event) {
    // if a scene has already been set, update its size
    if (scene())
        scene()->setSceneRect(0, 0, event->size().width(), event->size().height());
    QGraphicsView::resizeEvent(event);

    std::stringstream stream;
    stream << "Size: " << event->size().width() << " X " << event->size().height();
    emit printSizes(QString::fromStdString(stream.str()));
}
