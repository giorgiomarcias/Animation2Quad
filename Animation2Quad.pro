#-------------------------------------------------
#
# Project created by QtCreator 2013-02-09T13:33:58
#
#-------------------------------------------------

greaterThan(QT_MAJOR_VERSION, 4): cache()

QT       += core gui opengl concurrent

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets


TARGET = Animation2Quad
TEMPLATE = app

SOURCES +=  main.cpp\
            mainwindow.cpp \
            graphicsview.cpp \
            openglscene.cpp

HEADERS  += mainwindow.h \
            graphicsview.h \
            openglscene.h \
            meshes.h \
            spinbox.h \
            compute.h \
            animation.h \
            meshsequencefields.h \
            qanimationwrapper.h \
            qmeshsequencefieldswrapper.h \
            SmootherWrapper.h \
            doublespinbox.h

FORMS    += mainwindow.ui



DEFINES += NOT_NORMALIZE

# define paths
VCGLIBDIR = $$PWD/external/vcglib
GLEWDIR   = $$PWD/external/glew

# include paths
INCLUDEPATH += $$VCGLIBDIR
INCLUDEPATH += $$VCGLIBDIR/eigenlib
INCLUDEPATH += $$GLEWDIR/include

# compile glew
DEFINES += GLEW_STATIC
SOURCES += $$GLEWDIR/src/glew.c

# compile vcg
SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackball.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackmode.cpp



# include libigl and CoMISo
DEFINES += INCLUDE_TEMPLATES
unix|win32: LIBS += -L$$PWD/external/build_CoMISo/ -lCoMISo
win32: {
    DEFINES += _USE_MATH_DEFINES
    DEFINES += NOMINMAX
    QMAKE_CXXFLAGS += _SCL_SECURE_NO_DEPRECAT
    LIBS += -L$$PWD/external/libigl/external/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/lib -lopenblas.dll.a.lib
}
unix: LIBS += -lBLAS
macos: LIBS += -framework Accelerate

INCLUDEPATH += $$PWD/external/libigl/external
DEPENDPATH += $$PWD/external/libigl/external
INCLUDEPATH += $$PWD/external/libigl/external/CoMISo/ext/gmm-4.2/include
DEPENDPATH += $$PWD/external/libigl/external/CoMISo/ext/gmm-4.2/include

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/external/build_CoMISo/CoMISo.lib
else:unix|win32-g++: PRE_TARGETDEPS += $$PWD/external/build_CoMISo/libCoMISo.a

INCLUDEPATH += $$PWD/external/libigl/include
DEPENDPATH += $$PWD/external/libigl/include
